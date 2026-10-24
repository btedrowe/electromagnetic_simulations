#include "cpml.h"

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

pml_profile_t create_profile(int n, double offset, double dt, const pml_params_t *p) {
	pml_profile_t pr = { 0 };
	pr.n = n;
	pr.b = malloc(n*sizeof(float));
	pr.a = malloc(n*sizeof(float));
	pr.ik1 = malloc(n*sizeof(float));

	if (!pr.b || !pr.a || !pr.ik1) {
		fprintf(stderr, "Allocation failed\n");
		exit(1);
	}

	for (int i = 0; i < n; i++) {
		double pos = i + offset;
		double depth = fmax(p->npml - pos, pos - (n - 1 - p->npml));
		depth = fmin(fmax(depth, 0.0), (double)p->npml) / p->npml;

		double dm = pow(depth, p->m);
		double sigma = p->sigma_max * dm;
		double kappa = 1.0 + (p->kappa_max - 1.0) * dm;
		double alpha = p->alpha_max * (1.0 - depth);

		double b = exp(-(sigma / kappa + alpha) * dt);
		double a = 0.0;
		if (sigma > 0.0) {
			a = sigma * (b - 1.0) / (kappa * (sigma + kappa * alpha));
		}

		pr.b[i] = (float)b;
		pr.a[i] = (float)a;
		pr.ik1[i] = (float)(1.0 / kappa - 1.0);
	}

	return pr;
}

void destroy_profile(pml_profile_t *pr) {
	free(pr->b);
	free(pr->a);
	free(pr->ik1);
}

static void add_term(
    cpml_slab_t *out, int *count, const grid_t *g,
    float *tgt, const float *src, int ax, int sign,
    const pml_profile_t *pr, int fwd, float k
) {
    const int n[3] = { g->nx, g->ny, g->nz };
    const size_t stride[3] = { g->sx, g->sy, 1 };

    int vlo[3], vhi[3];
    for (int d = 0; d < 3; d++) {
        vlo[d] = fwd ? 0 : 1;
        vhi[d] = fwd ? n[d] - 1 : n[d];
    }

    int L = vlo[ax];
    while (L < vhi[ax] && pr->a[L] != 0.0f) { L++; }
    int H = vhi[ax];
    while (H > vlo[ax] && pr->a[H - 1] != 0.0f) { H--; }

    const int s0[2] = { vlo[ax], H }, s1[2] = { L, vhi[ax] };
    for (int s = 0; s < 2; s++) {
        if (s1[s] <= s0[s]) { continue; }

        cpml_slab_t *sl = &out[(*count)++];
        sl->tgt = tgt; sl->src = src; sl->ax = ax; sl->fwd = fwd;
        sl->stride = stride[ax]; sl->c = sign * k;
        sl->b = pr->b; sl->a = pr->a; sl->ik1 = pr->ik1;

        for (int d = 0; d < 3; d++) {
            sl->lo[d] = vlo[d];
            sl->hi[d] = vhi[d];
        }
        sl->lo[ax] = s0[s]; sl->hi[ax] = s1[s];

        size_t cnt = (size_t)(sl->hi[0] - sl->lo[0]) * \
                     (sl->hi[1] - sl->lo[1]) * \
                     (sl->hi[2] - sl->lo[2]);
        sl->psi = calloc(cnt, sizeof(float));

        if (!sl->psi) {
            fprintf(stderr, "Allocation failed\n");
            exit(1);
        }
    }
}

cpml_t create_cpml(const grid_t *g, float k, double dt, const pml_params_t *pp) {
    cpml_t p = { 0 };
    
    const int n[3] = { g->nx, g->ny, g->nz };
    for (int i = 0; i < 3; i++) {
        p.Hprof[i] = create_profile(n[i], 0.5, dt, pp);
        p.Eprof[i] = create_profile(n[i], 0.0, dt, pp);
    }

    #define H_(t, s, ax, sg) add_term(p.sH, &p.nH, g, g->t, g->s, ax, sg, &p.Hprof[ax], 1, k)
    #define E_(t, s, ax, sg) add_term(p.sE, &p.nE, g, g->t, g->s, ax, sg, &p.Eprof[ax], 0, k)
    H_(Hx, Ey, 2, +1); H_(Hx, Ez, 1, -1); H_(Hy, Ez, 0, +1);
    H_(Hy, Ex, 2, -1); H_(Hz, Ex, 1, +1); H_(Hz, Ey, 0, -1);
    E_(Ex, Hz, 1, +1); E_(Ex, Hy, 2, -1); E_(Ey, Hx, 2, +1);
    E_(Ey, Hz, 0, -1); E_(Ez, Hy, 0, +1); E_(Ez, Hx, 1, -1);
    #undef H_
    #undef E_

    return p;
}

void destroy_cpml(cpml_t *p) {
    for (int i = 0; i < p->nH; i++) {
        free(p->sH[i].psi);
    }
    for (int i = 0; i < p->nE; i++) {
        free(p->sE[i].psi);
    }
    for (int i = 0; i < 3; i++) {
        destroy_profile(&p->Hprof[i]);
        destroy_profile(&p->Eprof[i]);
    }
}

static void cpml_apply(const cpml_slab_t *s, size_t sx, size_t sy) {
    const int ax = s->ax, fwd = s->fwd;
    const int nj = s->hi[1] - s->lo[1], nk = s->hi[2] - s->lo[2];
    const size_t st = s->stride;
    const float c = s->c;
    float *restrict tgt = s->tgt;
    const float *restrict src = s->src;
    float *restrict psi = s->psi;

    #pragma omp parallel for schedule(static)
    for (int i = s->lo[0]; i < s->hi[0]; i++) {
        for (int j = s->lo[1]; j < s->hi[1]; j++) {
            for (int z = s->lo[2]; z < s->hi[2]; z++) {
                size_t cell = i*sx + j*sy + z;
                size_t q = ((size_t)(i - s->lo[0]) * nj + (j - s->lo[1])) * nk + (z - s->lo[2]);
                int ia = (ax == 0) ? i : (ax == 1) ? j : z;

                float d = fwd ? src[cell + st] - src[cell]
                              : src[cell] - src[cell - st];
                psi[q] = s->b[ia] * psi[q] + s->a[ia] * d;
                tgt[cell] += c * (s->ik1[ia] * d + psi[q]);
            }
        }
    }
}

void cpml_H(grid_t *g, cpml_t *p) {
    for (int i = 0; i < p->nH; i++) {
        cpml_apply(&p->sH[i], g->sx, g->sy);
    }

}

void cpml_E(grid_t *g, cpml_t *p) {
    for (int i = 0; i < p->nE; i++) {
        cpml_apply(&p->sE[i], g->sx, g->sy);
    }
}