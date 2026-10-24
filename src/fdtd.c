#include "fdtd.h"

#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <stdio.h>

grid_t create_grid(int N) {
    grid_t g = { 0 };

    g.nx = N; g.ny = N; g.nz = N;
    g.nzp = (g.nz + 15) & ~15;
    g.sx = (size_t)g.ny * g.nzp;
    g.sy = g.nzp;

    size_t bytes = (size_t)g.nx * g.ny * g.nzp * sizeof(float);
    float **f[6] = { &g.Ex, &g.Ey, &g.Ez, &g.Hx, &g.Hy, &g.Hz };

    for (int i = 0; i < 6; i++) {
      *f[i] = aligned_alloc(64, bytes);
      if (!*f[i]) {
        fprintf(stderr, "Allocation failed\n");
        exit(1);
      }

      memset(*f[i], 0, bytes);
    }

    return g;
}

void destroy_grid(grid_t *g) {
    free(g->Ex);
    free(g->Ey);
    free(g->Ez);
    free(g->Hx);
    free(g->Hy);
    free(g->Hz);
}

void update_H(grid_t *g, float k) {
    const size_t sx = g->sx, sy = g->sy;
    float *restrict Ex = g->Ex, *restrict Ey = g->Ey, *restrict Ez = g->Ez;
    float *restrict Hx = g->Hx, *restrict Hy = g->Hy, *restrict Hz = g->Hz;

    #pragma omp parallel for schedule(static)
    for (int i = 0; i < g->nx - 1; i++) {
		for (int j = 0; j < g->ny - 1; j++) {
			size_t base = i*sx + j*sy;
			#pragma omp simd
			for (int z = 0; z < g->nz - 1; z++) {
				size_t c = base + z;

				Hx[c] += k * ((Ey[c+1] - Ey[c]) - (Ez[c+sy] - Ez[c]));
				Hy[c] += k * ((Ez[c+sx] - Ez[c]) - (Ex[c+1] - Ex[c]));
				Hz[c] += k * ((Ex[c+sy] - Ex[c]) - (Ey[c+sx] - Ey[c]));
			}
		}
    }
}

void update_E(grid_t *g, float k) {
    const size_t sx = g->sx, sy = g->sy;
    float *restrict Ex = g->Ex, *restrict Ey = g->Ey, *restrict Ez = g->Ez;
    float *restrict Hx = g->Hx, *restrict Hy = g->Hy, *restrict Hz = g->Hz;

    #pragma omp parallel for schedule(static)
    for (int i = 1; i < g->nx; i++) {
		for (int j = 1; j < g->ny; j++) {
			size_t base = i*sx + j*sy;
			#pragma omp simd
			for (int z = 1; z < g->nz; z++) {
				size_t c = base + z;

				Ex[c] += k * ((Hz[c] - Hz[c-sy]) - (Hy[c] - Hy[c-1]));
				Ey[c] += k * ((Hx[c] - Hx[c-1]) - (Hz[c] - Hz[c-sx]));
				Ez[c] += k * ((Hy[c] - Hy[c-sx]) - (Hx[c] - Hx[c-sy]));
			}
		}
    }
}
