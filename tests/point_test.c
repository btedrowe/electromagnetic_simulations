#include "fdtd.h"
#include "cpml.h"

#include <stdio.h>
#include <math.h>
#include <string.h>

#define C0 299792458.0
#define MU0 1.25663706212e-6
#define EPS0 (1.0 / (MU0 * C0 * C0))
#define ETA0 (MU0 * C0)

#define NPML 10
#define EXTENT 1

typedef struct {
	int	N, steps;
	double	dx, Rs;
	const char *out;
} cfg_t;

static int parse_args(int argc, char **argv, cfg_t *c) {
	*c = (cfg_t) {
		.N = 64,
		.steps = 3000,
		.dx = 5e-3,
		.Rs = 50.0,
		.out = "point_ez.bin"
	};
	
	for (int i = 1; i < argc; i++) {
		char *eq = strchr(argv[i], '=');
		if (!eq) {
			fprintf(stderr, "bag arg: %s (expected key=value)\n", argv[i]);
			return 1;
		}
		*eq = '\0';
		const char *key = argv[i], *val = eq + 1;
		
		if (!strcmp(key, "N")) c->N = atoi(val);
		else if (!strcmp(key, "steps")) c->steps = atoi(val);
		else if (!strcmp(key, "dx")) c->dx = atof(val);
		else if (!strcmp(key, "Rs")) c->Rs = atof(val);
		else if (!strcmp(key, "out")) c->out = val;
		else {
			fprintf(stderr, "unknown key: %s\n", key);
			return 1;
		}
	}

	int need = EXTENT + 2 * (NPML + 1);
	if (c->N < need) {
		fprintf(stderr, "N=%d too small, need >= %d\n", c->N, need);
		return 1;
	}

	return 0;
}

int main(int argc, char **argv) {
	cfg_t cfg = { 0 };
	if (parse_args(argc, argv, &cfg)) {
		return 1;
	}
	
	float S = 0.90 / sqrtf(3);
	float dt = S * cfg.dx / C0;
	float k = dt / cfg.dx;

	int t0 = 25;
	float tau = 8.0;

	int f = cfg.N / 2;

	pml_params_t pp = {
		.npml = NPML,
		.m = 3,
		.sigma_max = 0.8 * (3 + 1) / cfg.dx,
		.kappa_max = 5.0,
		.alpha_max = 0.1 / (ETA0 * cfg.dx),
		.eps0 = EPS0
	};

	grid_t g = create_grid(cfg.N);
	cpml_t pml = create_cpml(&g, k, (double)dt, &pp);

	for (int i = 0; i < cfg.steps; i++) {
		update_H(&g, k);
		cpml_H(&g, &pml);
		update_E(&g, k);
		cpml_E(&g, &pml);

		float u = (float)(i - t0) / tau;

		g.Ez[f*g.sx + f*g.sy + f] += -2.0f * u * expf(-u*u);
	}

	FILE *fp = fopen("ez_c.bin", "wb");
	for (int i = 0; i < cfg.N; i++) {
		for (int j = 0; j < cfg.N; j++) {
			fwrite(&g.Ez[i*g.sx + j*g.sy + f], sizeof(float), 1, fp);
		}
	}
	fclose(fp);

	destroy_cpml(&pml);
	destroy_grid(&g);

	return 0;
}
