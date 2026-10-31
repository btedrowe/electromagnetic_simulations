#include "fdtd.h"
#include "cpml.h"
#include "dipole.h"

#include <stdio.h>
#include <math.h>
#include <string.h>

#define C0 299792458.0
#define MU0 1.25663706212e-6
#define EPS0 (1.0 / (MU0 * C0 * C0))
#define ETA0 (MU0 * C0)

typedef struct {
	int	N, half, steps;
	double	dx, Rs;
	const char *out;
} cfg_t;

static int parse_args(int argc, char **argv, cfg_t *c) {
	*c = (cfg_t) {
		.N = 64,
		.half = 9,
		.steps = 3000,
		.dx = 5e-3,
		.Rs = 50.0,
		.out = "dipole_vi.bin"
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
		else if (!strcmp(key, "half")) c->half = atoi(val);
		else if (!strcmp(key, "steps")) c->steps = atoi(val);
		else if (!strcmp(key, "dx")) c->dx = atof(val);
		else if (!strcmp(key, "Rs")) c->Rs = atof(val);
		else if (!strcmp(key, "out")) c->out = val;
		else {
			fprintf(stderr, "unknown key: %s\n", key);
			return 1;
		}
	}
	
	int need = 2 * c->half + 1 + 2 * (10 + 8);
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

	double	S = 0.9f / sqrtf(3.0f);
	double	dt = S * cfg.dx / C0;
	float	k = (float)(C0 * dt / cfg.dx);
	
	double	tau = 10.0 * cfg.dx / C0;
	double	t0 = 4.0 * tau;
	int	c = cfg.N / 2;

	pml_params_t pp = {
		.npml = 10,
		.m = 3,
		.sigma_max = 0.8 * (3 + 1) / (ETA0 * cfg.dx),
		.kappa_max = 5.0,
		.alpha_max = 0.1 / (ETA0 * cfg.dx),
		.eps0 = EPS0
	};

	grid_t	g = create_grid(cfg.N);
	cpml_t	pml = create_cpml(&g, k, dt, &pp);
	dipole_t ant = create_dipole(&g, c, c, c, cfg.half, dt, cfg.Rs, cfg.dx, EPS0);
	float	Vprev = 0.0f;
	
	double *rec = malloc((size_t)cfg.steps * 3 * sizeof(double));
	if (!rec) {
		fprintf(stderr, "Allocation failed\n");
		return 1;
	}

	for (int i = 0; i < cfg.steps; i++) {
		update_H(&g, k);
		cpml_H(&g, &pml);
		
		dipole_pre_E(&g, &ant);
		update_E(&g, k);
		cpml_E(&g, &pml);

		double time = (i + 0.5) * dt;
		double u = (time - t0) / tau;
		double Vs = exp(-u * u);
		double V = dipole_post_E(&g, &ant, (float)Vs);
		
		rec[3*i + 0] = time;
		rec[3*i + 1] = Vs;
		rec[3*i + 2] = 0.5f * (Vprev + V);
		Vprev = V;
	}
	
	FILE *fp = fopen(cfg.out, "wb");
	if (!fp) {
		perror("fopen");
		return 1;
	}
	size_t n = (size_t)cfg.steps * 3;
	if (fwrite(rec, sizeof(double), n, fp) != n) {
		fprintf(stderr, "Write failed\n");
		return 1;
	}
	fclose(fp);

	free(rec);
	destroy_dipole(&ant);
	destroy_cpml(&pml);
	destroy_grid(&g);

	return 0;
}
