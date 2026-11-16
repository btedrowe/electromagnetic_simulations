#include "fdtd.h"
#include "cpml.h"
#include "units.h"
#include "cli.h"

#include <stdio.h>
#include <math.h>
#include <string.h>
#include <limits.h>
#include <float.h>

#define NPML 10
#define MARGIN 10
#define EXTENT 1

typedef struct {
	int	N, steps;
	double	dx, Rs;
	const char *out;
} cfg_t;

static const opt_t OPTS[] = {
	OPT(cfg_t, N,		1, INT_MAX,	"grid cells per side"),
	OPT(cfg_t, steps,	1, INT_MAX,	"time steps"),
	OPT(cfg_t, dx,		DBL_MIN, 1.0,	"cell size, m"),
	OPT(cfg_t, Rs,		DBL_MIN, 1e6,	"source resistance, ohms"),
	OPT(cfg_t, out,		0, 0,		"output file name")
};
#define NOPTS (sizeof OPTS / sizeof OPTS[0])

int main(int argc, char **argv) {
	cfg_t cfg = { 
		.N = 64,
		.steps = 3000,
		.dx = 5e-3,
		.Rs = 50.0,
		.out = "point_ez.bin"
	};
	
	int rc = cli_parse(argc, argv, &cfg, OPTS, NOPTS);
	if (rc) {
		return rc == 2 ? 0 : 1;
	}
	
	int need = EXTENT + 2 * (NPML + MARGIN);
	if (cfg.N < need) {
		fprintf(stderr, "N=%d too small, need >= %d\n", cfg.N, need);
		return 1;
	}
	
	float S = 0.90 / sqrtf(3);
	float dt = S * cfg.dx / C0;
	float k = (float)(C0 * dt / cfg.dx);

	int t0 = 25;
	float tau = 8.0;

	int f = cfg.N / 2;

	pml_params_t pp = {
		.npml = NPML,
		.m = 3,
		.sigma_max = 0.8 * (3 + 1) / (ETA0 * cfg.dx),
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

	FILE *fp = fopen(cfg.out, "wb");
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
