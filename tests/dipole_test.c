#include "fdtd.h"
#include "cpml.h"
#include "dipole.h"
#include "units.h"
#include "cli.h"

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <limits.h>
#include <float.h>
#include <string.h>

#define NPML 10
#define MARGIN 8

typedef struct {
	int	N, half, steps;
	double	dx, Rs, fmax;
	const char *out;
} cfg_t;

static const opt_t OPTS[] = {
	OPT(cfg_t, N,		1, INT_MAX,	"grid cells per side"),
	OPT(cfg_t, half,	1, INT_MAX,	"dipole arm length, cells"),
	OPT(cfg_t, steps,	1, INT_MAX,	"time steps"),
	OPT(cfg_t, dx,		DBL_MIN, 1.0,	"cell size, m"),
	OPT(cfg_t, Rs,		DBL_MIN, 1e6,	"source resistance, ohms"),
	OPT(cfg_t, fmax,	DBL_MIN, 1e12,	"highest frequency of interest, Hz"),
	OPT(cfg_t, out,		0, 0,		"output file name")
};
#define NOPTS (sizeof OPTS / sizeof OPTS[0])

int main(int argc, char **argv) {
	cfg_t cfg = { 
		.N = 64,
		.half = 9,
		.steps = 3000,
		.dx = 5e-3,
		.Rs = 50.0,
		.fmax = 6e9,
		.out = "dipole_vi.bin"
	};
	
	int rc = cli_parse(argc, argv, &cfg, OPTS, NOPTS);
	if (rc) {
		return rc == 2 ? 0 : 1;
	}
	
	double cpw = C0 / (cfg.fmax * cfg.dx);
	if (cpw < 10) {
		fprintf(stderr, "dx too coarse: %0.1f cells/lambda\n", cpw);
		return 1;
	} else if (cpw < 20) {
		fprintf(stderr, "warning: only %0.1f cells/lambda at fmax\n", cpw);
	}
	
	int need = 2 * cfg.half + 1 + 2 * (NPML + MARGIN);
	if (cfg.N < need) {
		fprintf(stderr, "N=%d too small, need >= %d\n", cfg.N, need);
		return 1;
	}
	
	int status = 1;

	double	S = 0.9 / sqrt(3.0);
	double	dt = S * cfg.dx / C0;
	float	k = (float)(C0 * dt / cfg.dx);
	
	double	tau = 0.5 / cfg.fmax;
	double	t0 = 4.0 * tau;
	int	c = cfg.N / 2;

	pml_params_t pp = {
		.npml = NPML,
		.m = 3,
		.sigma_max = 0.8 * (3 + 1) / (ETA0 * cfg.dx),
		.kappa_max = 5.0,
		.alpha_max = 0.1 / (ETA0 * cfg.dx),
		.eps0 = EPS0
	};

	grid_t	g = create_grid(cfg.N);
	cpml_t	pml = create_cpml(&g, k, dt, &pp);
	dipole_t ant = create_dipole(&g, c, c, c, cfg.half, dt, cfg.Rs, cfg.dx, EPS0);
	double	Vprev = 0.0;
	
	double *rec = malloc((size_t)cfg.steps * 3 * sizeof(double));
	if (!rec) {
		fprintf(stderr, "Allocation failed\n");
		goto done;
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
		
		if (!isfinite(V)) {
			fprintf(stderr, "non-finite V at step %d\n", i);
			status = 2;
			break;
		}
		
		rec[3*i + 0] = time;
		rec[3*i + 1] = Vs;
		rec[3*i + 2] = 0.5f * (Vprev + V);
		Vprev = V;
	}
	
	FILE *fp = fopen(cfg.out, "wb");
	if (!fp) {
		perror(cfg.out);
		goto done;
	}
	size_t n = (size_t)cfg.steps * 3;
	if (fwrite(rec, sizeof(double), n, fp) != n) {
		fprintf(stderr, "Write failed\n");
		goto done;
	}
	
	status = 0;
	
	done:
		if (fp) {
			fclose(fp);
		}
		free(rec);
		destroy_dipole(&ant); destroy_cpml(&pml); destroy_grid(&g);
		return status;
}
