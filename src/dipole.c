#include "dipole.h"
#include <stdio.h>

dipole_t create_dipole(const grid_t *g, int cx, int cy, int cz, int half, \
		       float dt, float Rs, double dx, double eps0) {
	dipole_t d = { 0 };
	d.wire = malloc(2*(size_t)half * sizeof(size_t));
	if (!d.wire) {
		fprintf(stderr, "Allocation failed\n");
		exit(1);
	}
	
	size_t base = (size_t)cx * g->sx + (size_t)cy * g->sy;
	for (int m = 1; m <= half; m++) {
		d.wire[d.nwire++] = base + (cz + m);
		d.wire[d.nwire++] = base + (cz - m);
	}
	
	d.gap = base + cz;
	d.beta = (float)(dt / (2.0 * Rs * eps0 * dx));
	d.dx = (float)dx;
	
	return d;
}

void destroy_dipole(dipole_t *d) {
	free(d->wire);
	d->wire = NULL;
}

void dipole_pre_E(const grid_t *g, dipole_t *d) {
	d->ez_old = g->Ez[d->gap];
}

float dipole_post_E(grid_t *g, dipole_t *d, float Vs) {
	for (int n = 0; n < d->nwire; n++) {
		g->Ez[d->wire[n]] = 0.0f;
	}
	
	float Es = Vs / d->dx;
	float *e = &g->Ez[d->gap];
	*e = (*e - d->beta * d->ez_old + 2.0f * d->beta * Es) / (1.0f + d->beta);
	return *e * d->dx;
}
