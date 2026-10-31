#ifndef DIPOLE_H
#define DIPOLE_H
#include "fdtd.h"

typedef struct {
	size_t	*wire;
	int	nwire;
	size_t	gap;
	float	beta;
	float	ez_old;
	double dx;
} dipole_t;

dipole_t create_dipole(const grid_t *g, int cx, int cy, int cz, int half, \
		       float dt, float Rs, double dx, double eps0);
void destroy_dipole(dipole_t *d);
void dipole_pre_E(const grid_t *g, dipole_t *d);
float dipole_post_E(grid_t *g, dipole_t *d, float Vs);

#endif
