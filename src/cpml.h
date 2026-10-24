#ifndef CPML_H
#define CPML_H

#include <stddef.h>
#include "fdtd.h"

typedef struct {
	int npml, m;
	double sigma_max, kappa_max, alpha_max;
} pml_params_t;

typedef struct {
	int n;
	float *b, *a, *ik1;
} pml_profile_t;

typedef struct {
    float *tgt;
    const float *src;
    int ax, fwd;
    int lo[3], hi[3];
    size_t stride;
    float c;
    const float *b, *a, *ik1;
    float *psi;
} cpml_slab_t;

typedef struct {
    pml_profile_t Hprof[3], Eprof[3];
    cpml_slab_t sH[12], sE[12];
    int nH, nE;
} cpml_t;

pml_profile_t create_profile(int n, double offset, double dt, const pml_params_t *p);
void destroy_profile(pml_profile_t *pr);

cpml_t create_cpml(const grid_t *g, float k, double dt, const pml_params_t *pp);
void destroy_cpml(cpml_t *p);

void cpml_H(grid_t *g, cpml_t *p);
void cpml_E(grid_t *g, cpml_t *p);

#endif