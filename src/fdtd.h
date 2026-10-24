#ifndef FDTD_H
#define FDTD_H

#include <stdlib.h>

typedef struct {
    int nx, ny, nz, nzp;
    size_t sx, sy;
    float *Ex, *Ey, *Ez, *Hx, *Hy, *Hz;
} grid_t;

grid_t create_grid(int N);
void destroy_grid(grid_t *g);

void update_H(grid_t *g, float k);
void update_E(grid_t *g, float k);

#endif