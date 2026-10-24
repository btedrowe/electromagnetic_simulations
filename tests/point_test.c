#include "fdtd.h"
#include "cpml.h"

#include <stdio.h>
#include <math.h>

int main() {
    int time_steps = 90;
    int N = 50;

    float dx = 10.0;
    float c0 = 1.0;
    float S = 0.90 / sqrtf(3);
    float dt = S * dx / c0;
    float k = dt / dx;

    int t0 = 25;
    float tau = 8.0;

    int f = N / 2;

    pml_params_t pp = {
	    .npml = 10,
	    .m = 3,
	    .kappa_max = 5.0,
	    .alpha_max = 0.1,
	    .sigma_max = 0.8 * (3 + 1) / dx
	};

    grid_t g = create_grid(N);
    cpml_t pml = create_cpml(&g, k, (double)dt, &pp);

    for (int i = 0; i < time_steps; i++) {
		update_H(&g, k);
        cpml_H(&g, &pml);
		update_E(&g, k);
        cpml_E(&g, &pml);

		float u = (float)(i - t0) / tau;

		g.Ez[f*g.sx + f*g.sy + f] += -2.0f * u * expf(-u*u);
    }
    
    FILE *fp = fopen("ez_c.bin", "wb");
    for (int i = 0; i < N; i++) {
		for (int j = 0; j < N; j++) {
			fwrite(&g.Ez[i*g.sx + j*g.sy + f], sizeof(float), 1, fp);
		}
    }
    fclose(fp);

    destroy_cpml(&pml);
    destroy_grid(&g);
    
    return 0;
}