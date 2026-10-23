import numpy as np

N = 50
time_steps = 90

c0, dx = 1.0, 10
S = 0.90 / np.sqrt(3)
dt = S * dx / c0
k = dt / dx

fx = fy = fz = N // 2
t0, tau = 25, 8

def get_data():
    Ex, Ey, Ez, Hx, Hy, Hz = (np.zeros((N, N, N)) for _ in range(6))

    for t in range(time_steps):
        Hx[:-1, :-1, :-1] += k * ((Ey[:-1, :-1, 1:] - Ey[:-1, :-1, :-1]) - (Ez[:-1, 1:, :-1] - Ez[:-1, :-1, :-1]))
        Hy[:-1, :-1, :-1] += k * ((Ez[1:, :-1, :-1] - Ez[:-1, :-1, :-1]) - (Ex[:-1, :-1, 1:] - Ex[:-1, :-1, :-1]))
        Hz[:-1, :-1, :-1] += k * ((Ex[:-1, 1:, :-1] - Ex[:-1, :-1, :-1]) - (Ey[1:, :-1, :-1] - Ey[:-1, :-1, :-1]))

        Ex[1:, 1:, 1:] += k * ((Hz[1:, 1:, 1:] - Hz[1:, :-1, 1:]) - (Hy[1:, 1:, 1:] - Hy[1:, 1:, :-1]))
        Ey[1:, 1:, 1:] += k * ((Hx[1:, 1:, 1:] - Hx[1:, 1:, :-1]) - (Hz[1:, 1:, 1:] - Hz[:-1, 1:, 1:]))
        Ez[1:, 1:, 1:] += k * ((Hy[1:, 1:, 1:] - Hy[:-1, 1:, 1:]) - (Hx[1:, 1:, 1:] - Hx[1:, :-1, 1:]))

        u = (t - t0) / tau
        Ez[fx, fy, fz] += -2 * u * np.exp(-u**2)

    slice_idx = N // 2
    data = Ez[:, :, slice_idx].T
    vmax = np.abs(data).max()

    return slice_idx, data, vmax
