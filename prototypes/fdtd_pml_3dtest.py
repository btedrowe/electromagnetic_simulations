import numpy as np

N = 50
npml = 10
time_steps = 90

c0, dx = 1.0, 10
dt = 0.9 * dx / (c0 * np.sqrt(3))

m, kappa_max, alpha_max = 3, 5.0, 0.1
sigma_max = 0.8 * (m + 1) / dx

def profile(offset, axis):
    pos = np.arange(N) + offset

    depth = np.clip(
        np.maximum(
            npml - pos,
            pos - (N - 1 - npml)
        ),
        0,
        npml
    ) / npml
    sigma = sigma_max * depth**m
    kappa = 1 + (kappa_max - 1) * depth**m
    alpha = alpha_max * (1 - depth)

    b = np.exp(-(sigma / kappa + alpha) * dt)
    a = np.zeros(N)

    nz = sigma > 0
    a[nz] = sigma[nz] * (b[nz] - 1) / (kappa[nz] * (sigma[nz] + kappa[nz] * alpha[nz]))
    shape = [1, 1, 1]; shape[axis] = N

    return b.reshape(shape), a.reshape(shape), (1 / kappa).reshape(shape)

Hp = [profile(0.5, ax) for ax in range(3)]
Ep = [profile(0.0, ax) for ax in range(3)]

def dfwd(F, ax):
    out = np.zeros_like(F)
    sl = [slice(None)] * 3; sl[ax] = slice(0, -1)
    out[tuple(sl)] = np.diff(F, axis=ax)
    return out

def dbwd(F, ax):
    out = np.zeros_like(F)
    sl = [slice(None)] * 3; sl[ax] = slice(1, None)
    out[tuple(sl)] = np.diff(F, axis=ax)
    return out

psi = {}
def T(F, name, ax, prof, diff):
    b, a, invk = prof[ax]
    d = diff(F, ax) / dx
    p = psi.setdefault((name, ax), np.zeros_like(F))
    p *= b
    p += a * d
    return d * invk + p

def get_data():
    Ex, Ey, Ez, Hx, Hy, Hz = (np.zeros((N, N, N)) for _ in range(6))

    fx = fy = fz = N // 2
    t0, tau = 25, 8

    for t in range(time_steps):
        Hx += dt * (T(Ey, 'Ey', 2, Hp, dfwd) - T(Ez, 'Ez', 1, Hp, dfwd))
        Hy += dt * (T(Ez, 'Ez', 0, Hp, dfwd) - T(Ex, 'Ex', 2, Hp, dfwd))
        Hz += dt * (T(Ex, 'Ex', 1, Hp, dfwd) - T(Ey, 'Ey', 0, Hp, dfwd))

        Ex += dt * (T(Hz, 'Hz', 1, Ep, dbwd) - T(Hy, 'Hy', 2, Ep, dbwd))
        Ey += dt * (T(Hx, 'Hx', 2, Ep, dbwd) - T(Hz, 'Hz', 0, Ep, dbwd))
        Ez += dt * (T(Hy, 'Hy', 0, Ep, dbwd) - T(Hx, 'Hx', 1, Ep, dbwd))

        u = (t - t0) / tau
        Ez[fx, fy, fz] += -2 * u * np.exp(-u**2)

    slice_idx = N // 2
    data = Ez[:, :, slice_idx].T
    vmax = np.abs(data).max()

    return slice_idx, data, vmax