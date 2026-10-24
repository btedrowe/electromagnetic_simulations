import numpy as np

N = 50
npml = 10
time_steps = 90

c0, dx = 1.0, 10
dt = 0.9 * dx / (c0 * np.sqrt(3))
k = dt / dx

m, kappa_max, alpha_max = 3, 5.0, 0.1
sigma_max = 0.8 * (m + 1) / dx

def profile(offset):
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

    return b, a, 1 / kappa - 1

Hprof, Eprof = profile(0.5), profile(0.0)

def sl(ax, i0, i1):
    s = [slice(None)] * 3
    s[ax] = slice(i0, i1)

    return tuple(s)

def bshape(v, ax):
    s = [1, 1, 1]
    s[ax] = len(v)
    return v.reshape(s).astype(np.float32)

class Term:
    def __init__(self, target, src, ax, sign, prof, fwd):
        self.t, self.s, self.ax, self.fwd = target, src, ax, fwd
        self.c = sign * k

        self.full = (0, N - 1) if fwd else (1, N)
        ranges = [(0, npml), (N - npml, N - 1)] if fwd else [(1, npml), (N - npml, N)]
        b, a, ik1 = prof
        self.slabs = []

        for i0, i1 in ranges:
            shape = list(src.shape); shape[ax] = i1 - i0
            self.slabs.append((
                i0, i1,
                bshape(b[i0:i1], ax),
                bshape(a[i0:i1], ax),
                bshape(ik1[i0:i1], ax),
                np.zeros(shape, np.float32)
            ))

    def diff(self, i0, i1):
        ax, F = self.ax, self.s
        if self.fwd:
            return F[sl(ax, i0 + 1, i1 + 1)] - F[sl(ax, i0, i1)]
        return F[sl(ax, i0, i1)] - F[sl(ax, i0 - 1, i1 - 1)]

    def step(self):
        ax = self.ax
        self.t[sl(ax, *self.full)] += self.c * self.diff(*self.full)
        for i0, i1, b, a, ik1, psi in self.slabs:
            d = self.diff(i0, i1) / dx
            psi *= b
            psi += a * d
            self.t[sl(ax, i0, i1)] += self.c * dx * (ik1 * d + psi)


def get_data():
    Ex, Ey, Ez, Hx, Hy, Hz = (np.zeros((N, N, N)) for _ in range(6))

    TH = [
        Term(Hx, Ey, 2, +1, Hprof, True), Term(Hx, Ez, 1, -1, Hprof, True),
        Term(Hy, Ez, 0, +1, Hprof, True), Term(Hy, Ex, 2, -1, Hprof, True),
        Term(Hz, Ex, 1, +1, Hprof, True), Term(Hz, Ey, 0, -1, Hprof, True)
    ]
    TE = [
        Term(Ex, Hz, 1, +1, Eprof, False), Term(Ex, Hy, 2, -1, Eprof, False),
        Term(Ey, Hx, 2, +1, Eprof, False), Term(Ey, Hz, 0, -1, Eprof, False),
        Term(Ez, Hy, 0, +1, Eprof, False), Term(Ez, Hx, 1, -1, Eprof, False)
    ]
    fx = fy = fz = N // 2
    t0, tau = 25, 8

    for t in range(time_steps):
        for term in TH: term.step()
        for term in TE: term.step()

        u = (t - t0) / tau
        Ez[fx, fy, fz] += -2 * u * np.exp(-u**2)

    slice_idx = N // 2
    data = Ez[:, :, slice_idx].T
    vmax = np.abs(data).max()

    return slice_idx, data, vmax