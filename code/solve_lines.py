"""Step 2: find all lines in the pencil f - lam*Q by total-degree homotopy continuation.

Chart: in a random real basis g of R^4, S = span(g0 + a g2 + b g3, g1 + c g2 + d g3).
Unknowns z = (a, b, c, d, lam). Equations: the 5 coefficients of the binary quartic
(f - lam Q)(u*(g0 + a g2 + b g3) + r*(g1 + c g2 + d g3)).
"""
import numpy as np, sympy as sp, pickle, time
from pathlib import Path
HERE = Path(__file__).resolve().parent
BUILD = HERE / "build"
BUILD.mkdir(exist_ok=True)


rng = np.random.default_rng(5)
D = pickle.load(open(BUILD / 'pencil.pkl', 'rb'))
X = sp.symbols('x0:4')
a, b, c, d, lam, u, r = sp.symbols('a b c d lam u r')
Z = [a, b, c, d, lam]

G = np.linalg.qr(rng.standard_normal((4, 4)))[0]  # random orthogonal basis (columns)
vecL = [sp.Matrix(G[:, 0]) + a * sp.Matrix(G[:, 2]) + b * sp.Matrix(G[:, 3]),
        sp.Matrix(G[:, 1]) + c * sp.Matrix(G[:, 2]) + d * sp.Matrix(G[:, 3])]
pt = u * vecL[0] + r * vecL[1]
sub = {X[i]: pt[i] for i in range(4)}
fr = sp.Poly(sp.expand(D['f'].subs(sub)), u, r)
Qr = sp.Poly(sp.expand(D['Q'].subs(sub)), u, r)
eqs = [sp.expand(fr.coeff_monomial(u**(4 - k) * r**k) - lam * Qr.coeff_monomial(u**(4 - k) * r**k))
       for k in range(5)]
T = sp.lambdify(Z, eqs, 'numpy')
JT = sp.lambdify(Z, [[sp.diff(e, v) for v in Z] for e in eqs], 'numpy')


def Tval(z):
    return np.array(T(*z.T), dtype=complex).T  # (N,5)


def Jval(z):
    J = JT(*z.T)
    N = z.shape[0]
    out = np.empty((N, 5, 5), dtype=complex)
    for i in range(5):
        for j in range(5):
            out[:, i, j] = J[i][j]
    return out



def ssolve(A, B):
    try:
        return np.linalg.solve(A, B)
    except np.linalg.LinAlgError:
        return np.linalg.pinv(A) @ B

deg = 5
gamma = None


def H(z, t):
    S = z**deg - 1
    return (1 - t)[:, None] * gamma * S + t[:, None] * Tval(z)


def Hz(z, t):
    JS = np.zeros((z.shape[0], 5, 5), dtype=complex)
    idx = np.arange(5)
    JS[:, idx, idx] = deg * z**(deg - 1)
    return (1 - t)[:, None, None] * gamma * JS + t[:, None, None] * Jval(z)


def Ht(z, t):
    return Tval(z) - gamma * (z**deg - 1)


def dzdt(z, t):
    return -ssolve(Hz(z, t), Ht(z, t)[..., None])[..., 0]


def track(gam):
    global gamma
    gamma = gam
    roots = np.exp(2j * np.pi * np.arange(deg) / deg)
    import itertools
    start = np.array(list(itertools.product(roots, repeat=5)))
    N = len(start)
    z = start.copy()
    t = np.zeros(N)
    h = np.full(N, 0.01)
    active = np.ones(N, bool)
    diverged = np.zeros(N, bool)
    t0 = time.time()
    it = 0
    while active.any() and it < 20000:
        it += 1
        idx = np.where(active)[0]
        zi, ti, hi = z[idx], t[idx], np.minimum(h[idx], 1 - t[idx])
        # RK4 predictor
        k1 = dzdt(zi, ti)
        k2 = dzdt(zi + hi[:, None] / 2 * k1, ti + hi / 2)
        k3 = dzdt(zi + hi[:, None] / 2 * k2, ti + hi / 2)
        k4 = dzdt(zi + hi[:, None] * k3, ti + hi)
        zp = zi + hi[:, None] / 6 * (k1 + 2 * k2 + 2 * k3 + k4)
        tn = ti + hi
        ok = np.ones(len(idx), bool)
        for _ in range(3):
            dz = ssolve(Hz(zp, tn), H(zp, tn)[..., None])[..., 0]
            zp = zp - dz
        nd = np.linalg.norm(dz, axis=1)
        ok = nd < 1e-8 * (1 + np.linalg.norm(zp, axis=1))
        # accept
        acc = idx[ok]
        z[acc] = zp[ok]
        t[acc] = tn[ok]
        h[acc] = np.minimum(h[acc] * 1.6, 0.05)
        rej = idx[~ok]
        h[rej] *= 0.5
        big = np.linalg.norm(z, axis=1) > 1e6
        diverged |= big
        active = (t < 1 - 1e-14) & ~diverged & (h > 1e-12)
        if it % 200 == 0:
            print(it, 'active', active.sum(), 'done', (t >= 1 - 1e-14).sum(), 'div', diverged.sum(), f'{time.time()-t0:.0f}s', flush=True)
    
    
    return z[(t >= 1 - 1e-14) & ~diverged]

allfin = []
sols = []
for run in range(6):
    fin = track(np.exp(2j*np.pi*rng.random()))
    for _ in range(5):
        fin = fin - ssolve(Jval(fin), Tval(fin)[..., None])[..., 0]
    res = np.linalg.norm(Tval(fin), axis=1)
    fin = fin[(res < 1e-8) & (np.linalg.norm(fin, axis=1) < 1e5)]
    for p in fin:
        if all(np.linalg.norm(p - q) > 1e-6 * (1 + np.linalg.norm(p)) for q in sols):
            sols.append(p)
    print('run', run, 'total distinct', len(sols), flush=True)
    if len(sols) >= 320: break
sols = np.array(sols)
real = sols[np.abs(sols.imag).max(axis=1) < 1e-8].real
print('distinct finite solutions:', len(sols), ' real:', len(real))
pickle.dump({'G': G, 'sols': sols, 'real': real}, open(BUILD / 'lines.pkl', 'wb'))
