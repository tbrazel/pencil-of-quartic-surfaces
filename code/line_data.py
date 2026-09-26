"""Step 3: per-line data (normal derivatives F1, F2, Q|_L, determinant functional, index sign) with sign-consistency checks."""
import numpy as np, sympy as sp, pickle, json, itertools
from pathlib import Path
HERE = Path(__file__).resolve().parent
BUILD = HERE / "build"
BUILD.mkdir(exist_ok=True)


D = pickle.load(open(BUILD / 'pencil.pkl', 'rb'))
S = pickle.load(open(BUILD / 'lines.pkl', 'rb'))
X = sp.symbols('x0:4')
fN = sp.lambdify(X, D['f'], 'numpy')
QN = sp.lambdify(X, D['Q'], 'numpy')
G = S['G']
rng = np.random.default_rng(1)


def binary_coeffs(func, A, B, deg=4):
    """coeffs (k=0..deg) of func(u*A + r*B) in basis u^{deg-k} r^k, via interpolation."""
    ts = np.linspace(-1.3, 1.7, deg + 1)
    vals = np.array([func(*(A + t * B)) for t in ts])  # u=1, r=t
    V = np.vander(ts, deg + 1, increasing=True)
    return np.linalg.solve(V, vals)


def jet(F, P):
    """For F homogeneous quartic, basis P (cols p0,p1,n1,n2): return
    (F|_L coeffs, F1, F2) where F(x0 p0 + x1 p1 + x2 n1 + x3 n2) = F|_L + x2 F1 + x3 F2 + O(2)."""
    p0, p1, n1, n2 = P.T
    base = binary_coeffs(F, p0, p1)
    eps = 1e-5
    F1 = (binary_coeffs(lambda *x: F(*x), p0 + eps * n1, p1) - binary_coeffs(F, p0 - eps * n1, p1)) / (2 * eps)
    # the above perturbs only the u-part; do it properly: derivative of F(u p0 + r p1 + s n) in s at s=0
    def dF(n):
        def g(*x):
            return 0
        cs = []
        ts = np.linspace(-1.3, 1.7, 4)
        vals = []
        for t in ts:
            pt = p0 + t * p1
            vals.append((F(*(pt + eps * n)) - F(*(pt - eps * n))) / (2 * eps))
        V = np.vander(ts, 4, increasing=True)
        return np.linalg.solve(V, np.array(vals))  # cubic coeffs, basis u^{3-k} r^k
    return base, dF(n1), dF(n2)


def detfunctional(F1, F2):
    """c with det[g, uF1, uF2, rF1, rF2] = sum c_k g_k (basis u^{4-k} r^k)."""
    uF1 = np.r_[F1, 0]; rF1 = np.r_[0, F1]
    uF2 = np.r_[F2, 0]; rF2 = np.r_[0, F2]
    c = np.zeros(5)
    for k in range(5):
        e = np.zeros(5); e[k] = 1
        c[k] = np.linalg.det(np.column_stack([e, uF1, uF2, rF1, rF2]))
    return c


out = []
for sol in S['real']:
    a, b, c_, d, lam = sol
    v0 = G[:, 0] + a * G[:, 2] + b * G[:, 3]
    v1 = G[:, 1] + c_ * G[:, 2] + d * G[:, 3]
    F = lambda *x: fN(*x) - lam * QN(*x)
    # --- sign in the (non-centered) solving chart, via the determinant lemma ---
    def sigmaF_coeffs(z):
        aa, bb, cc, dd = z
        w0 = G[:, 0] + aa * G[:, 2] + bb * G[:, 3]
        w1 = G[:, 1] + cc * G[:, 2] + dd * G[:, 3]
        return binary_coeffs(F, w0, w1)
    z0 = np.array([a, b, c_, d]); h = 1e-6
    cols = [binary_coeffs(QN, v0, v1)]
    for j in range(4):
        e = np.zeros(4); e[j] = h
        cols.append((sigmaF_coeffs(z0 + e) - sigmaF_coeffs(z0 - e)) / (2 * h))
    sign_chart = np.sign(np.linalg.det(np.column_stack(cols)))
    # --- affine data: L meets chart x3 = 1 ---
    A2 = np.array([v0, v1])
    # point with x3=1 closest to origin, direction with x3=0
    dirh = v0 * v1[3] - v1 * v0[3]
    dirv = dirh[:3] / np.linalg.norm(dirh[:3])
    # some point with x3 = 1
    pth = v0 / v0[3] if abs(v0[3]) > abs(v1[3]) else v1 / v1[3]
    p = pth[:3]
    p = p - np.dot(p, dirv) * dirv
    n1 = np.cross(dirv, [0.3, 0.5, 0.8]); n1 /= np.linalg.norm(n1)
    n2 = np.cross(dirv, n1)
    P = np.column_stack([np.r_[p, 1], np.r_[dirv, 0], np.r_[n1, 0], np.r_[n2, 0]])
    base, F1, F2 = jet(F, P)
    q = binary_coeffs(QN, P[:, 0], P[:, 1])
    cvec = detfunctional(F1, F2)
    delta = cvec @ q
    # consistency checks: F vanishes on L; sign under random re-bases
    assert np.abs(base).max() < 1e-6 * (1 + np.abs(q).max()), base
    signs = []
    for _ in range(6):
        M = rng.standard_normal((2, 2)); N_ = rng.standard_normal((2, 2))
        Lb = P[:, :2] @ M
        Nb = P[:, 2:] @ N_ + P[:, :2] @ rng.standard_normal((2, 2))
        Pb = np.column_stack([Lb, Nb])
        _, G1, G2 = jet(F, Pb)
        qb = binary_coeffs(QN, Pb[:, 0], Pb[:, 1])
        signs.append(np.sign(detfunctional(G1, G2) @ qb))
    out.append(dict(lam=float(lam), p=p.tolist(), dir=dirv.tolist(), n1=n1.tolist(), n2=n2.tolist(),
                    F1=F1.tolist(), F2=F2.tolist(), q=q.tolist(), c=cvec.tolist(), delta=float(delta),
                    sign=int(np.sign(delta)), sign_chart=int(sign_chart), rebase_signs=[int(s) for s in signs],
                    dist=float(np.linalg.norm(p))))

out.sort(key=lambda o: o['lam'])
for o in out:
    print(f"lam={o['lam']:+.4f} dist={o['dist']:.3f} sign={o['sign']:+d} chart={o['sign_chart']:+d} rebases={o['rebase_signs']}")
print('signature (sum of signs) =', sum(o['sign'] for o in out))
json.dump(out, open(BUILD / 'lines.json', 'w'))
