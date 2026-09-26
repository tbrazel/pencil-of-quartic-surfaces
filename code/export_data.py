"""Step 4: export everything the site needs to code/data/pencil.json
(affine coefficients of f and Q in the chart w = 1, lambda range, per-line data with resultants)."""
import json, pickle
from pathlib import Path
import numpy as np, sympy as sp

HERE = Path(__file__).resolve().parent
BUILD = HERE / "build"
D = pickle.load(open(BUILD / "pencil.pkl", "rb"))
lines = json.load(open(BUILD / "lines.json"))
X = sp.symbols("x0:4")
fN = sp.lambdify(X, D["f"], "numpy"); QN = sp.lambdify(X, D["Q"], "numpy")

# lambda range: real members are level sets of f/Q (Q > 0 on RP^3)
rng = np.random.default_rng(0)
P = rng.standard_normal((400000, 4)); P /= np.linalg.norm(P, axis=1)[:, None]
R = fN(*P.T) / QN(*P.T)

mons = [(i, j, k) for i in range(5) for j in range(5) for k in range(5) if i + j + k <= 4]
def affine(expr):
    Pl = sp.Poly(expr.subs(X[3], 1), X[0], X[1], X[2])
    return [float(f"{float(Pl.coeff_monomial(X[0]**i * X[1]**j * X[2]**k)):.10g}") for (i, j, k) in mons]

t = sp.symbols("t")
for i, L in enumerate(lines):
    L["res"] = float(sp.resultant(sum(sp.Float(v, 30) * t**k for k, v in enumerate(L["F1"])),
                                  sum(sp.Float(v, 30) * t**k for k, v in enumerate(L["F2"])), t))
    L["name"] = f"L{i+1}"
    for key in ["p", "dir", "n1", "n2", "F1", "F2", "q", "c"]:
        L[key] = [float(f"{v:.10g}") for v in L[key]]
    for key in ["sign_chart", "rebase_signs", "dist"]:
        L.pop(key, None)

out = {"mons": mons, "f": affine(D["f"]), "Q": affine(D["Q"]),
       "lamMin": float(R.min()), "lamMax": float(R.max()), "ballR": 3.3, "lines": lines}
(HERE / "data").mkdir(exist_ok=True)
json.dump(out, open(HERE / "data" / "pencil.json", "w"), separators=(",", ":"))
print(f"wrote data/pencil.json: {len(lines)} real lines, lambda in [{R.min():.4f}, {R.max():.4f}]")
