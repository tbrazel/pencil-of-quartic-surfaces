"""Step 1: choose the pencil f - lambda*Q (random Kostlan f; Q = (x^2+y^2+z^2+w^2)^2 + 0.12 g, positive on RP^3)."""
from pathlib import Path
HERE = Path(__file__).resolve().parent
BUILD = HERE / "build"
BUILD.mkdir(exist_ok=True)
import numpy as np, sympy as sp, pickle, itertools
from math import factorial
rng=np.random.default_rng(12)
X=sp.symbols('x0:4')   # homogeneous coords; affine chart x3=1 (x3 = w)
mons=[m for m in itertools.product(range(5),repeat=4) if sum(m)==4]
def kostlan():
    return {m: rng.standard_normal()*np.sqrt(factorial(4)/np.prod([factorial(e) for e in m])) for m in mons}
def topoly(C): return sum(sp.Float(c)*sp.prod([X[i]**m[i] for i in range(4)]) for m,c in C.items())
g=kostlan(); f=kostlan()
s2=sum(x**2 for x in X)
Qp=sp.expand(s2**2+0.12*topoly(g))
fp=sp.expand(topoly(f))
# positivity of Q on S^3
Qn=sp.lambdify(X,Qp,'numpy')
P=rng.standard_normal((200000,4)); P/=np.linalg.norm(P,axis=1)[:,None]
print('min Q on S^3:',Qn(*P.T).min())
pickle.dump({'f':fp,'Q':Qp},open(BUILD / 'pencil.pkl','wb'))
