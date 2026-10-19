"""Beam deflection by double integration: EI y'' = M(x), y positive UP.
Constants come from continuity (slope, deflection) between sections and the support conditions."""
import numpy as np
import sympy as sp
from .core import _analyze, _ev, x, tex


def deflection(beam, EIv=1, at=None):
    _, c = _analyze(beam)
    segs, sup = c["segs"], c["sup"]
    n, EI = len(segs), sp.Symbol("EI", positive=True)
    C = sp.symbols(f"C_1:{2*n+1}")
    thC = [sp.integrate(M, x) / EI + C[2 * k] for k, (_, _, _, M) in enumerate(segs)]
    yC = [sp.integrate(thC[k], x) + C[2 * k + 1] for k in range(n)]
    eqs, conds = [], []
    for k in range(n - 1):
        b = segs[k][1]
        eqs += [(thC[k] - thC[k + 1]).subs(x, b), (yC[k] - yC[k + 1]).subs(x, b)]
        conds.append(rf"x={tex(b)}:\ \theta_{{{k+1}}}=\theta_{{{k+2}}},\ y_{{{k+1}}}=y_{{{k+2}}}")
    for t in sup:
        p = t["x"]
        k = next(i for i, sg in enumerate(segs) if sg[0] <= p <= sg[1])
        eqs.append(yC[k].subs(x, p)); conds.append(rf"y({tex(p)})=0")
        if t["type"] == "fixed":
            eqs.append(thC[k].subs(x, p)); conds.append(rf"\theta({tex(p)})=0")
    sol = sp.solve(eqs, C, dict=True)[0]
    th = [sp.expand(e.subs(sol)) for e in thC]
    y = [sp.expand(e.subs(sol)) for e in yC]
    steps = [
        dict(title="Governing equation", eqs=[r"EI\,\dfrac{d^2y}{dx^2}=M(x),\quad \theta=\dfrac{dy}{dx},\ \ y\ \text{positive upward}"]),
        dict(title="Integrate M(x) in each section", eqs=[
            q for k in range(n) for q in (rf"\theta_{{{k+1}}}={tex(thC[k])}", rf"y_{{{k+1}}}={tex(yC[k])}")]),
        dict(title="Boundary and continuity conditions", eqs=conds),
        dict(title="Solve for the constants", eqs=[tex(sp.Eq(k_, sp.simplify(sol[k_]))) for k_ in C]),
        dict(title="Slope and deflection equations", eqs=[
            rf"{tex(segs[k][0])}\le x\le {tex(segs[k][1])}:\ \theta={tex(th[k])},\ \ y={tex(y[k])}" for k in range(n)])]
    thn = [e.subs(EI, sp.nsimplify(EIv)) for e in th]
    yn = [e.subs(EI, sp.nsimplify(EIv)) for e in y]
    X, Y, T = [], [], []
    cands, atr = [], None
    for k, (x0, x1, _, _) in enumerate(segs):
        xs = np.linspace(float(x0), float(x1), 60)
        X += xs.tolist(); Y += _ev(yn[k], xs); T += _ev(thn[k], xs)
        for r in {x0, x1} | {r for r in sp.solve(thn[k], x) if r.is_real and x0 < r < x1}:
            cands.append((float(r), float(yn[k].subs(x, r))))
        if at is not None and x0 <= sp.nsimplify(at) <= x1 and atr is None:
            atr = dict(x=float(at), y=float(yn[k].subs(x, sp.nsimplify(at))), theta=float(thn[k].subs(x, sp.nsimplify(at))))
    if at is not None and atr is None:
        raise ValueError("x is outside the beam")
    lo, hi = min(cands, key=lambda q: q[1]), max(cands, key=lambda q: q[1])
    return dict(steps=steps, plot=dict(x=X, y=Y, theta=T), at=atr, EI=float(EIv),
                extremes=dict(y_min=dict(x=lo[0], value=lo[1]), y_max=dict(x=hi[0], value=hi[1])))
