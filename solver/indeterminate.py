"""Statically indeterminate beams by the compatibility (double-integration) method.
Unknowns: every support reaction (vertical force, plus a moment at fixed supports) and the 2n integration
constants. Equations: 2 equilibrium + continuity of slope/deflection between sections + y=0 at every support
(+ slope=0 at fixed supports). Conventions as in core.py; y positive up."""
import numpy as np
import sympy as sp
from .core import _parse, _ev, _w, x, s, tex


def indeterminate(beam, EIv=1):
    L, sup, loads = _parse(beam, determinate=False)
    EI = sp.Symbol("EI", positive=True)
    nm = [chr(65 + k) for k in range(len(sup))]
    Rs = [sp.Symbol(f"R_{n}") for n in nm]
    Ms = {k: sp.Symbol(f"M_{nm[k]}") for k, t in enumerate(sup) if t["type"] == "fixed"}
    SF = SM0 = SC = 0
    for l in loads:
        if l["k"] == "p":
            SF += l["F"]; SM0 += l["F"] * l["a"]
        elif l["k"] == "m":
            SC += l["C"]
        else:
            w = _w(l)
            SF += sp.integrate(w, (s, l["a"], l["b"])); SM0 += sp.integrate(s * w, (s, l["a"], l["b"]))
    eq0 = [sp.Eq(sum(Rs), SF), sp.Eq(sum(R * t["x"] for R, t in zip(Rs, sup)) + sum(Ms.values()), SM0 + SC)]
    bps = sorted({sp.Integer(0), L, *[t["x"] for t in sup], *[l["a"] for l in loads], *[l["b"] for l in loads if l["k"] == "d"]})
    segs = []
    for x0, x1 in zip(bps, bps[1:]):
        V = sum(R for R, t in zip(Rs, sup) if t["x"] <= x0)
        M = sum(R * (x - t["x"]) for R, t in zip(Rs, sup) if t["x"] <= x0) - sum(Mk for k, Mk in Ms.items() if sup[k]["x"] <= x0)
        for l in loads:
            if l["a"] > x0:
                continue
            if l["k"] == "p":
                V -= l["F"]; M -= l["F"] * (x - l["a"])
            elif l["k"] == "m":
                M += l["C"]
            else:
                u = l["b"] if l["b"] <= x0 else x
                V -= sp.integrate(_w(l), (s, l["a"], u)); M -= sp.integrate(_w(l) * (x - s), (s, l["a"], u))
        segs.append((x0, x1, V, M))
    n = len(segs)
    C = sp.symbols(f"C_1:{2*n+1}")
    thC = [sp.integrate(sg[3], x) / EI + C[2 * k] for k, sg in enumerate(segs)]
    yC = [sp.integrate(thC[k], x) + C[2 * k + 1] for k in range(n)]
    eqs, conds = [], []
    for k in range(n - 1):
        b = segs[k][1]
        eqs += [(thC[k] - thC[k + 1]).subs(x, b), (yC[k] - yC[k + 1]).subs(x, b)]
        conds.append(rf"x={tex(b)}:\ \theta_{{{k+1}}}=\theta_{{{k+2}}},\ y_{{{k+1}}}=y_{{{k+2}}}")
    for t, name in zip(sup, nm):
        p = t["x"]
        k = next(i for i, sg in enumerate(segs) if sg[0] <= p <= sg[1])
        eqs.append(yC[k].subs(x, p)); conds.append(rf"\text{{support }}{name}:\ y({tex(p)})=0")
        if t["type"] == "fixed":
            eqs.append(thC[k].subs(x, p)); conds.append(rf"\text{{support }}{name}:\ \theta({tex(p)})=0")
    unk = Rs + list(Ms.values()) + list(C)
    sol = sp.solve(eq0 + eqs, unk, dict=True)
    if not sol or any(u not in sol[0] for u in unk):
        raise ValueError("The beam is unstable or the supports are inconsistent")
    sol = sol[0]
    Vf = [sp.expand(sg[2].subs(sol)) for sg in segs]
    Mf = [sp.expand(sg[3].subs(sol)) for sg in segs]
    thf = [sp.expand(e.subs(sol)) for e in thC]
    yf = [sp.expand(e.subs(sol)) for e in yC]
    r = len(Rs) + len(Ms)
    reds = list(Rs) + list(Ms.values())
    steps = [
        dict(title="Degree of indeterminacy", eqs=[rf"r={r}\ \text{{unknown reactions}},\ 2\ \text{{equilibrium equations}}\ \Rightarrow\ \text{{degree}}={r-2}"]),
        dict(title="Equilibrium (all reactions unknown)", eqs=[tex(e) for e in eq0]),
        dict(title="Bending moment in each section, in terms of the reactions", eqs=[
            rf"{tex(sg[0])}\le x\le {tex(sg[1])}:\ M={tex(sp.expand(sg[3]))}" for sg in segs]),
        dict(title="Integrate EI y'' = M(x) twice in each section", eqs=[
            r"EI\,y''=M(x)\ \Rightarrow\ \theta_k=\int \frac{M}{EI}dx+C_{2k-1},\quad y_k=\int\theta_k\,dx+C_{2k}"]),
        dict(title="Compatibility conditions", eqs=conds),
        dict(title="Solve all unknowns together", eqs=[tex(sp.Eq(u, sp.simplify(sol[u]))) for u in reds]),
        dict(title="Final shear, moment, slope and deflection", eqs=[
            rf"{tex(segs[k][0])}\le x\le {tex(segs[k][1])}:\ V={tex(Vf[k])},\ M={tex(Mf[k])}" for k in range(n)]
            + [rf"\ \ \ \ \ \theta={tex(thf[k])},\ y={tex(yf[k])}" for k in range(n)])]
    e = sp.nsimplify(EIv)
    X, Vp, Mp, Yp, cm, cy = [], [], [], [], [], []
    for k, (x0, x1, _, _) in enumerate(segs):
        xs = np.linspace(float(x0), float(x1), 60)
        yn, thn = yf[k].subs(EI, e), thf[k].subs(EI, e)
        X += xs.tolist(); Vp += _ev(Vf[k], xs); Mp += _ev(Mf[k], xs); Yp += _ev(yn, xs)
        for q in {x0, x1} | {q for q in sp.solve(Vf[k], x) if q.is_real and x0 < q < x1}:
            cm.append((float(q), float(Mf[k].subs(x, q))))
        for q in {x0, x1} | {q for q in sp.solve(thn, x) if q.is_real and x0 < q < x1}:
            cy.append((float(q), float(yn.subs(x, q))))
    mx, mn = max(cm, key=lambda c: c[1]), min(cm, key=lambda c: c[1])
    ylo, yhi = min(cy, key=lambda c: c[1]), max(cy, key=lambda c: c[1])
    return dict(reactions={str(u): float(sol[u]) for u in reds}, steps=steps, plot=dict(x=X, V=Vp, M=Mp, y=Yp),
                extremes=dict(M_max=dict(x=mx[0], value=mx[1]), M_min=dict(x=mn[0], value=mn[1]),
                              y_min=dict(x=ylo[0], value=ylo[1]), y_max=dict(x=yhi[0], value=yhi[1])))
