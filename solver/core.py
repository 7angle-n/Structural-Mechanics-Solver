"""Statically determinate beam solver (SFD/BMD) with symbolic steps.
Units are free (examples: kN, m). Conventions: applied loads and distributed
loads positive DOWN, reactions positive UP, applied couples positive CLOCKWISE,
fixed-support moment M_A positive COUNTER-CLOCKWISE, sagging M positive."""
import numpy as np
import sympy as sp

x, s = sp.symbols("x s")
RA, RB, MA = sp.symbols("R_A R_B M_A")
Q = lambda v: sp.nsimplify(v, rational=True)
tex = sp.latex


def _w(l):
    return l["w1"] + (l["w2"] - l["w1"]) * (s - l["a"]) / (l["b"] - l["a"])


def _parse(beam, determinate=True):
    L = Q(beam["length"])
    if L <= 0:
        raise ValueError("Beam length must be positive")

    def inside(v):
        if not 0 <= v <= L:
            raise ValueError(f"Position {v} is outside the beam")

    sup = sorted(({**t, "x": Q(t["x"])} for t in beam["supports"]), key=lambda t: t["x"])
    for t in sup:
        inside(t["x"])
    kinds = [t["type"] for t in sup]
    ok_fixed = kinds == ["fixed"] and sup[0]["x"] in (0, L)
    ok_two = len(sup) == 2 and "fixed" not in kinds and sup[0]["x"] != sup[1]["x"]
    if determinate and not (ok_fixed or ok_two):
        raise ValueError("Use one fixed support at an end, or two pin/roller supports (determinate beams only)")
    loads = []
    for l in beam["loads"]:
        k = l["type"]
        if k in ("point", "moment"):
            a, v = Q(l["x"]), Q(l["magnitude"])
            inside(a)
            loads.append(dict(k="p", a=a, F=v) if k == "point" else dict(k="m", a=a, C=v))
        elif k in ("udl", "uvl"):
            a, b = Q(l["x_start"]), Q(l["x_end"])
            inside(a), inside(b)
            if b <= a:
                raise ValueError("Distributed load needs x_end > x_start")
            w1, w2 = (Q(l["w"]),) * 2 if k == "udl" else (Q(l["w_start"]), Q(l["w_end"]))
            loads.append(dict(k="d", a=a, b=b, w1=w1, w2=w2))
        else:
            raise ValueError(f"Unknown load type {k}")
    return L, sup, loads


def _ev(e, xs):
    f = sp.lambdify(x, e, "numpy")
    return np.broadcast_to(np.asarray(f(xs), dtype=float), xs.shape).tolist()


def _analyze(beam):
    L, sup, loads = _parse(beam)
    steps = [dict(title="Sign convention", eqs=[
        r"\text{Applied loads: down}+;\ \text{reactions: up}+;\ \text{couples: clockwise}+;\ \text{sagging } M+"])]

    # 1. equivalent point loads
    eq = []
    for i, l in enumerate(loads):
        if l["k"] == "p":
            l["Fr"], l["M0"] = l["F"], l["F"] * l["a"]
        elif l["k"] == "d":
            w = _w(l)
            l["Fr"] = sp.integrate(w, (s, l["a"], l["b"]))
            l["M0"] = sp.integrate(s * w, (s, l["a"], l["b"]))
            if l["Fr"] != 0:
                eq.append(rf"W_{{{i+1}}}={tex(l['Fr'])},\ \text{{acting at }}\bar x={tex(sp.nsimplify(l['M0']/l['Fr']))}")
        else:
            l["Fr"] = l["M0"] = 0
    if eq:
        steps.append(dict(title="Equivalent point loads of distributed loads", eqs=eq))
    SF = sum(l["Fr"] for l in loads)
    SM0 = sum(l["M0"] for l in loads)
    SC = sum(l["C"] for l in loads if l["k"] == "m")

    # 2. reactions
    fixed = sup[0]["type"] == "fixed"
    if fixed:
        p = sup[0]["x"]
        Ry, Ma = SF, SM0 - p * SF + SC
        rx, vals = [(p, Ry)], {"R_A": Ry, "M_A": Ma}
        steps.append(dict(title="Equilibrium: find the support reactions", eqs=[
            r"\sum F_y=0:\ " + tex(sp.Eq(RA, SF)),
            r"\sum M_A=0:\ " + tex(sp.Eq(MA, Ma))]))
    else:
        a, b = sup[0]["x"], sup[1]["x"]
        Rb = (SM0 - a * SF + SC) / (b - a)
        Ra = SF - Rb
        rx, vals = [(a, Ra), (b, Rb)], {"R_A": Ra, "R_B": Rb}
        steps.append(dict(title="Equilibrium: find the support reactions", eqs=[
            rf"\sum M_A=0:\ R_B\,({tex(b-a)})={tex(SM0 - a*SF + SC)}\ \Rightarrow\ " + tex(sp.Eq(RB, Rb)),
            r"\sum F_y=0:\ R_A+R_B=" + tex(SF) + r"\ \Rightarrow\ " + tex(sp.Eq(RA, Ra))]))

    # 3. section equations (cut from the left)
    pts = {sp.Integer(0), L, *[t["x"] for t in sup], *[l["a"] for l in loads],
           *[l["b"] for l in loads if l["k"] == "d"]}
    bps = sorted(pts)
    segs = []
    for x0, x1 in zip(bps, bps[1:]):
        V = sum(r for p_, r in rx if p_ <= x0)
        M = sum(r * (x - p_) for p_, r in rx if p_ <= x0)
        if fixed and rx[0][0] <= x0:
            M -= Ma
        for l in loads:
            if l["a"] > x0:
                continue
            if l["k"] == "p":
                V -= l["F"]; M -= l["F"] * (x - l["a"])
            elif l["k"] == "m":
                M += l["C"]
            else:
                u = l["b"] if l["b"] <= x0 else x
                V -= sp.integrate(_w(l), (s, l["a"], u))
                M -= sp.integrate(_w(l) * (x - s), (s, l["a"], u))
        segs.append((x0, x1, sp.expand(V), sp.expand(M)))
    sections = []
    for k, (x0, x1, V, M) in enumerate(segs, 1):
        sections.append(dict(range=rf"{tex(x0)}\le x\le {tex(x1)}", V=tex(V), M=tex(M)))
        steps.append(dict(title=f"Section {k}: cut at x, {x0} ≤ x ≤ {x1}", eqs=[
            rf"V(x)={tex(V)}", rf"M(x)={tex(M)}",
            rf"\text{{check: }}\frac{{dM}}{{dx}}={tex(sp.expand(sp.diff(M, x)))}=V(x)"]))

    # 4. numeric results
    X, Vp, Mp = [0.0], [0.0], [0.0]
    for x0, x1, V, M in segs:
        xs = np.linspace(float(x0), float(x1), 60)
        X += xs.tolist(); Vp += _ev(V, xs); Mp += _ev(M, xs)
    X.append(float(L)); Vp.append(0.0); Mp.append(0.0)

    cands, table = [], []
    for x0, x1, V, M in segs:
        for r in {x0, x1} | {r for r in sp.solve(V, x) if r.is_real and x0 < r < x1}:
            cands.append((float(r), float(M.subs(x, r))))
        table.append(dict(x=float(x0), V=float(V.subs(x, x0)), M=float(M.subs(x, x0))))
    lx, lV, lM = segs[-1][1], segs[-1][2], segs[-1][3]
    table.append(dict(x=float(lx), V=float(lV.subs(x, lx)), M=float(lM.subs(x, lx))))
    mx, mn = max(cands, key=lambda c: c[1]), min(cands, key=lambda c: c[1])
    res = dict(
        reactions={k: float(v) for k, v in vals.items()},
        extremes=dict(M_max=dict(x=mx[0], value=mx[1]), M_min=dict(x=mn[0], value=mn[1]),
                      V_max_abs=max(abs(v) for v in Vp)),
        steps=steps, sections=sections, table=table, plot=dict(x=X, V=Vp, M=Mp))

    return res, dict(L=L, sup=sup, segs=segs, fixed=fixed)


def solve(beam):
    return _analyze(beam)[0]
