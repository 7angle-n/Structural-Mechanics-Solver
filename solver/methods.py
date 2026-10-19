"""Deflection at one point by several classical methods (determinate beams). y positive UP, slope CCW +.
Conjugate-beam and moment-area are offered for simply supported beams and cantilevers; the other methods work for any determinate beam."""
import sympy as sp
from .core import _analyze, x, s, tex, Q
from .deflection import deflection

g6 = lambda v: f"{v:.6g}"


def _cases(sup, L):
    k = [t["type"] for t in sup]
    if len(sup) == 2 and "fixed" not in k and sup[0]["x"] == 0 and sup[1]["x"] == L:
        return "ss"
    if k == ["fixed"] and sup[0]["x"] == 0:
        return "cl"
    if k == ["fixed"] and sup[0]["x"] == L:
        return "cr"
    return None


def _pw(segs, f, lo, hi):
    tot = 0
    for x0, x1, V, M in segs:
        a, b = max(x0, lo), min(x1, hi)
        if a < b:
            tot += sp.integrate(sp.expand(f(sp.sympify(M).subs(x, s), s)), (s, a, b))
    return sp.simplify(tot)


def _pair(A, B, EI):
    pts = sorted({*(q[0] for q in A), *(q[1] for q in A), *(q[0] for q in B), *(q[1] for q in B)})
    tot = 0
    for u0, u1 in zip(pts, pts[1:]):
        m = (u0 + u1) / 2
        fa = next(q[3] for q in A if q[0] <= m <= q[1])
        fb = next(q[3] for q in B if q[0] <= m <= q[1])
        tot += sp.integrate(sp.expand(fa * fb), (x, u0, u1))
    return sp.simplify(tot / EI)


def methods(beam, EIv=1, at=None):
    if at is None:
        raise ValueError("Enter the position x where the deflection is wanted")
    _, c = _analyze(beam)
    segs, sup, L = c["segs"], c["sup"], c["L"]
    a = Q(at)
    if not 0 <= a <= L:
        raise ValueError("x is outside the beam")
    EI, e = sp.Symbol("EI", positive=True), sp.nsimplify(EIv)
    num = lambda v: float(sp.sympify(v).subs(EI, e))
    case, out = _cases(sup, L), []
    d = deflection(beam, EIv, at)
    out.append(dict(name="Double integration", steps=d["steps"], y=d["at"]["y"], theta=d["at"]["theta"]))
    f1, f0 = (lambda M, t: M / EI), None
    kw = lambda lo, hi, fn: _pw(segs, fn, lo, hi)
    if case:
        if case == "ss":
            Ra = kw(0, L, lambda M, t: M / EI * (L - t)) / L
            Rb = kw(0, L, f1) - Ra
            Ms = Ra * a - kw(0, a, lambda M, t: M / EI * (a - t))
            Vs = Ra - kw(0, a, f1)
            y, th = -Ms, -Vs
            cj = [r"\text{Conjugate beam: same span, same (pin/roller) supports, loaded with } M(x)/EI",
                  rf"R^*_A=\frac1L\int_0^L\frac{{M}}{{EI}}(L-s)\,ds={tex(Ra)},\quad R^*_B={tex(Rb)}",
                  rf"M^*(x)=R^*_Ax-\int_0^x\frac{{M}}{{EI}}(x-s)\,ds={tex(Ms)},\quad V^*(x)={tex(Vs)}"]
        elif case == "cl":
            y, th = kw(0, a, lambda M, t: M / EI * (a - t)), kw(0, a, f1)
            cj = [r"\text{Fixed end} \to \text{free end, free end} \to \text{fixed end; load the conjugate beam with } M/EI",
                  rf"M^*={tex(-y)},\quad V^*={tex(-th)}"]
        else:
            y, th = kw(a, L, lambda M, t: M / EI * (t - a)), -kw(a, L, f1)
            cj = [r"\text{Fixed end} \to \text{free end, free end} \to \text{fixed end; load the conjugate beam with } M/EI",
                  rf"M^*={tex(-y)},\quad V^*={tex(-th)}"]
        cj.append(rf"y=-M^*={tex(sp.simplify(y))},\quad \theta=-V^*={tex(sp.simplify(th))}")
        out.append(dict(name="Conjugate beam", steps=[dict(title="Conjugate beam method", eqs=cj)], y=num(y), theta=num(th)))
        if case == "ss":
            tBA = kw(0, L, lambda M, t: M / EI * (L - t))
            thA = -tBA / L
            y2 = thA * a + kw(0, a, lambda M, t: M / EI * (a - t))
            th2 = thA + kw(0, a, f1)
            ma = [rf"t_{{B/A}}=\int_0^L\frac{{M}}{{EI}}(L-s)\,ds={tex(tBA)}\quad(\text{{moment of the }}M/EI\text{{ area about B}})",
                  rf"\theta_A=-\frac{{t_{{B/A}}}}{{L}}={tex(thA)}\quad(y_B=0)",
                  rf"t_{{x/A}}=\int_0^x\frac{{M}}{{EI}}(x-s)\,ds,\quad y=\theta_Ax+t_{{x/A}}={tex(sp.simplify(y2))}",
                  rf"\theta=\theta_A+\int_0^x\frac{{M}}{{EI}}ds={tex(sp.simplify(th2))}"]
        elif case == "cl":
            y2, th2 = y, th
            ma = [r"\text{Tangent at the fixed end A is horizontal: } \theta_A=0,\ y_A=0",
                  rf"y=t_{{x/A}}=\int_0^x\frac{{M}}{{EI}}(x-s)\,ds={tex(sp.simplify(y2))}",
                  rf"\theta=\theta_A+\int_0^x\frac{{M}}{{EI}}ds={tex(sp.simplify(th2))}"]
        else:
            y2, th2 = y, th
            ma = [r"\text{Tangent at the fixed end B is horizontal: } \theta_B=0,\ y_B=0",
                  rf"y=\int_x^L\frac{{M}}{{EI}}(s-x)\,ds={tex(sp.simplify(y2))}",
                  rf"\theta=-\int_x^L\frac{{M}}{{EI}}ds={tex(sp.simplify(th2))}"]
        out.append(dict(name="Moment-area", steps=[dict(title="Moment-area theorems", eqs=ma)], y=num(y2), theta=num(th2)))
    else:
        note = "Available for simply supported beams (supports at both ends) and cantilevers. Use another method for overhangs."
        out += [dict(name="Conjugate beam", note=note), dict(name="Moment-area", note=note)]
    # superposition: one load at a time
    ys, ths, eq = 0.0, 0.0, []
    for i, l in enumerate(beam["loads"], 1):
        di = deflection({**beam, "loads": [l]}, EIv, at)["at"]
        ys += di["y"]; ths += di["theta"]
        eq.append(rf"\text{{Load {i} ({l['type']}) alone: }}y={g6(di['y'])},\ \theta={g6(di['theta'])}")
    eq.append(rf"\text{{Sum: }}y={g6(ys)},\ \theta={g6(ths)}")
    out.append(dict(name="Superposition", steps=[dict(title="Each load acting alone, then added", eqs=eq)], y=ys, theta=ths))
    # unit load / Castigliano
    b1 = {**beam, "loads": [dict(type="point", x=at, magnitude=1)]}
    b2 = {**beam, "loads": [dict(type="moment", x=at, magnitude=1)]}
    s1, s2 = _analyze(b1)[1]["segs"], _analyze(b2)[1]["segs"]
    dd, rr = _pair(segs, s1, EI), _pair(segs, s2, EI)
    out.append(dict(name="Unit load (Castigliano)", y=num(-dd), theta=num(-rr), steps=[dict(title="Virtual work with a unit load at the point", eqs=[
        r"M(x,P)=M(x)+P\,m(x)\ \Rightarrow\ \delta=\frac{\partial U}{\partial P}\Big|_{P=0}=\int\frac{M\,m}{EI}\,dx",
        *[rf"{tex(q[0])}\le x\le {tex(q[1])}:\ m_{{\text{{unit load}}}}={tex(sp.expand(q[3]))}" for q in s1],
        rf"\delta_{{down}}=\int\frac{{M\,m}}{{EI}}dx={tex(dd)}\ \Rightarrow\ y=-\delta={tex(-dd)}",
        *[rf"{tex(q[0])}\le x\le {tex(q[1])}:\ m_{{\text{{unit couple}}}}={tex(sp.expand(q[3]))}" for q in s2],
        rf"\varphi_{{cw}}=\int\frac{{M\,m_\theta}}{{EI}}dx={tex(rr)}\ \Rightarrow\ \theta=-\varphi={tex(-rr)}"])]))
    return dict(results=out, at=float(at), EI=float(EIv))
