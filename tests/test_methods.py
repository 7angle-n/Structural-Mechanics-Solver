import os, sys
import pytest
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from solver import methods

SS = [{"type": "pin", "x": 0}, {"type": "roller", "x": 10}]


def vals(r):
    return {m["name"]: m for m in r["results"]}


def test_ss_udl_all_methods_agree():
    r = methods(dict(length=10, supports=SS, loads=[{"type": "udl", "x_start": 0, "x_end": 10, "w": 4}]), 2, 5)
    for m in r["results"]:
        assert m["y"] == pytest.approx(-5 * 4 * 10 ** 4 / (384 * 2)), m["name"]
    assert vals(r)["Conjugate beam"]["theta"] == pytest.approx(0, abs=1e-9)


def test_cantilever_left_and_right():
    for sx, px in ((0, 4), (4, 0)):
        r = methods(dict(length=4, supports=[{"type": "fixed", "x": sx}], loads=[{"type": "point", "x": px, "magnitude": 10}]), 1, px)
        for m in r["results"]:
            assert m["y"] == pytest.approx(-10 * 64 / 3), (sx, m["name"])


def test_overhang_available_methods_agree():
    r = methods(dict(length=8, supports=[{"type": "pin", "x": 0}, {"type": "roller", "x": 6}],
                     loads=[{"type": "point", "x": 8, "magnitude": 10}]), 3, 8)
    ok = [m for m in r["results"] if "y" in m]
    assert {m["name"] for m in ok} == {"Double integration", "Superposition", "Unit load (Castigliano)"}
    assert all(m["y"] == pytest.approx(ok[0]["y"]) for m in ok)
    assert all(m["theta"] == pytest.approx(ok[0]["theta"]) for m in ok)
