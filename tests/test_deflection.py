import os, sys
import pytest
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from solver import deflection


def test_cantilever_tip_load():
    d = deflection(dict(length=4, supports=[{"type": "fixed", "x": 0}], loads=[{"type": "point", "x": 4, "magnitude": 10}]), 1, 4)
    assert d["at"]["y"] == pytest.approx(-10 * 64 / 3)
    assert d["at"]["theta"] == pytest.approx(-10 * 16 / 2)


def test_simply_supported_udl_midspan():
    d = deflection(dict(length=10, supports=[{"type": "pin", "x": 0}, {"type": "roller", "x": 10}],
                        loads=[{"type": "udl", "x_start": 0, "x_end": 10, "w": 4}]), 2)
    assert d["extremes"]["y_min"]["value"] == pytest.approx(-5 * 4 * 10 ** 4 / (384 * 2))
    assert d["extremes"]["y_min"]["x"] == pytest.approx(5)


def test_simply_supported_point_load_off_centre():
    d = deflection(dict(length=10, supports=[{"type": "pin", "x": 0}, {"type": "roller", "x": 10}],
                        loads=[{"type": "point", "x": 3, "magnitude": 10}]), 1, 3)
    a, b, P, L = 3, 7, 10, 10
    assert d["at"]["y"] == pytest.approx(-P * a * a * b * b / (3 * L))
