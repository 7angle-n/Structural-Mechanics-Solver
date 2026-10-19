import os, sys
import pytest
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from solver import indeterminate

U = lambda a, b, w=4: {"type": "udl", "x_start": a, "x_end": b, "w": w}


def test_propped_cantilever_udl():
    r = indeterminate(dict(length=10, supports=[{"type": "fixed", "x": 0}, {"type": "roller", "x": 10}], loads=[U(0, 10)]))
    assert r["reactions"]["R_A"] == pytest.approx(25) and r["reactions"]["R_B"] == pytest.approx(15)
    assert r["reactions"]["M_A"] == pytest.approx(50)


def test_fixed_fixed_udl():
    r = indeterminate(dict(length=10, supports=[{"type": "fixed", "x": 0}, {"type": "fixed", "x": 10}], loads=[U(0, 10)]))
    assert r["reactions"]["R_A"] == pytest.approx(20)
    assert r["reactions"]["M_A"] == pytest.approx(4 * 100 / 12)
    assert r["extremes"]["y_min"]["value"] == pytest.approx(-4 * 10 ** 4 / 384)


def test_two_span_continuous():
    r = indeterminate(dict(length=20, supports=[{"type": "pin", "x": 0}, {"type": "roller", "x": 10}, {"type": "roller", "x": 20}], loads=[U(0, 20)]))
    assert r["reactions"]["R_A"] == pytest.approx(15)
    assert r["reactions"]["R_B"] == pytest.approx(50)
    assert r["reactions"]["R_C"] == pytest.approx(15)
    assert r["extremes"]["M_min"]["value"] == pytest.approx(-50)
