import os, sys
import pytest
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from solver import solve

SS = [{"type": "pin", "x": 0}, {"type": "roller", "x": 10}]


def test_simply_supported_point_load():
    r = solve(dict(length=10, supports=SS, loads=[{"type": "point", "x": 5, "magnitude": 10}]))
    assert r["reactions"]["R_A"] == pytest.approx(5)
    assert r["extremes"]["M_max"]["value"] == pytest.approx(25)


def test_simply_supported_udl():
    r = solve(dict(length=10, supports=SS, loads=[{"type": "udl", "x_start": 0, "x_end": 10, "w": 4}]))
    assert r["reactions"]["R_B"] == pytest.approx(20)
    assert r["extremes"]["M_max"]["value"] == pytest.approx(4 * 100 / 8)


def test_cantilever_udl():
    r = solve(dict(length=4, supports=[{"type": "fixed", "x": 0}], loads=[{"type": "udl", "x_start": 0, "x_end": 4, "w": 5}]))
    assert r["reactions"]["R_A"] == pytest.approx(20)
    assert r["reactions"]["M_A"] == pytest.approx(40)
    assert r["extremes"]["M_min"]["value"] == pytest.approx(-40)


def test_overhang():
    r = solve(dict(length=8, supports=[{"type": "pin", "x": 0}, {"type": "roller", "x": 6}],
                   loads=[{"type": "point", "x": 8, "magnitude": 10}]))
    assert r["reactions"]["R_B"] == pytest.approx(80 / 6)
    assert r["extremes"]["M_min"]["value"] == pytest.approx(-20)


def test_triangular_load():
    r = solve(dict(length=6, supports=[{"type": "pin", "x": 0}, {"type": "roller", "x": 6}],
                   loads=[{"type": "uvl", "x_start": 0, "x_end": 6, "w_start": 0, "w_end": 6}]))
    assert r["reactions"]["R_A"] == pytest.approx(6)
    assert r["extremes"]["M_max"]["value"] == pytest.approx(6 * 36 / (9 * 3 ** 0.5))


def test_applied_moment():
    r = solve(dict(length=10, supports=SS, loads=[{"type": "moment", "x": 5, "magnitude": 20}]))
    assert r["reactions"]["R_A"] == pytest.approx(-2)
    assert r["reactions"]["R_B"] == pytest.approx(2)


def test_rejects_indeterminate():
    with pytest.raises(ValueError):
        solve(dict(length=10, supports=SS + [{"type": "roller", "x": 5}], loads=[]))
