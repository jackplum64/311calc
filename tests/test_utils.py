import math
import pytest
from utils import deg2rad, rad2deg

def test_deg2rad():
    # Test conversion for several angles
    assert math.isclose(deg2rad(180), math.pi, rel_tol=1e-9)
    assert math.isclose(deg2rad(90), math.pi / 2, rel_tol=1e-9)
    assert math.isclose(deg2rad(0), 0.0, rel_tol=1e-9)
    # Test with a negative angle
    assert math.isclose(deg2rad(-45), -math.pi / 4, rel_tol=1e-9)

def test_rad2deg():
    # Test conversion for several angles
    assert math.isclose(rad2deg(math.pi), 180.0, rel_tol=1e-9)
    assert math.isclose(rad2deg(math.pi / 2), 90.0, rel_tol=1e-9)
    assert math.isclose(rad2deg(0), 0.0, rel_tol=1e-9)
    # Test with a negative angle
    assert math.isclose(rad2deg(-math.pi / 4), -45.0, rel_tol=1e-9)
