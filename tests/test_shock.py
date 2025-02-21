import math
import pytest
from shock import (
    mach_angle,
    normal_shock_rho,
    normal_shock_P,
    normal_shock_T,
    normal_shock_P0,
    pitot_tube,
)


class TestMachAngle:
    def test_sample(self):
        ...

class TestNormalShockRho:
    def test_mach_from_ratio(self):
        assert math.isclose(normal_shock_rho(rho1=1, rho2= 2.66666666), 2.0, rel_tol=1e-5)
        

class TestNormalShockP:
    def test_mach_from_ratio(self):
        assert math.isclose(normal_shock_P(P1=1, P2=4.5), 2.0, rel_tol=1e-5)

class TestNormalShockT:
    def test_mach_from_ratio(self):
        assert math.isclose(normal_shock_T(T1=1, T2=1.6875), 2.0, rel_tol=1e-5)

class TestNormalShockP0:
    def sample_test(self):
        ...

class TestPitotTube:
    def sample_test(self):
        ...

