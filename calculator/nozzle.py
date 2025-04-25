import math
from typing import Optional
from scipy.optimize import fsolve

def area_mach_relation(
    *,
    gamma: Optional[float] = None,     # – (ratio of specific heats)
    M: Optional[float] = None,         # – (Mach number)
    A: Optional[float] = None,         # m² (area)
    A_star: Optional[float] = None,    # m² (throat area)
    area_ratio: Optional[float] = None # – (A/A*)
) -> float:
    """
    Solve the isentropic area–Mach relation:

        (A/A*)² = 1/M² * [2/(γ+1) * (1 + (γ–1)/2 * M²)]^((γ+1)/(γ–1))

    Exactly one of gamma, M, area_ratio, A, or A_star must be None.
    Returns whichever one is missing.
    """
    # determine A/A* if given
    if area_ratio is not None:
        R = area_ratio
    elif A is not None and A_star is not None:
        if A_star == 0.0:
            raise ValueError("A_star cannot be zero")
        R = A / A_star
    else:
        R = None

    # compute ratio directly if needed
    if R is None:
        if gamma is None or M is None:
            raise ValueError("gamma and M required to compute ratio")
        if M <= 0.0:
            raise ValueError("Mach number must be positive")
        t = (2.0 / (gamma + 1.0)) * (1.0 + (gamma - 1.0) / 2.0 * M**2)
        e = (gamma + 1.0) / (gamma - 1.0)
        ratio = math.sqrt((1.0 / M**2) * t**e)
        if A is not None:
            return A / ratio      # solve for A_star
        if A_star is not None:
            return ratio * A_star # solve for A
        return ratio             # return A/A*

    # solve for gamma if missing
    if gamma is None:
        if M is None:
            raise ValueError("M required to solve for gamma")
        def f_gamma(g: float) -> float:
            t = (2.0 / (g + 1.0)) * (1.0 + (g - 1.0) / 2.0 * M**2)
            e = (g + 1.0) / (g - 1.0)
            return math.sqrt((1.0 / M**2) * t**e) - R
        gamma_guess = 1.4
        sol = fsolve(f_gamma, gamma_guess)
        return float(sol[0])

    # solve for M if missing
    if M is None:
        def f_M(m: float) -> float:
            t = (2.0 / (gamma + 1.0)) * (1.0 + (gamma - 1.0) / 2.0 * m**2)
            e = (gamma + 1.0) / (gamma - 1.0)
            return math.sqrt((1.0 / m**2) * t**e) - R
        M_guess = 2.0
        sol = fsolve(f_M, M_guess)
        return float(sol[0])

    # compute A or A_star if one is missing
    if A is None:
        if A_star is None:
            raise ValueError("A_star required to compute A")
        return R * A_star
    if A_star is None:
        return A / R

    # all provided, return ratio
    return R


def main():
    print('nozzle.py was run but there is no code in main()')



if __name__ == '__main__':
    main()
