import math
from .utils import (rad2deg, deg2rad)


def prandtl_meyer(gamma=None, mach=None):
    """
    Compute the Prandtl-Meyer function (in degrees) for a given Mach number and specific heat ratio (gamma).

    The Prandtl-Meyer function is defined as:
    
        ν(M) = sqrt((γ+1)/(γ-1)) * arctan( sqrt((γ-1)/(γ+1) * (M^2 - 1) ) ) - arctan( sqrt(M^2 - 1) )
    
    Parameters:
        gamma (float): Specific heat ratio (dimensionless). Must be greater than 1.
        mach (float): Mach number (dimensionless). Must be greater than or equal to 1.
        
    Returns:
        float: Prandtl-Meyer angle in degrees.

    Raises:
        ValueError: If gamma is not greater than 1 or if mach is less than 1.
    """
    
    # Validate inputs
    if gamma is None or mach is None:
        raise ValueError("Both 'gamma' and 'mach' must be provided.")
    if gamma <= 1:
        raise ValueError("The specific heat ratio 'gamma' must be greater than 1.")
    if mach < 1:
        raise ValueError("The Mach number 'mach' must be greater than or equal to 1.")
    
    # Compute the intermediate term sqrt(M^2 - 1) [dimensionless]
    sqrt_mach_term = math.sqrt(mach**2 - 1)
    
    # Compute the term inside the first arctan function
    sqrt_ratio = math.sqrt((gamma - 1) / (gamma + 1) * (mach**2 - 1))
    
    # Compute the first and second arctan values (in radians)
    first_term = math.sqrt((gamma + 1) / (gamma - 1)) * math.atan(sqrt_ratio)
    second_term = math.atan(sqrt_mach_term)
    
    # Return the Prandtl-Meyer angle in degrees
    return rad2deg(first_term - second_term)


def prandtl_meyer_difference(M1=None, M2=None, gamma=None):
    """
    Compute the difference in Prandtl-Meyer angles (in degrees) between two Mach numbers.

    The difference is calculated as:
    
        θ = ν(M₂) - ν(M₁)
    
    where ν(M) is the Prandtl-Meyer function evaluated at Mach number M.

    Parameters:
        M1 (float): The first Mach number (dimensionless), must be greater than or equal to 1.
        M2 (float): The second Mach number (dimensionless), must be greater than or equal to 1.
        gamma (float): Specific heat ratio (dimensionless). Default is 1.4.

    Returns:
        float: Difference in Prandtl-Meyer angle (in degrees).

    Raises:
        ValueError: If any of the Mach numbers or gamma do not meet the input requirements.
    """
    # Compute the Prandtl-Meyer angles for both Mach numbers
    angle_M1 = prandtl_meyer(gamma, M1)
    angle_M2 = prandtl_meyer(gamma, M2)
    
    # Return the difference in angles (in degrees)
    return angle_M2 - angle_M1


def expansion_T2_T1(gamma=None, M1=None, M2=None, T1=None, T2=None, T2_T1_ratio=None):
    """
    Solve the expansion temperature ratio relation

        T2/T1 = (1 + [(gamma - 1) / 2] * M1^2) / (1 + [(gamma - 1) / 2] * M2^2)

    The function may be used in two modes:
      1. Ratio mode: If T2_T1_ratio is provided (non-None), it is used and the individual temperatures T1 and T2 are ignored.
      2. Component mode: If both T1 and T2 are provided, the effective ratio is computed as T2/T1.
         In either case, if both T1 and T2 are unavailable, the function returns the computed ratio.
         
    Exactly one of the four independent quantities (T2_T1_ratio (or T2/T1), gamma, M1, or M2) must be unknown.
    The other parameters must be provided.

    Parameters:
      gamma (float): Specific heat ratio (None if unknown)
      M1 (float): Upstream Mach number (None if unknown)
      M2 (float): Downstream Mach number (None if unknown)
      T1 (float): Upstream temperature (None if unknown)
      T2 (float): Downstream temperature (None if unknown)
      T2_T1_ratio (float): Temperature ratio T2/T1 (None if unknown)

    Returns:
      float: The computed missing variable. When T2_T1_ratio is solved for, the result is the ratio;
             otherwise, the computed gamma, M1, or M2 is returned.

    Raises:
      ValueError: if not exactly one independent variable among the effective temperature ratio,
                  gamma, M1, and M2 is None or if a math error occurs.
    """
    # Determine the effective temperature ratio.
    if T2_T1_ratio is not None:
        eff_ratio = T2_T1_ratio
    elif (T1 is not None or T2 is not None):
        if T1 is None or T2 is None:
            raise ValueError("Both T1 and T2 must be provided or neither, if T2_T1_ratio is not given.")
        if T1 == 0:
            raise ValueError("T1 is zero, cannot compute temperature ratio.")
        eff_ratio = T2 / T1
    else:
        eff_ratio = None

    # Count the number of unknowns among the independent quantities.
    unknowns = [eff_ratio is None, gamma is None, M1 is None, M2 is None]
    if sum(unknowns) != 1:
        raise ValueError("Exactly one variable among the effective T2/T1 ratio, gamma, M1, and M2 must be None.")

    # Case 1: Compute T2_T1 ratio when it is unknown.
    if eff_ratio is None:
        # gamma, M1, and M2 are provided.
        a = (gamma - 1) / 2
        return (1 + a * M1**2) / (1 + a * M2**2)

    # For the following cases, gamma must be known.
    a = (gamma - 1) / 2

    # Case 2: Compute gamma when it is unknown.
    if gamma is None:
        # Solve for gamma from:
        # eff_ratio = (1 + ((gamma - 1)/2)*M1^2) / (1 + ((gamma - 1)/2)*M2^2)
        # Rearranged: ((gamma - 1)/2) = (eff_ratio - 1) / (M1^2 - eff_ratio * M2^2)
        denom = M1**2 - eff_ratio * M2**2
        if denom == 0:
            raise ValueError("Denominator is zero when computing gamma.")
        a_val = (eff_ratio - 1) / denom
        return 2 * a_val + 1

    # Case 3: Compute M1 when it is unknown.
    if M1 is None:
        # Solve: 1 + a*M1^2 = eff_ratio * (1 + a*M2^2)
        numerator = eff_ratio * (1 + a * M2**2) - 1
        if a == 0:
            raise ValueError("gamma equals 1, leading to division by zero.")
        if numerator < 0:
            raise ValueError("Computed M1^2 is negative.")
        return math.sqrt(numerator / a)

    # Case 4: Compute M2 when it is unknown.
    if M2 is None:
        # Solve: 1 + a*M2^2 = (1 + a*M1^2) / eff_ratio
        numerator = (1 + a * M1**2) / eff_ratio - 1
        if a == 0:
            raise ValueError("gamma equals 1, leading to division by zero.")
        if numerator < 0:
            raise ValueError("Computed M2^2 is negative.")
        return math.sqrt(numerator / a)


def expansion_P2_P1(gamma=None, M1=None, M2=None, P1=None, P2=None, P2_P1_ratio=None):
    """
    Solve the expansion pressure ratio relation

        P2/P1 = [ (1 + [(gamma - 1) / 2] * M1^2) / (1 + [(gamma - 1) / 2] * M2^2) ]^(gamma/(gamma - 1))

    The function may be used in two modes:
      1. Ratio mode: If P2_P1_ratio is provided (non-None), it is used and the individual pressures P1 and P2 are ignored.
      2. Component mode: If both P1 and P2 are provided, the effective ratio is computed as P2/P1.
         In either case, if both P1 and P2 are unavailable, the function returns the computed ratio.
         
    Exactly one of the four independent quantities (P2_P1_ratio (or P2/P1), gamma, M1, or M2) must be unknown.
    The other parameters must be provided.

    Parameters:
      gamma (float): Specific heat ratio (None if unknown)
      M1 (float): Upstream Mach number (None if unknown)
      M2 (float): Downstream Mach number (None if unknown)
      P1 (float): Upstream pressure (None if unknown)
      P2 (float): Downstream pressure (None if unknown)
      P2_P1_ratio (float): Pressure ratio P2/P1 (None if unknown)

    Returns:
      float: The computed missing variable. When P2_P1_ratio is solved for, the result is the ratio;
             otherwise, the computed gamma, M1, or M2 is returned.

    Raises:
      ValueError: if not exactly one variable among the effective P2/P1 ratio, gamma, M1, and M2 is None
                  or if a math error occurs.
    """
    # Determine the effective pressure ratio.
    if P2_P1_ratio is not None:
        eff_ratio = P2_P1_ratio
    elif (P1 is not None or P2 is not None):
        if P1 is None or P2 is None:
            raise ValueError("Both P1 and P2 must be provided or neither, if P2_P1_ratio is not given.")
        if P1 == 0:
            raise ValueError("P1 is zero, cannot compute pressure ratio.")
        eff_ratio = P2 / P1
    else:
        eff_ratio = None

    # Count the number of unknowns among the independent quantities.
    unknowns = [eff_ratio is None, gamma is None, M1 is None, M2 is None]
    if sum(unknowns) != 1:
        raise ValueError("Exactly one variable among the effective P2/P1 ratio, gamma, M1, and M2 must be None.")

    # Case 1: Compute P2/P1 ratio when it is unknown.
    if eff_ratio is None:
        # gamma, M1, and M2 are provided.
        a = (gamma - 1) / 2
        R = (1 + a * M1**2) / (1 + a * M2**2)
        return R ** (gamma / (gamma - 1))

    # For the following cases, gamma is known.
    a = (gamma - 1) / 2

    # Case 2: Compute gamma when it is unknown.
    if gamma is None:
        # The equation to solve is:
        # eff_ratio = [ (1 + ((gamma - 1)/2)*M1^2) / (1 + ((gamma - 1)/2)*M2^2) ]^(gamma/(gamma-1))
        # This is a transcendental equation in gamma. Use the bisection method to solve for gamma.
        def f(g):
            if g <= 1:
                return float('inf')
            a_val = (g - 1) / 2
            R = (1 + a_val * M1**2) / (1 + a_val * M2**2)
            return R ** (g / (g - 1)) - eff_ratio
        low = 1 + 1e-6
        high = 50.0  # Upper bound selected based on typical values of gamma.
        if f(low) * f(high) > 0:
            raise ValueError("Cannot find a valid root for gamma in the interval.")
        tol = 1e-10
        while high - low > tol:
            mid = (low + high) / 2
            if f(low) * f(mid) <= 0:
                high = mid
            else:
                low = mid
        return (low + high) / 2

    # Case 3: Compute M1 when it is unknown.
    if M1 is None:
        # Starting from:
        # eff_ratio = [ (1 + a*M1^2) / (1 + a*M2^2) ]^(gamma/(gamma-1))
        # Taking both sides to the power (gamma-1)/gamma:
        # eff_ratio^((gamma-1)/gamma) = (1 + a*M1^2) / (1 + a*M2^2)
        base = eff_ratio ** ((gamma - 1) / gamma)
        numerator = base * (1 + a * M2**2) - 1
        if a == 0:
            raise ValueError("gamma equals 1, leading to division by zero.")
        if numerator < 0:
            raise ValueError("Computed M1^2 is negative.")
        return math.sqrt(numerator / a)

    # Case 4: Compute M2 when it is unknown.
    if M2 is None:
        base = eff_ratio ** ((gamma - 1) / gamma)
        numerator = (1 + a * M1**2) / base - 1
        if a == 0:
            raise ValueError("gamma equals 1, leading to division by zero.")
        if numerator < 0:
            raise ValueError("Computed M2^2 is negative.")
        return math.sqrt(numerator / a)




def main():
    print('prandtly_meyer.py was run but there is no code in main()')



if __name__ == '__main__':
    main()