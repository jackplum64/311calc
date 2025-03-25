import math

    

def mach1_mach2(M1=None, M2=None):
    """
    Solve for the missing variable in the relation

        M2^2 = (1 + ((gamma-1)/2) * M1^2) / (gamma * M1^2 - ((gamma-1)/2))

    Exactly one of M1 or M2 must be None.

    Parameters:
      M1 (float): Mach number before the shock (None if unknown)
      M2 (float): Mach number after the shock (None if unknown)

    Returns:
      float: The computed M2 (if M1 was provided) or computed M1 (if M2 was provided).

    Raises:
      ValueError: if not exactly one variable is None, if division by zero occurs,
                  or if a math error occurs (e.g. negative value under a square root).
    """
    # Use a common value for the ratio of specific heats (air)
    gamma = 1.4

    # Ensure exactly one of M1 or M2 is unknown.
    if (M1 is None and M2 is None) or (M1 is not None and M2 is not None):
        raise ValueError("Exactly one of M1 or M2 must be None.")

    # If M1 is provided, solve for M2.
    if M1 is not None:
        denominator = gamma * M1**2 - (gamma - 1) / 2
        if denominator == 0:
            raise ValueError("Division by zero encountered when computing M2.")
        value = (1 + ((gamma - 1) / 2) * M1**2) / denominator
        if value < 0:
            raise ValueError("Computed M2^2 is negative, cannot take square root.")
        return math.sqrt(value)
    
    # If M2 is provided, solve for M1.
    else:
        denominator = gamma * M2**2 - (gamma - 1) / 2
        if denominator == 0:
            raise ValueError("Division by zero encountered when computing M1.")
        M1_squared = (1 + ((gamma - 1) / 2) * M2**2) / denominator
        if M1_squared < 0:
            raise ValueError("Computed M1^2 is negative, cannot take square root.")
        return math.sqrt(M1_squared)


def normal_shock_rho(M1=None, rho1=None, rho2=None, rho2_rho1_ratio=None):
    """
    Solve for the missing variable in the normal shock density relation:
    
       rho2 / rho1 = (6 * M1^2) / (5 + M1^2)
       
    The missing variable can be any one of M1, rho1, rho2, or the density ratio.
    The ratio may be provided explicitly via rho2_rho1_ratio.
    
    Parameters:
      M1 (float): Upstream Mach number (None if unknown)
      rho1 (float): Upstream density (None if unknown)
      rho2 (float): Downstream density (None if unknown)
      rho2_rho1_ratio (float): Density ratio (rho2/rho1) (None if unknown)
      
    Returns:
      float: The computed value of the missing variable.
      
    Raises:
      ValueError: if the input combination is invalid or a math error occurs.
    """
    # If M1 is given, compute the theoretical density ratio.
    theo_ratio = (6 * M1**2) / (5 + M1**2) if M1 is not None else None

    # Determine the provided density ratio.
    if rho2_rho1_ratio is not None:
        prop_ratio = rho2_rho1_ratio
    elif (rho1 is not None) and (rho2 is not None):
        prop_ratio = rho2 / rho1
    else:
        prop_ratio = None

    # --- Cases ---
    # (1) Only M1 is provided (both rho1 and rho2 and the ratio are missing): return ratio.
    if (M1 is not None) and (rho1 is None and rho2 is None and rho2_rho1_ratio is None):
        return theo_ratio

    # (2) M1 is missing but the density ratio is provided (either explicitly or via both densities):
    elif (M1 is None) and (prop_ratio is not None):
        if (6 - prop_ratio) == 0:
            raise ValueError("Division by zero encountered while solving for M1.")
        val = (5 * prop_ratio) / (6 - prop_ratio)
        if val < 0:
            raise ValueError("No real solution for M1 (negative value under square root).")
        return math.sqrt(val)

    # (3) M1 and rho1 are provided; compute missing rho2.
    elif (M1 is not None) and (rho1 is not None) and (rho2 is None):
        return rho1 * theo_ratio

    # (4) M1 and rho2 are provided; compute missing rho1.
    elif (M1 is not None) and (rho2 is not None) and (rho1 is None):
        if theo_ratio == 0:
            raise ValueError("Division by zero encountered while solving for rho1.")
        return rho2 / theo_ratio

    # (5) M1 is missing and both densities are provided; compute M1.
    elif (M1 is None) and (rho1 is not None) and (rho2 is not None):
        prop_ratio = rho2 / rho1
        if (6 - prop_ratio) == 0:
            raise ValueError("Division by zero encountered while solving for M1.")
        val = (5 * prop_ratio) / (6 - prop_ratio)
        if val < 0:
            raise ValueError("No real solution for M1 (negative value under square root).")
        return math.sqrt(val)

    else:
        raise ValueError("Invalid combination of inputs for normal_shock_rho.")


def normal_shock_P(M1=None, P1=None, P2=None, P2_P1_ratio=None):
    """
    Solve for the missing variable in the normal shock pressure relation:
    
       P2 / P1 = (7*M1^2 - 1) / 6
       
    The missing variable can be any one of M1, P1, P2, or the pressure ratio.
    The ratio may be provided explicitly via P2_P1_ratio.
    
    Parameters:
      M1 (float): Upstream Mach number (None if unknown)
      P1 (float): Upstream pressure (None if unknown)
      P2 (float): Downstream pressure (None if unknown)
      P2_P1_ratio (float): Pressure ratio (P2/P1) (None if unknown)
      
    Returns:
      float: The computed value of the missing variable.
      
    Raises:
      ValueError: if the input combination is invalid or a math error occurs.
    """
    # If M1 is given, compute the theoretical pressure ratio.
    theo_ratio = (7 * M1**2 - 1) / 6 if M1 is not None else None

    # Determine the provided pressure ratio.
    if P2_P1_ratio is not None:
        prop_ratio = P2_P1_ratio
    elif (P1 is not None) and (P2 is not None):
        prop_ratio = P2 / P1
    else:
        prop_ratio = None

    # --- Cases ---
    # (1) Only M1 is provided: return the computed ratio.
    if (M1 is not None) and (P1 is None and P2 is None and P2_P1_ratio is None):
        return theo_ratio

    # (2) M1 is missing but the pressure ratio is provided: solve for M1.
    elif (M1 is None) and (prop_ratio is not None):
        val = (6 * prop_ratio + 1) / 7
        if val < 0:
            raise ValueError("No real solution for M1 (negative value under square root).")
        return math.sqrt(val)

    # (3) M1 and P1 provided; compute missing P2.
    elif (M1 is not None) and (P1 is not None) and (P2 is None):
        return P1 * theo_ratio

    # (4) M1 and P2 provided; compute missing P1.
    elif (M1 is not None) and (P2 is not None) and (P1 is None):
        if theo_ratio == 0:
            raise ValueError("Division by zero encountered while solving for P1.")
        return P2 / theo_ratio

    # (5) M1 missing and both pressures provided; compute M1.
    elif (M1 is None) and (P1 is not None) and (P2 is not None):
        prop_ratio = P2 / P1
        val = (6 * prop_ratio + 1) / 7
        if val < 0:
            raise ValueError("No real solution for M1 (negative value under square root).")
        return math.sqrt(val)

    else:
        raise ValueError("Invalid combination of inputs for normal_shock_P.")


def normal_shock_T(M1=None, T1=None, T2=None, T2_T1_ratio=None):
    """
    Solve for the missing variable in the normal shock temperature relation:
    
       T2 / T1 = ((7*M1^2 - 1) * (5 + M1^2)) / (36 * M1^2)
       
    The missing variable can be any one of M1, T1, T2, or the temperature ratio.
    The ratio may be provided explicitly via T2_T1_ratio.
    
    Parameters:
      M1 (float): Upstream Mach number (None if unknown)
      T1 (float): Upstream temperature (None if unknown)
      T2 (float): Downstream temperature (None if unknown)
      T2_T1_ratio (float): Temperature ratio (T2/T1) (None if unknown)
      
    Returns:
      float: The computed value of the missing variable.
      
    Raises:
      ValueError: if the input combination is invalid or a math error occurs.
    """
    if M1 is not None:
        if M1 == 0:
            raise ValueError("M1 cannot be zero.")
        theo_ratio = ((7 * M1**2 - 1) * (5 + M1**2)) / (36 * M1**2)
    else:
        theo_ratio = None

    # Determine the provided temperature ratio.
    if T2_T1_ratio is not None:
        prop_ratio = T2_T1_ratio
    elif (T1 is not None) and (T2 is not None):
        prop_ratio = T2 / T1
    else:
        prop_ratio = None

    # --- Cases ---
    # (1) Only M1 is provided: return computed temperature ratio.
    if (M1 is not None) and (T1 is None and T2 is None and T2_T1_ratio is None):
        return theo_ratio

    # (2) M1 is missing but the temperature ratio is provided: solve for M1.
    elif (M1 is None) and (prop_ratio is not None):
        # The relation can be rearranged to a quadratic in X = M1^2:
        #   7*X^2 + (34 - 36*prop_ratio)*X - 5 = 0.
        a = 7
        b = 34 - 36 * prop_ratio
        c = -5
        disc = b**2 - 4 * a * c
        if disc < 0:
            raise ValueError("No real solution for M1 (discriminant is negative).")
        X1 = (-b + math.sqrt(disc)) / (2 * a)
        X2 = (-b - math.sqrt(disc)) / (2 * a)
        X = None
        if X1 > 0 and X2 > 0:
            X = min(X1, X2)
        elif X1 > 0:
            X = X1
        elif X2 > 0:
            X = X2
        else:
            raise ValueError("No positive solution for M1.")
        return math.sqrt(X)

    # (3) M1 and T1 provided; compute missing T2.
    elif (M1 is not None) and (T1 is not None) and (T2 is None):
        return T1 * theo_ratio

    # (4) M1 and T2 provided; compute missing T1.
    elif (M1 is not None) and (T2 is not None) and (T1 is None):
        if theo_ratio == 0:
            raise ValueError("Division by zero encountered while solving for T1.")
        return T2 / theo_ratio

    # (5) M1 is missing and both T1 and T2 are provided; compute M1.
    elif (M1 is None) and (T1 is not None) and (T2 is not None):
        prop_ratio = T2 / T1
        a = 7
        b = 34 - 36 * prop_ratio
        c = -5
        disc = b**2 - 4 * a * c
        if disc < 0:
            raise ValueError("No real solution for M1 (discriminant is negative).")
        X1 = (-b + math.sqrt(disc)) / (2 * a)
        X2 = (-b - math.sqrt(disc)) / (2 * a)
        X = None
        if X1 > 0 and X2 > 0:
            X = min(X1, X2)
        elif X1 > 0:
            X = X1
        elif X2 > 0:
            X = X2
        else:
            raise ValueError("No positive solution for M1.")
        return math.sqrt(X)

    else:
        raise ValueError("Invalid combination of inputs for normal_shock_T.")


def normal_shock_P0(P1=None, P2=None, T1=None, T2=None, P01=None, P02=None,
                     P2_P1_ratio=None, T1_T2_ratio=None, P02_P01_ratio=None):
    """
    Solve for the missing variable in the normal shock total pressure relation:
    
       P02 / P01 = (P2 / P1) * (T1 / T2)^(7/2)
       
    The missing variable can be any one of P1, P2, T1, T2, P01, or P02,
    or the corresponding ratio (static, temperature, or total pressure).
    The ratios may be provided explicitly via P2_P1_ratio, T1_T2_ratio, or P02_P01_ratio.
    
    Parameters:
      P1 (float): Upstream static pressure (None if unknown)
      P2 (float): Downstream static pressure (None if unknown)
      T1 (float): Upstream temperature (None if unknown)
      T2 (float): Downstream temperature (None if unknown)
      P01 (float): Upstream total pressure (None if unknown)
      P02 (float): Downstream total pressure (None if unknown)
      P2_P1_ratio (float): Static pressure ratio (P2/P1) (None if unknown)
      T1_T2_ratio (float): Temperature ratio (T1/T2) (None if unknown)
      P02_P01_ratio (float): Total pressure ratio (P02/P01) (None if unknown)
      
    Returns:
      float: The computed value of the missing variable.
      
    Raises:
      ValueError: if the input combination is invalid or a math error occurs.
    """
    # Determine the static pressure ratio.
    if P2_P1_ratio is not None:
        static_ratio = P2_P1_ratio
    elif (P1 is not None) and (P2 is not None):
        static_ratio = P2 / P1
    else:
        static_ratio = None

    # Determine the temperature ratio (T1/T2).
    if T1_T2_ratio is not None:
        temp_ratio = T1_T2_ratio
    elif (T1 is not None) and (T2 is not None):
        temp_ratio = T1 / T2
    else:
        temp_ratio = None

    # Determine the total pressure ratio.
    if P02_P01_ratio is not None:
        total_ratio = P02_P01_ratio
    elif (P01 is not None) and (P02 is not None):
        total_ratio = P02 / P01
    else:
        total_ratio = None

    # --- Cases ---
    # (1) Missing P02: require P01, static ratio, and temperature ratio.
    if (P02 is None) and (P01 is not None) and (static_ratio is not None) and (temp_ratio is not None):
        return P01 * static_ratio * (temp_ratio)**(7/2)

    # (2) Missing P01.
    elif (P01 is None) and (P02 is not None) and (static_ratio is not None) and (temp_ratio is not None):
        factor = static_ratio * (temp_ratio)**(7/2)
        if factor == 0:
            raise ValueError("Division by zero encountered while solving for P01.")
        return P02 / factor

    # (3) Missing P2 (with P1 provided).
    elif (P2 is None) and (P1 is not None) and (P01 is not None) and (P02 is not None) and (temp_ratio is not None):
        # Compute the implied static ratio from the total relation.
        computed_static = (P02 / P01) / (temp_ratio)**(7/2)
        return P1 * computed_static

    # (4) Missing P1 (with P2 provided).
    elif (P1 is None) and (P2 is not None) and (P01 is not None) and (P02 is not None) and (temp_ratio is not None):
        computed_static = (P02 / P01) / (temp_ratio)**(7/2)
        if computed_static == 0:
            raise ValueError("Division by zero encountered while solving for P1.")
        return P2 / computed_static

    # (5) Missing T2 (with T1 provided).
    elif (T2 is None) and (T1 is not None) and (P1 is not None) and (P2 is not None) and (P01 is not None) and (P02 is not None):
        computed_static = P2 / P1
        if computed_static == 0:
            raise ValueError("Division by zero encountered while solving for T2.")
        # Invert the relation: (T1/T2)^(7/2) = (P02/P01)/(P2/P1)
        ratio_val = (P02 / P01) / computed_static
        return T1 / (ratio_val**(2/7))

    # (6) Missing T1 (with T2 provided).
    elif (T1 is None) and (T2 is not None) and (P1 is not None) and (P2 is not None) and (P01 is not None) and (P02 is not None):
        computed_static = P2 / P1
        if computed_static == 0:
            raise ValueError("Division by zero encountered while solving for T1.")
        ratio_val = (P02 / P01) / computed_static
        return T2 * (ratio_val**(2/7))

    # (7) If both P1 and P2 are missing (static group missing) but total and temperature groups are provided.
    elif (P1 is None and P2 is None) and (P01 is not None) and (P02 is not None) and (temp_ratio is not None):
        return (P02 / P01) / (temp_ratio)**(7/2)

    # (8) If both T1 and T2 are missing (temperature group missing) but static and total groups are provided.
    elif (T1 is None and T2 is None) and (P1 is not None) and (P2 is not None) and (P01 is not None) and (P02 is not None):
        return ((P02 / P01) / (P2 / P1))**(2/7)

    # (9) If both P01 and P02 are missing (total pressure group missing) but static and temperature groups are provided.
    # This branch expects individual pressures and temperatures.
    elif (P01 is None and P02 is None) and (P1 is not None) and (P2 is not None) and (T1 is not None) and (T2 is not None):
        return (P2 / P1) * (T1 / T2)**(7/2)

    # (10) New branch: If total pressure ratio is the missing variable and only the ratios are provided.
    elif (P01 is None and P02 is None and P02_P01_ratio is None) and (static_ratio is not None and temp_ratio is not None):
        return static_ratio * (temp_ratio)**(7/2)

    else:
        raise ValueError("Invalid combination of inputs for normal_shock_P0.")


def pitot_tube(M=None, P1=None, P02=None):
    """
    Solve for the missing variable in the pitot tube relation:
    
       P02 / P1 = (6^6 * M^7) / (5^(7/2) * (7*M^2 - 1)^(5/2))
       
    Exactly one of M, P1, or P02 must be None.
    
    Parameters:
      M (float): Mach number (None if unknown)
      P1 (float): Static pressure (None if unknown)
      P02 (float): Pitot (total) pressure (None if unknown)
      
    Returns:
      float: The computed value of the missing variable.
      
    Raises:
      ValueError: if not exactly one variable is None or if a math error occurs.
    """
    variables = {"M": M, "P1": P1, "P02": P02}
    missing = [name for name, value in variables.items() if value is None]
    if len(missing) != 1:
        raise ValueError("Exactly one variable must be None (the one to solve for).")
    
    if P02 is None:
        # Solve for P02: P02 = P1 * (6^6 * M^7) / (5^(7/2) * (7*M^2 - 1)^(5/2))
        if M is None or P1 is None:
            raise ValueError("Insufficient variables to solve for P02.")
        return P1 * (6**6 * M**7) / (5**(7/2) * (7 * M**2 - 1)**(5/2))
    
    elif P1 is None:
        # Solve for P1: P1 = P02 / [ (6^6 * M^7) / (5^(7/2) * (7*M^2 - 1)^(5/2)) ]
        if M is None or P02 is None:
            raise ValueError("Insufficient variables to solve for P1.")
        factor = (6**6 * M**7) / (5**(7/2) * (7 * M**2 - 1)**(5/2))
        if factor == 0:
            raise ValueError("Division by zero encountered while solving for P1.")
        return P02 / factor
    
    elif M is None:
        # Solve for M numerically from:
        # (P02/P1) = (6^6 * M^7) / (5^(7/2) * (7*M^2 - 1)^(5/2))
        if P1 is None or P02 is None:
            raise ValueError("Insufficient variables to solve for M.")
        target = P02 / P1

        def f(m):
            # Stay in the valid domain: 7*m^2 - 1 must be positive.
            if 7 * m**2 - 1 <= 0:
                return -target
            return (6**6 * m**7) / (5**(7/2) * (7 * m**2 - 1)**(5/2)) - target

        # The function F(M) has a minimum at m = sqrt(0.5). The MATLAB solution uses the principal value,
        # which corresponds to the upper branch (m >= sqrt(0.5)).
        m_min = math.sqrt(0.5)
        # Check if the target is below the minimum possible value.
        if abs(f(m_min)) < 1e-6:
            return m_min
        if f(m_min) > 0:
            raise ValueError("No solution exists: target is below the minimum value of the function.")
        
        lower = m_min
        upper = lower * 10
        # Expand the upper bound until f(upper) > 0
        while f(upper) < 0:
            upper *= 2
            if upper > 1e6:
                raise ValueError("Unable to bracket the root for M.")
        tol = 1e-6
        for _ in range(100):
            mid = (lower + upper) / 2
            f_mid = f(mid)
            if abs(f_mid) < tol:
                return mid
            if f(lower) * f_mid < 0:
                upper = mid
            else:
                lower = mid
        raise ValueError("Root finding did not converge for M.")


def main():
    print('normal_shock.py was run but there is no code in main()')



if __name__ == '__main__':
    main()
