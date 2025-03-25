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


def normal_shock_rho(M1=None, rho1=None, rho2=None):
    """
    Solve for the missing variable in the normal shock density relation:
    
       rho2 / rho1 = (6 * M1^2) / (5 + M1^2)
       
    Exactly one of M1, rho1, or rho2 must be None.
    
    Parameters:
      M1 (float): Upstream Mach number (None if unknown)
      rho1 (float): Upstream density (None if unknown)
      rho2 (float): Downstream density (None if unknown)
      
    Returns:
      float: The computed value of the missing variable.
      
    Raises:
      ValueError: if not exactly one variable is None or if a math error occurs.
    """
    variables = {"M1": M1, "rho1": rho1, "rho2": rho2}
    missing = [name for name, value in variables.items() if value is None]
    if len(missing) != 1:
        raise ValueError("Exactly one variable must be None (the one to solve for).")
    
    if rho2 is None:
        # Solve for rho2: rho2 = rho1 * (6*M1^2) / (5 + M1^2)
        if M1 is None or rho1 is None:
            raise ValueError("Insufficient variables to solve for rho2.")
        return rho1 * (6 * M1**2) / (5 + M1**2)
    elif rho1 is None:
        # Solve for rho1: rho1 = rho2 / ((6*M1^2) / (5 + M1^2))
        if M1 is None or rho2 is None:
            raise ValueError("Insufficient variables to solve for rho1.")
        ratio = (6 * M1**2) / (5 + M1**2)
        if ratio == 0:
            raise ValueError("Division by zero encountered while solving for rho1.")
        return rho2 / ratio
    elif M1 is None:
        # Solve for M1 from: rho2/rho1 = (6*M1^2) / (5 + M1^2)
        if rho1 is None or rho2 is None:
            raise ValueError("Insufficient variables to solve for M1.")
        R = rho2 / rho1
        if (6 - R) == 0:
            raise ValueError("Division by zero encountered while solving for M1.")
        val = (5 * R) / (6 - R)
        if val < 0:
            raise ValueError("No real solution for M1 (negative value under square root).")
        return math.sqrt(val)


def normal_shock_P(M1=None, P1=None, P2=None):
    """
    Solve for the missing variable in the normal shock pressure relation:
    
       P2 / P1 = (7*M1^2 - 1) / 6
       
    Exactly one of M1, P1, or P2 must be None.
    
    Parameters:
      M1 (float): Upstream Mach number (None if unknown)
      P1 (float): Upstream pressure (None if unknown)
      P2 (float): Downstream pressure (None if unknown)
      
    Returns:
      float: The computed value of the missing variable.
      
    Raises:
      ValueError: if not exactly one variable is None or if a math error occurs.
    """
    variables = {"M1": M1, "P1": P1, "P2": P2}
    missing = [name for name, value in variables.items() if value is None]
    if len(missing) != 1:
        raise ValueError("Exactly one variable must be None (the one to solve for).")
    
    if P2 is None:
        # Solve for P2: P2 = P1 * ((7*M1^2 - 1)/6)
        if M1 is None or P1 is None:
            raise ValueError("Insufficient variables to solve for P2.")
        return P1 * ((7 * M1**2 - 1) / 6)
    elif P1 is None:
        # Solve for P1: P1 = P2 / ((7*M1^2 - 1)/6)
        if M1 is None or P2 is None:
            raise ValueError("Insufficient variables to solve for P1.")
        factor = (7 * M1**2 - 1) / 6
        if factor == 0:
            raise ValueError("Division by zero encountered while solving for P1.")
        return P2 / factor
    elif M1 is None:
        # Solve for M1: (7*M1^2 - 1)/6 = P2/P1  ->  M1^2 = (6*(P2/P1) + 1)/7
        if P1 is None or P2 is None:
            raise ValueError("Insufficient variables to solve for M1.")
        val = (6 * (P2 / P1) + 1) / 7
        if val < 0:
            raise ValueError("No real solution for M1 (negative value under square root).")
        return math.sqrt(val)


def normal_shock_T(M1=None, T1=None, T2=None):
    """
    Solve for the missing variable in the normal shock temperature relation:
    
       T2 / T1 = ((7*M1^2 - 1) * (5 + M1^2)) / (36 * M1^2)
       
    Exactly one of M1, T1, or T2 must be None.
    
    Parameters:
      M1 (float): Upstream Mach number (None if unknown)
      T1 (float): Upstream temperature (None if unknown)
      T2 (float): Downstream temperature (None if unknown)
      
    Returns:
      float: The computed value of the missing variable.
      
    Raises:
      ValueError: if not exactly one variable is None or if a math error occurs.
    """
    variables = {"M1": M1, "T1": T1, "T2": T2}
    missing = [name for name, value in variables.items() if value is None]
    if len(missing) != 1:
        raise ValueError("Exactly one variable must be None (the one to solve for).")
    
    if M1 is not None:
        if M1 == 0:
            raise ValueError("M1 cannot be zero.")
        expr = ((7 * M1**2 - 1) * (5 + M1**2)) / (36 * M1**2)
    
    if T2 is None:
        # Solve for T2: T2 = T1 * expression
        if M1 is None or T1 is None:
            raise ValueError("Insufficient variables to solve for T2.")
        return T1 * expr
    elif T1 is None:
        # Solve for T1: T1 = T2 / expression
        if M1 is None or T2 is None:
            raise ValueError("Insufficient variables to solve for T1.")
        if expr == 0:
            raise ValueError("Division by zero encountered while solving for T1.")
        return T2 / expr
    elif M1 is None:
        # Solve for M1 from: T2/T1 = ((7*M1^2 - 1)*(5 + M1^2))/(36*M1^2)
        # Let r = T2/T1. Then: 7*M1^4 + (34 - 36*r)*M1^2 - 5 = 0, with X = M1^2.
        r = T2 / T1
        a = 7
        b = 34 - 36 * r
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


def normal_shock_P0(P1=None, P2=None, T1=None, T2=None, P01=None, P02=None):
    """
    Solve for the missing variable in the normal shock total pressure relation:
    
       P02 / P01 = (P2 / P1) * (T1 / T2)^(7/2)
       
    Exactly one of P1, P2, T1, T2, P01, or P02 must be None.
    
    Parameters:
      P1 (float): Upstream static pressure (None if unknown)
      P2 (float): Downstream static pressure (None if unknown)
      T1 (float): Upstream temperature (None if unknown)
      T2 (float): Downstream temperature (None if unknown)
      P01 (float): Upstream total pressure (None if unknown)
      P02 (float): Downstream total pressure (None if unknown)
      
    Returns:
      float: The computed value of the missing variable.
      
    Raises:
      ValueError: if not exactly one variable is None or if a math error occurs.
    """
    variables = {"P1": P1, "P2": P2, "T1": T1, "T2": T2, "P01": P01, "P02": P02}
    missing = [name for name, value in variables.items() if value is None]
    if len(missing) != 1:
        raise ValueError("Exactly one variable must be None (the one to solve for).")
    
    # Relation: P02/P01 = (P2/P1) * (T1/T2)^(7/2)
    if P02 is None:
        if P01 is None or P1 is None or P2 is None or T1 is None or T2 is None:
            raise ValueError("Insufficient variables to solve for P02.")
        return P01 * (P2 / P1) * (T1 / T2)**(7/2)
    elif P01 is None:
        if P02 is None or P1 is None or P2 is None or T1 is None or T2 is None:
            raise ValueError("Insufficient variables to solve for P01.")
        factor = (P2 / P1) * (T1 / T2)**(7/2)
        if factor == 0:
            raise ValueError("Division by zero encountered while solving for P01.")
        return P02 / factor
    elif P2 is None:
        if P1 is None or P01 is None or P02 is None or T1 is None or T2 is None:
            raise ValueError("Insufficient variables to solve for P2.")
        factor = (T1 / T2)**(7/2)
        return P1 * (P02 / P01) / factor
    elif P1 is None:
        if P2 is None or P01 is None or P02 is None or T1 is None or T2 is None:
            raise ValueError("Insufficient variables to solve for P1.")
        factor = (T1 / T2)**(7/2)
        if (P02 / P01) * factor == 0:
            raise ValueError("Division by zero encountered while solving for P1.")
        return P2 / ((P02 / P01) * factor)
    elif T2 is None:
        if T1 is None or P1 is None or P2 is None or P01 is None or P02 is None:
            raise ValueError("Insufficient variables to solve for T2.")
        # Solve for T2: (T1/T2)^(7/2) = (P02/P01) / (P2/P1)
        ratio = (P02 / P01) / (P2 / P1)
        return T1 / (ratio**(2/7))
    elif T1 is None:
        if T2 is None or P1 is None or P2 is None or P01 is None or P02 is None:
            raise ValueError("Insufficient variables to solve for T1.")
        ratio = (P02 / P01) / (P2 / P1)
        return T2 * (ratio**(2/7))


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
