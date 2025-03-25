import math
import scipy.optimize as opt


def mach_angle(angle=None, mach=None):
    """ 
    Solve for the missing variable in the relation

        mu = asin(1 / M)

    Exactly one of angle or mach must be None.

    Parameters:
      angle (float): Mach angle in degrees (None if unknown)
      mach (float): Mach number (None if unknown)

    Returns:
      float: The computed mach angle (in degrees) or mach number.

    Raises:
      ValueError: if not exactly one variable is None or if a math error occurs.
    """
    if angle is None:
        if mach is None:
            raise ValueError("Insufficient variables to solve for angle.")
        return math.degrees(math.asin(1 / mach))
    elif mach is None:
        if angle is None:
            raise ValueError("Insufficient variables to solve for mach.")
        # math.sin expects angle in radians.
        return 1 / math.sin(math.radians(angle))
    

def mach_beta_relation(Mn1=None, M1=None, beta=None):
    """
    Solve for the missing variable in the relation

        Mn1 = M1 * sin(beta)

    Exactly one of Mn1, M1, or beta must be None.

    Parameters:
      Mn1 (float): The normal component of Mach number (None if unknown)
      M1 (float): The free-stream Mach number (None if unknown)
      beta (float): The shock angle in degrees (None if unknown)

    Returns:
      float: The computed value of the missing variable. If beta is computed,
             the result is in degrees.

    Raises:
      ValueError: if not exactly one variable is None or if a math error occurs.
    """
    # Count the number of None values among the parameters.
    none_count = sum(x is None for x in [Mn1, M1, beta])
    if none_count != 1:
        raise ValueError("Exactly one variable must be None.")

    if Mn1 is None:
        # Compute Mn1 = M1 * sin(beta)
        try:
            return M1 * math.sin(math.radians(beta))
        except Exception as e:
            raise ValueError(f"Math error when computing Mn1: {e}")
    elif M1 is None:
        # Compute M1 = Mn1 / sin(beta)
        try:
            sin_beta = math.sin(math.radians(beta))
            if sin_beta == 0:
                raise ValueError("sin(beta) is zero, cannot divide by zero.")
            return Mn1 / sin_beta
        except Exception as e:
            raise ValueError(f"Math error when computing M1: {e}")
    elif beta is None:
        # Compute beta = arcsin(Mn1 / M1)
        try:
            ratio = Mn1 / M1
            if not -1 <= ratio <= 1:
                raise ValueError("The ratio Mn1/M1 must be between -1 and 1.")
            return math.degrees(math.asin(ratio))
        except Exception as e:
            raise ValueError(f"Math error when computing beta: {e}")


def theta_beta_mach(theta=None, beta=None, M1=None, gamma=None):
    """
    Solve for the missing variable in the relation

        tan(theta) = 2 * cot(beta) * [(M1^2 * sin^2(beta) - 1) / (M1^2 * (gamma + cos(2*beta)) + 2)]

    Exactly one of theta, beta, M1, or gamma must be None.

    Parameters:
      theta (float): Flow deflection angle in degrees (None if unknown)
      beta (float): Shock angle in degrees (None if unknown)
      M1 (float): Free-stream Mach number (None if unknown)
      gamma (float): Specific heat ratio (None if unknown)

    Returns:
      float: The computed value of the missing variable. If theta or beta is computed,
             the result is in degrees.

    Raises:
      ValueError: if not exactly one variable is None or if a math error occurs.
    """
    # Count the number of None values among the parameters.
    none_count = sum(x is None for x in [theta, beta, M1, gamma])
    if none_count != 1:
        raise ValueError("Exactly one variable must be None.")

    # Case 1: Compute theta when it is unknown.
    if theta is None:
        try:
            beta_rad = math.radians(beta)
            term = (M1**2 * (math.sin(beta_rad))**2 - 1) / (M1**2 * (gamma + math.cos(2 * beta_rad)) + 2)
            tan_theta = 2 * (math.cos(beta_rad) / math.sin(beta_rad)) * term
            return math.degrees(math.atan(tan_theta))
        except Exception as e:
            raise ValueError(f"Math error when computing theta: {e}")

    # Case 2: Compute beta when it is unknown.
    elif beta is None:
        try:
            theta_rad = math.radians(theta)
            # Define the function f(beta) = 2*cot(beta)*[(M1^2*sin^2(beta)-1)/(M1^2*(gamma+cos(2*beta))+2)] - tan(theta)
            def f(beta_rad: float) -> float:
                # Avoid division by zero.
                if math.sin(beta_rad) == 0:
                    return float('inf')
                cot_beta = math.cos(beta_rad) / math.sin(beta_rad)
                return 2 * cot_beta * ((M1**2 * math.sin(beta_rad)**2 - 1) /
                                         (M1**2 * (gamma + math.cos(2 * beta_rad)) + 2)) - math.tan(theta_rad)
            # For a valid oblique shock, beta must be greater than the Mach angle.
            beta_lower = math.asin(1 / M1)
            beta_upper = math.pi / 2 - 1e-6  # Slightly less than 90° to avoid singularity.
            if f(beta_lower) * f(beta_upper) > 0:
                raise ValueError("Cannot find a valid root for beta in the interval.")
            # Use the bisection method to solve for beta.
            tol = 1e-10
            while beta_upper - beta_lower > tol:
                beta_mid = (beta_lower + beta_upper) / 2
                if f(beta_lower) * f(beta_mid) <= 0:
                    beta_upper = beta_mid
                else:
                    beta_lower = beta_mid
            beta_rad_sol = (beta_lower + beta_upper) / 2
            return math.degrees(beta_rad_sol)
        except Exception as e:
            raise ValueError(f"Math error when computing beta: {e}")

    # Case 3: Compute M1 when it is unknown.
    elif M1 is None:
        try:
            beta_rad = math.radians(beta)
            tan_theta = math.tan(math.radians(theta))
            cot_beta = math.cos(beta_rad) / math.sin(beta_rad)
            # Rearranged relation yields:
            # M1^2 = [-2*cot(beta) - 2*tan(theta)] / [tan(theta)*(gamma + cos(2*beta)) - 2*cot(beta)*sin^2(beta)]
            numerator = -2 * cot_beta - 2 * tan_theta
            denominator = tan_theta * (gamma + math.cos(2 * beta_rad)) - 2 * cot_beta * (math.sin(beta_rad)**2)
            if denominator == 0:
                raise ValueError("Denominator is zero when computing M1.")
            M1_sq = numerator / denominator
            if M1_sq <= 0:
                raise ValueError("Computed M1^2 is non-positive.")
            return math.sqrt(M1_sq)
        except Exception as e:
            raise ValueError(f"Math error when computing M1: {e}")

    # Case 4: Compute gamma when it is unknown.
    elif gamma is None:
        try:
            beta_rad = math.radians(beta)
            tan_theta = math.tan(math.radians(theta))
            cot_beta = math.cos(beta_rad) / math.sin(beta_rad)
            # Starting from the relation and solving for gamma:
            # tan(theta) * (M1^2*(gamma + cos(2*beta)) + 2) = 2*cot(beta)*(M1^2*sin^2(beta) - 1)
            # Rearranged to yield:
            # gamma = [2*M1^2*sin^2(beta)*cot(beta) - 2*cot(beta) - M1^2*cos(2*beta)*tan(theta) - 2*tan(theta)] / (M1^2*tan(theta))
            numerator = 2 * M1**2 * (math.sin(beta_rad)**2) * cot_beta - 2 * cot_beta - M1**2 * math.cos(2 * beta_rad) * tan_theta - 2 * tan_theta
            denominator = M1**2 * tan_theta
            if denominator == 0:
                raise ValueError("Denominator is zero when computing gamma.")
            return numerator / denominator
        except Exception as e:
            raise ValueError(f"Math error when computing gamma: {e}")
        

def solve_beta(theta: float, M1: float, gamma: float, branch: str = 'weak') -> float:
    """
    Solve for the shock angle beta using the theta-beta-Mach relation for the specified shock branch.

        tan(theta) = 2 * cot(beta) * [(M1^2 * sin^2(beta) - 1) / (M1^2 * (gamma + cos(2*beta)) + 2)]

    Parameters:
      theta (float): Flow deflection angle in degrees.
      M1 (float): Free-stream Mach number.
      gamma (float): Specific heat ratio.
      branch (str): 'weak' for the weak shock solution or 'strong' for the strong shock solution.

    Returns:
      float: Shock angle beta in degrees.

    Raises:
      ValueError: if an invalid branch is provided or a math error occurs.
    """
    # Define the shock relation function f(beta) = 0.
    def shock_eq(beta: float) -> float:
        beta_rad = math.radians(beta)
        theta_rad = math.radians(theta)
        term = (M1**2 * (math.sin(beta_rad))**2 - 1) / (M1**2 * (gamma + math.cos(2 * beta_rad)) + 2)
        return 2 * (math.cos(beta_rad) / math.sin(beta_rad)) * term - math.tan(theta_rad)

    # Choose the initial guess based on the desired shock branch.
    if branch == 'weak':
        beta_guess = theta + 1.0  # A small increment above theta for the weak solution.
    elif branch == 'strong':
        beta_guess = 80.0  # An initial guess near 90° for the strong solution.
    else:
        raise ValueError("Invalid branch specified. Use 'weak' or 'strong'.")

    beta_solution = opt.fsolve(shock_eq, beta_guess)
    return beta_solution[0]
        

def machn1_machn2_relation(Mn1=None, Mn2=None, gamma=None):
    """
    Solve for the missing variable in the relation

        Mn2^2 = (1 + ((gamma - 1)/2) * Mn1^2) / (gamma * Mn1^2 - ((gamma - 1)/2))

    Exactly one of Mn1, Mn2, or gamma must be None.

    Parameters:
      Mn1 (float): The upstream normal Mach number (None if unknown)
      Mn2 (float): The downstream normal Mach number (None if unknown)
      gamma (float): The specific heat ratio (None if unknown)

    Returns:
      float: The computed value of the missing variable. For Mn1 and Mn2,
             the result is positive.

    Raises:
      ValueError: if not exactly one variable is None, if a math error occurs,
                  or if the computed value is non-physical.
    """
    # Count the number of None values among the parameters.
    none_count = sum(x is None for x in [Mn1, Mn2, gamma])
    if none_count != 1:
        raise ValueError("Exactly one variable must be None.")

    # Define a helper constant function: a = (gamma - 1) / 2
    try:
        # Case 1: Compute Mn1 when it is unknown.
        if Mn1 is None:
            # Mn2 and gamma must be provided.
            a = (gamma - 1) / 2
            denominator = gamma * (Mn2**2) - a
            if denominator == 0:
                raise ValueError("Denominator is zero when computing Mn1.")
            # Rearranged relation:
            # Mn2^2 = (1 + a * Mn1^2) / (gamma * Mn1^2 - a)
            # Let x = Mn1^2, then: x = (1 + a * Mn2^2) / (gamma * Mn2^2 - a)
            x = (1 + a * (Mn2**2)) / (gamma * (Mn2**2) - a)
            if x <= 0:
                raise ValueError("Computed Mn1^2 is non-positive.")
            return math.sqrt(x)

        # Case 2: Compute Mn2 when it is unknown.
        elif Mn2 is None:
            # Mn1 and gamma must be provided.
            a = (gamma - 1) / 2
            denominator = gamma * (Mn1**2) - a
            if denominator == 0:
                raise ValueError("Denominator is zero when computing Mn2.")
            expression = (1 + a * (Mn1**2)) / denominator
            if expression < 0:
                raise ValueError("Computed Mn2^2 is negative.")
            return math.sqrt(expression)

        # Case 3: Compute gamma when it is unknown.
        elif gamma is None:
            numerator = 2 - Mn1**2 - Mn2**2
            denominator = 2 * Mn2**2 * Mn1**2 - Mn2**2 - Mn1**2
            if denominator == 0:
                raise ValueError("Denominator is zero when computing gamma.")
            return numerator / denominator

    except Exception as e:
        raise ValueError(f"Math error when computing the missing variable: {e}")


def oblique_rho2_rho1(rho2=None, rho1=None, Mn1=None, gamma=None):
    """
    Solve for the missing variable in the relation

        rho2 / rho1 = ((gamma + 1) * Mn1^2) / (2 + (gamma - 1) * Mn1^2)

    Exactly one of rho2, rho1, Mn1, or gamma must be None.

    Parameters:
      rho2 (float): Downstream density (None if unknown)
      rho1 (float): Upstream density (None if unknown)
      Mn1 (float): Upstream normal Mach number (None if unknown)
      gamma (float): Specific heat ratio (None if unknown)

    Returns:
      float: The computed value of the missing variable. In the case of rho2 or rho1,
             the result is a density; for Mn1, a Mach number; and for gamma, a dimensionless value.

    Raises:
      ValueError: if not exactly one variable is None, if a math error occurs, or if the computed value is non-physical.
    """
    # Count the number of None values among the parameters.
    none_count = sum(x is None for x in [rho2, rho1, Mn1, gamma])
    if none_count != 1:
        raise ValueError("Exactly one variable must be None.")

    try:
        # If rho2/rho1 ratio is the unknown, compute the ratio first.
        if rho2 is None:
            # rho2 = (rho2/rho1) * rho1, where the ratio is given by:
            ratio = ((gamma + 1) * Mn1**2) / (2 + (gamma - 1) * Mn1**2)
            return ratio * rho1

        elif rho1 is None:
            # rho1 = rho2 / (rho2/rho1)
            ratio = ((gamma + 1) * Mn1**2) / (2 + (gamma - 1) * Mn1**2)
            if ratio == 0:
                raise ValueError("Density ratio computed as zero, cannot divide by zero.")
            return rho2 / ratio

        elif Mn1 is None:
            # Solve for Mn1 from: (rho2/rho1) = ((gamma + 1) * Mn1^2) / (2 + (gamma - 1) * Mn1^2)
            # Let R = rho2/rho1. Then:
            # R * (2 + (gamma - 1)*Mn1^2) = (gamma + 1)*Mn1^2
            # => 2*R = Mn1^2 * [(gamma + 1) - R*(gamma - 1)]
            # => Mn1^2 = (2 * R) / [(gamma + 1) - R*(gamma - 1)]
            R = rho2 / rho1
            denominator = (gamma + 1) - R * (gamma - 1)
            if denominator == 0:
                raise ValueError("Denominator is zero when computing Mn1.")
            Mn1_sq = (2 * R) / denominator
            if Mn1_sq <= 0:
                raise ValueError("Computed Mn1^2 is non-positive.")
            return math.sqrt(Mn1_sq)

        elif gamma is None:
            # Solve for gamma from: R = ((gamma + 1)*Mn1^2) / (2 + (gamma - 1)*Mn1^2)
            # Let R = rho2/rho1. Multiply both sides by the denominator:
            # R*(2 + (gamma - 1)*Mn1^2) = (gamma + 1)*Mn1^2
            # => 2*R + R*(gamma - 1)*Mn1^2 = (gamma + 1)*Mn1^2
            # Rearranging:
            # R*(gamma - 1)*Mn1^2 - (gamma + 1)*Mn1^2 = -2*R
            # Factor Mn1^2: Mn1^2 * [gamma*(R - 1) - (R + 1)] = -2*R
            # => gamma*(R - 1) = (R + 1) - (2*R)/Mn1^2
            # => gamma = [(R + 1) - (2*R)/Mn1^2] / (R - 1)
            R = rho2 / rho1
            if R == 1:
                raise ValueError("rho2/rho1 equals 1; cannot solve for gamma uniquely.")
            return ((R + 1) - (2 * R) / (Mn1**2)) / (R - 1)

    except Exception as e:
        raise ValueError(f"Math error when computing the missing variable: {e}")


def oblique_P2_P1(P1=None, P2=None, gamma=None, Mn1=None):
    """
    Solve for the missing variable in the relation

        P2 / P1 = 1 + [2 * gamma / (gamma + 1)] * (Mn1^2 - 1)

    Exactly one of P1, P2, gamma, or Mn1 must be None.

    Parameters:
      P1 (float): Upstream pressure (None if unknown)
      P2 (float): Downstream pressure (None if unknown)
      gamma (float): Specific heat ratio (None if unknown)
      Mn1 (float): Upstream normal Mach number (None if unknown)

    Returns:
      float: The computed value of the missing variable. For P1 or P2, the result is a pressure;
             for Mn1, a Mach number; and for gamma, a dimensionless value.

    Raises:
      ValueError: if not exactly one variable is None, if a math error occurs, or if the computed value is non-physical.
    """
    # Count the number of None values among the parameters.
    none_count = sum(x is None for x in [P1, P2, gamma, Mn1])
    if none_count != 1:
        raise ValueError("Exactly one variable must be None.")

    try:
        # If P2/P1 ratio is the unknown, compute the ratio first.
        if P2 is None:
            # P2 = P1 * (1 + [2*gamma/(gamma+1)]*(Mn1^2 - 1))
            return P1 * (1 + (2 * gamma / (gamma + 1)) * (Mn1**2 - 1))

        elif P1 is None:
            # P1 = P2 / (1 + [2*gamma/(gamma+1)]*(Mn1^2 - 1))
            factor = 1 + (2 * gamma / (gamma + 1)) * (Mn1**2 - 1)
            if factor == 0:
                raise ValueError("Computed factor is zero when computing P1.")
            return P2 / factor

        elif Mn1 is None:
            # Solve for Mn1 from:
            # P2/P1 = 1 + [2*gamma/(gamma+1)]*(Mn1^2 - 1)
            # Let R = P2/P1. Then:
            # R - 1 = [2*gamma/(gamma+1)]*(Mn1^2 - 1)
            # => Mn1^2 = ((R - 1) * (gamma + 1))/(2*gamma) + 1
            R = P2 / P1
            Mn1_sq = ((R - 1) * (gamma + 1)) / (2 * gamma) + 1
            if Mn1_sq <= 0:
                raise ValueError("Computed Mn1^2 is non-positive.")
            return math.sqrt(Mn1_sq)

        elif gamma is None:
            # Solve for gamma from:
            # P2/P1 = 1 + [2*gamma/(gamma+1)]*(Mn1^2 - 1)
            # Let R = P2/P1. Multiply both sides by (gamma+1):
            # R*(gamma + 1) = (gamma + 1) + 2*gamma*(Mn1^2 - 1)
            # => R*gamma + R = gamma + 1 + 2*gamma*(Mn1^2 - 1)
            # Rearranging gamma terms:
            # gamma*(R - 1 - 2*(Mn1^2 - 1)) = 1 - R
            # => gamma = (1 - R) / (R - 1 - 2*(Mn1^2 - 1))
            denominator = (P2 / P1) - 1 - 2 * (Mn1**2 - 1)
            if denominator == 0:
                raise ValueError("Denominator is zero when computing gamma.")
            return (1 - (P2 / P1)) / denominator

    except Exception as e:
        raise ValueError(f"Math error when computing the missing variable: {e}")

def temperature_ratio(T1=None, T2=None, P1=None, P2=None, rho1=None, rho2=None,
                      T_ratio=None, P_ratio=None, rho_ratio=None):
    """
    Solve for the missing variable in the relation

        T2 / T1 = (P2 / P1) * (rho1 / rho2)

    The function may be used in two modes:
      1. Ratio mode: If T_ratio, P_ratio, or rho_ratio is provided (non-None), the corresponding ratio
         is used and the individual temperatures, pressures, or densities are ignored.
      2. Component mode: If both T1 and T2 (or both P1 and P2, or both rho1 and rho2) are provided,
         the corresponding ratio is computed.

    If both the individual values are unknown (i.e. T1 and T2 for temperature, etc.), the function
    returns the computed ratio.

    Exactly one of the three independent quantities (T_ratio, P_ratio, or rho_ratio) must be unknown
    (i.e. determined by the relation). The other two must be provided either directly or via their components.

    Parameters:
      T1 (float): Upstream temperature (None if unknown)
      T2 (float): Downstream temperature (None if unknown)
      P1 (float): Upstream pressure (None if unknown)
      P2 (float): Downstream pressure (None if unknown)
      rho1 (float): Upstream density (None if unknown)
      rho2 (float): Downstream density (None if unknown)
      T_ratio (float): Temperature ratio T2/T1 (None if unknown)
      P_ratio (float): Pressure ratio P2/P1 (None if unknown)
      rho_ratio (float): Density ratio rho1/rho2 (None if unknown)

    Returns:
      float: The computed missing ratio (T_ratio, P_ratio, or rho_ratio). In the case where T1 and T2 (or P1/P2 or rho1/rho2)
             are both unknown, the function returns the computed ratio.

    Raises:
      ValueError: if not exactly one independent ratio is unknown or if a math error occurs.
    """
    # Determine effective ratios from provided values.
    # For temperature: use T_ratio if provided; otherwise, if both T1 and T2 are provided, compute T2/T1.
    if T_ratio is not None:
        eff_T = T_ratio
    elif (T1 is not None and T2 is not None):
        if T1 == 0:
            raise ValueError("T1 is zero, cannot compute temperature ratio.")
        eff_T = T2 / T1
    else:
        eff_T = None

    # For pressure: use P_ratio if provided; otherwise, if both P1 and P2 are provided, compute P2/P1.
    if P_ratio is not None:
        eff_P = P_ratio
    elif (P1 is not None and P2 is not None):
        if P1 == 0:
            raise ValueError("P1 is zero, cannot compute pressure ratio.")
        eff_P = P2 / P1
    else:
        eff_P = None

    # For density: use rho_ratio if provided; otherwise, if both rho1 and rho2 are provided, compute rho1/rho2.
    # Note: The relation uses rho1/rho2 (not rho2/rho1).
    if rho_ratio is not None:
        eff_rho = rho_ratio
    elif (rho1 is not None and rho2 is not None):
        if rho2 == 0:
            raise ValueError("rho2 is zero, cannot compute density ratio.")
        eff_rho = rho1 / rho2
    else:
        eff_rho = None

    # Count how many of the effective ratios are None.
    none_count = sum(x is None for x in [eff_T, eff_P, eff_rho])
    if none_count != 1:
        raise ValueError("Exactly one independent ratio (T_ratio, P_ratio, or rho_ratio) must be unknown, with the other two provided or computable.")

    try:
        # Solve for the missing ratio using the relation: T2/T1 = (P2/P1) * (rho1/rho2)
        if eff_T is None:
            # Compute T_ratio = P_ratio * (rho_ratio)
            return eff_P * eff_rho
        elif eff_P is None:
            # Compute P_ratio = (T2/T1) / (rho1/rho2)
            if eff_rho == 0:
                raise ValueError("Density ratio is zero, cannot compute P_ratio.")
            return eff_T / eff_rho
        elif eff_rho is None:
            # Compute rho_ratio = (T2/T1) / (P2/P1)
            if eff_P == 0:
                raise ValueError("Pressure ratio is zero, cannot compute rho_ratio.")
            return eff_T / eff_P
    except Exception as e:
        raise ValueError(f"Math error when computing the missing ratio: {e}")
    


def phi_beta_theta(phi=None, beta=None, theta=None):
    """
    Solve for the missing variable in the relation

        phi = beta - theta

    Exactly one of phi, beta, or theta must be None.

    Parameters:
      phi (float): The angle difference in degrees (None if unknown)
      beta (float): The shock angle in degrees (None if unknown)
      theta (float): The flow deflection angle in degrees (None if unknown)

    Returns:
      float: The computed value of the missing variable in degrees.

    Raises:
      ValueError: if not exactly one variable is None or if a math error occurs.
    """
    # Count the number of None values among the parameters.
    none_count = sum(x is None for x in [phi, beta, theta])
    if none_count != 1:
        raise ValueError("Exactly one variable must be None.")

    try:
        if phi is None:
            # Compute phi = beta - theta
            return beta - theta
        elif beta is None:
            # Compute beta = phi + theta
            return phi + theta
        elif theta is None:
            # Compute theta = beta - phi
            return beta - phi
    except Exception as e:
        raise ValueError(f"Math error when computing the missing variable: {e}")
    

def M2_Mn2_relation(M2=None, Mn2=None, beta=None, theta=None):
    """
    Solve for the missing variable in the relation

        M2 = Mn2 / sin(beta - theta)

    Exactly one of M2, Mn2, beta, or theta must be None.

    Parameters:
      M2 (float): Downstream Mach number (None if unknown)
      Mn2 (float): Downstream normal Mach number (None if unknown)
      beta (float): Shock angle in degrees (None if unknown)
      theta (float): Flow deflection angle in degrees (None if unknown)

    Returns:
      float: The computed value of the missing variable. In the case of beta or theta,
             the result is in degrees.

    Raises:
      ValueError: if not exactly one variable is None, if a math error occurs,
                  or if the computed value is non-physical.
    """
    # Count the number of None values among the parameters.
    none_count = sum(x is None for x in [M2, Mn2, beta, theta])
    if none_count != 1:
        raise ValueError("Exactly one variable must be None.")

    # Case 1: Compute M2 when it is unknown.
    if M2 is None:
        try:
            delta = beta - theta
            delta_rad = math.radians(delta)
            sin_delta = math.sin(delta_rad)
            if sin_delta == 0:
                raise ValueError("sin(beta - theta) is zero, cannot compute M2.")
            return Mn2 / sin_delta
        except Exception as e:
            raise ValueError(f"Math error when computing M2: {e}")

    # Case 2: Compute Mn2 when it is unknown.
    elif Mn2 is None:
        try:
            return M2 * math.sin(math.radians(beta - theta))
        except Exception as e:
            raise ValueError(f"Math error when computing Mn2: {e}")

    # Case 3: Compute beta when it is unknown.
    elif beta is None:
        try:
            if M2 == 0:
                raise ValueError("M2 is zero, cannot compute arcsin(Mn2/M2).")
            ratio = Mn2 / M2
            if ratio < -1 or ratio > 1:
                raise ValueError("Invalid ratio Mn2/M2 for arcsin; it must be between -1 and 1.")
            angle_diff = math.degrees(math.asin(ratio))
            return theta + angle_diff
        except Exception as e:
            raise ValueError(f"Math error when computing beta: {e}")

    # Case 4: Compute theta when it is unknown.
    elif theta is None:
        try:
            if M2 == 0:
                raise ValueError("M2 is zero, cannot compute arcsin(Mn2/M2).")
            ratio = Mn2 / M2
            if ratio < -1 or ratio > 1:
                raise ValueError("Invalid ratio Mn2/M2 for arcsin; it must be between -1 and 1.")
            angle_diff = math.degrees(math.asin(ratio))
            return beta - angle_diff
        except Exception as e:
            raise ValueError(f"Math error when computing theta: {e}")









def main():
    print('oblique_shock.py was run but there is no code in main()')



if __name__ == '__main__':
    main()
