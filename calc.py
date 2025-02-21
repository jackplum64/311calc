import math

from isentropic import (isentropic_T_T0, isentropic_mach_T_T0, 
                       isentropic_P_P0, isentropic_mach_P_P0,
                       isentropic_rho_rho0, isentropic_mach_rho_rho0)

from shock import (mach_angle, mach1_mach2, normal_shock_P,
                   normal_shock_P0, normal_shock_T, normal_shock_rho, pitot_tube)

from utils import pa2atm, atm2pa, rad2deg, deg2rad


def bernoulli_eqn(P=None, P0=None, rho=None, U=None):
    """
    Solve for the missing variable in the relation:

        P0 = P + 0.5 * rho * U^2

    Exactly one of P, P0, rho, or U must be None.

    Parameters:
      P (float): Static Pressure (None if unknown)
      P0 (float): Total Pressure (None if unknown)
      rho (float): Static density (None if unknown)
      U (float): Flow velocity (None if unknown)
      
    Returns:
      float: The computed value of the missing variable.
      
    Raises:
      ValueError: if not exactly one variable is None or if a math error occurs.
    """
    variables = {"P": P, "P0": P0, "rho": rho, "U": U}
    missing = [name for name, value in variables.items() if value is None]
    if len(missing) != 1:
        raise ValueError("Exactly one variable must be None (the one to solve for).")

    if P0 is None:
        # Solve for P0: P0 = P + 0.5 * rho * U^2
        if P is None or rho is None or U is None:
            raise ValueError("Insufficient variables to solve for P0.")
        return P + 0.5 * rho * U**2
    elif P is None:
        # Solve for P: P = P0 - 0.5 * rho * U^2
        if P0 is None or rho is None or U is None:
            raise ValueError("Insufficient variables to solve for P.")
        return P0 - 0.5 * rho * U**2
    elif rho is None:
        # Solve for rho: rho = 2*(P0 - P) / U^2
        if P0 is None or P is None or U is None:
            raise ValueError("Insufficient variables to solve for rho.")
        if U == 0:
            raise ValueError("U cannot be zero when solving for rho (division by zero).")
        rho_value = 2 * (P0 - P) / U**2
        if rho_value < 0:
            raise ValueError("Computed density is negative. Check the input values for consistency.")
        return rho_value
    elif U is None:
        # Solve for U: U = sqrt(2*(P0 - P)/rho)
        if P0 is None or P is None or rho is None:
            raise ValueError("Insufficient variables to solve for U.")
        if rho == 0:
            raise ValueError("rho cannot be zero when solving for U (division by zero).")
        val = 2 * (P0 - P) / rho
        if val < 0:
            raise ValueError("No real solution for U (negative value under square root).")
        return math.sqrt(val)


def mach_eqn(M=None, U=None, gamma=None, R=None, T=None):
    """
    Solve for the missing variable in the Mach relation:

        M = U / sqrt(gamma * R * T)

    Exactly one of M, U, gamma, R, or T must be None.

    Parameters:
      M (float): Mach number (None if unknown)
      U (float): Flow velocity (None if unknown)
      gamma (float): Specific heat ratio (None if unknown)
      R (float): Specific gas constant (None if unknown)
      T (float): Static temperature (None if unknown)

    Returns:
      float: The computed value of the missing variable.

    Raises:
      ValueError: if not exactly one variable is None or if a math error occurs.
    """
    variables = {"M": M, "U": U, "gamma": gamma, "R": R, "T": T}
    missing = [name for name, value in variables.items() if value is None]
    if len(missing) != 1:
        raise ValueError("Exactly one variable must be None (the one to solve for).")
    
    if M is None:
        # Solve for M: M = U / sqrt(gamma * R * T)
        if U is None or gamma is None or R is None or T is None:
            raise ValueError("Insufficient variables to solve for M.")
        if gamma <= 0 or R <= 0 or T <= 0:
            raise ValueError("gamma, R, and T must be positive.")
        return U / math.sqrt(gamma * R * T)
    elif U is None:
        # Solve for U: U = M * sqrt(gamma * R * T)
        if M is None or gamma is None or R is None or T is None:
            raise ValueError("Insufficient variables to solve for U.")
        if gamma <= 0 or R <= 0 or T <= 0:
            raise ValueError("gamma, R, and T must be positive.")
        return M * math.sqrt(gamma * R * T)
    elif gamma is None:
        # Solve for gamma: gamma = (U / M)^2 / (R * T)
        if M is None or U is None or R is None or T is None:
            raise ValueError("Insufficient variables to solve for gamma.")
        if M == 0:
            raise ValueError("M cannot be zero when solving for gamma.")
        if R <= 0 or T <= 0:
            raise ValueError("R and T must be positive.")
        return (U / M)**2 / (R * T)
    elif R is None:
        # Solve for R: R = (U / M)^2 / (gamma * T)
        if M is None or U is None or gamma is None or T is None:
            raise ValueError("Insufficient variables to solve for R.")
        if M == 0:
            raise ValueError("M cannot be zero when solving for R.")
        if gamma <= 0 or T <= 0:
            raise ValueError("gamma and T must be positive.")
        return (U / M)**2 / (gamma * T)
    elif T is None:
        # Solve for T: T = U^2 / (M^2 * gamma * R)
        if M is None or U is None or gamma is None or R is None:
            raise ValueError("Insufficient variables to solve for T.")
        if M == 0:
            raise ValueError("M cannot be zero when solving for T.")
        if gamma <= 0 or R <= 0:
            raise ValueError("gamma and R must be positive.")
        return U**2 / (M**2 * gamma * R)


def velocity_T0_T(U=None, Cp=None, T0=None, T=None):
    """
    Solve for the missing variable in the relation:

        U^2 = 2 * Cp * (T0 - T)

    Exactly one of U, Cp, T0, or T must be None.

    Parameters:
      U (float): Flow velocity (None if unknown)
      Cp (float): Specific heat at constant pressure (None if unknown)
      T0 (float): Total temperature (None if unknown)
      T (float): Static temperature (None if unknown)

    Returns:
      float: The computed value of the missing variable.

    Raises:
      ValueError: if not exactly one variable is None or if a math error occurs.
    """
    variables = {"U": U, "Cp": Cp, "T0": T0, "T": T}
    missing = [name for name, value in variables.items() if value is None]
    if len(missing) != 1:
        raise ValueError("Exactly one variable must be None (the one to solve for).")
    
    if U is None:
        # U = sqrt(2*Cp*(T0-T))
        if Cp is None or T0 is None or T is None:
            raise ValueError("Insufficient variables to solve for U.")
        value = 2 * Cp * (T0 - T)
        if value < 0:
            raise ValueError("No real solution for U (negative value under square root).")
        return math.sqrt(value)
    elif Cp is None:
        # Cp = U^2 / (2*(T0-T))
        if U is None or T0 is None or T is None:
            raise ValueError("Insufficient variables to solve for Cp.")
        denominator = 2 * (T0 - T)
        if denominator == 0:
            raise ValueError("T0-T cannot be zero when solving for Cp.")
        return U**2 / denominator
    elif T0 is None:
        # T0 = T + U^2/(2*Cp)
        if U is None or Cp is None or T is None:
            raise ValueError("Insufficient variables to solve for T0.")
        if Cp == 0:
            raise ValueError("Cp cannot be zero when solving for T0.")
        return T + U**2 / (2 * Cp)
    elif T is None:
        # T = T0 - U^2/(2*Cp)
        if U is None or Cp is None or T0 is None:
            raise ValueError("Insufficient variables to solve for T.")
        if Cp == 0:
            raise ValueError("Cp cannot be zero when solving for T.")
        return T0 - U**2 / (2 * Cp)


def Cv_rho_T(rho0=None, rho=None, T0=None, T=None, Cv=None, R=None):
    """
    Solve for the missing variable in the relation:
    
        (rho0 / rho) = (T0 / T)^(Cv / R)
    
    Exactly one of rho0, rho, T0, T, Cv, or R must be None.
    
    Parameters:
      rho0 (float): Total density (None if unknown)
      rho (float): Static density (None if unknown)
      T0 (float): Total temperature (None if unknown)
      T (float): Static temperature (None if unknown)
      Cv (float): Specific heat at constant volume (None if unknown)
      R (float): Specific gas constant (None if unknown)
      
    Returns:
      float: The computed value of the missing variable.
      
    Raises:
      ValueError: if not exactly one variable is None or if math errors occur.
    """
    variables = {"rho0": rho0, "rho": rho, "T0": T0, "T": T, "Cv": Cv, "R": R}
    missing = [name for name, value in variables.items() if value is None]
    if len(missing) != 1:
        raise ValueError("Exactly one variable must be None (the one to solve for).")
    
    if rho0 is None:
        if rho is None or T0 is None or T is None or Cv is None or R is None:
            raise ValueError("Insufficient variables to solve for rho0.")
        return rho * (T0 / T) ** (Cv / R)
    elif rho is None:
        if rho0 is None or T0 is None or T is None or Cv is None or R is None:
            raise ValueError("Insufficient variables to solve for rho.")
        return rho0 / (T0 / T) ** (Cv / R)
    elif T0 is None:
        if rho0 is None or rho is None or T is None or Cv is None or R is None:
            raise ValueError("Insufficient variables to solve for T0.")
        return T * (rho0 / rho) ** (R / Cv)
    elif T is None:
        if rho0 is None or rho is None or T0 is None or Cv is None or R is None:
            raise ValueError("Insufficient variables to solve for T.")
        return T0 / (rho0 / rho) ** (R / Cv)
    elif Cv is None:
        if rho0 is None or rho is None or T0 is None or T is None or R is None:
            raise ValueError("Insufficient variables to solve for Cv.")
        # Ensure arguments for logarithms are positive.
        if (rho0 / rho) <= 0 or (T0 / T) <= 0:
            raise ValueError("Temperatures and densities must be positive for logarithm.")
        return R * math.log(rho0 / rho) / math.log(T0 / T)
    elif R is None:
        if rho0 is None or rho is None or T0 is None or T is None or Cv is None:
            raise ValueError("Insufficient variables to solve for R.")
        # Ensure arguments for logarithms are positive.
        if (rho0 / rho) <= 0 or (T0 / T) <= 0:
            raise ValueError("Temperatures and densities must be positive for logarithm.")
        return Cv * math.log(T0 / T) / math.log(rho0 / rho)
    

def compute_entropies(Cp_metric, R_metric, T1, T2, P1, P2, s1_initial=0.0):
    """
    Compute the entropy at two states given the relation:
    
        Δs = Cp_metric * ln(T2/T1) - R_metric * ln(P2/P1)
        s₂ = s₁ + Δs
    
    Parameters:
      Cp_metric (float): Specific heat at constant pressure (in consistent units).
      R_metric (float): Specific gas constant (in consistent units).
      T1 (float): Temperature at state 1 (must be positive).
      T2 (float): Temperature at state 2 (must be positive).
      P1 (float): Pressure at state 1 (must be positive).
      P2 (float): Pressure at state 2 (must be positive).
      s1_initial (float, optional): Reference entropy at state 1 (default is 0.0).
      
    Returns:
      tuple: A tuple (s1, s2) where:
             s1 is the entropy at state 1 (s1_initial),
             s2 is the entropy at state 2.
             
    Raises:
      ValueError: If any temperature or pressure is non-positive.
    """
    if T1 <= 0 or T2 <= 0:
        raise ValueError("Temperatures must be positive.")
    if P1 <= 0 or P2 <= 0:
        raise ValueError("Pressures must be positive.")
    
    # Compute the change in entropy using the given relation.
    delta_s = Cp_metric * math.log(T2 / T1) - R_metric * math.log(P2 / P1)
    
    s1 = s1_initial
    s2 = s1 + delta_s
    return s1, s2


def Cp_P_T(P0=None, P=None, T0=None, T=None, Cp=None, R=None):
    """
    Solve for the missing variable in the relation:
    
        (P0 / P) = (T0 / T)^(Cp / R)
    
    Exactly one of P0, P, T0, T, Cp, or R must be None.
    
    Parameters:
      P0 (float): Total pressure (None if unknown)
      P (float): Static pressure (None if unknown)
      T0 (float): Total temperature (None if unknown)
      T (float): Static temperature (None if unknown)
      Cp (float): Specific heat at constant pressure (None if unknown)
      R (float): Specific gas constant (None if unknown)
      
    Returns:
      float: The computed value of the missing variable.
      
    Raises:
      ValueError: if not exactly one variable is None or if math errors occur.
    """
    variables = {"P0": P0, "P": P, "T0": T0, "T": T, "Cp": Cp, "R": R}
    missing = [name for name, value in variables.items() if value is None]
    if len(missing) != 1:
        raise ValueError("Exactly one variable must be None (the one to solve for).")
    
    if P0 is None:
        if P is None or T0 is None or T is None or Cp is None or R is None:
            raise ValueError("Insufficient variables to solve for P0.")
        return P * (T0 / T) ** (Cp / R)
    elif P is None:
        if P0 is None or T0 is None or T is None or Cp is None or R is None:
            raise ValueError("Insufficient variables to solve for P.")
        return P0 / (T0 / T) ** (Cp / R)
    elif T0 is None:
        if P0 is None or P is None or T is None or Cp is None or R is None:
            raise ValueError("Insufficient variables to solve for T0.")
        return T * (P0 / P) ** (R / Cp)
    elif T is None:
        if P0 is None or P is None or T0 is None or Cp is None or R is None:
            raise ValueError("Insufficient variables to solve for T.")
        return T0 / (P0 / P) ** (R / Cp)
    elif Cp is None:
        if P0 is None or P is None or T0 is None or T is None or R is None:
            raise ValueError("Insufficient variables to solve for Cp.")
        # Ensure arguments for logarithms are positive.
        if (P0 / P) <= 0 or (T0 / T) <= 0:
            raise ValueError("Temperatures and pressures must be positive for logarithm.")
        return R * math.log(P0 / P) / math.log(T0 / T)
    elif R is None:
        if P0 is None or P is None or T0 is None or T is None or Cp is None:
            raise ValueError("Insufficient variables to solve for R.")
        # Ensure arguments for logarithms are positive.
        if (P0 / P) <= 0 or (T0 / T) <= 0:
            raise ValueError("Temperatures and pressures must be positive for logarithm.")
        return Cp * math.log(T0 / T) / math.log(P0 / P)
    




def main():
    print(f'\n----- AerE 311 Calculator -----')

    

    #isentropic_T_T0(T0=None, T=None, gamma=None, M=None)
    #isentropic_P_P0(P0=None, P=None, T0=None, T=None, Cp=None, R=None)
    #isentropic_rho_rho0(rho0=None, rho=None, T0=None, T=None, Cv=None, R=None)
    #isentropic_mach_T_T0(M=None, T0=None, T=None)
    #isentropic_mach_P_P0(M=None, P0=None, P=None)
    #isentropic_mach_rho_rho0(M=None, rho0=None, rho=None)

    #bernoulli_eqn(P=None, P0=None, rho=None, U=None)
    #mach_eqn(M=None, U=None, gamma=None, R=None, T=None)
    #velocity_T0_T(U=None, Cp=None, T0=None, T=None)
    #Cv_rho_T(rho0=None, rho=None, T0=None, T=None, Cv=None, R=None)
    #compute_entropies(Cp_metric, R_metric, T1, T2, P1, P2, s1_initial=0.0)
    #Cp_P_T(P0=None, P=None, T0=None, T=None, Cp=None, R=None)

    #mach_angle(angle=None, mach=None)
    #mach1_mach2(M1=None, M2=None)
    #normal_shock_rho(M1=None, rho1=None, rho2=None)
    #normal_shock_P(M1=None, P1=None, P2=None)
    #normal_shock_T(M1=None, T1=None, T2=None)
    #normal_shock_P0(P1=None, P2=None, T1=None, T2=None, P01=None, P02=None)
    #pitot_tube(M=None, P1=None, P02=None)
    








    
    



if __name__ == '__main__':
    main()
