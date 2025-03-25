# In isentropic flow, both rho0 and T0 are constnat throughout the flow
import math

def isentropic_T_T0(T0=None, T=None, gamma=None, M=None):
    """
    Solve for the missing variable in the isentropic relation:
    
        T0 / T = 1 + ((gamma - 1) / 2) * M^2
        
    Exactly one of T0, T, gamma, or M must be None.
    
    Parameters:
      T0 (float): Total temperature (None if unknown)
      T (float): Static temperature (None if unknown)
      gamma (float): Specific heat ratio (None if unknown)
      M (float): Mach number (None if unknown)
      
    Returns:
      float: The computed value of the missing variable.
      
    Raises:
      ValueError: if not exactly one variable is None or if a math error occurs.
    """
    variables = {"T0": T0, "T": T, "gamma": gamma, "M": M}
    missing = [name for name, value in variables.items() if value is None]
    if len(missing) != 1:
        raise ValueError("Exactly one variable must be None (the one to solve for).")
    
    if T0 is None:
        # Given T, gamma, and M: T0 = T * (1 + ((gamma-1)/2)*M^2)
        if T is None or gamma is None or M is None:
            raise ValueError("Insufficient variables to solve for T0.")
        T0 = T * (1 + ((gamma - 1) / 2) * M**2)
        return T0
    elif T is None:
        # Given T0, gamma, and M: T = T0 / (1 + ((gamma-1)/2)*M^2)
        if T0 is None or gamma is None or M is None:
            raise ValueError("Insufficient variables to solve for T.")
        T = T0 / (1 + ((gamma - 1) / 2) * M**2)
        return T
    elif gamma is None:
        # Given T0, T, and M: gamma = 1 + (2*(T0/T - 1))/M^2
        if T0 is None or T is None or M is None:
            raise ValueError("Insufficient variables to solve for gamma.")
        if M == 0:
            raise ValueError("M cannot be zero when solving for gamma (division by zero).")
        gamma = 1 + (2 * (T0 / T - 1)) / (M**2)
        return gamma
    elif M is None:
        # Given T0, T, and gamma: M = sqrt( 2*(T0/T - 1)/(gamma-1) )
        if T0 is None or T is None or gamma is None:
            raise ValueError("Insufficient variables to solve for M.")
        if gamma == 1:
            raise ValueError("gamma must not equal 1 when solving for M (division by zero).")
        val = 2 * (T0 / T - 1) / (gamma - 1)
        if val < 0:
            raise ValueError("No real solution for M (negative under square root).")
        M = math.sqrt(val)
        return M

def isentropic_P_P0(P0=None, P=None, T0=None, T=None, Cp=None, R=None):
    """
    Solve for the missing variable in the isentropic relation:
    
        P0 / P = (T0 / T)^(Cp / R)
        
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
        P0 = P * (T0 / T)**(Cp / R)
        return P0
    elif P is None:
        if P0 is None or T0 is None or T is None or Cp is None or R is None:
            raise ValueError("Insufficient variables to solve for P.")
        P = P0 / (T0 / T)**(Cp / R)
        return P
    elif T0 is None:
        if P0 is None or P is None or T is None or Cp is None or R is None:
            raise ValueError("Insufficient variables to solve for T0.")
        T0 = T * (P0 / P)**(R / Cp)
        return T0
    elif T is None:
        if P0 is None or P is None or T0 is None or Cp is None or R is None:
            raise ValueError("Insufficient variables to solve for T.")
        T = T0 / (P0 / P)**(R / Cp)
        return T
    elif Cp is None:
        if P0 is None or P is None or T0 is None or T is None or R is None:
            raise ValueError("Insufficient variables to solve for Cp.")
        if T0 / T <= 0 or P0 / P <= 0:
            raise ValueError("Temperatures and pressures must be positive for logarithm.")
        Cp = R * math.log(P0 / P) / math.log(T0 / T)
        return Cp
    elif R is None:
        if P0 is None or P is None or T0 is None or T is None or Cp is None:
            raise ValueError("Insufficient variables to solve for R.")
        if T0 / T <= 0 or P0 / P <= 0:
            raise ValueError("Temperatures and pressures must be positive for logarithm.")
        R = Cp * math.log(T0 / T) / math.log(P0 / P)
        return R

def isentropic_rho_rho0(rho0=None, rho=None, T0=None, T=None, Cv=None, R=None):
    """
    Solve for the missing variable in the isentropic relation:
    
        rho0 / rho = (T0 / T)^(Cv / R)
        
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
        rho0 = rho * (T0 / T)**(Cv / R)
        return rho0
    elif rho is None:
        if rho0 is None or T0 is None or T is None or Cv is None or R is None:
            raise ValueError("Insufficient variables to solve for rho.")
        rho = rho0 / (T0 / T)**(Cv / R)
        return rho
    elif T0 is None:
        if rho0 is None or rho is None or T is None or Cv is None or R is None:
            raise ValueError("Insufficient variables to solve for T0.")
        T0 = T * (rho0 / rho)**(R / Cv)
        return T0
    elif T is None:
        if rho0 is None or rho is None or T0 is None or Cv is None or R is None:
            raise ValueError("Insufficient variables to solve for T.")
        T = T0 / (rho0 / rho)**(R / Cv)
        return T
    elif Cv is None:
        if rho0 is None or rho is None or T0 is None or T is None or R is None:
            raise ValueError("Insufficient variables to solve for Cv.")
        if T0 / T <= 0 or rho0 / rho <= 0:
            raise ValueError("Temperatures and densities must be positive for logarithm.")
        Cv = R * math.log(rho0 / rho) / math.log(T0 / T)
        return Cv
    elif R is None:
        if rho0 is None or rho is None or T0 is None or T is None or Cv is None:
            raise ValueError("Insufficient variables to solve for R.")
        if T0 / T <= 0 or rho0 / rho <= 0:
            raise ValueError("Temperatures and densities must be positive for logarithm.")
        R = Cv * math.log(T0 / T) / math.log(rho0 / rho)
        return R

def isentropic_mach_T_T0(M=None, T0=None, T=None):
    """
    Solve for the missing variable in the relation:
    
        M = sqrt(5 * (T0 / T - 1))
        
    Exactly one of M, T0, or T must be None.
    
    Parameters:
      M (float): Mach number (None if unknown)
      T0 (float): Total temperature (None if unknown)
      T (float): Static temperature (None if unknown)
      
    Returns:
      float: The computed value of the missing variable.
      
    Raises:
      ValueError: if not exactly one variable is None or if a math error occurs.
    """
    variables = {"M": M, "T0": T0, "T": T}
    missing = [name for name, value in variables.items() if value is None]
    if len(missing) != 1:
        raise ValueError("Exactly one variable must be None (the one to solve for).")
    
    if M is None:
        if T0 is None or T is None:
            raise ValueError("Insufficient variables to solve for M.")
        val = 5 * (T0 / T - 1)
        if val < 0:
            raise ValueError("No real solution for M (negative under square root).")
        M = math.sqrt(val)
        return M
    elif T0 is None:
        if M is None or T is None:
            raise ValueError("Insufficient variables to solve for T0.")
        # Rearranged from: M^2 = 5*(T0/T - 1)  =>  T0 = T*(1 + M^2/5)
        T0 = T * (1 + M**2 / 5)
        return T0
    elif T is None:
        if M is None or T0 is None:
            raise ValueError("Insufficient variables to solve for T.")
        T = T0 / (1 + M**2 / 5)
        return T

def isentropic_mach_P_P0(M=None, P0=None, P=None):
    """
    Solve for the missing variable in the relation:
    
        M = sqrt(5 * ((P0 / P)^(2/7) - 1))
        
    Exactly one of M, P0, or P must be None.
    
    Parameters:
      M (float): Mach number (None if unknown)
      P0 (float): Total pressure (None if unknown)
      P (float): Static pressure (None if unknown)
      
    Returns:
      float: The computed value of the missing variable.
      
    Raises:
      ValueError: if not exactly one variable is None or if a math error occurs.
    """
    variables = {"M": M, "P0": P0, "P": P}
    missing = [name for name, value in variables.items() if value is None]
    if len(missing) != 1:
        raise ValueError("Exactly one variable must be None (the one to solve for).")
    
    if M is None:
        if P0 is None or P is None:
            raise ValueError("Insufficient variables to solve for M.")
        val = 5 * ((P0 / P)**(2/7) - 1)
        if val < 0:
            raise ValueError("No real solution for M (negative under square root).")
        M = math.sqrt(val)
        return M
    elif P0 is None:
        if M is None or P is None:
            raise ValueError("Insufficient variables to solve for P0.")
        # From: M^2 = 5*((P0/P)^(2/7) - 1)  =>  (P0/P)^(2/7) = M^2/5 + 1
        P0 = P * (M**2 / 5 + 1)**(7/2)
        return P0
    elif P is None:
        if M is None or P0 is None:
            raise ValueError("Insufficient variables to solve for P.")
        P = P0 / (M**2 / 5 + 1)**(7/2)
        return P

def isentropic_mach_rho_rho0(M=None, rho0=None, rho=None):
    """
    Solve for the missing variable in the relation:
    
        M = sqrt(5 * ((rho0 / rho)^(2/5) - 1))
        
    Exactly one of M, rho0, or rho must be None.
    
    Parameters:
      M (float): Mach number (None if unknown)
      rho0 (float): Total density (None if unknown)
      rho (float): Static density (None if unknown)
      
    Returns:
      float: The computed value of the missing variable.
      
    Raises:
      ValueError: if not exactly one variable is None or if a math error occurs.
    """
    variables = {"M": M, "rho0": rho0, "rho": rho}
    missing = [name for name, value in variables.items() if value is None]
    if len(missing) != 1:
        raise ValueError("Exactly one variable must be None (the one to solve for).")
    
    if M is None:
        if rho0 is None or rho is None:
            raise ValueError("Insufficient variables to solve for M.")
        val = 5 * ((rho0 / rho)**(2/5) - 1)
        if val < 0:
            raise ValueError("No real solution for M (negative under square root).")
        M = math.sqrt(val)
        return M
    elif rho0 is None:
        if M is None or rho is None:
            raise ValueError("Insufficient variables to solve for rho0.")
        # From: M^2 = 5*((rho0/rho)^(2/5) - 1)  =>  (rho0/rho)^(2/5) = M^2/5 + 1
        rho0 = rho * (M**2 / 5 + 1)**(5/2)
        return rho0
    elif rho is None:
        if M is None or rho0 is None:
            raise ValueError("Insufficient variables to solve for rho.")
        rho = rho0 / (M**2 / 5 + 1)**(5/2)
        return rho


def main():
    print('isentropic.py was run but there is no code in main()')


if __name__ == '__main__':
    main()