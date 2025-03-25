import math
from utils import (rad2deg, deg2rad)


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



def main():
    print('prandtly_meyer.py was run but there is no code in main()')



if __name__ == '__main__':
    main()