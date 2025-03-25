from isentropic import (isentropic_T_T0, isentropic_mach_T_T0, 
                       isentropic_P_P0, isentropic_mach_P_P0,
                       isentropic_rho_rho0, isentropic_mach_rho_rho0)
from normal_shock import (mach1_mach2, normal_shock_P,
                   normal_shock_P0, normal_shock_T, normal_shock_rho, pitot_tube)
from oblique_shock import (mach_angle, mach_beta_relation, theta_beta_mach,
                               oblique_rho2_rho1, oblique_P2_P1, temperature_ratio,
                               phi_beta_theta, M2_Mn2_relation)
from utils import (pa2atm, atm2pa, rad2deg, deg2rad, r2k, k2r)
from general import (bernoulli_eqn, mach_eqn, velocity_T0_T,
                     Cv_rho_T, compute_entropies, Cp_P_T)


def main():
    print(f'\n----- AerE 311 Calculator -----\n')




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

    #mach1_mach2(M1=None, M2=None)
    #normal_shock_rho(M1=None, rho1=None, rho2=None)
    #normal_shock_P(M1=None, P1=None, P2=None)
    #normal_shock_T(M1=None, T1=None, T2=None)
    #normal_shock_P0(P1=None, P2=None, T1=None, T2=None, P01=None, P02=None)
    #pitot_tube(M=None, P1=None, P02=None)

    ### || theta_beta_mach is deprecated, use solve_beta.||  ###
    #theta_beta_mach(theta=None, beta=None, M1=None, gamma=None)
    #solve_beta(theta: float, M1: float, gamma: float, branch: str = 'weak') -> float

    #mach_angle(angle=None, mach=None)oblique_rho2_rho1, oblique_P2_P1, temperature_ratio
    #oblique_rho2_rho1(rho2=None, rho1=None, Mn1=None, gamma=None)
    #oblique_P2_P1(P1=None, P2=None, gamma=None, Mn1=None)
    #temperature_ratio(T1=None, T2=None, P1=None, P2=None, rho1=None, rho2=None, T_ratio=None, P_ratio=None, rho_ratio=None)
    #phi_beta_theta(phi=None, beta=None, theta=None)
    #M2_Mn2_relation(M2=None, Mn2=None, beta=None, theta=None)

    #deg2rad(angle)
    #rad2deg(angle)
    #pa2atm(pressure)
    #atm2pa(pressure)
    #r2k(temperature)
    #k2r(temperature)
    








    
    



if __name__ == '__main__':
    main()
