from ..calculator.calc import *
from ..calculator.isentropic import *
from ..calculator.normal_shock import *
from ..calculator.utils import *


def hw4_1():
    R = 287.05
    gamma = 1.4
    Cv = (5 / 2) * R
    Cp = (7/2) * R
    
    T1 = 288
    M1 = 2.500
    P1 = 1 * 101325

    U = mach_eqn(M=M1, gamma=gamma, R=R, T=T1)
    T0 = isentropic_T_T0(T=T1, gamma=gamma, M=M1)
    P0 = isentropic_P_P0(P=P1, T0=T0, T=T1, Cp=Cp, R=R)
    rho1 = bernoulli_eqn(P=P1, P0=P0, U=U)


    P2 = normal_shock_P(M1=M1, P1=P1)
    T2 = normal_shock_T(M1=M1, T1=T1)
    rho2 = normal_shock_rho(M1=M1, rho1=rho1)
    M2 = mach1_mach2(M1=M1)

    P02 = normal_shock_P0(P1=P1, P2=P2, T1=T1, T2=T2, P01=P0)

    T02 = isentropic_T_T0(T=T2, gamma=gamma, M=M2)

    s1, s2 = compute_entropies(Cp, R, T1, T2, P1, P2, s1_initial=0.0)

    print(s1, s2)

    print("The properties just downstream of the shock are\n")
    print(f"p2 = {P2 / 101325:.2f} atm")
    print(f"T2 = {T2:.0f} K")
    print(f"ρ2 = {rho2:.3f} kg/m3")
    print(f"M2 = {M2:.3f}")
    print(f"p02 = {P02 / 101325:.2f} atm")
    print(f"T02 = {T02:.0f} K")


def hw4_2():
    R = 287.05
    gamma = 1.4
    Cv = (5 / 2) * R
    Cp = (7/2) * R

    P1 = 1 * 101325
    P2 = 10.33 * 101325
    T2 = 1340 * 5 / 9 # R to K

    M1 = normal_shock_P(P1=P1, P2=P2)
    T1 = normal_shock_T(M1=M1, T2=T2)

    M2 = mach1_mach2(M1=M1)

    T02 = isentropic_T_T0(T=T2, gamma=gamma, M=M2)
    P02 = isentropic_P_P0(P=P2, T0=T02, T=T2, Cp=Cp, R=R)

    print("\nResults for the second case:")
    print(f"The Mach number along upstream of the wave is {M1:.4f}")
    print(f"The temperature along upstream of the wave is {(T1 * 1.8):.3f} °R")
    print(f"The total temperature along downstream of the wave is {(T02 * 1.8):.3f} °R")
    print(f"The total pressure along downstream of the wave is {(P02 / 101325):.1f} atm")


def hw4_3():
    R = 287.05
    gamma = 1.4
    Cv = (5 / 2) * R
    Cp = (7/2) * R

    P02 = 1.325 * 101325
    P1 = 0.1 * 101325

    M = pitot_tube(P1=P1, P02=P02)
    print(f'The Mach number in the tube is {M:.3f}')


def main():
    hw4_1()
    hw4_2()
    hw4_3()


if __name__ == '__main__':
    main()