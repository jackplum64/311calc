from isentropic import (isentropic_T_T0, isentropic_mach_T_T0, 
                       isentropic_P_P0, isentropic_mach_P_P0,
                       isentropic_rho_rho0, isentropic_mach_rho_rho0)


class State:
    def __init__(self, T=None, T0=None, P=None, P0=None, rho=None, rho0=None, velocity=None, mach=None):
        self.T = T
        self.T0 = T0
        self.P = P
        self.P0 = P0
        self.rho = rho
        self.rho0 = rho0
        self.velocity = velocity
        self.mach = mach

        self.R = 287.052874 # For dry air: J / (kg * K)
        self.Cv = 717.632185 # (5/2) * R
        self.Cp = 1004.68506 # (7/2) * R

    def solve_state_isentropic(self):
        """
        Solve for all state variables using isentropic relations

        Raises:
            ValueError: if state is not solvable
        """

        groups = [
            {
                "func": self.isentropic_mach_T_T0,
                "mapping": {"M": "mach", "T": "T", "T0": "T0"},
            },
            {
                "func": self.isentropic_mach_P_P0,
                "mapping": {"M": "mach", "P": "P", "P0": "P0"},
            },
            {
                "func": self.isentropic_mach_rho_rho0,
                "mapping": {"M": "mach", "rho": "rho", "rho0": "rho0"},
            },
            {
                "func": self.isentropic_T_T0,
                "mapping": {"M": "mach", "gamma": "gamma", "T": "T", "T0": "T0"},
            },
            {
                "func": self.isentropic_P_P0,
                "mapping": {"P0": "P0", "P": "P", "T0": "T0", "T": "T", "Cp": "Cp", "R": "R"},
            },
            {
                "func": self.isentropic_rho_rho0,
                "mapping": {"rho0": "rho0", "rho": "rho", "T0": "T0", "T": "T", "Cv": "Cv", "R": "R"},
            },
        ]

        solved_any = False
        for group in groups:
            func = group["func"]
            mapping = group["mapping"]
            # Build arguments from the current state values.
            args = {param: getattr(self, attr) for param, attr in mapping.items()}
            # Identify which variable is missing.
            missing = [(param, attr) for param, attr in mapping.items()
                    if getattr(self, attr) is None]
            if len(missing) == 1:
                # Exactly one variable is missing, so solve for it.
                param, attr = missing[0]
                setattr(self, attr, func(**args))
                solved_any = True

        








def main():
    print('state.py was run but there is no code in main()')


if __name__ == '__main__':
    main()