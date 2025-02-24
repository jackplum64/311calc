import math

def deg2rad(angle_deg: float) -> float:
    return angle_deg * math.pi / 180


def rad2deg(angle_rad: float) -> float:
    return angle_rad * 180 / math.pi


def atm2pa(pres_atm: float) -> float:
    return pres_atm * 101325


def pa2atm(pres_pa: float) -> float:
    return pres_pa / 101325

def k2r(temp_K: float) -> float:
    return temp_K * 1.8

def r2k(temp_R: float) -> float:
    return temp_R / 1.8


def main():
    print('utils.py was run but there is no code in main()')


if __name__ == '__main__':
    main()