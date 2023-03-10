# Step 0: Add NRPy's directory to the path
# https://stackoverflow.com/questions/16780014/import-file-from-parent-directory
import os,sys
nrpy_dir_path = os.path.join("..")
if nrpy_dir_path not in sys.path:
    sys.path.append(nrpy_dir_path)

# Step 0.a: Import the NRPy+ core modules and set the reference metric to Cartesian
import sympy as sp               # SymPy: The Python computer algebra package upon which NRPy+ depends
from sympy import Rational as rl
import NRPy_param_funcs as par   # NRPy+: Parameter interface
import indexedexp as ixp         # NRPy+: Symbolic indexed expression (e.g., tensors, vectors, etc.) support
import reference_metric as rfm   # NRPy+: Reference metric support
import Min_Max_and_Piecewise_Expressions as noif

par.set_parval_from_str("grid::DIM", 3)
DIM = par.parval_from_str("grid::DIM")

par.set_parval_from_str("reference_metric::CoordSystem","Cartesian")
rfm.reference_metric()


def balsara1(x, bound = 0.0):
    vU = ixp.zerorank1()
    BU = ixp.zerorank1()

    rho_left = rl(1.0)
    rho_right = rl(1,8)
    rho = noif.coord_less_bound(x, bound)*rho_left \
        + noif.coord_geq_bound(x,bound)*rho_right

    press_left = rl(1.0)
    press_right = rl(1,10)
    press = noif.coord_less_bound(x, bound)*press_left \
          + noif.coord_geq_bound(x,bound)*press_right

    BU[0] = rl(1,2)

    By_left = rl(1.0)
    By_right = rl(-1.0)
    BU[1] = noif.coord_less_bound(x, bound)*By_left \
          + noif.coord_geq_bound(x,bound)*By_right

    return rho, press, vU, BU

def balsara2(x, bound = 0.0):
    vU = ixp.zerorank1()
    BU = ixp.zerorank1()

    rho = rl(1.0)

    press_left = rl(30.0)
    press_right = rl(1.0)
    press = noif.coord_less_bound(x, bound)*press_left \
          + noif.coord_geq_bound(x,bound)*press_right

    BU[0] = rl(1,2)

    By_left = rl(6.0)
    By_right = rl(7,10)
    BU[1] = BU[2] = noif.coord_less_bound(x, bound)*By_left \
                  + noif.coord_geq_bound(x,bound)*By_right

    return rho, press, vU, BU

def balsara3(x, bound = 0.0):
    vU = ixp.zerorank1()
    BU = ixp.zerorank1()

    rho = rl(1.0)

    press_left = rl(1000.0)
    press_right = rl(1,10)
    press = noif.coord_less_bound(x, bound)*press_left \
          + noif.coord_geq_bound(x,bound)*press_right

    BU[0] = rl(10.)

    By_left = rl(7.0)
    By_right = rl(7,10)
    BU[1] = BU[2] = noif.coord_less_bound(x, bound)*By_left \
                  + noif.coord_geq_bound(x,bound)*By_right

    return rho, press, vU, BU


def balsara4(x, bound = 0.0):
    vU = ixp.zerorank1()
    BU = ixp.zerorank1()

    rho = rl(1.0)

    press = rl(1,10)

    vx_left = rl(999,1000)
    vx_right = rl(-999,1000)
    vU[0] = noif.coord_less_bound(x, bound)*vx_left \
          + noif.coord_geq_bound(x,bound)*vx_right

    BU[0] = rl(10.)

    By_left = rl(7.0)
    By_right = rl(-7.0)
    BU[1] = BU[2] = noif.coord_less_bound(x, bound)*By_left \
                  + noif.coord_geq_bound(x,bound)*By_right

    return rho, press, vU, BU

def balsara5(x, bound = 0.0):
    vU = ixp.zerorank1()
    BU = ixp.zerorank1()

    rho_left = rl(108,100)
    rho_right = rl(1.0)
    rho = noif.coord_less_bound(x, bound)*rho_left \
        + noif.coord_geq_bound(x,bound)*rho_right

    press_left = rl(95, 100)
    press_right = rl(1.0)
    press = noif.coord_less_bound(x, bound)*press_left \
          + noif.coord_geq_bound(x,bound)*press_right

    vx_left = rl(4,10)
    vx_right = rl(-45,100)
    vU[0] = noif.coord_less_bound(x, bound)*vx_left \
          + noif.coord_geq_bound(x,bound)*vx_right

    vy_left = rl(3,10)
    vy_right = rl(-2,10)
    vU[1] = noif.coord_less_bound(x, bound)*vy_left \
          + noif.coord_geq_bound(x,bound)*vy_right

    vU[2] = rl(2,10)

    BU[0] = rl(2.)

    By_left = rl(3,10)
    By_right = rl(-7,10)
    BU[1] = noif.coord_less_bound(x, bound)*By_left \
          + noif.coord_geq_bound(x,bound)*By_right

    Bz_left = rl(3,10)
    Bz_right = rl(5,10)
    BU[2] = noif.coord_less_bound(x, bound)*Bz_left \
          + noif.coord_geq_bound(x,bound)*Bz_right

    return rho, press, vU, BU


def BtoA_piecewise_constant_Cart_Flat(x, y, z, BU):
    AD = ixp.zerorank1()
    AD[0] = z*BU[1]
    AD[1] = x*BU[2]
    AD[2] = y*BU[0]

    return AD


def hydro1(x, bound = 0.5):
    vU = ixp.zerorank1()
    BU = ixp.zerorank1()

    rho_left = rl(10.)
    rho_right = rl(1.0)
    rho = noif.coord_less_bound(x, bound)*rho_left \
        + noif.coord_geq_bound(x,bound)*rho_right

    press_left = rl(40, 3)
    press_right = rl(0.0)
    press = noif.coord_less_bound(x, bound)*press_left \
          + noif.coord_geq_bound(x,bound)*press_right

    return rho, press, vU, BU


def hydro2(x, bound = 0.5):
    vU = ixp.zerorank1()
    BU = ixp.zerorank1()

    rho = rl(1.0)

    press_left = rl(1000.0)
    press_right = rl(1, 100)
    press = noif.coord_less_bound(x, bound)*press_left \
          + noif.coord_geq_bound(x,bound)*press_right

    return rho, press, vU, BU