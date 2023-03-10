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


def cylindrical_explosion(r, r_in=0.8, r_out=1.0):
    vU = ixp.zerorank1()
    BU = ixp.zerorank1()

    r_out_minus_r_in = r_out - r_in
    r_out_minus_r = r_out - r
    r_minus_r_in = r - r_in

    rho_in = rl(1,100)
    rho_out = rl(1,10000)
    rho_mid = sp.exp( (r_out_minus_r*sp.ln(rho_in)  +  r_minus_r_in*sp.ln(rho_out) ) / (r_out_minus_r_in) )

    rho = noif.coord_leq_bound(r, r_in)*rho_in\
         +noif.coord_greater_bound(r, r_in)*noif.coord_less_bound(r, r_out)*rho_mid\
         +noif.coord_geq_bound(r, r_out)*rho_out


    press_in = rl(1.0)
    press_out = rl(3,100000)
    press_mid = sp.exp( (r_out_minus_r*sp.ln(press_in)  +  r_minus_r_in*sp.ln(press_out) ) / (r_out_minus_r_in) )

    press = noif.coord_leq_bound(r, r_in)*press_in\
           +noif.coord_greater_bound(r, r_in)*noif.coord_less_bound(r, r_out)*press_mid\
           +noif.coord_geq_bound(r, r_out)*press_out

    BU[0] = rl(1,10)

    return rho, press, vU, BU

def magnetic_rotor(r, r_in=0.1, Omega=9.95, cartx=rfm.Cartx, carty=rfm.Carty):
    cart_list = [cartx, carty]

    vU_cyl = ixp.zerorank1()
    BU = ixp.zerorank1()

    rho_in = rl(10.0)
    rho_out = rl(1.0)

    rho = noif.coord_leq_bound(r, r_in)*rho_in\
         +noif.coord_greater_bound(r, r_in)*rho_out

    press = rl(1.0)

    par.set_parval_from_str("reference_metric::CoordSystem","Cylindrical")
    rfm.reference_metric()

    vU_Cyl = ixp.zerorank1()
    vU_Cyl[1] = noif.coord_leq_bound(r, r_in)*Omega

    Jac_dUCart_dDrfmUD, Jac_dUrfm_dDCartUD = rfm.compute_Jacobian_and_inverseJacobian_tofrom_Cartesian()
    vU = rfm.basis_transform_vectorU_from_rfmbasis_to_Cartesian(Jac_dUCart_dDrfmUD, vU_Cyl)

    for i in range(DIM):
        for j in range(DIM):
            vU[i] = vU[i].subs(rfm.xx[j], rfm.Cart_to_xx[j])

    par.set_parval_from_str("reference_metric::CoordSystem","Cartesian")
    rfm.reference_metric()

    for i in range(DIM):
        for j in range(2):
            vU[i] = vU[i].subs(rfm.Cart[j], cart_list[j])

    BU[0] = rl(1.0)

    return rho, press, vU, BU

def loop_advection(r, R_loop=rl(3,10), A_loop=rl(1/1000), vz_equals_zero=False):
    vU = ixp.zerorank1()
    AD = ixp.zerorank1()

    rho = rl(1.0)

    press = rl(3.0)

    vx_in = rl(1,12)
    vy_in = rl(1,24)

    vU[0] = noif.coord_less_bound(r, R_loop)*vx_in
    vU[1] = noif.coord_less_bound(r, R_loop)*vy_in
    if not vz_equals_zero:
        vU[2] = noif.coord_less_bound(r, R_loop)*vy_in

    AD[2] = noif.max_noif(0, A_loop*(R_loop - r))

    return rho, press, vU, AD

