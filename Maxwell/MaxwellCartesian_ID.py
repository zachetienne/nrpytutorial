import NRPy_param_funcs as par
import indexedexp as ixp
import grid as gri
import finite_difference as fin
from outputC import *

def MaxwellCartesian_ID():
    DIM = par.parval_from_str("grid::DIM")

    x,y,z = gri.register_gridfunctions("AUX",["x","y","z"])
    gammaDD = ixp.register_gridfunctions_for_single_rank2("AUX","gammaDD", "sym12") # The AUX or EVOL designation is *not*
                                                                                    # used in diagnostic modules.

    # Step 1: Declare free parameters intrinsic to these initial data
    amp = par.Cparameters("REAL",thismodule,"amp")
    lam = par.Cparameters("REAL",thismodule,"lam")

    # Step 2: Set the initial data
    global AD_ID,ED_ID,psi_ID
    AD_ID = ixp.zerorank1()

    ED_ID = ixp.zerorank1()
    EU_ID = ixp.zerorank1()
    # Set the coordinate transformations:
    radial = sp.sqrt(x*x + y*y + z*z)
    polar = sp.atan2(sp.sqrt(x*x + y*y),z)
    EU_phi = 8*amp*radial*sp.sin(polar)*lam*lam*sp.exp(lam*radial*radial)
    EU_ID[0] = (y * EU_phi)/sp.sqrt(x*x + y*y)
    EU_ID[1] = (y * EU_phi)/sp.sqrt(x*x + y*y)
    # The z component (2)is zero. 
    for i in range(DIM):
        for j in range(DIM):
            ED_ID[i] += gammaDD[i][j] * EU_ID[j]

    psi_ID = sp.sympify(0)
    Gamma_ID = sp.sympify(0)