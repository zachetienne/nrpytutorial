import NRPy_param_funcs as par
import indexedexp as ixp
import grid as gri
import finite_difference as fin
from outputC import *

def MaxwellCartesian_ID():
    DIM = par.parval_from_str("grid::DIM")

    x,y,z = gri.register_gridfunctions("AUX",["x","y","z"])
    gammaDD = ixp.register_gridfunctions_for_single_rank2("AUX","gammaDD", "sym01") # The AUX or EVOL designation is *not*
                                                                                    # used in diagnostic modules.

    # Step 1: Declare free parameters intrinsic to these initial data
    amp = par.Cparameters("REAL",__name__,"amp") # __name__ = "MaxwellCartesian_ID", this module's name
    lam = par.Cparameters("REAL",__name__,"lam") # __name__ = "MaxwellCartesian_ID", this module's name

    # Step 2: Set the initial data
    system = par.parval_from_str("System_to_use")
    if system is "System_I" or "System_II":
        global AidD,EidD,psi_ID
        AidD = ixp.zerorank1()

        EidD = ixp.zerorank1()
        EidU = ixp.zerorank1()
        # Set the coordinate transformations:
        radial = sp.sqrt(x*x + y*y + z*z)
        polar = sp.atan2(sp.sqrt(x*x + y*y),z)
        EU_phi = 8*amp*radial*sp.sin(polar)*lam*lam*sp.exp(-lam*radial*radial)
        EidU[0] = -(y * EU_phi)/sp.sqrt(x*x + y*y)
        EidU[1] = (x * EU_phi)/sp.sqrt(x*x + y*y)
        # The z component (2)is zero. 
        for i in range(DIM):
            for j in range(DIM):
                EidD[i] += gammaDD[i][j] * EidU[j]

        psi_ID = sp.sympify(0)
        if system is "System_II":
            global Gamma_ID
            Gamma_ID = sp.sympify(0)
    else:
        print("Invalid choice of system: System_to_use must be either System_I or System_II")
