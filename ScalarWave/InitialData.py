# Generating C code for plane wave initial
#  data for the scalar wave equation in
#  ***Cartesian*** coordinates, in up to
#  *three* spatial dimensions
#
# Authors: Zachariah B. Etienne
#          zachetie **at** gmail **dot* com
#          Thiago Assumpcao
#          assumpcaothiago **at** gmail **dot* com
#
# License: BSD 2-Clause

# COMPLETE DOCUMENTATION (JUPYTER NOTEBOOKS):
# START PAGE (start here!):  ../NRPy+_Tutorial.ipynb
# THIS MODULE: ../Tutorial-Scalarwave.ipynb

# Step P1: Import needed NRPy+ core modules:
import NRPy_param_funcs as par                # NRPy+: Parameter interface
import sympy as sp                            # SymPy: The Python computer algebra package upon which NRPy+ depends
import reference_metric as rfm                # NRPy+: Reference metric support
import sys                                    # Standard Python module for multiplatform OS-level functions
from ScalarWave.CommonParams import wavespeed # NRPy+: Common parameters for all ScalarWave modules (defines wavespeed)

# The name of this module ("InitialData") is given by __name__:
thismodule = __name__

# Set up spherically-symmetric Gaussian initial data
def SphericalGaussian(CoordSystem="Cartesian",default_time=0,default_sigma=3):
    # Step 1: Set parameters for the wave
    DIM = par.parval_from_str("grid::DIM")

    # Step 2: Set up Cartesian coordinates in terms of the native CoordSystem we have chosen.
    #         E.g., if CoordSystem="Cartesian", then xx_to_Cart = [xx[0],xx[1],xx[2]]
    #         or if CoordSystem="Spherical", then xx_to_Cart = [xx[0]*sp.sin(xx[1])*sp.cos(xx[2]),
    #                                                       xx[0]*sp.sin(xx[1])*sp.sin(xx[2]),
    #                                                       xx[0]*sp.cos(xx[1])]
    par.set_parval_from_str("reference_metric::CoordSystem",CoordSystem)
    rfm.reference_metric() # Must call this function to specify rfm.xx_to_Cart
    xx_to_Cart = rfm.xx_to_Cart

    # Step 3: Declare free parameters intrinsic to these initial data
    time  = par.Cparameters("REAL", thismodule, "time",  default_time)
    sigma = par.Cparameters("REAL", thismodule, "sigma", default_sigma)

    # Step 4: Compute r
    r = sp.sympify(0)
    for i in range(DIM):
        r += xx_to_Cart[i]**2
    r = sp.sqrt(r)

    # Step 5: Set initial data for uu and vv, where vv_ID = \partial_t uu_ID.
    global uu_ID, vv_ID
    # uu_ID = (r - wavespeed*time)/r * sp.exp(- (r - wavespeed*time)**2 / (2*sigma**2) )
    uu_ID = (+(r - wavespeed * time) / r * sp.exp(- (r - wavespeed * time) ** 2 / (2 * sigma ** 2))
             +(r + wavespeed * time) / r * sp.exp(- (r + wavespeed * time) ** 2 / (2 * sigma ** 2))) + sp.sympify(1)
    vv_ID = sp.diff(uu_ID, time)

# Set up monochromatic plane-wave initial data
def PlaneWave(CoordSystem="Cartesian",default_time=0,default_k0=1,default_k1=1,default_k2=1):
    # Step 1: Set parameters defined in other modules
    DIM = par.parval_from_str("grid::DIM")

    # Step 2: Set up Cartesian coordinates in terms of the native CoordSystem we have chosen.
    #         E.g., if CoordSystem="Cartesian", then xx_to_Cart = [xx[0],xx[1],xx[2]]
    #         or if CoordSystem="Spherical", then xx_to_Cart = [xx[0]*sp.sin(xx[1])*sp.cos(xx[2]),
    #                                                       xx[0]*sp.sin(xx[1])*sp.sin(xx[2]),
    #                                                       xx[0]*sp.cos(xx[1])]
    par.set_parval_from_str("reference_metric::CoordSystem",CoordSystem)
    rfm.reference_metric()
    xx_to_Cart = rfm.xx_to_Cart

    # Step 3: Declare free parameters intrinsic to these initial data
    time = par.Cparameters("REAL", thismodule, "time", default_time)
    kk   = par.Cparameters("REAL", thismodule, ["kk0", "kk1", "kk2"], [default_k0,default_k1,default_k2])

    # Step 4: Normalize the k vector
    kk_norm_factor = sp.sqrt(kk[0] ** 2 + kk[1] ** 2 + kk[2] ** 2)

    # Step 5: Compute k_norm.x
    dot_product = sp.sympify(0)
    for i in range(DIM):
        dot_product += kk[i] * xx_to_Cart[i]
    dot_product /= kk_norm_factor

    # Step 6: Set initial data for uu and vv, where vv_ID = \partial_t uu_ID.
    global uu_ID, vv_ID
    uu_ID = sp.sin(dot_product - wavespeed * time) + 2
    vv_ID = sp.diff(uu_ID, time)

# Initial data driver routine
def InitialData(Type="PlaneWave",CoordSystem="Cartesian",
                default_time=0,
                default_k0=1,default_k1=1,default_k2=1,
                default_sigma=3):
    if Type=="PlaneWave":
        PlaneWave(CoordSystem=CoordSystem, default_time=default_time,
                  default_k0=default_k0,default_k1=default_k1,default_k2=default_k2)
    elif Type=="SphericalGaussian":
        SphericalGaussian(CoordSystem=CoordSystem, default_time=default_time,   default_sigma=default_sigma)
    else:
        print("Error: ScalarWave initial data Type="+str(Type)+" not supported.")
        sys.exit(1)
