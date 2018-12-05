# Generating C code for plane wave initial
#  data for the scalar wave equation in
#  ***Cartesian*** coordinates, in up to
#  *three* spatial dimensions
#
# Author: Zachariah B. Etienne
#         zachetie **at** gmail **dot* com
#
# License: BSD 2-Clause

# COMPLETE DOCUMENTATION (JUPYTER NOTEBOOKS):
# START PAGE (start here!):  ../NRPy+_Tutorial.ipynb
# THIS MODULE: ../Tutorial-Scalarwave.ipynb

# Step P1: Import needed NRPy+ core modules:
#from outputC import * # Needed for lhrh() named tuple
import grid as gri
import NRPy_param_funcs as par
import sympy as sp
import ScalarWave.ScalarWave_RHSs as swrhs

thismodule = __name__

def InitialData_PlaneWave():
    # Step 1: Set parameters defined in other modules
    wavespeed = swrhs.wavespeed
    DIM = par.parval_from_str("grid::DIM")
    xx = gri.xx

    # Step 2: Declare free parameters intrinsic to these initial data
    time = par.Cparameters("REAL", thismodule, "time")
    kk = par.Cparameters("REAL", thismodule, ["kk0", "kk1", "kk2"])

    #  Step 3: Normalize the k vector
    kk_norm = sp.sqrt(kk[0] ** 2 + kk[1] ** 2 + kk[2] ** 2)

    # Step 4: Compute k.x
    dot_product = sp.sympify(0)
    for i in range(DIM):
        dot_product += xx[i] * kk[i]
    dot_product /= kk_norm

    # Step 5: Set initial data for uu and vv, where vv_ID = \partial_t uu_ID.
    global uu_ID,vv_ID
    uu_ID = sp.sin(dot_product - wavespeed * time)+2
    vv_ID = sp.diff(uu_ID, time)