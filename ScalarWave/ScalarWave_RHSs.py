# Generating C code for the right-hand-side
#  of the scalar wave equation, in
#  ***Cartesian*** coordinates, in
#  arbitrary spatial dimension
#  (up to four spatial dimensions
#   supported at the moment)
#
# Author: Zachariah B. Etienne
#         zachetie **at** gmail **dot* com
#
# License: BSD 2-Clause

# COMPLETE DOCUMENTATION (JUPYTER NOTEBOOKS):
# START PAGE (start here!):  ../NRPy+_Tutorial.ipynb
# THIS MODULE: ../Tutorial-Scalarwave.ipynb

# Depends on: finite_difference.py, which itself depends on grid.py.
#             Everything depends on outputC.py.
#
# \partial_t^2 u  =  c^2 nabla^2 u
#
# where u is a function of space and time
# and nabla is the Laplacian differential operator.
#
# We rewrite as set of two first-order-in-time PDEs,
# defining \partial_t u = v.
#
# Then the above becomes
#
# \partial_t u = v
# \partial_t v = c^2 nabla^2 u
#
# where u and v are functions of space and time,
# and nabla is the Laplaciaion differential operator.

# Step P1: Import needed NRPy+ core modules:
from outputC import * # Needed for lhrh() named tuple
import indexedexp as ixp
import grid as gri
import finite_difference as fin

# Step P2: Define the C parameter wavespeed. The `wavespeed`
#          variable is a proper SymPy variable, so it can be
#          used in below expressions. In the C code, it acts
#          just like a usual parameter, whose value is
#          specified in the parameter file.
# The name of this module ("ScalarWave") is given by __name__:
thismodule = __name__
global wavespeed
wavespeed = par.Cparameters("REAL", thismodule, "wavespeed")

# Define the main ScalarWave() function, which outputs C code
def ScalarWave_RHSs():
    # Step 1: Get the spatial dimension, defined in the
    #         NRPy+ "grid" module.
    DIM = par.parval_from_str("DIM")

    # Step 2: Register gridfunctions that are needed as input
    #         to the scalar wave RHS expressions.
    uu, vv = gri.register_gridfunctions("EVOL",["uu","vv"])

    # Step 3: Declare the rank-2 indexed expression \partial_{ij} u,
    #         which is symmetric about interchange of indices i and j
    #         Derivative variables like these must have an underscore
    #         in them, so the finite difference module can parse the
    #         variable name properly.
    uu_dDD = ixp.declarerank2("uu_dDD","sym01")

    # Step 4: Specify RHSs as global variables,
    #         to enable access outside this
    #         function (e.g., for C code output)
    global uu_rhs,vv_rhs

    # Step 5: Define right-hand sides for the evolution.
    uu_rhs = vv
    vv_rhs = 0
    for i in range(DIM):
        vv_rhs += wavespeed*wavespeed*uu_dDD[i][i]


    # # Step 5: Generate C code for scalarwave evolution equations,
    # #         print output to the screen (standard out, or stdout).
    # fin.FD_outputC("stdout",
    #                [lhrh(lhs=gri.gfaccess("out_gfs","UU_rhs"),rhs=uu_rhs),
    #                 lhrh(lhs=gri.gfaccess("out_gfs","VV_rhs"),rhs=vv_rhs)])