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
# and nabla is the Laplacian differential operator.

# Step P1: Import needed NRPy+ core modules:
import NRPy_param_funcs as par                # NRPy+: Parameter interface
import grid as gri                            # NRPy+: Functionality for handling numerical grids
import indexedexp as ixp                      # NRPy+: Symbolic indexed expression (e.g., tensors, vectors, etc.) support
from ScalarWave.CommonParams import wavespeed # NRPy+: Common parameters for all ScalarWave modules (defines wavespeed)

# The name of this module ("ScalarWave") is given by __name__:
thismodule = __name__

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


    # # Step 6: Generate C code for scalarwave evolution equations,
    # #         print output to the screen (standard out, or stdout).
    # fin.FD_outputC("stdout",
    #                [lhrh(lhs=gri.gfaccess("out_gfs","UU_rhs"),rhs=uu_rhs),
    #                 lhrh(lhs=gri.gfaccess("out_gfs","VV_rhs"),rhs=vv_rhs)])
