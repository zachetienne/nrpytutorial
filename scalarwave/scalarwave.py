# scalarwave.py:
#  Outputs C code for solving the simple,
#  Cartesian scalar wave equation in
#  arbitrary dimension:
#
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

import NRPy_param_funcs as par
import indexedexp as ixp
import grid as gri
import finite_difference as fin
from outputC import * # Needed for lhrh() named tuple

# The name of this module ("scalarwave") is given by __name__:
thismodule = __name__

# Define the main scalarwave() function, which outputs C code
def scalarwave():
    # Step 0: Get the spatial dimension, as it was set in 
    #         the grid module.
    DIM = par.parval_from_str("DIM")

    # Step 1: Register gridfunctions that are needed as input.
    uu, vv = gri.register_gridfunctions("EVOL",["uu","vv"])

    # Step 2: Declare the rank-2 indexed expression \partial_{ij} u,
    #         which is symmetric about interchange of indices i and j
    #         Derivative variables like these must have an underscore
    #         in them, so the finite difference module can parse the
    #         variable name properly.
    uu_dDD = ixp.declarerank2("uu_dDD","sym12")

    # Step 3: Define the C parameter wavespeed. The `wavespeed`
    #         variable is a proper SymPy variable, so it can be
    #         used in below expressions. In the C code, it acts
    #         just like a usual parameter, whose value is 
    #         specified in the parameter file.
    wavespeed = par.Cparameters("REAL",thismodule,"wavespeed")

    # Step 4: Define right-hand sides for the evolution.
    uu_rhs = vv
    vv_rhs = 0
    for i in range(DIM):
        vv_rhs += wavespeed*wavespeed*uu_dDD[i][i]

    # Step 5: Generate C code for scalarwave evolution equations,
    #         print output to the screen (standard out, or stdout).
    fin.FD_outputC("stdout",
                   [lhrh(lhs=gri.gfaccess("out_gfs","UU_rhs"),rhs=uu_rhs),
                    lhrh(lhs=gri.gfaccess("out_gfs","VV_rhs"),rhs=vv_rhs)])
