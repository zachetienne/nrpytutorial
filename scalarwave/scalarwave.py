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
import tensor as ten
import grid as gri
import finite_difference as fin
from outputC import *

import sympy as sp

# import finite_difference as fin

# Step 1: Initialize free parameters for this module:
# modulename here will be set to "scalarwave", based on the filename.
thismodule = __name__

# Step 3: Define the main scalarwave() function, which outputs C code
def scalarwave():
    # Step 3a: Register gridfunctions that are needed as input.
    DIM = par.parval_from_str("DIM")
    uu, vv = gri.register_gridfunctions("EVOL",["uu","vv"])
    # Step 3b: Declare the rank-2 indexed expression \partial_{ij} u,
    #          which is symmetric about interchange of indices i and j
    #          Derivative variables like these must have an underscore
    #          in them, so the finite difference module can parse the
    #          variable name properly.
    uu_dDD = ten.declarerank2("uu_dDD","sym12")
    hDD = ten.register_gridfunctions_for_single_rank2("AUX","hDD","sym12")
    hDD_dDD = ten.declarerank4("hDD_dDD", "sym12_sym34")
    hDD_dKOD = ten.declarerank3("hDD_dKOD", "sym12")
    wavespeed = par.Cparameters("REAL",thismodule,"wavespeed")
    uu_rhs = vv
    vv_rhs = 0
    for i in range(DIM):
        vv_rhs += wavespeed*wavespeed*uu_dDD[i][i]
        # for j in range(DIM):
        #     for k in range(DIM):
        #         vv_rhs += hDD_dKOD[i][j][k]
        #         for l in range(DIM):
        #             vv_rhs += hDD_dDD[i][j][k][l]

    fin.FD_outputC_ToFile(str(thismodule)+"/NRPy_codegen/scalarwave_RHS.h",
                          [lhrh(lhs=gri.gfaccess("out_gfs","UU_rhs"),rhs=uu_rhs),
                           lhrh(lhs=gri.gfaccess("out_gfs","VV_rhs"),rhs=vv_rhs)])
