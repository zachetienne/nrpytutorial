# finite_difference.py:
#  As documented in the NRPy+ tutorial notebook:
#    Tutorial-Finite_Difference_Derivatives.ipynb ,
#  This module generates C kernels for numerically
#   solving PDEs with finite differences.
#
# Depends primarily on: outputC.py and grid.py.

# Author: Zachariah B. Etienne
#         zachetie **at** gmail **dot* com

from outputC import parse_outCparams_string, outC_function_dict  # NRPy+: Core C code output module
import NRPy_param_funcs as par   # NRPy+: parameter interface
import sympy as sp               # SymPy: The Python computer algebra package upon which NRPy+ depends
import grid as gri               # NRPy+: Functions having to do with numerical grids
import os, sys                   # Standard Python module for multiplatform OS-level functions
from finite_difference_helpers import extract_from_list_of_deriv_vars__base_gfs_and_deriv_ops_lists
from finite_difference_helpers import generate_list_of_deriv_vars_from_lhrh_sympyexpr_list
from finite_difference_helpers import read_gfs_from_memory, FDparams, construct_Ccode

# Step 1: Initialize free parameters for this module:
modulename = __name__
# Centered finite difference accuracy order
par.initialize_param(par.glb_param("int",  modulename, "FD_CENTDERIVS_ORDER",          4))
par.initialize_param(par.glb_param("bool", modulename, "FD_functions_enable",      False))
par.initialize_param(par.glb_param("int",  modulename, "FD_KO_ORDER__CENTDERIVS_PLUS", 2))

def FD_outputC(filename,sympyexpr_list, params="", upwindcontrolvec=""):
    outCparams = parse_outCparams_string(params)

    # Step 0.a:
    # In case sympyexpr_list is a single sympy expression,
    #     convert it to a list with just one element.
    #     This enables the rest of the routine to assume
    #     sympyexpr_list is indeed a list.
    if type(sympyexpr_list) is not list:
        sympyexpr_list = [sympyexpr_list]

    # Step 0.b:
    # finite_difference.py takes control over outCparams.includebraces here,
    #     which is necessary because outputC() is called twice:
    #     first for the reads from main memory and finite difference
    #     stencil expressions, and second for the SymPy expressions and
    #     writes to main memory.
    # If outCparams.includebraces==True, then it will close off the braces
    #     after the finite difference stencil expressions and start new ones
    #     for the SymPy expressions and writes to main memory, resulting
    #     in a non-functioning C code.
    # To get around this issue, we create braces around the entire
    #     string of C output from this function, only if
    #     outCparams.includebraces==True.
    # See Step 5 for open and close braces
    if outCparams.includebraces == "True":
        indent = "   "
    else:
        indent = ""

    # Step 0.c: FDparams named tuple stores parameters used in the finite-difference codegen
    FDparams.SIMD_enable         = outCparams.SIMD_enable
    FDparams.PRECISION           = par.parval_from_str("PRECISION")
    FDparams.FD_CD_order         = par.parval_from_str("FD_CENTDERIVS_ORDER")
    FDparams.FD_functions_enable = par.parval_from_str("FD_functions_enable")
    FDparams.DIM                 = par.parval_from_str("DIM")
    FDparams.MemAllocStyle       = par.parval_from_str("MemAllocStyle")
    FDparams.upwindcontrolvec    = upwindcontrolvec
    FDparams.fullindent          = indent + outCparams.preindent
    FDparams.outCparams          = params

    # Step 1: Generate from list of SymPy expressions in the form
    #     [lhrh(lhs=var, rhs=expr),lhrh(...),...]
    #     all derivative expressions, which we will process next.
    list_of_deriv_vars = generate_list_of_deriv_vars_from_lhrh_sympyexpr_list(sympyexpr_list, FDparams)

    # Step 2a: Extract from list_of_deriv_vars a list of base gridfunctions
    #         and a list of derivative operators. Usually takes list of SymPy
    #         symbols as input, but could just take strings, as this function
    #         does only string manipulations.
    # Example:
    # >>> extract_from_list_of_deriv_vars__base_gfs_and_deriv_ops_lists(["aDD_dD012","aDD_dKOD012","vetU_dKOD21","hDD_dDD0112"])
    # (['aDD01', 'aDD01', 'vetU2', 'hDD01'], ['dD2', 'dKOD2', 'dKOD1', 'dDD12'])
    list_of_base_gridfunction_names_in_derivs, list_of_deriv_operators = \
        extract_from_list_of_deriv_vars__base_gfs_and_deriv_ops_lists(list_of_deriv_vars)

    # Step 2b:
    # Next, check each base gridfunction to determine whether
    #     it is indeed registered as a gridfunction.
    #     If not, exit with error.
    for basegf in list_of_base_gridfunction_names_in_derivs:
        is_gf = False
        for gf in gri.glb_gridfcs_list:
            if basegf == str(gf.name):
                is_gf = True
        if not is_gf:
            print("Error: Attempting to take the derivative of "+basegf+", which is not a registered gridfunction.")
            print("       Make sure your gridfunction name does not have any underscores in it!")
            sys.exit(1)

    # Step 2c:
    # Check each derivative operator to make sure it is
    #     supported. If not, error out.
    for i in range(len(list_of_deriv_operators)):
        found_derivID = False
        for derivID in ["dD","dupD","ddnD","dKOD"]:
            if derivID in list_of_deriv_operators[i]:
                found_derivID = True
        if not found_derivID:
            print("Error: Valid derivative operator in "+list_of_deriv_operators[i]+" not found.")
            sys.exit(1)

    # Step 3:
    # Evaluate the finite difference stencil for each
    #     derivative operator, being careful not to
    #     needlessly recompute.
    # Note: Each finite difference stencil consists
    #     of two parts:
    #     1) The coefficient, and
    #     2) The index corresponding to the coefficient.
    #     The former is stored as a rational number, and
    #     the latter as a simple string, such that e.g.,
    #     in 3D, the empty string corresponds to (i,j,k),
    #     the string "ip1" corresponds to (i+1,j,k),
    #     the string "ip1kp1" corresponds to (i+1,j,k+1),
    #     etc.
    fdcoeffs = [[] for i in range(len(list_of_deriv_operators))]
    fdstencl = [[[] for i in range(4)] for j in range(len(list_of_deriv_operators))]
    for i in range(len(list_of_deriv_operators)):
        fdcoeffs[i], fdstencl[i] = compute_fdcoeffs_fdstencl(list_of_deriv_operators[i])

    # Step 4: Create C code to read gridfunctions from memory
    read_from_memory_Ccode = read_gfs_from_memory(list_of_base_gridfunction_names_in_derivs, fdstencl, sympyexpr_list,
                                                  FDparams)

    # Step 5: construct C code.
    Coutput = ""
    if outCparams.includebraces == "True":
        Coutput = outCparams.preindent + "{\n"
    Coutput = construct_Ccode(sympyexpr_list, list_of_deriv_vars,
                           list_of_base_gridfunction_names_in_derivs, list_of_deriv_operators,
                           fdcoeffs, fdstencl, read_from_memory_Ccode, FDparams, Coutput)
    if outCparams.includebraces == "True":
        Coutput += outCparams.preindent+"}"

    # Step 6: Output the C code in desired format: stdout, string, or file.
    if filename == "stdout":
        print(Coutput)
    elif filename == "returnstring":
        return Coutput
    else:
        # Output to the file specified by outCfilename
        with open(filename, outCparams.outCfileaccess) as file:
            file.write(Coutput)
        successstr = ""
        if outCparams.outCfileaccess == "a":
            successstr = "Appended "
        elif outCparams.outCfileaccess == "w":
            successstr = "Wrote "
        print(successstr + "to file \"" + filename + "\"")

def output_finite_difference_functions_h(path=os.path.join(".")):
    with open(os.path.join(path, "finite_difference_functions.h"), "w") as file:
        file.write("""
#ifndef __FD_FUNCTIONS_H__
#define __FD_FUNCTIONS_H__
#include "math.h"
#include "stdio.h"
#include "stdlib.h"
""")
        UNUSED   = "__attribute__((unused))"
        NOINLINE = "__attribute__((noinline))"
        if par.parval_from_str("grid::GridFuncMemAccess") == "ETK":
            UNUSED   = "CCTK_ATTRIBUTE_UNUSED"
            NOINLINE = "CCTK_ATTRIBUTE_NOINLINE"
        file.write("#define _UNUSED   " + UNUSED   + "\n")
        file.write("#define _NOINLINE " + NOINLINE + "\n")

        for key, item in outC_function_dict.items():
            if "__FD_OPERATOR_FUNC__" in item:
                file.write(item.replace("const REAL_SIMD_ARRAY _NegativeOne_ =",
                                        "const REAL_SIMD_ARRAY "+UNUSED+" _NegativeOne_ =")) # Many of the NegativeOne's get optimized away in the SIMD postprocessing step. No need for all the warnings

        file.write("#endif // #ifndef __FD_FUNCTIONS_H__\n")

#######################################################
#  FINITE-DIFFERENCE COEFFICIENT ALGORITHM

#  Define the to-be-inverted matrix, A.
#  We define A row-by-row, according to the prescription
#  derived in notes/notes.pdf, via the following pattern
#  that applies for arbitrary order.
#
#  As an example, consider a 5-point finite difference
#  stencil (4th-order accurate), where we wish to compute
#  some derivative at the center point.
#
#  Then A is given by:
#
#  -2^0  -1^0  1  1^0   2^0
#  -2^1  -1^1  0  1^1   2^1
#  -2^2  -1^2  0  1^2   2^2
#  -2^3  -1^3  0  1^3   2^3
#  -2^4  -1^4  0  1^4   2^4
#
#  Then right-multiplying A^{-1}
#  by (1 0 0 0 0)^T will yield 0th deriv. stencil
#  by (0 1 0 0 0)^T will yield 1st deriv. stencil
#  by (0 0 1 0 0)^T will yield 2nd deriv. stencil
#  etc.
#
#  Next suppose we want an upwinded, 4th-order accurate
#  stencil. For this case, A is given by:
#
#  -1^0  1  1^0   2^0   3^0
#  -1^1  0  1^1   2^1   3^1
#  -1^2  0  1^2   2^2   3^2
#  -1^3  0  1^3   2^3   3^3
#  -1^4  0  1^4   2^4   3^4
#
#  ... and similarly for the downwinded derivative.
#
#  Finally, let's consider a 3rd-order accurate
#  stencil. This would correspond to an in-place
#  upwind stencil with stencil radius of 2 gridpoints,
#  where other, centered derivatives are 4th-order
#  accurate. For this case, A is given by:
#
#  -1^0  1  1^0   2^0
#  -1^1  0  1^1   2^1
#  -1^2  0  1^2   2^2
#  -1^3  0  1^3   2^3
#  -1^4  0  1^4   2^4
#
#  ... and similarly for the downwinded derivative.
#
#  The general pattern is as follows:
#
#  1) The top row is all 1's,
#  2) If the second row has N elements (N must be odd),
#  .... then the radius of the stencil is rs = (N-1)/2
#  .... and the j'th row e_j = j-rs-1. For example,
#  .... for 4th order, we have rs = 2
#  .... j  | element
#  .... 1  | -2
#  .... 2  | -1
#  .... 3  |  0
#  .... 4  |  1
#  .... 5  |  2
#  3) The L'th row, L>2 will be the same as the second
#  .... row, but with each element e_j -> e_j^(L-1)
#  A1 is used later to validate the inverted
#  matrix.

def compute_fdcoeffs_fdstencl(derivstring,FDORDER=-1):
    # Step 0: Set finite differencing order, stencil size, and up/downwinding
    if FDORDER == -1:
        FDORDER = par.parval_from_str("FD_CENTDERIVS_ORDER")
        if "dKOD" in derivstring:
            FDORDER += par.parval_from_str("FD_KO_ORDER__CENTDERIVS_PLUS")

    STENCILSIZE = FDORDER+1
    UPDOWNWIND_stencil_shift = 0
    # dup/dnD = single-point-offset upwind/downwinding.
    if "dupD" in derivstring:
        UPDOWNWIND_stencil_shift =  1
    elif "ddnD" in derivstring:
        UPDOWNWIND_stencil_shift = -1
    # dfullup/dnD = full upwind/downwinding.
    elif "dfullupD" in derivstring:
        UPDOWNWIND_stencil_shift =  int(FDORDER/2)
    elif "dfulldnD" in derivstring:
        UPDOWNWIND_stencil_shift = -int(FDORDER/2)

    # Step 1: Set up matrix based on the stencil size (FDORDER+1).
    #         See documentation above for details on how this
    #         matrix is set up.
    M = sp.zeros(STENCILSIZE,STENCILSIZE)
    for i in range(STENCILSIZE):
        for j in range(STENCILSIZE):
            if i == 0:
                M[(i,j)] = 1 # Setting n^0 = 1 for all n, including n=0, because this matches the pattern
            else:
                dist_from_xeq0_col = j - sp.Rational((STENCILSIZE - 1),2) + UPDOWNWIND_stencil_shift
                if dist_from_xeq0_col==0:
                    M[(i,j)] = 0
                else:
                    M[(i, j)] = dist_from_xeq0_col**(i)
    Minv = sp.zeros(STENCILSIZE,STENCILSIZE)
    Minv = M**(-1)

    # Step 2:
    #     Based on the input derivative string,
    #     pick out the relevant row of the matrix
    #     inverse, as outlined in the detailed code
    #     comments prior to this function definition.
    derivtype = "FirstDeriv"
    matrixrow = 1
    if "DDD" in derivstring:
        print("Error: Only derivatives up to second order currently supported.")
        print("       Feel free to contribute to NRPy+ to extend its functionality!")
        sys.exit(1)
    elif "DD" in derivstring:

        if derivstring[len(derivstring)-1] == derivstring[len(derivstring)-2]:
            # Assuming i==j, we call \partial_i \partial_j gf an "unmixed" second derivative,
            #     or more simply, just "SecondDeriv":
            derivtype = "SecondDeriv"
            matrixrow = 2
        else:
            # Assuming i!=j, we call \partial_i \partial_j gf a MIXED second derivative,
            #     which is computed using a composite of first derivative operations.
            derivtype = "MixedSecondDeriv"
    elif "dKOD" in derivstring:
        derivtype = "KreissOligerDeriv"
        matrixrow = STENCILSIZE - 1
    else:
        # Up/downwinded and first derivs are all of "FirstDeriv" type
        pass

    # Step 3:
    #     Set finite difference coefficients
    #     and stencil points corresponding to
    #     each finite difference coefficient.
    fdcoeffs = []
    fdstencl = []
    if derivtype != "MixedSecondDeriv":
        for i in range(STENCILSIZE):
            idx4 = [0, 0, 0, 0]
            # First compute finite difference coefficient.
            fdcoeff = sp.factorial(matrixrow)*Minv[(i,matrixrow)]
            # Do not store fdcoeff or fdstencil if
            # finite difference coefficient is zero.
            if fdcoeff != 0:
                fdcoeffs.append(fdcoeff)
                if derivtype == "KreissOligerDeriv":
                    fdcoeffs[i] *= (-1)**(sp.Rational((STENCILSIZE+1),2))/2**matrixrow

                # Next store finite difference stencil point
                # corresponding to coefficient.
                gridpt_posn = i - int((STENCILSIZE-1)/2) + UPDOWNWIND_stencil_shift
                if gridpt_posn != 0:
                    dirn = int(derivstring[len(derivstring)-1])
                    idx4[dirn] = gridpt_posn
                fdstencl.append(idx4)
    else:
        # Mixed second derivative finite difference coeffs
        #     consist of products of first deriv coeffs,
        #     defined in first Minv matrix row.
        for i in range(STENCILSIZE):
            for j in range(STENCILSIZE):
                idx4 = [0, 0, 0, 0]

                # First compute finite difference coefficient.
                fdcoeff = (sp.factorial(matrixrow)*Minv[(i,matrixrow)]) * \
                          (sp.factorial(matrixrow)*Minv[(j,matrixrow)])

                # Do not store fdcoeff or fdstencil if
                # finite difference coefficient is zero.
                if fdcoeff != 0:
                    fdcoeffs.append(fdcoeff)

                    # Next store finite difference stencil point
                    # corresponding to coefficient.
                    gridpt_posn1 = i - int((STENCILSIZE - 1) / 2)
                    gridpt_posn2 = j - int((STENCILSIZE - 1) / 2)
                    dirn1 = int(derivstring[len(derivstring) - 1])
                    dirn2 = int(derivstring[len(derivstring) - 2])
                    idx4[dirn1] = gridpt_posn1
                    idx4[dirn2] = gridpt_posn2
                    fdstencl.append(idx4)
    return fdcoeffs, fdstencl
