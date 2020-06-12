# finite_difference.py:
#  As documented in the NRPy+ tutorial notebook:
#    Tutorial-Finite_Difference_Derivatives.ipynb ,
#  This module generates C kernels for numerically
#   solving PDEs with finite differences.
#
# Depends primarily on: outputC.py and grid.py.

# Author: Zachariah B. Etienne
#         zachetie **at** gmail **dot* com

from outputC import parse_outCparams_string,superfast_uniq,outputC # NRPy+: Core C code output module
import NRPy_param_funcs as par   # NRPy+: parameter interface
import sympy as sp               # SymPy: The Python computer algebra package upon which NRPy+ depends
import grid as gri               # NRPy+: Functions having to do with numerical grids
import sys                       # Standard Python module for multiplatform OS-level functions
from finite_difference_helpers import extract_from_list_of_deriv_vars__base_gfs_and_deriv_ops_lists
from finite_difference_helpers import generate_list_of_deriv_vars_from_lhrh_sympyexpr_list
from finite_difference_helpers import read_gfs_from_memory, type__var, FDparams, varsuffix

# Step 1: Initialize free parameters for this module:
modulename = __name__
# Centered finite difference accuracy order
par.initialize_param(par.glb_param("int", modulename, "FD_CENTDERIVS_ORDER",  4))
par.initialize_param(par.glb_param("int", modulename, "FD_KO_ORDER__CENTDERIVS_PLUS", 2))

def FD_outputC(filename,sympyexpr_list, params="", upwindcontrolvec=""):
    outCparams = parse_outCparams_string(params)

    FDparams.SIMD_enable      = outCparams.SIMD_enable
    FDparams.PRECISION        = par.parval_from_str("PRECISION")
    FDparams.DIM              = par.parval_from_str("DIM")
    FDparams.MemAllocStyle    = par.parval_from_str("MemAllocStyle")
    FDparams.upwindcontrolvec = upwindcontrolvec

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
    # See Step 6 for corresponding end brace.
    if outCparams.includebraces == "True":
        Coutput = outCparams.preindent+"{\n"
        indent = "   "
    else:
        Coutput = ""
        indent = ""

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

    # Step 2d (Upwinded derivatives algorithm, part 1):
    # If an upwinding control vector is specified, determine
    #    which of the elements of the vector will be required.
    #    This ensures that those elements are read from memory.
    # For example, if a symmetry axis is specified,
    #     upwind derivatives with respect to only
    #     two of the three dimensions are used. Here
    #     we find all directions used for upwinding.
    if upwindcontrolvec != "":
        upwind_directions_unsorted_withdups = []
        for deriv_op in list_of_deriv_operators:
            if "dupD" in deriv_op:
                if deriv_op[len(deriv_op)-1].isdigit():
                    dirn = int(deriv_op[len(deriv_op)-1])
                    upwind_directions_unsorted_withdups.append(dirn)
                else:
                    print("Error: Derivative operator "+deriv_op+" does not contain a direction")
                    sys.exit(1)
        upwind_directions = []
        if len(upwind_directions_unsorted_withdups)>0:
            upwind_directions = superfast_uniq(upwind_directions_unsorted_withdups)
            upwind_directions = sorted(upwind_directions,key=sp.default_sort_key)

    # Step 3:
    # Evaluate the finite difference stencil for each
    #     derivative operator,
    # TODO: being careful not to needlessly recompute.
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

    # Step 5: Output C code. C code consists of three parts
    #         a) Read gridfunctions from memory at needed pts.
    #         b) Perform arithmetic needed for input expressions
    #            provided in sympyexpr_list[].rhs and associated
    #            finite differences.
    #         c) Write output to gridfunctions specified in
    #            sympyexpr_list[].lhs.
    def indent_Ccode(Ccode):
        Ccodesplit = Ccode.splitlines()
        outstring = ""
        for i in range(len(Ccodesplit)):
            outstring += outCparams.preindent+indent+Ccodesplit[i]+'\n'
        return outstring

    # Step 5a: Read gridfunctions from memory at needed pts.
    # *** No need to do anything here; already set in
    #     string "read_from_memory_Ccode". ***

    # FIXME: Update these code comments:
    # Step 5b: Perform arithmetic needed for finite differences
    #          associated with input expressions provided in
    #          sympyexpr_list[].rhs.
    #          Note that exprs and lhsvarnames contain
    #          i)  finite difference expressions (constructed
    #              in steps above) and associated variable names,
    #              and
    #          ii) Input expressions sympyexpr_list[], which
    #              in general depend on finite difference
    #              variables.
    exprs       = []
    lhsvarnames = []
    # Step 5b.i: Output finite difference expressions to
    #            Coutput string
    for i in range(len(list_of_deriv_vars)):
        exprs.append(sp.sympify(0)) # Append a new element to the list of derivative expressions.
        lhsvarnames.append(type__var(list_of_deriv_vars[i], FDparams))
        var = list_of_base_gridfunction_names_in_derivs[i]
        for j in range(len(fdcoeffs[i])):
            varname = str(var)+varsuffix(fdstencl[i][j], FDparams)
            exprs[i] += fdcoeffs[i][j]*sp.sympify(varname)

        # Multiply each expression by the appropriate power
        #   of 1/dx[i]
        invdx = []
        for d in range(par.parval_from_str("DIM")):
            invdx.append(sp.sympify("invdx"+str(d)))
        # First-order or Kreiss-Oliger derivatives:
        if (len(list_of_deriv_operators[i]) == 5 and "dKOD" in list_of_deriv_operators[i]) or \
           (len(list_of_deriv_operators[i]) == 3 and "dD" in list_of_deriv_operators[i]) or \
           (len(list_of_deriv_operators[i]) == 5 and ("dupD" in list_of_deriv_operators[i] or "ddnD" in list_of_deriv_operators[i])):
            dirn = int(list_of_deriv_operators[i][len(list_of_deriv_operators[i])-1])
            exprs[i] *= invdx[dirn]
        # Second-order derivs:
        elif len(list_of_deriv_operators[i]) == 5 and "dDD" in list_of_deriv_operators[i]:
            dirn1 = int(list_of_deriv_operators[i][len(list_of_deriv_operators[i]) - 2])
            dirn2 = int(list_of_deriv_operators[i][len(list_of_deriv_operators[i]) - 1])
            exprs[i] *= invdx[dirn1]*invdx[dirn2]
        else:
            print("Error: was unable to parse derivative operator: ",list_of_deriv_operators[i])
            sys.exit(1)
    # Step 5b.ii: If upwind control vector is specified,
    #             add upwind control vectors to the
    #             derivative expression list, so its
    #             needed elements are read from memory.
    if upwindcontrolvec != "":
        for i in range(len(upwind_directions)):
            exprs.append(upwindcontrolvec[upwind_directions[i]])
            lhsvarnames.append(type__var("UpwindControlVectorU"+str(upwind_directions[i]),FDparams))

    # Step 5b.iii: Output useful code comment regarding
    #              which step we are on. *At most* this
    #              is a 3-step process:
    #           1. Read from memory & compute FD stencils,
    #           2. Perform upwinding, and
    #           3. Evaluate remaining expressions+write
    #              results to main memory.
    NRPy_FD_StepNumber = 1
    NRPy_FD__Number_of_Steps = 1
    if len(read_from_memory_Ccode) > 0:
        NRPy_FD__Number_of_Steps += 1
    if upwindcontrolvec != "" and len(upwind_directions) > 0:
        NRPy_FD__Number_of_Steps += 1

    if len(read_from_memory_Ccode) > 0:
        Coutput += indent_Ccode("/* \n * NRPy+ Finite Difference Code Generation, Step "
                                + str(NRPy_FD_StepNumber) + " of " + str(NRPy_FD__Number_of_Steps)+
                                ": Read from main memory and compute finite difference stencils:\n */\n")
        NRPy_FD_StepNumber = NRPy_FD_StepNumber + 1
        # Prefix chosen CSE variables with "FD", for the finite difference coefficients:
        Coutput += indent_Ccode(outputC(exprs,lhsvarnames,"returnstring",params=params + ",CSE_varprefix=FDPart1,includebraces=False,CSE_preprocess=True,SIMD_find_more_subs=True",
                                        prestring=read_from_memory_Ccode))

    # Step 5b.iv: Implement control-vector upwinding algorithm.
    if upwindcontrolvec != "":
        if len(upwind_directions) > 0:
            Coutput += indent_Ccode("/* \n * NRPy+ Finite Difference Code Generation, Step "
                                    + str(NRPy_FD_StepNumber) + " of " + str(NRPy_FD__Number_of_Steps) +
                                    ": Implement upwinding algorithm:\n */\n")
            NRPy_FD_StepNumber = NRPy_FD_StepNumber + 1
            if FDparams.SIMD_enable == "True":
                Coutput += """
const double tmp_upwind_Integer_1 = 1.000000000000000000000000000000000;
const REAL_SIMD_ARRAY upwind_Integer_1 = ConstSIMD(tmp_upwind_Integer_1);
const double tmp_upwind_Integer_0 = 0.000000000000000000000000000000000;
const REAL_SIMD_ARRAY upwind_Integer_0 = ConstSIMD(tmp_upwind_Integer_0);
"""
            for dirn in upwind_directions:
                Coutput += indent_Ccode(type__var("UpWind" + str(dirn),FDparams) +
                                        " = UPWIND_ALG(UpwindControlVectorU" + str(dirn) + ");\n")
        upwindU = [sp.sympify(0) for i in range(par.parval_from_str("DIM"))]
        for dirn in upwind_directions:
            upwindU[dirn] = sp.sympify("UpWind"+str(dirn))
        upwind_expr_list, var_list = [], []
        for i in range(len(list_of_deriv_vars)):
            if len(list_of_deriv_operators[i]) == 5 and ("dupD" in list_of_deriv_operators[i]):
                var_dupD = sp.sympify("UpwindAlgInput"+str(list_of_deriv_vars[i]))
                var_ddnD = sp.sympify("UpwindAlgInput"+str(list_of_deriv_vars[i]).replace("_dupD","_ddnD"))
                upwind_dirn = int(list_of_deriv_operators[i][len(list_of_deriv_operators[i])-1])
                upwind_expr = upwindU[upwind_dirn]*(var_dupD - var_ddnD) + var_ddnD
                upwind_expr_list.append(upwind_expr)
                var_list.append(type__var(str(list_of_deriv_vars[i]),FDparams,AddPrefix_for_UpDownWindVars=False))
        # For convenience, we require type__var() above to
        # prefix up/downwinded variables with "UpwindAlgInput".
        # Here we do not wish to have this prefix.
        Coutput += indent_Ccode(outputC(upwind_expr_list,var_list,
                                                "returnstring",params=params + ",CSE_varprefix=FDPart2,includebraces=False"))

    # Step 5b.v: Add input RHS & LHS expressions from
    #             sympyexpr_list[]
    Coutput += indent_Ccode("/* \n * NRPy+ Finite Difference Code Generation, Step "
                            + str(NRPy_FD_StepNumber) + " of " + str(NRPy_FD__Number_of_Steps) +
                            ": Evaluate SymPy expressions and write to main memory:\n */\n")
    exprs       = []
    lhsvarnames = []
    for i in range(len(sympyexpr_list)):
        exprs.append(sympyexpr_list[i].rhs)
        if FDparams.SIMD_enable == "True":
            lhsvarnames.append("const REAL_SIMD_ARRAY __RHS_exp_"+str(i))
        else:
            lhsvarnames.append(sympyexpr_list[i].lhs)

    # Step 5c: Write output to gridfunctions specified in
    #          sympyexpr_list[].lhs.
    write_to_mem_string = ""
    if FDparams.SIMD_enable == "True":
        for i in range(len(sympyexpr_list)):
            write_to_mem_string += "WriteSIMD(&"+sympyexpr_list[i].lhs+", __RHS_exp_"+str(i)+");\n"
    Coutput += indent_Ccode(outputC(exprs,lhsvarnames,"returnstring", params = params+",CSE_varprefix=FDPart3,includebraces=False,preindent=0", prestring="",poststring=write_to_mem_string))

    # Step 6: Add consistent indentation to the output end brace.
    #         See Step 0.b for corresponding start brace.
    if outCparams.includebraces == "True":
        Coutput += outCparams.preindent+"}\n"

    # Step 7: Output the C code in desired format: stdout, string, or file.
    if filename == "stdout":
        print(Coutput)
    elif filename == "returnstring":
        return Coutput+'\n'
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
#    print(gri.glb_gridfcs_list[1].name,list_of_points_read_from_memory[1])


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
    UPDOWNWIND = 0
    if "dupD" in derivstring:
        UPDOWNWIND =  1
    elif "ddnD" in derivstring:
        UPDOWNWIND = -1

    # Step 1: Set up matrix based on the stencil size (FDORDER+1).
    #         See documentation above for details on how this
    #         matrix is set up.
    M = sp.zeros(STENCILSIZE,STENCILSIZE)
    for i in range(STENCILSIZE):
        for j in range(STENCILSIZE):
            if i == 0:
                M[(i,j)] = 1 # Setting n^0 = 1 for all n, including n=0, because this matches the pattern
            else:
                dist_from_xeq0_col = j - sp.Rational((STENCILSIZE - 1),2) + UPDOWNWIND
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
                gridpt_posn = i - int((STENCILSIZE-1)/2) + UPDOWNWIND
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
    return fdcoeffs,fdstencl
