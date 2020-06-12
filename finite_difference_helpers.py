# finite_difference.py:
#  This module provides supporting functions for
#    finite_difference.py, which is documented in
#    the NRPy+ tutorial notebook:
#   Tutorial-Finite_Difference_Derivatives.ipynb ,
#
# Depends primarily on: outputC.py and grid.py.

# Author: Zachariah B. Etienne
#         zachetie **at** gmail **dot* com
from outputC import parse_outCparams_string,superfast_uniq,outputC # NRPy+: Core C code output module
import sympy as sp               # SymPy: The Python computer algebra package upon which NRPy+ depends
import grid as gri               # NRPy+: Functions having to do with numerical grids
import sys                       # Standard Python module for multiplatform OS-level functions

def generate_list_of_deriv_vars_from_lhrh_sympyexpr_list(sympyexpr_list,upwindcontrolvec):
    """
    Generate from list of SymPy expressions in the form
    [lhrh(lhs=var, rhs=expr),lhrh(...),...]
    all derivative expressions.
    :param sympyexpr_list <- list of SymPy expressions in the form [lhrh(lhs=var, rhs=expr),lhrh(...),...]:
    :return list of derivative variables; creating _ddnD in case upwinding is enabled with control vector:
    >>> from outputC import lhrh
    >>> import indexedexp as ixp
    >>> import grid as gri
    >>> import NRPy_param_funcs as par
    >>> from finite_difference_helpers import generate_list_of_deriv_vars_from_lhrh_sympyexpr_list
    >>> aDD     = ixp.register_gridfunctions_for_single_rank2("EVOL","aDD","sym01")
    >>> aDD_dDD = ixp.declarerank4("aDD_dDD","sym01_sym23")
    >>> aDD_dupD = ixp.declarerank3("aDD_dupD","sym01")
    >>> betaU   = ixp.register_gridfunctions_for_single_rank1("EVOL","betaU")
    >>> a0,a1,b,c = par.Cparameters("REAL",__name__,["a0","a1","b","c"],1)
    >>> upwindcontrolvec="betaU"
    >>> exprlist = [lhrh(lhs=a0,rhs=b*aDD[1][0] + b*aDD_dDD[2][1][2][1] + c*aDD_dDD[0][1][1][0]), \
                    lhrh(lhs=a1,rhs=aDD_dDD[1][0][0][1] + c*aDD_dupD[0][2][1]*betaU[1])]
    >>> generate_list_of_deriv_vars_from_lhrh_sympyexpr_list(exprlist,upwindcontrolvec)
    [aDD_dDD0101, aDD_dDD1212, aDD_ddnD021, aDD_dupD021]
    """
    # Step 1a:
    # Create a list of free symbols in the sympy expr list
    #     that are registered neither as gridfunctions nor
    #     as C parameters. These *must* be derivatives,
    #     so we call the list "list_of_deriv_vars"
    list_of_deriv_vars_with_duplicates = []
    for expr in sympyexpr_list:
        for var in expr.rhs.free_symbols:
            vartype = gri.variable_type(var)
            if vartype == "other":
                # vartype=="other" should ONLY refer to derivatives, so
                #    if "_dD" or variants do not appear in a variable classified
                #    neither as a gridfunction nor a Cparameter, then error out.
                if ("_dD"   in str(var)) or \
                   ("_dKOD" in str(var)) or \
                   ("_dupD" in str(var)) or \
                   ("_ddnD" in str(var)):
                    list_of_deriv_vars_with_duplicates.append(var)
                else:
                    print("Error: Unregistered variable \""+str(var)+"\" in SymPy expression for "+expr.lhs)
                    print("All variables in SymPy expressions passed to FD_outputC() must be registered")
                    print("in NRPy+ as either a gridfunction or Cparameter, by calling")
                    print(str(var)+" = register_gridfunctions...() (in ixp/grid) if \""+str(var)+"\" is a gridfunction, or")
                    print(str(var)+" = Cparameters() (in par) otherwise (e.g., if it is a free parameter set at C runtime).")
                    sys.exit(1)
    list_of_deriv_vars = superfast_uniq(list_of_deriv_vars_with_duplicates)

    # Upwinding with respect to a control vector: algorithm description.
    #   To enable, set the FD_outputC()'s fourth function argument to the
    #   desired control vector. In BSSN, the betaU vector controls the upwinding.
    #   See https://arxiv.org/pdf/gr-qc/0206072.pdf for motivation and
    #   https://arxiv.org/pdf/gr-qc/0109032.pdf for implementation details,
    #   at second order. Note that the BSSN shift vector behaves like a *negative*
    #   velocity. See http://www.damtp.cam.ac.uk/user/naweb/ii/advection/advection.php
    #   for a very basic example motivating this choice.

    # Step 1b: For each variable with suffix _dupD, append to
    #          the list_of_deriv_vars the corresponding _ddnD.
    #          Both are required for control-vector upwinding. See
    #          the above print() block for further documentation
    #          on upwinding--both motivation and implementation
    #          details.
    if upwindcontrolvec != "":
        for var in list_of_deriv_vars:
            if "_dupD" in str(var):
                list_of_deriv_vars.append(sp.sympify(str(var).replace("_dupD","_ddnD")))

    # Finally, sort the list_of_deriv_vars. This ensures
    #     consistency in the C code output, and might even be
    #     tuned to reduce cache misses.
    #     Thanks to Aaron Meurer for this nice one-liner!
    return sorted(list_of_deriv_vars,key=sp.default_sort_key)

def extract_from_list_of_deriv_vars__base_gfs_and_deriv_ops_lists(list_of_deriv_vars):
    """ Extract from list_of_deriv_vars a list of base gridfunctions
        and a list of derivative operators.
    :param list_of_deriv_vars:
    :return list_of_base_gridfunctions,, list_of_deriv_operators:
    >>> from finite_difference_helpers import extract_from_list_of_deriv_vars__base_gfs_and_deriv_ops_lists
    >>> extract_from_list_of_deriv_vars__base_gfs_and_deriv_ops_lists(["aDD_dD012","aDD_dKOD012","vetU_dKOD21","hDD_dDD0112"])
    (['aDD01', 'aDD01', 'vetU2', 'hDD01'], ['dD2', 'dKOD2', 'dKOD1', 'dDD12'])
    """
    deriv__base_gridfunction_name = []
    deriv__operator = []
    # Step 2a:
    # For each var in "list_of_deriv_vars", determine the
    #     base gridfunction name and derivative operator.
    for var in list_of_deriv_vars:
        # Step 2a.1: Check that the number of integers appearing
        #            in the suffix of a variable name matches the
        #            number of U's + D's in the variable name:
        varstr = str(var)
        num_UDs = 0
        for i in range(len(varstr)):
            if varstr[i] == 'D' or varstr[i] == 'U':
                num_UDs += 1
        num_digits = 0
        i = len(varstr) - 1
        while varstr[i].isdigit():
            num_digits += 1
            i -= 1
        if num_UDs != num_digits:
            print("Error: " + varstr + " has " + str(num_UDs) + " U's and D's, but ")
            print(str(num_digits) + " integers at the end. These must be equal.")
            print("Please rename your gridfunction.")
            sys.exit(1)
        # Step 2a.2: Based on the variable name, find the rank of
        #            the underlying gridfunction of which we're
        #            trying to take the derivative.
        rank = 0  # rank = "number of juxtaposed U's and D's before the underscore in a derivative expression"
        underscore_position = -1
        for i in range(len(varstr) - 1, -1, -1):
            if underscore_position > 0 and (varstr[i] == "U" or varstr[i] == "D"):
                rank += 1
            if varstr[i] == "_":
                underscore_position = i

        # Step 2a.3: Based on the variable name, find the order
        #            of the derivative we're trying to take.
        deriv_order = 0  # deriv_order = "number of D's after the underscore in a derivative expression"
        for i in range(underscore_position + 1, len(varstr)):
            if (varstr[i] == "D"):
                deriv_order += 1

        # Step 2a.4: Based on derivative order and rank,
        #            store the base gridfunction name in
        #            deriv__base_gridfunction_name[]
        deriv__base_gridfunction_name.append(varstr[0:underscore_position] +
                                             varstr[len(varstr) - deriv_order - rank:len(varstr) - deriv_order])
        deriv__operator.append(varstr[underscore_position + 1:len(varstr) - deriv_order - rank] +
                               varstr[len(varstr) - deriv_order:len(varstr)])
    return deriv__base_gridfunction_name, deriv__operator
