# finite_difference.py:
#  This module provides supporting functions for
#    finite_difference.py, which is documented in
#    the NRPy+ tutorial notebook:
#   Tutorial-Finite_Difference_Derivatives.ipynb ,
#
# Depends primarily on: outputC.py and grid.py.

# Author: Zachariah B. Etienne
#         zachetie **at** gmail **dot* com
from outputC import superfast_uniq, outputC, outC_function_dict, add_to_Cfunction_dict  # NRPy+: Core C code output module
from suffixes import getsuffix
import NRPy_param_funcs as par      # NRPy+: parameter interface
import sympy as sp                  # SymPy: The Python computer algebra package upon which NRPy+ depends
import grid as gri                  # NRPy+: Functions having to do with numerical grids
import sys                          # Standard Python module for multiplatform OS-level functions
from collections import namedtuple  # Standard Python: Enable namedtuple data type
from fstr import f

FDparams = namedtuple('FDparams', 'PRECISION FD_CD_order enable_FD_functions enable_SIMD DIM MemAllocStyle upwindcontrolvec fullindent outCparams')

#########################################
# STEP 1: EXTRACT DERIVATIVES TO COMPUTE
#         FROM LIST OF SYMPY EXPRESSIONS

def generate_list_of_deriv_vars_from_lhrh_sympyexpr_list(sympyexpr_list,FDparams):
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
    >>> from finite_difference_helpers import generate_list_of_deriv_vars_from_lhrh_sympyexpr_list,FDparams
    >>> aDD     = ixp.register_gridfunctions_for_single_rank2("EVOL","aDD","sym01")
    >>> aDD_dDD = ixp.declarerank4("aDD_dDD","sym01_sym23")
    >>> aDD_dupD = ixp.declarerank3("aDD_dupD","sym01")
    >>> betaU   = ixp.register_gridfunctions_for_single_rank1("EVOL","betaU")
    >>> a0,a1,b,c = par.Cparameters("REAL",__name__,["a0","a1","b","c"],1)
    >>> FDparams.upwindcontrolvec=betaU
    >>> exprlist = [lhrh(lhs=a0,rhs=b*aDD[1][0] + b*aDD_dDD[2][1][2][1] + c*aDD_dDD[0][1][1][0]), \
                    lhrh(lhs=a1,rhs=aDD_dDD[1][0][0][1] + c*aDD_dupD[0][2][1]*betaU[1])]
    >>> generate_list_of_deriv_vars_from_lhrh_sympyexpr_list(exprlist,FDparams)
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
                    raise Exception()
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
    if FDparams.upwindcontrolvec != "":
        for var in list_of_deriv_vars:
            if "_dupD" in str(var):
                list_of_deriv_vars.append(sp.sympify(str(var).replace("_dupD", "_ddnD")))

    # Finally, sort the list_of_deriv_vars. This ensures
    #     consistency in the C code output, and might even be
    #     tuned to reduce cache misses.
    #     Thanks to Aaron Meurer for this nice one-liner!
    return sorted(list_of_deriv_vars, key=sp.default_sort_key)
########################################

########################################
# STEP 2: EXTRACT DERIVATIVE INFORMATION
#         FROM LIST OF SYMPY EXPRESSIONS

def extract_from_list_of_deriv_vars__base_gfs_and_deriv_ops_lists(list_of_deriv_vars):
    """ Extract from list_of_deriv_vars a list of base gridfunctions
        and a list of derivative operators.
    :param list_of_deriv_vars:
    :return list_of_base_gridfunctions, list_of_deriv_operators:
    >>> from finite_difference_helpers import extract_from_list_of_deriv_vars__base_gfs_and_deriv_ops_lists
    >>> extract_from_list_of_deriv_vars__base_gfs_and_deriv_ops_lists(["aDD_dD012","aDD_dKOD012","vetU_dKOD21","hDD_dDD0112"])
    (['aDD01', 'aDD01', 'vetU2', 'hDD01'], ['dD2', 'dKOD2', 'dKOD1', 'dDD12'])
    """
    list_of_base_gridfunction_names_in_derivs = []
    list_of_deriv_operators = []
    # Step 2a:
    # For each var in "list_of_deriv_vars", determine the
    #     base gridfunction name and derivative operator.
    for var in list_of_deriv_vars:
        # Step 2a.1: Check that the number of integers appearing
        #            in the suffix of a variable name matches the
        #            number of U's + D's in the variable name:
        varstr = str(var)
        # Count number of "D" and "U" in the varstr.
        num_UDs = varstr.count('D') + varstr.count('U')
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
        #            list_of_base_gridfunction_names_in_derivs[]
        list_of_base_gridfunction_names_in_derivs.append(varstr[0:underscore_position] +
                                             varstr[len(varstr) - deriv_order - rank:len(varstr) - deriv_order])
        list_of_deriv_operators.append(varstr[underscore_position + 1:len(varstr) - deriv_order - rank] +
                               varstr[len(varstr) - deriv_order:len(varstr)])
    return list_of_base_gridfunction_names_in_derivs, list_of_deriv_operators
########################################

########################################
# STEP 4: GENERATE C CODE FOR READING
#         NEEDED GRIDPOINTS FROM MEMORY

from operator import itemgetter

def type__var(in_var,FDparams, AddPrefix_for_UpDownWindVars=True,is_const=True):
    """ Outputs [type] [variable name]; e.g.,
    "const double variable"

    :param in_var: Variable name
    :param FDparams: Parameters used in the finite-difference codegen
    :param AddPrefix_for_UpDownWindVars: Boolean -- add a prefix to up/downwind variables?
    :return: Return [type] [variable name]
    >>> from finite_difference_helpers import type__var, FDparams
    >>> FDparams.enable_SIMD = "True"
    >>> type__var("aDD00",FDparams)
    \'const REAL_SIMD_ARRAY aDD00\'

    >>> from finite_difference_helpers import type__var, FDparams
    >>> FDparams.enable_SIMD = "False"
    >>> FDparams.PRECISION = "double"
    >>> type__var("variable",FDparams)
    \'const double variable\'
    """
    varname = str(in_var)
    # Disable prefixing upwinded and downwinded variables
    # if the upwind control vector algorithm is disabled.
    if FDparams.upwindcontrolvec == "":
        AddPrefix_for_UpDownWindVars = False
    if AddPrefix_for_UpDownWindVars:
        if "_dupD" in varname:  # Variables suffixed with "_dupD" are set
            #                    to be the "pure" upwinded derivative,
            #                    before the upwinding algorithm has been
            #                    applied. However, when they are used
            #                    in the RHS expressions, it is assumed
            #                    that the up. algorithm has been applied.
            #                    To ensure consistency we rename all
            #                    _dupD suffixed variables as
            #                    _dupDPUREUPWIND, and use them as input
            #                    into the upwinding algorithm. The output
            #                    will be the original _dupD variable.
            varname = "UpwindAlgInput"+varname
        if "_ddnD" in varname:  # For consistency with _dupD
            varname = "UpwindAlgInput"+varname
    if FDparams.enable_SIMD == "True":
        return "const REAL_SIMD_ARRAY " + varname
    vartype = FDparams.PRECISION
    if is_const:
        return "const " + vartype + " " + varname
    else:
        return vartype + " " + varname

def read_from_memory_Ccode_onept(gfname,idx, FDparams, idxs=set()):
    """

    :param gfname: gridfunction name; a string
    :param idx: Grid index relative to (i0,i1,i2,i3); e.g., "0,1,2,3" -> (i0,i1+1,i2+2,i3+3); later indices ignored for DIM<4
    :param FDparams: Parameters used in the finite-difference codegen
    :param idxs: list of indexes used by the code
    :return: C code string for reading in this gridfunction at point idx from memory
    >>> import indexedexp as ixp
    >>> from finite_difference_helpers import FDparams, read_from_memory_Ccode_onept
    >>> FDparams.DIM = 3
    >>> FDparams.upwindcontrolvec = ""
    >>> FDparams.enable_SIMD = "True"
    >>> vetU = ixp.register_gridfunctions_for_single_rank1("EVOL","vetU",FDparams.DIM)
    >>> read_from_memory_Ccode_onept("vetU0","0,1,-2,300",FDparams)
    \'const REAL_SIMD_ARRAY vetU0_i0_i1p1_i2m2 = ReadSIMD(&in_gfs[IDX4S(VETU0GF, i0,i1+1,i2-2)]);\\n\'
    """
    idxsplit = idx.split(',')
    idx4 = [int(idxsplit[0]),int(idxsplit[1]),int(idxsplit[2]),int(idxsplit[3])]
    gf_array_name = "in_gfs" # Default array name.
    if gri.ET_driver == "CarpetX":
        ijklstring = ijk_carpetx(idx4)
        assert FDparams.DIM == 3
    else:
        ijklstring = ijkl_string(idx4, FDparams)
    gfaccess_str = gri.gfaccess(gf_array_name,gfname,ijklstring,"USE")
    idxs.add(",".join([str(ii) for ii in idx4]))
    if FDparams.enable_SIMD == "True":
        retstring = type__var(gfname, FDparams) + varsuffix(gfname, idx4, FDparams) + " = ReadSIMD(&" + gfaccess_str + ");"
    elif gri.find_gftype(gfname) == "SCALAR_TMP":
        retstring = ""
    else:
        assert varsuffix(gfname, idx4, FDparams) is not None
        assert gfaccess_str is not None
        retstring = type__var(gfname, FDparams) + varsuffix(gfname, idx4, FDparams) + " = " + gfaccess_str + ";"
    return retstring+"\n"

def ijk_carpetx(idx4):
    result = ""
    for i in range(len(idx4)):
        ff = idx4[i]
        if ff == 0:
            pass
        elif ff == 1:
            result += f("+p.DI[{i}]")
        elif ff == -1:
            result += f("-p.DI[{i}]")
        elif ff < 0:
            result += f("-{-ff}*p.DI[{i}]")
        else:
            result += f("+{ff}*p.DI[{i}]")
    return result

def ijkl_string(idx4, FDparams):
    """Generate string for reading gridfunction from specific location in memory
    if DIM==4:
        input: [i,j,k,l]
        output: "i0+i,i1+j,i2+k,i3+l"
    if DIM==3:
        input: [i,j,k,l]
        output: "i0+i,i1+j,i2+k"
    etc.
    :param idx4: An array of 4 integers, indicating a finite difference grid index relative to where the fd is being computed
    :param FDparams: Parameters used in the finite-difference codegen
    :return: DIM==3 input [i,j,k,l] -> output "i0+i,i1+j,i2+k"
    >>> from finite_difference_helpers import ijkl_string, FDparams
    >>> FDparams.DIM = 4
    >>> ijkl_string([-2,1,0,-1], FDparams)
    \'i0-2,i1+1,i2,i3-1\'

    >>> from finite_difference_helpers import ijkl_string, FDparams
    >>> FDparams.DIM = 3
    >>> ijkl_string([-2,-1,-1,-300], FDparams)
    \'i0-2,i1-1,i2-1\'
    """
    if False and gri.ET_driver == "CarpetX":
        assert FDparams.DIM == 3
        return ijk_carpetx(idx4)
    retstring = ""
    for i in range(FDparams.DIM):
        if i > 0:
            # Add a comma
            retstring += ","
        retstring += "i" + str(i) + "+" + str(idx4[i])
    return retstring.replace("+-", "-").replace("+0", "")

def varsuffix(name, idx4, FDparams):
    """Generate string for suffixing single point read in from memory
    Example: If a gridfunction is named hDD00, and we want to read from memory data at i0+1,i1,i2-1,
    we store the value of this gridfunction as hDD00_i0p1_i1_i2m1; this function provides the suffix.
    if DIM==3:
    input: [0,2,1,-100]
    output: "_i0_i1p2_i2p1"

    :param idx4: An array of 4 integers, indicating a finite difference grid index relative to where the fd is being computed
    :param FDparams: Parameters used in the finite-difference codegen
    :return: returns suffix to uniquely name a point of data for a gridfunction
    >>> from finite_difference_helpers import varsuffix, FDparams
    >>> FDparams.DIM=3
    >>> varsuffix("", [-2,0,-1,-300], FDparams)
    \'_i0m2_i1_i2m1\'
    """
    base_suffix = getsuffix(name)
    if idx4 == [0, 0, 0, 0]:
        return base_suffix
    return base_suffix + "_" + ijkl_string(idx4, FDparams).replace(",", "_").replace("+", "p").replace("-", "m")

def read_gfs_from_memory(list_of_base_gridfunction_names_in_derivs, fdstencl, sympyexpr_list, FDparams, idxs=None):
    # with open(list_of_base_gridfunction_names_in_derivs[0]+".txt","w") as file:
    #     file.write(str(list_of_base_gridfunction_names_in_derivs))
    #     file.write(str(fdstencl))
    #     file.write(str(sympyexpr_list))
    #     file.write(str(FDparams))
    """

    :param list_of_base_gridfunction_names_in_derivs:
    :param fdstencl:
    :param sympyexpr_list:
    :param FDparams:
    :return:
    >>> from outputC import lhrh
    >>> import indexedexp as ixp
    >>> import NRPy_param_funcs as par
    >>> from finite_difference_helpers import generate_list_of_deriv_vars_from_lhrh_sympyexpr_list,FDparams
    >>> from finite_difference_helpers import extract_from_list_of_deriv_vars__base_gfs_and_deriv_ops_lists
    >>> from finite_difference_helpers import read_gfs_from_memory
    >>> from finite_difference import compute_fdcoeffs_fdstencl
    >>> import grid as gri
    >>> gri.glb_gridfcs_list = []
    >>> hDD      = ixp.register_gridfunctions_for_single_rank2("EVOL","hDD","sym01")
    >>> hDD_dD   = ixp.declarerank3("hDD_dD","sym01")
    >>> hDD_dupD = ixp.declarerank3("hDD_dupD","sym01")
    >>> vU       = ixp.register_gridfunctions_for_single_rank1("EVOL","vU")
    >>> a0,a1,b,c = par.Cparameters("REAL",__name__,["a0","a1","b","c"],1)
    >>> par.set_parval_from_str("finite_difference::FD_CENTDERIVS_ORDER",2)
    >>> FDparams.DIM=3
    >>> FDparams.enable_SIMD="False"
    >>> FDparams.PRECISION="double"
    >>> FDparams.MemAllocStyle="012"
    >>> FDparams.upwindcontrolvec=vU
    >>> exprlist = [lhrh(lhs=a0,rhs=b*hDD[1][0] + c*hDD_dD[0][1][1]), \
                    lhrh(lhs=a1,rhs=c*hDD_dupD[0][2][1]*vU[1])]
    >>> list_of_deriv_vars = generate_list_of_deriv_vars_from_lhrh_sympyexpr_list(exprlist,FDparams)
    >>> list_of_base_gridfunction_names_in_derivs, list_of_deriv_operators = extract_from_list_of_deriv_vars__base_gfs_and_deriv_ops_lists(list_of_deriv_vars)
    >>> fdcoeffs = [[] for i in range(len(list_of_deriv_operators))]
    >>> fdstencl = [[[] for i in range(4)] for j in range(len(list_of_deriv_operators))]
    >>> for i in range(len(list_of_deriv_operators)): fdcoeffs[i], fdstencl[i] = compute_fdcoeffs_fdstencl(list_of_deriv_operators[i])
    >>> print(read_gfs_from_memory(list_of_base_gridfunction_names_in_derivs, fdstencl, exprlist, FDparams))
    const double hDD01_i0_i1m1_i2 = in_gfs[IDX4S(HDD01GF, i0,i1-1,i2)];
    const double hDD01 = in_gfs[IDX4S(HDD01GF, i0,i1,i2)];
    const double hDD01_i0_i1p1_i2 = in_gfs[IDX4S(HDD01GF, i0,i1+1,i2)];
    const double hDD02_i0_i1m2_i2 = in_gfs[IDX4S(HDD02GF, i0,i1-2,i2)];
    const double hDD02_i0_i1m1_i2 = in_gfs[IDX4S(HDD02GF, i0,i1-1,i2)];
    const double hDD02 = in_gfs[IDX4S(HDD02GF, i0,i1,i2)];
    const double hDD02_i0_i1p1_i2 = in_gfs[IDX4S(HDD02GF, i0,i1+1,i2)];
    const double hDD02_i0_i1p2_i2 = in_gfs[IDX4S(HDD02GF, i0,i1+2,i2)];
    const double vU1 = in_gfs[IDX4S(VU1GF, i0,i1,i2)];
    <BLANKLINE>
    """

    # Step 4a: Compile list of points to read from memory
    #          for each gridfunction i, based on list
    #          provided in fdstencil[i][].
    list_of_points_read_from_memory_with_duplicates = [[] for i in range(len(gri.glb_gridfcs_list))]
    for j in range(len(list_of_base_gridfunction_names_in_derivs)):
        derivgfname = list_of_base_gridfunction_names_in_derivs[j]
        # Next find the corresponding gridfunction index:
        for i in range(len(gri.glb_gridfcs_list)):
            gfname = gri.glb_gridfcs_list[i].name
            # If the gridfunction for the derivative matches, then
            #    add to the list of points read from memory:
            if derivgfname == gfname:
                for k in range(len(fdstencl[j])):
                    list_of_points_read_from_memory_with_duplicates[i].append(str(fdstencl[j][k][0]) + "," +
                                                                              str(fdstencl[j][k][1]) + "," +
                                                                              str(fdstencl[j][k][2]) + "," +
                                                                              str(fdstencl[j][k][3]))

    def get_symbols(arg):
        if hasattr(arg, "free_symbols"):
            return list(arg.free_symbols)
        elif type(arg) == str:
            return [arg]
        else:
            assert False

    # Step 4b: "Zeroth derivative" case:
    #     If gridfunction appears in expression not
    #     as derivative (i.e., by itself), it must
    #     be read from memory as well.
    for expr in range(len(sympyexpr_list)):
        for var in get_symbols(sympyexpr_list[expr].rhs) + get_symbols(sympyexpr_list[expr].lhs):
            vartype = gri.variable_type(var)
            if vartype == "gridfunction":
                for i in range(len(gri.glb_gridfcs_list)):
                    gfname = gri.glb_gridfcs_list[i].name
                    if gfname == str(var):
                        list_of_points_read_from_memory_with_duplicates[i].append("0,0,0,0")


    # Step 4c: Remove duplicates when reading from memory;
    #     do not needlessly read the same variable
    #     from memory twice.
    list_of_points_read_from_memory = [[] for i in range(len(gri.glb_gridfcs_list))]
    for i in range(len(gri.glb_gridfcs_list)):
        list_of_points_read_from_memory[i] = superfast_uniq(list_of_points_read_from_memory_with_duplicates[i])

    # Step 4d: Minimize cache misses:
    #      Sort the list of points read from
    #      main memory by how they are stored
    #      in memory.

    # Step 4d.i: Define a function that maps a gridpoint
    #     index (i,j,k,l) to a unique memory "address",
    #     which will correspond to the correct ordering
    #     of actual memory addresses.
    #
    #     Input: a list of 4 indices, e.g., (i,j,k,l)
    #            corresponding to a gridpoint's *spatial*
    #            index in memory (thus we support up to
    #            4D in space). If spatial dimension is
    #            less than 4D, then just set latter
    #            index/indices to zero. E.g., for 2D
    #            spatial indexing, set (i,j,0,0).
    #     Output: a single number, which when sorted
    #            will yield a unique "address" in memory
    #            such that consecutive addresses are
    #            consecutive in memory.
    def unique_idx(idx4,FDparams):
        # os and sz are set *just for the purposes of ensuring indices are ordered in memory*
        #    Do not modify the values of os and sz.
        os = 50  # offset
        sz = 100 # assumed size in each direction
        if FDparams.MemAllocStyle == "210":
            return str(int(idx4[0])+os + sz*( (int(idx4[1])+os) + sz*( (int(idx4[2])+os) + sz*( int(idx4[3])+os ) ) ))
        if FDparams.MemAllocStyle == "012":
            return str(int(idx4[3])+os + sz*( (int(idx4[2])+os) + sz*( (int(idx4[1])+os) + sz*( int(idx4[0])+os ) ) ))
        print("Error: MemAllocStyle = "+FDparams.MemAllocStyle+" unsupported.")
        sys.exit(1)

    # Step 4d.ii: For each gridfunction and
    #      point read from memory, call unique_idx,
    #      then sort according to memory "address"
    # Input: list_of_points_read_from_memory[gridfunction][point],
    #        gri.glb_gridfcs_list[gridfunction]
    # Output: 1) A list of points to be read from
    #            memory, sorted according to memory
    #            "address":
    #            sorted_list_of_points_read_from_memory[gridfunction][point]
    #        2) A list containing the gridfunction
    #           read at each point, with the number
    #           of elements corresponding exactly
    #           to the total number of points read
    #           from memory for all gridfunctions:
    #           read_from_memory_gf[]
    read_from_memory_gf    = []
    sorted_list_of_points_read_from_memory = [[] for i in range(len(gri.glb_gridfcs_list))]
    for gfidx in range(len(gri.glb_gridfcs_list)):
        # Continue only if reading at least one point of gfidx from memory.
        #     The sorting algorithm at the end of this code block is not
        #     well-defined (will throw an error) if no points of gfidx are
        #     read from memory.
        if len(list_of_points_read_from_memory[gfidx]) > 0:
            read_from_memory_index = []
            for idx in list_of_points_read_from_memory[gfidx]:
                read_from_memory_gf.append(gri.glb_gridfcs_list[gfidx])
                idxsplit = idx.split(',')
                idx4 = [int(idxsplit[0]),int(idxsplit[1]),int(idxsplit[2]),int(idxsplit[3])]
                read_from_memory_index.append(unique_idx(idx4, FDparams))
                # https://stackoverflow.com/questions/13668393/python-sorting-two-lists
                _unused_list, sorted_list_of_points_read_from_memory[gfidx] = \
                    [list(x) for x in zip(*sorted(zip(read_from_memory_index, list_of_points_read_from_memory[gfidx]),
                                                  key=itemgetter(0)))]
    # Step 4e: Create the full C code string
    #      for reading from memory:

    read_from_memory_Ccode = ""
    count = 0
    if idxs is None:
        idxs = set()
    for gfidx in range(len(sorted_list_of_points_read_from_memory)):
        for pt in range(len(sorted_list_of_points_read_from_memory[gfidx])):
            read_from_memory_Ccode += read_from_memory_Ccode_onept(read_from_memory_gf[count].name,
                                                                   sorted_list_of_points_read_from_memory[gfidx][pt],
                                                                   FDparams,idxs)
            count += 1
    return read_from_memory_Ccode
#################################

#################################
# STEP 5: C CODE OUTPUT ROUTINES
def construct_FD_exprs_as_SymPy_exprs(list_of_deriv_vars,
                                      list_of_base_gridfunction_names_in_derivs, list_of_deriv_operators,
                                      fdcoeffs, fdstencl):
    FDexprs = []
    FDlhsvarnames = []
    # Step 5.a.ii.A: Output finite difference expressions to
    #                Coutput string
    for i in range(len(list_of_deriv_vars)):
        FDexprs.append(sp.sympify(0))  # Append a new element to the list of derivative expressions.
        FDlhsvarnames.append(type__var(list_of_deriv_vars[i], FDparams))
        var = list_of_base_gridfunction_names_in_derivs[i]
        for j in range(len(fdcoeffs[i])):
            varname = str(var) + varsuffix(str(var),fdstencl[i][j], FDparams)
            FDexprs[i] += fdcoeffs[i][j] * sp.sympify(varname)

        # Multiply each expression by the appropriate power
        #   of 1/dx[i]
        invdx = []
        for d in range(FDparams.DIM):
            invdx.append(sp.sympify("invdx" + str(d)))
        # First-order or Kreiss-Oliger derivatives:
        if (len(list_of_deriv_operators[i]) == 5 and "dKOD" in list_of_deriv_operators[i]) or \
                (len(list_of_deriv_operators[i]) == 3 and "dD" in list_of_deriv_operators[i]) or \
                (len(list_of_deriv_operators[i]) == 5 and (
                        "dupD" in list_of_deriv_operators[i] or "ddnD" in list_of_deriv_operators[i])):
            dirn = int(list_of_deriv_operators[i][len(list_of_deriv_operators[i]) - 1])
            FDexprs[i] *= invdx[dirn]
        # Second-order derivs:
        elif len(list_of_deriv_operators[i]) == 5 and "dDD" in list_of_deriv_operators[i]:
            dirn1 = int(list_of_deriv_operators[i][len(list_of_deriv_operators[i]) - 2])
            dirn2 = int(list_of_deriv_operators[i][len(list_of_deriv_operators[i]) - 1])
            FDexprs[i] *= invdx[dirn1] * invdx[dirn2]
        else:
            print("Error: was unable to parse derivative operator: ", list_of_deriv_operators[i])
            sys.exit(1)
    return FDexprs, FDlhsvarnames

def find_which_op_idx(op, list_of_deriv_operators):
    for j in range(len(list_of_deriv_operators)):
        if op == list_of_deriv_operators[j]:
            return j
    print("Error: could not find operator "+str(op)+" in ",list_of_deriv_operators)
    sys.exit(1)


def add_FD_func_to_outC_function_dict(list_of_deriv_vars,
                                      list_of_base_gridfunction_names_in_derivs, list_of_deriv_operators,
                                      fdcoeffs, fdstencl):
    # Step 5.a.ii.A: First construct a list of all the unique finite difference functions
    list_of_uniq_deriv_operators = superfast_uniq(list_of_deriv_operators)
    c_type = "REAL"
    if par.parval_from_str("grid::GridFuncMemAccess") == "ETK":
        c_type = "CCTK_REAL"
    func_prefix = "order_"+str(FDparams.FD_CD_order)+"_"
    if FDparams.enable_SIMD == "True":
        c_type = "REAL_SIMD_ARRAY"
        func_prefix = "SIMD_"+func_prefix

    # Stores the needed calls to the functions we're adding to outC_function_dict:
    FDfunccall_list = []
    for op in list_of_uniq_deriv_operators:
        which_op_idx = find_which_op_idx(op, list_of_deriv_operators)

        rhs_expr = sp.sympify(0)
        for j in range(len(fdcoeffs[which_op_idx])):
            var = sp.sympify("f" + varsuffix("f",fdstencl[which_op_idx][j], FDparams))
            rhs_expr += fdcoeffs[which_op_idx][j] * var

        # Multiply each expression by the appropriate power
        #   of 1/dx[i]
        invdx = []
        used_invdx = [False, False, False, False]
        for d in range(FDparams.DIM):
            invdx.append(sp.sympify("invdx" + str(d)))
        # First-order or Kreiss-Oliger derivatives:
        if ( (len(op) == 5 and "dKOD" in op) or
             (len(op) == 3 and   "dD" in op) or
             (len(op) == 5 and ("dupD" in op or "ddnD" in op)) ):
            dirn = int(op[len(op) - 1])
            rhs_expr *= invdx[dirn]
            used_invdx[dirn] = True
        # Second-order derivs:
        elif len(op) == 5 and "dDD" in op:
            dirn1 = int(op[len(op) - 2])
            dirn2 = int(op[len(op) - 1])
            used_invdx[dirn1] = used_invdx[dirn2] = True
            rhs_expr *= invdx[dirn1]*invdx[dirn2]
        else:
            print("Error: was unable to parse derivative operator: ", op)
            sys.exit(1)

        outfunc_params = ""
        for d in range(FDparams.DIM):
            if used_invdx[d]:
                outfunc_params += "const " + c_type + " invdx" + str(d) + ","

        for j in range(len(fdcoeffs[which_op_idx])):
            var = sp.sympify("f" + varsuffix("",fdstencl[which_op_idx][j], FDparams))
            outfunc_params += "const " + c_type + " " + str(var)
            if j != len(fdcoeffs[which_op_idx])-1:
                outfunc_params += ","

        for i in range(len(list_of_deriv_operators)):
            # print("comparing ",list_of_deriv_operators[i],op)
            if list_of_deriv_operators[i] == op:
                funccall = type__var(list_of_deriv_vars[i], FDparams) + " = " + func_prefix + "f_" + str(op) + "("
                for d in range(FDparams.DIM):
                    if used_invdx[d]:
                        funccall += "invdx" + str(d) + ","
                gfname = list_of_base_gridfunction_names_in_derivs[i]
                for j in range(len(fdcoeffs[which_op_idx])):
                    funccall += gfname + varsuffix(gfname, fdstencl[which_op_idx][j], FDparams)
                    if j != len(fdcoeffs[which_op_idx])-1:
                        funccall += ","
                funccall += ");"
                FDfunccall_list.append(funccall)

        # If the function already exists in the outC_function_dict, then do not add it; move to the next op.
        if func_prefix + "f_" + str(op) not in outC_function_dict:
            p = "preindent=1,enable_SIMD="+FDparams.enable_SIMD+",outCverbose=False,CSE_preprocess=True,includebraces=False"
            outFDstr = outputC(rhs_expr, "retval", "returnstring", params=p)
            outFDstr = outFDstr.replace("retval = ", "return ")
            add_to_Cfunction_dict(desc=" * (__FD_OPERATOR_FUNC__) Finite difference operator for "+str(op).replace("dDD", "second derivative: ").
                                  replace("dD", "first derivative: ").replace("dKOD", "Kreiss-Oliger derivative: ").
                                  replace("dupD", "upwinded derivative: ").replace("ddnD", "downwinded derivative: ") + " direction. In Cartesian coordinates, directions 0,1,2 correspond to x,y,z directions, respectively.",
                                  c_type="static " + c_type + " _NOINLINE _UNUSED",
                                  name=func_prefix+"f_" + str(op), enableCparameters=False,
                                  params=outfunc_params, preloop="", body=outFDstr)
    return FDfunccall_list

def construct_Ccode(sympyexpr_list, list_of_deriv_vars,
                    list_of_base_gridfunction_names_in_derivs,list_of_deriv_operators,
                    fdcoeffs, fdstencl, read_from_memory_Ccode, FDparams, Coutput):
    """
    C code is constructed in *up to* 3 parts:
     5.a) Read gridfunctions from memory at needed pts
          for finite differencing; compute finite-differencing
          stencils.
     5.b) Implement upwinding algorithm (if relevant)
     5.c) Evaluate SymPy expressions and write to main
          memory
    """
    # Failed Doctest. However, mathematically equivalent with Sympy 1.3
    # :param sympyexpr_list:
    # :param list_of_deriv_vars:
    # :param list_of_base_gridfunction_names_in_derivs:
    # :param list_of_deriv_operators:
    # :param fdcoeffs:
    # :param fdstencl:
    # :param read_from_memory_Ccode:
    # :param FDparams:
    # :param Coutput: The start of the Coutput string; this function's output will be pasted to a copy of Coutput
    # :return: Returns a C code string
    # >>> from outputC import lhrh
    # >>> import indexedexp as ixp
    # >>> import NRPy_param_funcs as par
    # >>> from finite_difference_helpers import generate_list_of_deriv_vars_from_lhrh_sympyexpr_list,FDparams
    # >>> from finite_difference_helpers import extract_from_list_of_deriv_vars__base_gfs_and_deriv_ops_lists
    # >>> from finite_difference_helpers import read_gfs_from_memory, construct_Ccode
    # >>> from finite_difference import compute_fdcoeffs_fdstencl
    # >>> import grid as gri
    # >>> gri.glb_gridfcs_list = []
    # >>> hDD      = ixp.register_gridfunctions_for_single_rank2("EVOL","hDD","sym01")
    # >>> hDD_dD   = ixp.declarerank3("hDD_dD","sym01")
    # >>> hDD_dupD = ixp.declarerank3("hDD_dupD","sym01")
    # >>> vU       = ixp.register_gridfunctions_for_single_rank1("EVOL","vU")
    # >>> a0,a1,b,c = par.Cparameters("REAL",__name__,["a0","a1","b","c"],1)
    # >>> par.set_parval_from_str("finite_difference::FD_CENTDERIVS_ORDER",2)
    # >>> FDparams.DIM=3
    # >>> FDparams.enable_SIMD="False"
    # >>> FDparams.enable_FD_functions=False
    # >>> FDparams.PRECISION="double"
    # >>> FDparams.MemAllocStyle="012"
    # >>> FDparams.upwindcontrolvec=vU
    # >>> FDparams.fullindent=""
    # >>> FDparams.outCparams="outCverbose=False"
    # >>> exprlist = [lhrh(lhs=a0,rhs=b*hDD[1][0] + c*hDD_dD[0][1][1]), \
    #                 lhrh(lhs=a1,rhs=c*hDD_dupD[0][2][1]*vU[1])]
    # >>> list_of_deriv_vars = generate_list_of_deriv_vars_from_lhrh_sympyexpr_list(exprlist,FDparams)
    # >>> list_of_base_gridfunction_names_in_derivs, list_of_deriv_operators = extract_from_list_of_deriv_vars__base_gfs_and_deriv_ops_lists(list_of_deriv_vars)
    # >>> fdcoeffs = [[] for i in range(len(list_of_deriv_operators))]
    # >>> fdstencl = [[[] for i in range(4)] for j in range(len(list_of_deriv_operators))]
    # >>> for i in range(len(list_of_deriv_operators)): fdcoeffs[i], fdstencl[i] = compute_fdcoeffs_fdstencl(list_of_deriv_operators[i])
    # >>> memread_Ccode = read_gfs_from_memory(list_of_base_gridfunction_names_in_derivs, fdstencl, exprlist, FDparams)
    # >>> print(construct_Ccode(exprlist, list_of_deriv_vars, \
    #           list_of_base_gridfunction_names_in_derivs, list_of_deriv_operators, \
    #           fdcoeffs, fdstencl, memread_Ccode, FDparams, ""))
    # /*
    #  * NRPy+ Finite Difference Code Generation, Step 1 of 3: Read from main memory and compute finite difference stencils:
    #  */
    # const double hDD01_i0_i1m1_i2 = in_gfs[IDX4S(HDD01GF, i0,i1-1,i2)];
    # const double hDD01 = in_gfs[IDX4S(HDD01GF, i0,i1,i2)];
    # const double hDD01_i0_i1p1_i2 = in_gfs[IDX4S(HDD01GF, i0,i1+1,i2)];
    # const double hDD02_i0_i1m2_i2 = in_gfs[IDX4S(HDD02GF, i0,i1-2,i2)];
    # const double hDD02_i0_i1m1_i2 = in_gfs[IDX4S(HDD02GF, i0,i1-1,i2)];
    # const double hDD02 = in_gfs[IDX4S(HDD02GF, i0,i1,i2)];
    # const double hDD02_i0_i1p1_i2 = in_gfs[IDX4S(HDD02GF, i0,i1+1,i2)];
    # const double hDD02_i0_i1p2_i2 = in_gfs[IDX4S(HDD02GF, i0,i1+2,i2)];
    # const double vU1 = in_gfs[IDX4S(VU1GF, i0,i1,i2)];
    # const double FDPart1_Rational_1_2 = 1.0/2.0;
    # const double FDPart1_Integer_2 = 2.0;
    # const double FDPart1_Rational_3_2 = 3.0/2.0;
    # const double hDD_dD011 = FDPart1_Rational_1_2*invdx1*(-hDD01_i0_i1m1_i2 + hDD01_i0_i1p1_i2);
    # const double UpwindAlgInputhDD_ddnD021 = invdx1*(-FDPart1_Integer_2*hDD02_i0_i1m1_i2 + FDPart1_Rational_1_2*hDD02_i0_i1m2_i2 + FDPart1_Rational_3_2*hDD02);
    # const double UpwindAlgInputhDD_dupD021 = invdx1*(FDPart1_Integer_2*hDD02_i0_i1p1_i2 - FDPart1_Rational_1_2*hDD02_i0_i1p2_i2 - FDPart1_Rational_3_2*hDD02);
    # const double UpwindControlVectorU1 = vU1;
    # /*
    #  * NRPy+ Finite Difference Code Generation, Step 2 of 3: Implement upwinding algorithm:
    #  */
    # const double UpWind1 = UPWIND_ALG(UpwindControlVectorU1);
    # const double hDD_dupD021 = UpWind1*(-UpwindAlgInputhDD_ddnD021 + UpwindAlgInputhDD_dupD021) + UpwindAlgInputhDD_ddnD021;
    # /*
    #  * NRPy+ Finite Difference Code Generation, Step 3 of 3: Evaluate SymPy expressions and write to main memory:
    #  */
    # a0 = b*hDD01 + c*hDD_dD011;
    # a1 = c*hDD_dupD021*vU1;
    # <BLANKLINE>

    def indent_Ccode(Ccode):
        Ccodesplit = Ccode.splitlines()
        outstring = ""
        for i in range(len(Ccodesplit)):
            if Ccodesplit[i] != "":
                if Ccodesplit[i].lstrip().startswith("#"):
                    # Remove all indentation from preprocessor statements (lines that start with "#")
                    outstring += Ccodesplit[i].lstrip() + '\n'
                else:
                    outstring += FDparams.fullindent + Ccodesplit[i] + '\n'
        return outstring.rstrip(" ")  # make sure to remove trailing whitespace!

    # Step 5.a.i: Read gridfunctions from memory at needed pts.
    # *** No need to do anything here; already set in
    #     string "read_from_memory_Ccode". ***

    # Step 5.a.ii: Perform arithmetic needed for finite differences
    #              associated with input expressions provided in
    #              sympyexpr_list[].rhs.
    #           Note that FDexprs and FDlhsvarnames contain
    #          A) Finite difference expressions (constructed
    #             in steps above) and associated variable names,
    #             and
    #          B) Input expressions sympyexpr_list[], which
    #             in general depend on finite difference
    #             variables.
    FDexprs = []
    FDlhsvarnames = []
    if not FDparams.enable_FD_functions:
        FDexprs, FDlhsvarnames = \
            construct_FD_exprs_as_SymPy_exprs(list_of_deriv_vars,
                                              list_of_base_gridfunction_names_in_derivs, list_of_deriv_operators,
                                              fdcoeffs, fdstencl)

    # Compute finite differences using function calls (instead of inlined calculations)?
    if FDparams.enable_FD_functions:
        # If so, add FD functions to outputC's outC_function_dict (C function dictionary),
        #   AND return the full set of needed calls to these functions (to funccall_list)
        funccall_list = \
            add_FD_func_to_outC_function_dict(list_of_deriv_vars,
                                              list_of_base_gridfunction_names_in_derivs, list_of_deriv_operators,
                                              fdcoeffs, fdstencl)

    # Step 5.b.i: (Upwinded derivatives algorithm, part 1):
    # If an upwinding control vector is specified, determine
    #    which of the elements of the vector will be required.
    #    This ensures that those elements are read from memory.
    # For example, if a symmetry axis is specified,
    #     upwind derivatives with respect to only
    #     two of the three dimensions are used. Here
    #     we find all directions used for upwinding.
    upwind_directions = []
    if FDparams.upwindcontrolvec != "":
        upwind_directions_unsorted_withdups = []
        for deriv_op in list_of_deriv_operators:
            if "dupD" in deriv_op:
                if deriv_op[len(deriv_op)-1].isdigit():
                    dirn = int(deriv_op[len(deriv_op)-1])
                    upwind_directions_unsorted_withdups.append(dirn)
                else:
                    print("Error: Derivative operator "+deriv_op+" does not contain a direction")
                    sys.exit(1)
        if len(upwind_directions_unsorted_withdups) > 0:
            upwind_directions = superfast_uniq(upwind_directions_unsorted_withdups)
            upwind_directions = sorted(upwind_directions,key=sp.default_sort_key)
        #   If upwind control vector is specified,
        #        add upwind control vectors to the
        #        derivative expression list, so its
        #        needed elements are read from memory.
        for dirn in upwind_directions:
            FDexprs.append(FDparams.upwindcontrolvec[dirn])
            FDlhsvarnames.append(type__var("UpwindControlVectorU" + str(dirn), FDparams))

    # Step 5.x: Output useful code comment regarding
    #           which step we are on. *At most* this
    #           is a 3-step process:
    #        1. Read from memory & compute FD stencils,
    #        2. Perform upwinding, and
    #        3. Evaluate remaining expressions+write
    #           results to main memory.
    NRPy_FD_StepNumber = 1
    NRPy_FD__Number_of_Steps = 1
    if len(read_from_memory_Ccode) > 0:
        NRPy_FD__Number_of_Steps += 1
    if FDparams.upwindcontrolvec != "" and len(upwind_directions) > 0:
        NRPy_FD__Number_of_Steps += 1

    if len(read_from_memory_Ccode) > 0:
        Coutput += indent_Ccode("/*\n * NRPy+ Finite Difference Code Generation, Step "
                                + str(NRPy_FD_StepNumber) + " of " + str(NRPy_FD__Number_of_Steps) +
                                ": Read from main memory and compute finite difference stencils:\n */\n")
        NRPy_FD_StepNumber = NRPy_FD_StepNumber + 1
        if FDparams.enable_FD_functions:
            # Compute finite differences using function calls (instead of inlined calculations)
            Coutput += indent_Ccode(read_from_memory_Ccode)
            for funccall in funccall_list:
                Coutput += indent_Ccode(funccall)
            if FDparams.upwindcontrolvec != "":
                # Compute finite differences using inlined calculations
                params = FDparams.outCparams
                # We choose the CSE temporary variable prefix "FDpart1" for the finite difference coefficients:
                params += ",CSE_varprefix=FDPart1,includebraces=False,CSE_preprocess=True,SIMD_find_more_subs=True"
                Coutput += indent_Ccode(outputC(FDexprs, FDlhsvarnames, "returnstring", params=params))

        else:
            # Compute finite differences using inlined calculations
            params = FDparams.outCparams.replace("preindent=1", "preindent=0")  # Remove an unnecessary indentation
            # We choose the CSE temporary variable prefix "FDpart1" for the finite difference coefficients:
            params += ",CSE_varprefix=FDPart1,includebraces=False,CSE_preprocess=True,SIMD_find_more_subs=True"
            Coutput += indent_Ccode(outputC(FDexprs, FDlhsvarnames, "returnstring",params=params,
                                            prestring=read_from_memory_Ccode))

    # Step 5.b.ii: Implement control-vector upwinding algorithm.
    if FDparams.upwindcontrolvec != "":
        if len(upwind_directions) > 0:
            Coutput += indent_Ccode("/*\n * NRPy+ Finite Difference Code Generation, Step "
                                    + str(NRPy_FD_StepNumber) + " of " + str(NRPy_FD__Number_of_Steps) +
                                    ": Implement upwinding algorithm:\n */\n")
            NRPy_FD_StepNumber = NRPy_FD_StepNumber + 1
            if FDparams.enable_SIMD == "True":
                for n in ["0", "1"]:
                    Coutput += indent_Ccode("const double tmp_upwind_Integer_"+n+" = "+n+".000000000000000000000000000000000;\n")
                    Coutput += indent_Ccode("const REAL_SIMD_ARRAY upwind_Integer_"+n+" = ConstSIMD(tmp_upwind_Integer_"+n+");\n")
            for dirn in upwind_directions:
                Coutput += indent_Ccode(type__var("UpWind" + str(dirn), FDparams) +
                                        " = UPWIND_ALG(UpwindControlVectorU" + str(dirn) + ");\n")
        upwindU = [sp.sympify(0) for i in range(FDparams.DIM)]
        for dirn in upwind_directions:
            upwindU[dirn] = sp.sympify("UpWind" + str(dirn))
        upwind_expr_list, var_list = [], []
        for i in range(len(list_of_deriv_vars)):
            if len(list_of_deriv_operators[i]) == 5 and ("dupD" in list_of_deriv_operators[i]):
                var_dupD = sp.sympify("UpwindAlgInput" + str(list_of_deriv_vars[i]))
                var_ddnD = sp.sympify("UpwindAlgInput" + str(list_of_deriv_vars[i]).replace("_dupD", "_ddnD"))
                upwind_dirn = int(list_of_deriv_operators[i][len(list_of_deriv_operators[i]) - 1])
                upwind_expr = upwindU[upwind_dirn] * (var_dupD - var_ddnD) + var_ddnD
                upwind_expr_list.append(upwind_expr)
                var_list.append(type__var(str(list_of_deriv_vars[i]), FDparams, AddPrefix_for_UpDownWindVars=False))
        # For convenience, we require type__var() above to
        # prefix up/downwinded variables with "UpwindAlgInput".
        # Here we do not wish to have this prefix.
        Coutput += indent_Ccode(outputC(upwind_expr_list, var_list,
                                        "returnstring", params=FDparams.outCparams + ",CSE_varprefix=FDPart2,includebraces=False"))

    # Step 5.c.i: Add input RHS & LHS expressions from
    #             sympyexpr_list[]
    Coutput += indent_Ccode("/*\n * NRPy+ Finite Difference Code Generation, Step "
                            + str(NRPy_FD_StepNumber) + " of " + str(NRPy_FD__Number_of_Steps) +
                            ": Evaluate SymPy expressions and write to main memory:\n */\n")
    exprs = []
    lhsvarnames = []
    for i in range(len(sympyexpr_list)):
        exprs.append(sympyexpr_list[i].rhs)
        if FDparams.enable_SIMD == "True":
            lhsvarnames.append("const REAL_SIMD_ARRAY __RHS_exp_" + str(i))
        else:
            lhsvarnames.append(sympyexpr_list[i].lhs)

    # Step 5.c.ii: Write output to gridfunctions specified in
    #              sympyexpr_list[].lhs.
    write_to_mem_string = ""
    if FDparams.enable_SIMD == "True":
        for i in range(len(sympyexpr_list)):
            write_to_mem_string += "WriteSIMD(&" + sympyexpr_list[i].lhs + ", __RHS_exp_" + str(i) + ");\n"

    # outputC requires as its second argument a list of strings.
    #   Sometimes when the lhs's are simple constants, but the inputs
    #   contain gridfunctions, it is necessary to convert the lhs's
    #   to strings:
    lhsvarnamestrings = []
    for lhs in lhsvarnames:
        lhsvarnamestrings.append(str(lhs))

    Coutput += indent_Ccode(outputC(exprs, lhsvarnamestrings, "returnstring",
                                    params=FDparams.outCparams + ",CSE_varprefix=FDPart3,includebraces=False,preindent=0",
                                    prestring="", poststring=write_to_mem_string))

    return Coutput
#################################

if __name__ == "__main__":
    import doctest
    sys.exit(doctest.testmod()[0])
