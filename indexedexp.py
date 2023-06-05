# indexedexp.py: functions related to indexed expressions,
# including e.g., tensors and pseudotensors:

# Step 1: Load needed modules
import NRPy_param_funcs as par   # NRPy+: Parameter interface
import grid as gri               # NRPy+: Functions having to do with numerical grids
import functional as func        # NRPy+: Python toolkit for functional programming
import sympy as sp               # SymPy: The Python computer algebra package upon which NRPy+ depends
import sys                       # Standard Python module for multiplatform OS-level functions
import re                        # Standard Python module for regular expressions

# =====================
# Things relating to grid function groups.
# Essentially, the group name is the basename
# for a rankN tensor.
index_group = {}
rev_index_group = {}


def get_group_name(gf_name):
    if gf_name in rev_index_group:
        return rev_index_group[gf_name]

    # Assume that if the name is not in the
    # rev_index_group, then it is in a group
    # with the same name as the gf
    assert gf_name in gri.glb_gridfcs_map(), "Not a valid grid function: '" + gf_name + "'"
    return gf_name


def get_all_group_names():
    group_names = {}
    for gf_group in index_group:
        group_names[gf_group] = 1
    for gf_name in gri.glb_gridfcs_map():
        if gf_name not in rev_index_group:
            group_names[gf_name] = 1
    return group_names


def get_gfnames_for_group(gf_group):
    if gf_group in index_group:
        return index_group[gf_group]
    assert gf_group in gri.glb_gridfcs_map(), "Not a valid grid group: '" + gf_group + "'"
    return {gf_group: 1}


def find_gftype_for_group(gf_group, fail_on_missing=True):
    for gf_name in get_gfnames_for_group(gf_group):
        return gri.find_gftype(gf_name, fail_on_missing=fail_on_missing)
    return None


def find_gfmodule_for_group(gf_group, fail_on_missing=True):
    for gf_name in get_gfnames_for_group(gf_group):
        return gri.find_gfmodule(gf_name, fail_on_missing=fail_on_missing)
    return None


def find_centering_for_group(gf_group, fail_on_missing=True):
    for gf_name in get_gfnames_for_group(gf_group):
        return gri.find_centering(gf_name)
    return None


# =====================

thismodule = __name__
par.initialize_param(par.glb_param("char", thismodule, "symmetry_axes", ""))


def declare_indexedexp(rank, symbol=None, symmetry=None, dimension=None, namefun=None):
    """ Generate an indexed expression of specified rank and dimension

        >>> ixp = declare_indexedexp(rank=2, symbol='M', dimension=3, symmetry='sym01')
        >>> assert func.pipe(ixp, lambda x: func.repeat(func.flatten, x, 1), set, len) == 6

        >>> ixp = declare_indexedexp(rank=3, symbol='M', dimension=3, symmetry='sym01')
        >>> assert len(set(func.repeat(func.flatten, ixp, 2))) == 18
        >>> ixp = declare_indexedexp(rank=3, symbol='M', dimension=3, symmetry='sym02')
        >>> assert len(set(func.repeat(func.flatten, ixp, 2))) == 18
        >>> ixp = declare_indexedexp(rank=3, symbol='M', dimension=3, symmetry='sym12')
        >>> assert len(set(func.repeat(func.flatten, ixp, 2))) == 18

        >>> ixp = declare_indexedexp(rank=3, symbol='M', dimension=3, symmetry='sym012')
        >>> assert len(set(func.repeat(func.flatten, ixp, 2))) == 10

        >>> ixp = declare_indexedexp(rank=4, symbol='M', dimension=3, symmetry='sym01')
        >>> assert len(set(func.repeat(func.flatten, ixp, 3))) == 54
        >>> ixp = declare_indexedexp(rank=4, symbol='M', dimension=3, symmetry='sym02')
        >>> assert len(set(func.repeat(func.flatten, ixp, 3))) == 54
        >>> ixp = declare_indexedexp(rank=4, symbol='M', dimension=3, symmetry='sym03')
        >>> assert len(set(func.repeat(func.flatten, ixp, 3))) == 54
        >>> ixp = declare_indexedexp(rank=4, symbol='M', dimension=3, symmetry='sym12')
        >>> assert len(set(func.repeat(func.flatten, ixp, 3))) == 54
        >>> ixp = declare_indexedexp(rank=4, symbol='M', dimension=3, symmetry='sym13')
        >>> assert len(set(func.repeat(func.flatten, ixp, 3))) == 54
        >>> ixp = declare_indexedexp(rank=4, symbol='M', dimension=3, symmetry='sym23')
        >>> assert len(set(func.repeat(func.flatten, ixp, 3))) == 54

        >>> ixp = declare_indexedexp(rank=4, symbol='M', dimension=3, symmetry='sym012')
        >>> assert len(set(func.repeat(func.flatten, ixp, 3))) == 30
        >>> ixp = declare_indexedexp(rank=4, symbol='M', dimension=3, symmetry='sym013')
        >>> assert len(set(func.repeat(func.flatten, ixp, 3))) == 30
        >>> ixp = declare_indexedexp(rank=4, symbol='M', dimension=3, symmetry='sym01_sym23')
        >>> assert len(set(func.repeat(func.flatten, ixp, 3))) == 36
        >>> ixp = declare_indexedexp(rank=4, symbol='M', dimension=3, symmetry='sym02_sym13')
        >>> assert len(set(func.repeat(func.flatten, ixp, 3))) == 36
        >>> ixp = declare_indexedexp(rank=4, symbol='M', dimension=3, symmetry='sym023')
        >>> assert len(set(func.repeat(func.flatten, ixp, 3))) == 30
        >>> ixp = declare_indexedexp(rank=4, symbol='M', dimension=3, symmetry='sym03_sym12')
        >>> assert len(set(func.repeat(func.flatten, ixp, 3))) == 36
        >>> ixp = declare_indexedexp(rank=4, symbol='M', dimension=3, symmetry='sym123')
        >>> assert len(set(func.repeat(func.flatten, ixp, 3))) == 30

        >>> ixp = declare_indexedexp(rank=4, symbol='M', dimension=3, symmetry='sym0123')
        >>> assert len(set(func.repeat(func.flatten, ixp, 3))) == 15

        >>> ixp = declare_indexedexp(rank=2, symbol='M', dimension=3, symmetry='anti01')
        >>> assert len(set(map(abs, func.repeat(func.flatten, ixp, 1))).difference({0})) == 3
        >>> ixp = declare_indexedexp(rank=3, symbol='M', dimension=3, symmetry='anti012')
        >>> assert len(set(map(abs, func.repeat(func.flatten, ixp, 2))).difference({0})) == 1
        >>> ixp = declare_indexedexp(rank=4, symbol='M', dimension=3, symmetry='anti0123')
        >>> assert len(set(map(abs, func.repeat(func.flatten, ixp, 3))).difference({0})) == 0
    """
    if not dimension or dimension == -1:
        dimension = par.parval_from_str('DIM')
    if symbol is not None:
        if not isinstance(symbol, str) or not re.match(r'[\w_]', symbol):
            raise ValueError('symbol must be an alphabetic string')
    if dimension is not None:
        if not isinstance(dimension, int) or dimension <= 0:
            raise ValueError('dimension must be a positive integer')
    indexedexp = _init(rank * [dimension], symbol, namefun=namefun)
    if symmetry: return symmetrize(rank, indexedexp, symmetry, dimension)
    return apply_symmetry_condition_to_derivatives(indexedexp)


def _init(shape, symbol, index=None, namefun=None):
    if isinstance(shape, int):
        shape = [shape]
    if not index: index = []
    if namefun is None:
        iterable = [sp.Symbol(symbol + ''.join(str(n) for n in index + [i]))
                    if symbol else sp.sympify(0) for i in range(shape[0])]
    else:
        iterable = namefun(symbol, index, shape)
    if len(shape) > 1:
        for i in range(shape[0]):
            iterable[i] = _init(shape[1:], symbol, index + [i], namefun=namefun)
    return iterable


def symmetrize(rank, indexedexp, symmetry, dimension):
    if rank == 1:
        if symmetry == 'nosym': return indexedexp
        raise Exception('cannot symmetrize indexed expression of rank 1')
    if rank == 2:
        indexedexp = symmetrize_rank2(indexedexp, symmetry, dimension)
    elif rank == 3:
        indexedexp = symmetrize_rank3(indexedexp, symmetry, dimension)
    elif rank == 4:
        indexedexp = symmetrize_rank4(indexedexp, symmetry, dimension)
    else:
        raise Exception('unsupported rank for indexed expression')
    return apply_symmetry_condition_to_derivatives(indexedexp)


def symmetrize_rank2(indexedexp, symmetry, dimension):
    for sym in symmetry.split('_'):
        sign = 1 if sym[:3] == 'sym' else -1
        for i, j in func.product(range(dimension), repeat=2):
            if sym[-2:] == '01':
                if j < i:
                    indexedexp[i][j] = sign * indexedexp[j][i]
                elif i == j and sign < 0:
                    indexedexp[i][j] = 0
            elif sym == 'nosym':
                pass
            else:
                raise Exception('unsupported symmetry option \'' + sym + '\'')
    return indexedexp


def symmetrize_rank3(indexedexp, symmetry, dimension):
    symmetry_, symmetry = symmetry, []
    for sym in symmetry_.split('_'):
        index = 3 if sym[:3] == 'sym' else 4
        if len(sym[index:]) == 3:
            prefix = sym[:index]
            symmetry.append(prefix + sym[index:(index + 2)])
            symmetry.append(prefix + sym[(index + 1):(index + 3)])
        else:
            symmetry.append(sym)
    for sym in (symmetry[k] for n in range(len(symmetry), 0, -1) for k in range(n)):
        sign = 1 if sym[:3] == 'sym' else -1
        for i, j, k in func.product(range(dimension), repeat=3):
            if sym[-2:] == '01':
                if j < i:
                    indexedexp[i][j][k] = sign * indexedexp[j][i][k]
                elif i == j and sign < 0:
                    indexedexp[i][j][k] = 0
            elif sym[-2:] == '02':
                if k < i:
                    indexedexp[i][j][k] = sign * indexedexp[k][j][i]
                elif i == k and sign < 0:
                    indexedexp[i][j][k] = 0
            elif sym[-2:] == '12':
                if k < j:
                    indexedexp[i][j][k] = sign * indexedexp[i][k][j]
                elif j == k and sign < 0:
                    indexedexp[i][j][k] = 0
            elif sym == 'nosym':
                pass
            else:
                raise Exception('unsupported symmetry option \'' + sym + '\'')
    return indexedexp


def symmetrize_rank4(indexedexp, symmetry, dimension):
    symmetry_, symmetry = symmetry, []
    for sym in symmetry_.split('_'):
        index = 3 if sym[:3] == 'sym' else 4
        if len(sym[index:]) in (3, 4):
            prefix = sym[:index]
            symmetry.append(prefix + sym[index:(index + 2)])
            symmetry.append(prefix + sym[(index + 1):(index + 3)])
            if len(sym[index:]) == 4:
                symmetry.append(prefix + sym[(index + 2):(index + 4)])
        else:
            symmetry.append(sym)
    for sym in (symmetry[k] for n in range(len(symmetry), 0, -1) for k in range(n)):
        sign = 1 if sym[:3] == 'sym' else -1
        for i, j, k, l in func.product(range(dimension), repeat=4):
            if sym[-2:] == '01':
                if j < i:
                    indexedexp[i][j][k][l] = sign * indexedexp[j][i][k][l]
                elif i == j and sign < 0:
                    indexedexp[i][j][k][l] = 0
            elif sym[-2:] == '02':
                if k < i:
                    indexedexp[i][j][k][l] = sign * indexedexp[k][j][i][l]
                elif i == k and sign < 0:
                    indexedexp[i][j][k][l] = 0
            elif sym[-2:] == '03':
                if l < i:
                    indexedexp[i][j][k][l] = sign * indexedexp[l][j][k][i]
                elif i == l and sign < 0:
                    indexedexp[i][j][k][l] = 0
            elif sym[-2:] == '12':
                if k < j:
                    indexedexp[i][j][k][l] = sign * indexedexp[i][k][j][l]
                elif j == k and sign < 0:
                    indexedexp[i][j][k][l] = 0
            elif sym[-2:] == '13':
                if l < j:
                    indexedexp[i][j][k][l] = sign * indexedexp[i][l][k][j]
                elif j == l and sign < 0:
                    indexedexp[i][j][k][l] = 0
            elif sym[-2:] == '23':
                if l < k:
                    indexedexp[i][j][k][l] = sign * indexedexp[i][j][l][k]
                elif k == l and sign < 0:
                    indexedexp[i][j][k][l] = 0
            elif sym == 'nosym':
                pass
            else:
                raise Exception('unsupported symmetry option \'' + sym + '\'')
    return indexedexp


def zerorank1(DIM=-1):
    return declare_indexedexp(rank=1, dimension=DIM)


def zerorank2(DIM=-1):
    return declare_indexedexp(rank=2, dimension=DIM)


def zerorank3(DIM=-1):
    return declare_indexedexp(rank=3, dimension=DIM)


def zerorank4(DIM=-1):
    return declare_indexedexp(rank=4, dimension=DIM)


def apply_symmetry_condition_to_derivatives(IDX_OBJ):
    symmetry_axes = par.parval_from_str("indexedexp::symmetry_axes")
    if symmetry_axes == "":
        return IDX_OBJ

    rank = 1
    if isinstance(IDX_OBJ[0], list):
        if not isinstance(IDX_OBJ[0][0], list):
            rank = 2
        elif not isinstance(IDX_OBJ[0][0][0], list):
            rank = 3
        elif not isinstance(IDX_OBJ[0][0][0][0], list):
            rank = 4
        else:
            print("Error: could not figure out rank for ", IDX_OBJ)
            sys.exit(1)

    def does_IDXOBJ_perform_derivative_across_symmetry_axis(idxobj_str):
        if "_d" in idxobj_str:
            # First we find the order of the derivative:
            deriv_order = 0
            for i in range(len(idxobj_str) - 1):
                if idxobj_str[i] == "_" and idxobj_str[i + 1] == "d":
                    # The order of the derivative is given by the number of D's in a row after the _d:
                    for k in range(i + 2, len(idxobj_str)):
                        if idxobj_str[k] == "D":
                            deriv_order = deriv_order + 1
            if deriv_order > 2:
                print("Error. Derivative order > 2 not supported. Found derivative order = " + str(deriv_order))
                sys.exit(1)
            end_idx_of_idxobj_str = len(idxobj_str) - 1
            for j in range(end_idx_of_idxobj_str, end_idx_of_idxobj_str - deriv_order, -1):
                if idxobj_str[j] in symmetry_axes:
                    return True
        return False

    if rank == 1:
        DIM = len(IDX_OBJ)
        for i0 in range(DIM):
            if does_IDXOBJ_perform_derivative_across_symmetry_axis(str(IDX_OBJ[i0])) == True:
                IDX_OBJ[i0] = sp.sympify(0)
    if rank == 2:
        DIM = len(IDX_OBJ[0])
        for i0 in range(DIM):
            for i1 in range(DIM):
                if does_IDXOBJ_perform_derivative_across_symmetry_axis(str(IDX_OBJ[i0][i1])) == True:
                    IDX_OBJ[i0][i1] = sp.sympify(0)
    if rank == 3:
        DIM = len(IDX_OBJ[0][0])
        for i0 in range(DIM):
            for i1 in range(DIM):
                for i2 in range(DIM):
                    if does_IDXOBJ_perform_derivative_across_symmetry_axis(str(IDX_OBJ[i0][i1][i2])) == True:
                        IDX_OBJ[i0][i1][i2] = sp.sympify(0)
    if rank == 4:
        DIM = len(IDX_OBJ[0][0][0])
        for i0 in range(DIM):
            for i1 in range(DIM):
                for i2 in range(DIM):
                    for i3 in range(DIM):
                        if does_IDXOBJ_perform_derivative_across_symmetry_axis(str(IDX_OBJ[i0][i1][i2][i3])) == True:
                            IDX_OBJ[i0][i1][i2][i3] = sp.sympify(0)
    return IDX_OBJ


def declarerank1(symbol, DIM=-1, namefun=None):
    return declare_indexedexp(rank=1, symbol=symbol, dimension=DIM, namefun=namefun)


def register_gridfunctions_for_single_rank1(gf_type, gf_basename, DIM=-1, f_infinity=0.0, wavespeed=1.0,
                                            external_module=None, centering=None, namefun=None):
    # Step 0: Verify the gridfunction basename is valid:
    gri.verify_gridfunction_basename_is_valid(gf_basename)

    # Step 1: Declare a list of SymPy variables,
    #         where IDX_OBJ_TMP[i] = gf_basename+str(i)
    IDX_OBJ_TMP = declarerank1(gf_basename, DIM, namefun=namefun)

    # Step 2: Register each gridfunction
    if DIM == -1:
        DIM = par.parval_from_str("DIM")
    gf_list = []
    for i in range(DIM):
        gf_list.append(str(IDX_OBJ_TMP[i]))
    gri.register_gridfunctions(gf_type, gf_list, rank=1, is_indexed=True, DIM=DIM, f_infinity=f_infinity,
                               wavespeed=wavespeed, external_module=external_module, centering=centering)

    # Step 3: Return array of SymPy variables
    return IDX_OBJ_TMP


def declarerank2(symbol, symmetry, DIM=-1, namefun=None):
    return declare_indexedexp(rank=2, symbol=symbol, symmetry=symmetry, dimension=DIM, namefun=namefun)


def register_gridfunctions_for_single_rank2(gf_type, gf_basename, symmetry_option, DIM=-1, f_infinity=0.0,
                                            wavespeed=1.0, external_module=None, centering=None, namefun=None):
    # Step 0: Verify the gridfunction basename is valid:
    gri.verify_gridfunction_basename_is_valid(gf_basename)

    # Step 1: Declare a list of lists of SymPy variables,
    #         where IDX_OBJ_TMP[i][j] = gf_basename+str(i)+str(j)
    IDX_OBJ_TMP = declarerank2(gf_basename, symmetry_option, DIM, namefun)

    # Step 2: register each gridfunction, being careful not
    #         not to store duplicates due to rank-2 symmetries.
    if DIM == -1:
        DIM = par.parval_from_str("DIM")
    # Register only unique gridfunctions. Otherwise
    # rank-2 symmetries might result in duplicates
    gf_list = []
    for i in range(DIM):
        for j in range(DIM):
            save = True
            for l in range(len(gf_list)):
                if gf_list[l] == str(IDX_OBJ_TMP[i][j]):
                    save = False
            if save == True:
                gf_list.append(str(IDX_OBJ_TMP[i][j]))

    gri.register_gridfunctions(gf_type, gf_list, rank=2, is_indexed=True, DIM=DIM,
                               f_infinity=f_infinity, wavespeed=wavespeed,
                               external_module=external_module, centering=centering)

    # Step 3: Return array of SymPy variables
    return IDX_OBJ_TMP


def make_gf_set(gf_set, rank, IDX_TMP_OBJ, DIM):
    """
    Called by register_gridfunctions_for_single_rankN to
    generate a list of gridfunction names.
    """
    if rank == 0:
        gf_set[str(IDX_TMP_OBJ)] = 1
        return
    for d in range(DIM):
        make_gf_set(gf_set, rank - 1, IDX_TMP_OBJ[d], DIM)


def add_index_group(basename, IDX_OBJ_TMP):
    if type(IDX_OBJ_TMP) == list:
        for row in IDX_OBJ_TMP:
            add_index_group(basename, row)
    else:
        s = str(IDX_OBJ_TMP)
        index_group[basename][s] = 1
        rev_index_group[s] = basename


def register_gridfunctions_for_single_rankN(rank, gf_type, gf_basename, symmetry_option="", DIM=-1, f_infinity=0.0,
                                            wavespeed=1.0, external_module=None, centering=None, namefun=None):
    if rank == 0:
        return gri.register_gridfunctions(gf_type, [gf_basename], rank=rank, is_indexed=False, DIM=DIM,
                                          f_infinity=f_infinity, wavespeed=wavespeed,
                                          external_module=external_module, centering=centering)

    # Step 0: Verify the gridfunction basename is valid:
    gri.verify_gridfunction_basename_is_valid(gf_basename)

    # Step 1: Declare a list of lists of SymPy variables,
    #         where IDX_OBJ_TMP[i][j] = gf_basename+str(i)+str(j)
    IDX_OBJ_TMP = declare_indexedexp(rank=rank, symbol=gf_basename, symmetry=symmetry_option, dimension=DIM,
                                     namefun=namefun)

    # Step 2: register each gridfunction, being careful not
    #         not to store duplicates due to rank-2 symmetries.
    if DIM == -1:
        DIM = par.parval_from_str("DIM")

    # Step 3: generate the list of grid function names
    gf_set = {}
    make_gf_set(gf_set, rank, IDX_OBJ_TMP, DIM)
    gf_list = list(gf_set.keys())

    gri.register_gridfunctions(gf_type, gf_list, rank=rank, is_indexed=True, DIM=DIM,
                               f_infinity=f_infinity, wavespeed=wavespeed,
                               external_module=external_module, centering=centering)

    assert gf_basename not in index_group, "Duplicate use of grid function " + gf_basename
    index_group[gf_basename] = {}
    add_index_group(gf_basename, IDX_OBJ_TMP)

    # Step 4: Return array of SymPy variables
    return IDX_OBJ_TMP


def declarerank3(symbol, symmetry, DIM=-1):
    return declare_indexedexp(rank=3, symbol=symbol, symmetry=symmetry, dimension=DIM)


def declarerank4(symbol, symmetry, DIM=-1):
    return declare_indexedexp(rank=4, symbol=symbol, symmetry=symmetry, dimension=DIM)


class NonInvertibleMatrixError(ZeroDivisionError):
    """ Matrix Not Invertible; Division By Zero """


# We use the following functions to evaluate 3-metric inverses
def symm_matrix_inverter2x2(a):
    # It is far more efficient to write out the matrix determinant and inverse by hand
    #   instead of using SymPy's built-in functions, since the matrix is symmetric.
    outDET = a[0][0] * a[1][1] - a[0][1] ** 2
    if outDET == 0: raise NonInvertibleMatrixError('matrix has determinant zero')

    outINV = [[sp.sympify(0) for i in range(2)] for j in range(2)]

    # First fill in the upper-triangle of the gPhysINV matrix...
    outINV[0][0] = a[1][1] / outDET
    outINV[0][1] = -a[0][1] / outDET
    outINV[1][1] = a[0][0] / outDET
    outINV[1][0] = outINV[0][1]
    return outINV, outDET


def symm_matrix_inverter3x3(a):
    # It is far more efficient to write out the matrix determinant and inverse by hand
    #   instead of using SymPy's built-in functions, since the matrix is symmetric.
    outDET = -a[0][2] ** 2 * a[1][1] + 2 * a[0][1] * a[0][2] * a[1][2] - \
             a[0][0] * a[1][2] ** 2 - a[0][1] ** 2 * a[2][2] + \
             a[0][0] * a[1][1] * a[2][2]
    if outDET == 0:
        # print(a)
        raise NonInvertibleMatrixError('matrix has determinant zero')

    outINV = [[sp.sympify(0) for i in range(3)] for j in range(3)]

    # First fill in the upper-triangle of the gPhysINV matrix...
    outINV[0][0] = (-a[1][2]**2              + a[1][1]*a[2][2])/outDET
    outINV[0][1] = (+a[0][2]*a[1][2] - a[0][1]*a[2][2])/outDET
    outINV[0][2] = (-a[0][2]*a[1][1] + a[0][1]*a[1][2])/outDET
    outINV[1][1] = (-a[0][2]**2              + a[0][0]*a[2][2])/outDET
    outINV[1][2] = (+a[0][1]*a[0][2] - a[0][0]*a[1][2])/outDET
    outINV[2][2] = (-a[0][1]**2              + a[0][0]*a[1][1])/outDET
    outINV[1][0] = outINV[0][1]
    outINV[2][0] = outINV[0][2]
    outINV[2][1] = outINV[1][2]
    return outINV, outDET


# Validation test code for symm_matrix_inverter4x4():
# import indexedexp as ixp
# R4DD = ixp.declarerank2("R4DD", "sym01", DIM=4)
#
# # Compute R4DD's inverse:
# R4DDinv, det = ixp.symm_matrix_inverter4x4(R4DD)
#
# # Next matrix multiply: IsUnit = R^{-1} R
# IsUnit = ixp.zerorank2(DIM=4)
# for i in range(4):
#     for j in range(4):
#         for k in range(4):
#             IsUnit[i][j] += R4DDinv[i][k] * R4DD[k][j]
#             # If you'd like to check R R^{-1} instead:
#             # IsUnit[i][j] += R4DD[i][k] * R4DDinv[k][j]
#
# # Next check, is IsUnit == Unit matrix?!
# from UnitTesting.assert_equal import check_zero
# for diag in range(4):
#     print(check_zero(IsUnit[diag][diag]-1))
# for offdiag_i in range(4):
#     for offdiag_j in range(4):
#         if offdiag_i != offdiag_j:
#             print(check_zero(IsUnit[offdiag_i][offdiag_j]))
# # ^^ all should output as True.
def symm_matrix_inverter4x4(a):
    # It is far more efficient to write out the matrix determinant and inverse by hand
    #   instead of using SymPy's built-in functions, since the matrix is symmetric.
    outDET = + a[0][2]*a[0][2]*a[1][3]*a[1][3] + a[0][3]*a[0][3]*a[1][2]*a[1][2] + a[0][1]*a[0][1]*a[2][3]*a[2][3] \
             - a[0][0]*a[1][3]*a[1][3]*a[2][2] - a[0][3]*a[0][3]*a[1][1]*a[2][2] - a[0][0]*a[1][1]*a[2][3]*a[2][3] \
             - 2*(+ a[0][1]*a[0][2]*a[1][3]*a[2][3] - a[0][0]*a[1][2]*a[1][3]*a[2][3]                              \
                  - a[0][3]*(- a[0][2]*a[1][2]*a[1][3] + a[0][1]*a[1][3]*a[2][2]                                   \
                             + a[0][2]*a[1][1]*a[2][3] - a[0][1]*a[1][2]*a[2][3]))                                 \
             - a[3][3] * (+ a[0][2]*a[0][2]*a[1][1] - a[0][1]*a[0][2]*a[1][2] - a[0][1]*a[0][2]*a[1][2]            \
                          + a[0][0]*a[1][2]*a[1][2] + a[0][1]*a[0][1]*a[2][2] - a[0][0]*a[1][1]*a[2][2])
    if outDET == 0: raise NonInvertibleMatrixError('matrix has determinant zero')

    outINV = [[sp.sympify(0) for i in range(4)] for j in range(4)]

    # First fill in the upper-triangle of the gPhysINV matrix...
    outINV[0][0] = (-a[1][3]*a[1][3]*a[2][2] + 2*a[1][2]*a[1][3]*a[2][3] - a[1][1]*a[2][3]*a[2][3] - a[1][2]*a[1][2]*a[3][3] + a[1][1]*a[2][2]*a[3][3])/outDET
    outINV[1][1] = (-a[0][3]*a[0][3]*a[2][2] + 2*a[0][2]*a[0][3]*a[2][3] - a[0][0]*a[2][3]*a[2][3] - a[0][2]*a[0][2]*a[3][3] + a[0][0]*a[2][2]*a[3][3])/outDET
    outINV[2][2] = (-a[0][3]*a[0][3]*a[1][1] + 2*a[0][1]*a[0][3]*a[1][3] - a[0][0]*a[1][3]*a[1][3] - a[0][1]*a[0][1]*a[3][3] + a[0][0]*a[1][1]*a[3][3])/outDET
    outINV[3][3] = (-a[0][2]*a[0][2]*a[1][1] + 2*a[0][1]*a[0][2]*a[1][2] - a[0][0]*a[1][2]*a[1][2] - a[0][1]*a[0][1]*a[2][2] + a[0][0]*a[1][1]*a[2][2])/outDET
    outINV[0][1] = (+a[0][3]*a[1][3]*a[2][2] -   a[0][3]*a[1][2]*a[2][3] - a[0][2]*a[1][3]*a[2][3] + a[0][1]*a[2][3]*a[2][3] + a[0][2]*a[1][2]*a[3][3] - a[0][1]*a[2][2]*a[3][3])/outDET
    outINV[0][2] = (-a[0][3]*a[1][2]*a[1][3] +   a[0][2]*a[1][3]*a[1][3] + a[0][3]*a[1][1]*a[2][3] - a[0][1]*a[1][3]*a[2][3] - a[0][2]*a[1][1]*a[3][3] + a[0][1]*a[1][2]*a[3][3])/outDET
    outINV[0][3] = (-a[0][2]*a[1][2]*a[1][3] +   a[0][1]*a[1][3]*a[2][2] + a[0][3]*a[1][2]*a[1][2] - a[0][3]*a[1][1]*a[2][2] + a[0][2]*a[1][1]*a[2][3] - a[0][1]*a[1][2]*a[2][3])/outDET
    outINV[1][2] = (+a[0][3]*a[0][3]*a[1][2] +   a[0][0]*a[1][3]*a[2][3] - a[0][3]*a[0][2]*a[1][3] - a[0][3]*a[0][1]*a[2][3] + a[0][1]*a[0][2]*a[3][3] - a[0][0]*a[1][2]*a[3][3])/outDET
    outINV[1][3] = (+a[0][2]*a[0][2]*a[1][3] +   a[0][1]*a[0][3]*a[2][2] - a[0][0]*a[1][3]*a[2][2] + a[0][0]*a[1][2]*a[2][3] - a[0][2]*a[0][3]*a[1][2] - a[0][2]*a[0][1]*a[2][3])/outDET
    outINV[2][3] = (+a[0][2]*a[0][3]*a[1][1] -   a[0][1]*a[0][3]*a[1][2] - a[0][1]*a[0][2]*a[1][3] + a[0][0]*a[1][2]*a[1][3] + a[0][1]*a[0][1]*a[2][3] - a[0][0]*a[1][1]*a[2][3])/outDET

    # Then we fill the lower triangle of the symmetric matrix
    outINV[1][0] = outINV[0][1]
    outINV[2][0] = outINV[0][2]
    outINV[2][1] = outINV[1][2]
    outINV[3][0] = outINV[0][3]
    outINV[3][1] = outINV[1][3]
    outINV[3][2] = outINV[2][3]

    return outINV, outDET


# SymPy's generic matrix inverter takes a long time to invert 3x3 matrices, so here we have an optimized version.
# We use the following functions to evaluate 3-metric inverses
def generic_matrix_inverter2x2(a):
    outDET = a[0][0] * a[1][1] - a[0][1] * a[1][0]
    if outDET == 0: raise NonInvertibleMatrixError('matrix has determinant zero')

    outINV = [[sp.sympify(0) for i in range(2)] for j in range(2)]

    outINV[0][0] = a[1][1] / outDET
    outINV[0][1] = -a[0][1] / outDET
    outINV[1][1] = a[0][0] / outDET
    outINV[1][0] = -a[1][0] / outDET
    return outINV, outDET


def generic_matrix_inverter3x3(a):
    outDET = -a[0][2]*a[1][1]*a[2][0] + a[0][1]*a[1][2]*a[2][0] + \
              a[0][2]*a[1][0]*a[2][1] - a[0][0]*a[1][2]*a[2][1] - \
              a[0][1]*a[1][0]*a[2][2] + a[0][0]*a[1][1]*a[2][2]
    if outDET == 0: raise NonInvertibleMatrixError('matrix has determinant zero')

    outINV = [[sp.sympify(0) for i in range(3)] for j in range(3)]

    outINV[0][0] = -a[1][2]*a[2][1] + a[1][1]*a[2][2]
    outINV[0][1] =  a[0][2]*a[2][1] - a[0][1]*a[2][2]
    outINV[0][2] = -a[0][2]*a[1][1] + a[0][1]*a[1][2]
    outINV[1][0] =  a[1][2]*a[2][0] - a[1][0]*a[2][2]
    outINV[1][1] = -a[0][2]*a[2][0] + a[0][0]*a[2][2]
    outINV[1][2] =  a[0][2]*a[1][0] - a[0][0]*a[1][2]
    outINV[2][0] = -a[1][1]*a[2][0] + a[1][0]*a[2][1]
    outINV[2][1] =  a[0][1]*a[2][0] - a[0][0]*a[2][1]
    outINV[2][2] = -a[0][1]*a[1][0] + a[0][0]*a[1][1]

    for i in range(3):
        for j in range(3):
            outINV[i][j] /= outDET

    return outINV, outDET


def generic_matrix_inverter4x4(a):
    # A = {{a00, a01, a02, a03},
    #      {a10, a11, a12, a13},
    #      {a20, a21, a22, a23},
    #      {a30, a31, a32, a33}}
    # A // MatrixForm
    # CForm[FullSimplify[Det[A]]] >>> t2.txt
    # cat t2.txt | sed "s/ //g" |sed "s/ //g;s/\([0-3]\)/[\1]/g"
    outDET = a[0][1]*a[1][3]*a[2][2]*a[3][0]-a[0][1]*a[1][2]*a[2][3]*a[3][0]-a[0][0]*a[1][3]*a[2][2]*a[3][1]+ \
             a[0][0]*a[1][2]*a[2][3]*a[3][1]-a[0][1]*a[1][3]*a[2][0]*a[3][2]+a[0][0]*a[1][3]*a[2][1]*a[3][2]+ \
             a[0][1]*a[1][0]*a[2][3]*a[3][2]-a[0][0]*a[1][1]*a[2][3]*a[3][2]+ \
        a[0][3]*(a[1][2]*a[2][1]*a[3][0]-a[1][1]*a[2][2]*a[3][0]-a[1][2]*a[2][0]*a[3][1]+a[1][0]*a[2][2]*a[3][1]+
                 a[1][1]*a[2][0]*a[3][2]-a[1][0]*a[2][1]*a[3][2])+ \
             (a[0][1]*a[1][2]*a[2][0]-a[0][0]*a[1][2]*a[2][1]-a[0][1]*a[1][0]*a[2][2]+a[0][0]*a[1][1]*a[2][2])*a[3][3]+\
        a[0][2]*(-(a[1][3]*a[2][1]*a[3][0])+a[1][1]*a[2][3]*a[3][0]+a[1][3]*a[2][0]*a[3][1]-a[1][0]*a[2][3]*a[3][1]-
                 a[1][1]*a[2][0]*a[3][3]+a[1][0]*a[2][1]*a[3][3])
    if outDET == 0: raise NonInvertibleMatrixError('matrix has determinant zero')

    outINV = [[sp.sympify(0) for i in range(4)] for j in range(4)]

    # CForm[FullSimplify[Inverse[A]*Det[A]]] >>> t.txt
    # cat t.txt | sed "s/,/\n/g;s/List(//g;s/))/)/g;s/)//g;s/(//g"|grep -v ^$|sed "s/ //g;s/\([0-3]\)/[\1]/g"| awk '{line[NR]=$0}END{count=1;for(i=0;i<4;i++) { for(j=0;j<4;j++) { printf "outINV[%d][%d] = %s\n", i,j,line[count];count++; }}}'
    outINV[0][0] = -a[1][3]*a[2][2]*a[3][1]+a[1][2]*a[2][3]*a[3][1]+a[1][3]*a[2][1]*a[3][2]-a[1][1]*a[2][3]*a[3][2]-a[1][2]*a[2][1]*a[3][3]+a[1][1]*a[2][2]*a[3][3]
    outINV[0][1] =  a[0][3]*a[2][2]*a[3][1]-a[0][2]*a[2][3]*a[3][1]-a[0][3]*a[2][1]*a[3][2]+a[0][1]*a[2][3]*a[3][2]+a[0][2]*a[2][1]*a[3][3]-a[0][1]*a[2][2]*a[3][3]
    outINV[0][2] = -a[0][3]*a[1][2]*a[3][1]+a[0][2]*a[1][3]*a[3][1]+a[0][3]*a[1][1]*a[3][2]-a[0][1]*a[1][3]*a[3][2]-a[0][2]*a[1][1]*a[3][3]+a[0][1]*a[1][2]*a[3][3]
    outINV[0][3] =  a[0][3]*a[1][2]*a[2][1]-a[0][2]*a[1][3]*a[2][1]-a[0][3]*a[1][1]*a[2][2]+a[0][1]*a[1][3]*a[2][2]+a[0][2]*a[1][1]*a[2][3]-a[0][1]*a[1][2]*a[2][3]
    outINV[1][0] =  a[1][3]*a[2][2]*a[3][0]-a[1][2]*a[2][3]*a[3][0]-a[1][3]*a[2][0]*a[3][2]+a[1][0]*a[2][3]*a[3][2]+a[1][2]*a[2][0]*a[3][3]-a[1][0]*a[2][2]*a[3][3]
    outINV[1][1] = -a[0][3]*a[2][2]*a[3][0]+a[0][2]*a[2][3]*a[3][0]+a[0][3]*a[2][0]*a[3][2]-a[0][0]*a[2][3]*a[3][2]-a[0][2]*a[2][0]*a[3][3]+a[0][0]*a[2][2]*a[3][3]
    outINV[1][2] =  a[0][3]*a[1][2]*a[3][0]-a[0][2]*a[1][3]*a[3][0]-a[0][3]*a[1][0]*a[3][2]+a[0][0]*a[1][3]*a[3][2]+a[0][2]*a[1][0]*a[3][3]-a[0][0]*a[1][2]*a[3][3]
    outINV[1][3] = -a[0][3]*a[1][2]*a[2][0]+a[0][2]*a[1][3]*a[2][0]+a[0][3]*a[1][0]*a[2][2]-a[0][0]*a[1][3]*a[2][2]-a[0][2]*a[1][0]*a[2][3]+a[0][0]*a[1][2]*a[2][3]
    outINV[2][0] = -a[1][3]*a[2][1]*a[3][0]+a[1][1]*a[2][3]*a[3][0]+a[1][3]*a[2][0]*a[3][1]-a[1][0]*a[2][3]*a[3][1]-a[1][1]*a[2][0]*a[3][3]+a[1][0]*a[2][1]*a[3][3]
    outINV[2][1] =  a[0][3]*a[2][1]*a[3][0]-a[0][1]*a[2][3]*a[3][0]-a[0][3]*a[2][0]*a[3][1]+a[0][0]*a[2][3]*a[3][1]+a[0][1]*a[2][0]*a[3][3]-a[0][0]*a[2][1]*a[3][3]
    outINV[2][2] = -a[0][3]*a[1][1]*a[3][0]+a[0][1]*a[1][3]*a[3][0]+a[0][3]*a[1][0]*a[3][1]-a[0][0]*a[1][3]*a[3][1]-a[0][1]*a[1][0]*a[3][3]+a[0][0]*a[1][1]*a[3][3]
    outINV[2][3] =  a[0][3]*a[1][1]*a[2][0]-a[0][1]*a[1][3]*a[2][0]-a[0][3]*a[1][0]*a[2][1]+a[0][0]*a[1][3]*a[2][1]+a[0][1]*a[1][0]*a[2][3]-a[0][0]*a[1][1]*a[2][3]
    outINV[3][0] =  a[1][2]*a[2][1]*a[3][0]-a[1][1]*a[2][2]*a[3][0]-a[1][2]*a[2][0]*a[3][1]+a[1][0]*a[2][2]*a[3][1]+a[1][1]*a[2][0]*a[3][2]-a[1][0]*a[2][1]*a[3][2]
    outINV[3][1] = -a[0][2]*a[2][1]*a[3][0]+a[0][1]*a[2][2]*a[3][0]+a[0][2]*a[2][0]*a[3][1]-a[0][0]*a[2][2]*a[3][1]-a[0][1]*a[2][0]*a[3][2]+a[0][0]*a[2][1]*a[3][2]
    outINV[3][2] =  a[0][2]*a[1][1]*a[3][0]-a[0][1]*a[1][2]*a[3][0]-a[0][2]*a[1][0]*a[3][1]+a[0][0]*a[1][2]*a[3][1]+a[0][1]*a[1][0]*a[3][2]-a[0][0]*a[1][1]*a[3][2]
    outINV[3][3] = -a[0][2]*a[1][1]*a[2][0]+a[0][1]*a[1][2]*a[2][0]+a[0][2]*a[1][0]*a[2][1]-a[0][0]*a[1][2]*a[2][1]-a[0][1]*a[1][0]*a[2][2]+a[0][0]*a[1][1]*a[2][2]

    for mu in range(4):
        for nu in range(4):
            outINV[mu][nu] /= outDET

    return outINV, outDET


# Define the rank-3 version of the Levi-Civita symbol.
def LeviCivitaSymbol_dim3_rank3():
    LeviCivitaSymbol = zerorank3(DIM=3)

    for i in range(3):
        for j in range(3):
            for k in range(3):
                # From https://codegolf.stackexchange.com/questions/160359/levi-civita-symbol :
                LeviCivitaSymbol[i][j][k] = (i - j) * (j - k) * (k - i) * sp.Rational(1, 2)
    return LeviCivitaSymbol


# Define the UUU rank-3 version of the Levi-Civita *tensor*; UUU divides by sqrtgammaDET
def LeviCivitaTensorUUU_dim3_rank3(sqrtgammaDET):
    # Here, we import the Levi-Civita tensor and compute the tensor with upper indices
    LeviCivitaSymbolDDD = LeviCivitaSymbol_dim3_rank3()
    LeviCivitaTensorUUU = zerorank3(DIM=3)
    for i in range(3):
        for j in range(3):
            for k in range(3):
                LeviCivitaTensorUUU[i][j][k] = LeviCivitaSymbolDDD[i][j][k] / sqrtgammaDET
    return LeviCivitaTensorUUU


# Define the DDD rank-3 version of the Levi-Civita *tensor*; DDD multiplies by sqrtgammaDET
def LeviCivitaTensorDDD_dim3_rank3(sqrtgammaDET):
    # Here, we import the Levi-Civita tensor and compute the tensor with lower indices
    LeviCivitaSymbolDDD = LeviCivitaSymbol_dim3_rank3()
    LeviCivitaTensorDDD = zerorank3(DIM=3)
    for i in range(3):
        for j in range(3):
            for k in range(3):
                LeviCivitaTensorDDD[i][j][k] = LeviCivitaSymbolDDD[i][j][k] * sqrtgammaDET
    return LeviCivitaTensorDDD


if __name__ == "__main__":
    import doctest

    sys.exit(doctest.testmod()[0])
