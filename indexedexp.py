# indexedexp.py: functions related to indexed expressions,
# including e.g., tensors and pseudotensors:

# Step 1: Load needed modules
import NRPy_param_funcs as par   # NRPy+: Parameter interface
import grid as gri               # NRPy+: Functions having to do with numerical grids
import sympy as sp               # SymPy: The Python computer algebra package upon which NRPy+ depends
import sys                       # Standard Python module for multiplatform OS-level functions
from itertools import product, chain

thismodule = __name__
par.initialize_param(par.glb_param("char", thismodule, "symmetry_axes",  ""))

def declare_indexedexp(rank, symbol=None, symmetry=None, dimension=None):
    """ Generate an indexed expression with specified rank

        >>> declare_indexedexp(rank=1, dimension=3)
        [0, 0, 0]
        >>> declare_indexedexp(rank=1, symbol="vU", dimension=3)
        [vU0, vU1, vU2]

        >>> declare_indexedexp(rank=2, dimension=3)
        [[0, 0, 0], [0, 0, 0], [0, 0, 0]]
        >>> declare_indexedexp(rank=2, symbol='gUU', dimension=3)
        [[gUU00, gUU01, gUU02], [gUU10, gUU11, gUU12], [gUU20, gUU21, gUU22]]
        >>> declare_indexedexp(rank=2, symbol='gUU', symmetry='sym01', dimension=3)
        [[gUU00, gUU01, gUU02], [gUU01, gUU11, gUU12], [gUU02, gUU12, gUU22]]

        >>> declare_indexedexp(rank=3, dimension=3)
        [[[0, 0, 0], [0, 0, 0], [0, 0, 0]], [[0, 0, 0], [0, 0, 0], [0, 0, 0]], \
[[0, 0, 0], [0, 0, 0], [0, 0, 0]]]
        >>> declare_indexedexp(rank=3, symbol='gUU', symmetry='sym01', dimension=3)
        [[[gUU000, gUU001, gUU002], [gUU010, gUU011, gUU012], [gUU020, gUU021, gUU022]], \
[[gUU010, gUU011, gUU012], [gUU110, gUU111, gUU112], [gUU120, gUU121, gUU122]], \
[[gUU020, gUU021, gUU022], [gUU120, gUU121, gUU122], [gUU220, gUU221, gUU222]]]
        >>> declare_indexedexp(rank=3, symbol='gUU', symmetry='sym12', dimension=3)
        [[[gUU000, gUU001, gUU002], [gUU001, gUU011, gUU012], [gUU002, gUU012, gUU022]], \
[[gUU100, gUU101, gUU102], [gUU101, gUU111, gUU112], [gUU102, gUU112, gUU122]], \
[[gUU200, gUU201, gUU202], [gUU201, gUU211, gUU212], [gUU202, gUU212, gUU222]]]
    """
    if not dimension or dimension == -1: dimension = par.parval_from_str('DIM')
    loop_index = ['str(%s)' % chr(97 + n) for n in range(rank)]
    indexing = ' + '.join(loop_index)
    interior = 'sp.sympify(0)' if not symbol \
        else 'sp.sympify(\'%s\' + %s)' % (symbol, indexing)
    indexedexp = '[' * rank + interior
    for i in range(1, rank + 1):
        indexedexp += ' for %s in range(%s)]' % (loop_index[rank - i][4:-1], dimension)
    indexedexp = eval(indexedexp)
    if symmetry: return symmetrize(rank, indexedexp, symmetry, dimension)
    return apply_symmetry_condition_to_derivatives(indexedexp)

def symmetrize(rank, indexedexp, symmetry, dimension):
    if rank == 1:
        raise Exception ('cannot symmetrize indexed expression of rank 1')
    if rank == 2:
        for sym in symmetry.split('_'):
            for i, j in product(range(dimension), repeat=2):
                if sym == 'sym01':
                    if j < i: indexedexp[i][j] = indexedexp[j][i]
                elif sym == 'nosym': pass
                else: raise Exception('unsupported symmetry option \'' + sym + '\'')
    elif rank == 3:
        symmetry = set(chain.from_iterable(('sym' + sym[3:5], 'sym' + sym[4:6]) 
            if len(sym[3:]) == 3 else (sym,) for sym in symmetry.split('_')))
        for sym in chain(symmetry, reversed(list(symmetry)[:-1])):
            for i, j, k in product(range(dimension), repeat=3):
                if sym == 'sym01':
                    if j < i: indexedexp[i][j][k] = indexedexp[j][i][k]
                elif sym == 'sym02':
                    if k < i: indexedexp[i][j][k] = indexedexp[k][j][i]
                elif sym == 'sym12':
                    if k < j: indexedexp[i][j][k] = indexedexp[i][k][j]
                elif sym == 'nosym': pass
                else: raise Exception('unsupported symmetry option \'' + sym + '\'')
    elif rank == 4:
        symmetry = set(chain.from_iterable(('sym' + sym[3:5], 'sym' + sym[4:6]) 
            if len(sym[3:]) == 3 else ('sym' + sym[3:5], 'sym' + sym[4:6], 'sym' + sym[5:7]) 
            if len(sym[3:]) == 4 else (sym,) for sym in symmetry.split('_')))
        for sym in chain(symmetry, reversed(list(symmetry)[:-1])):
            for i, j, k, l in product(range(dimension), repeat=4):
                if sym == 'sym01':
                    if j < i: indexedexp[i][j][k][l] = indexedexp[j][i][k][l]
                elif sym == 'sym02':
                    if k < i: indexedexp[i][j][k][l] = indexedexp[k][j][i][l]
                elif sym == 'sym03':
                    if l < i: indexedexp[i][j][k][l] = indexedexp[l][j][k][i]
                elif sym == 'sym12':
                    if k < j: indexedexp[i][j][k][l] = indexedexp[i][k][j][l]
                elif sym == 'sym13':
                    if l < j: indexedexp[i][j][k][l] = indexedexp[i][l][k][j]
                elif sym == 'sym23':
                    if l < k: indexedexp[i][j][k][l] = indexedexp[i][j][l][k]
                elif sym == 'nosym': pass
                else: raise Exception('unsupported symmetry option \'' + sym + '\'')
    else: raise Exception('unsupported rank for indexed expression')
    return apply_symmetry_condition_to_derivatives(indexedexp)

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
            print("Error: could not figure out rank for ",IDX_OBJ)
            sys.exit(1)

    def does_IDXOBJ_perform_derivative_across_symmetry_axis(idxobj_str):
        if "_d" in idxobj_str:
            # First we find the order of the derivative:
            deriv_order = 0
            for i in range(len(idxobj_str)-1):
                if idxobj_str[i] == "_" and idxobj_str[i+1]=="d":
                    # The order of the derivative is given by the number of D's in a row after the _d:
                    for k in range(i+2,len(idxobj_str)):
                        if idxobj_str[k] == "D":
                            deriv_order = deriv_order + 1
            if deriv_order > 2:
                print("Error. Derivative order > 2 not supported. Found derivative order = "+str(deriv_order))
                sys.exit(1)
            end_idx_of_idxobj_str = len(idxobj_str)-1
            for j in range(end_idx_of_idxobj_str,end_idx_of_idxobj_str-deriv_order,-1):
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

def declarerank1(symbol, DIM=-1):
    return declare_indexedexp(rank=1, symbol=symbol, dimension=DIM)

def register_gridfunctions_for_single_rank1(gf_type,gf_basename, DIM=-1):
    # Step 0: Verify the gridfunction basename is valid:
    gri.verify_gridfunction_basename_is_valid(gf_basename)

    # Step 1: Declare a list of SymPy variables,
    #         where IDX_OBJ_TMP[i] = gf_basename+str(i)
    IDX_OBJ_TMP = declarerank1(gf_basename, DIM)

    # Step 2: Register each gridfunction
    if DIM==-1:
        DIM = par.parval_from_str("DIM")
    gf_list = []
    for i in range(DIM):
        gf_list.append(str(IDX_OBJ_TMP[i]))
    gri.register_gridfunctions(gf_type, gf_list, rank=1, is_indexed=True, DIM=DIM)

    # Step 3: Return array of SymPy variables
    return IDX_OBJ_TMP

def declarerank2(symbol, symmetry, DIM=-1):
    return declare_indexedexp(rank=2, symbol=symbol, symmetry=symmetry, dimension=DIM)

def register_gridfunctions_for_single_rank2(gf_type,gf_basename, symmetry_option, DIM=-1):
    # Step 0: Verify the gridfunction basename is valid:
    gri.verify_gridfunction_basename_is_valid(gf_basename)

    # Step 1: Declare a list of lists of SymPy variables,
    #         where IDX_OBJ_TMP[i][j] = gf_basename+str(i)+str(j)
    IDX_OBJ_TMP = declarerank2(gf_basename,symmetry_option, DIM)

    # Step 2: register each gridfunction, being careful not
    #         not to store duplicates due to rank-2 symmetries.
    if DIM==-1:
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
    gri.register_gridfunctions(gf_type,gf_list,rank=2, is_indexed=True, DIM=DIM)

    # Step 3: Return array of SymPy variables
    return IDX_OBJ_TMP

def declarerank3(symbol, symmetry, DIM=-1):
    return declare_indexedexp(rank=3, symbol=symbol, symmetry=symmetry, dimension=DIM)

def declarerank4(symbol, symmetry, DIM=-1):
    return declare_indexedexp(rank=4, symbol=symbol, symmetry=symmetry, dimension=DIM)


# We use the following functions to evaluate 3-metric inverses
def symm_matrix_inverter3x3(a):
    # It is far more efficient to write out the matrix determinant and inverse by hand
    #   instead of using SymPy's built-in functions, since the matrix is symmetric.
    outDET = -a[0][2]**2*a[1][1] + 2*a[0][1]*a[0][2]*a[1][2] - \
                a[0][0]*a[1][2]**2 - a[0][1]**2*a[2][2] + \
                a[0][0]*a[1][1]*a[2][2]

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
def generic_matrix_inverter3x3(a):
    outDET = -a[0][2]*a[1][1]*a[2][0] + a[0][1]*a[1][2]*a[2][0] + \
              a[0][2]*a[1][0]*a[2][1] - a[0][0]*a[1][2]*a[2][1] - \
              a[0][1]*a[1][0]*a[2][2] + a[0][0]*a[1][1]*a[2][2]

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
                LeviCivitaSymbol[i][j][k] = (i - j) * (j - k) * (k - i) * sp.Rational(1,2)
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
