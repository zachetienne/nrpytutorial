# indexedexp.py: functions related to indexed expressions, 
# including e.g., tensors and pseudotensors:

# Step 1: Load needed modules
import NRPy_param_funcs as par
import grid as gri
import sympy as sp

thismodule = __name__
par.initialize_param(par.glb_param("char", thismodule, "symmetry_axes",  ""))

def zerorank1(DIM=-1):
    if DIM == -1:
        DIM = par.parval_from_str("DIM")
    return [sp.sympify(0) for i in range(DIM)]

def zerorank2(DIM=-1):
    if DIM == -1:
        DIM = par.parval_from_str("DIM")
    return [[sp.sympify(0) for i in range(DIM)] for j in range(DIM)]

def zerorank3(DIM=-1):
    if DIM == -1:
        DIM = par.parval_from_str("DIM")
    return [[[sp.sympify(0) for i in range(DIM)] for j in range(DIM)] for k in range(DIM)]

def zerorank4(DIM=-1):
    if DIM == -1:
        DIM = par.parval_from_str("DIM")
    return [[[[sp.sympify(0) for i in range(DIM)] for j in range(DIM)] for k in range(DIM)] for l in range(DIM)]

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
            exit(1)

    def does_IDXOBJ_perform_derivative_across_symmetry_axis(idxobj_str):
        returnval = False
        if "_d" in idxobj_str:
            # First we find the order of the derivative:
            deriv_order = 0
            underscore_position = -1000
            for i in range(len(idxobj_str)-1):
                if idxobj_str[i] == "_" and idxobj_str[i+1]=="d":
                    # The order of the derivative is given by the number of D's in a row after the _d:
                    for k in range(i+2,len(idxobj_str)):
                        if idxobj_str[k] == "D":
                            deriv_order = deriv_order + 1
            if deriv_order > 2:
                print("Error. Derivative order > 2 not supported. Found derivative order = "+str(deriv_order))
                exit(1)
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

def declarerank1(objname, DIM=-1):
    if DIM==-1:
        DIM = par.parval_from_str("DIM")
    IDX_OBJ_TMP = [sp.sympify(objname + str(i)) for i in range(DIM)]
    return apply_symmetry_condition_to_derivatives(IDX_OBJ_TMP)

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
    gri.register_gridfunctions(gf_type, gf_list, is_indexed=True)

    # Step 3: Return array of SymPy variables
    return IDX_OBJ_TMP

def declarerank2(objname, symmetry_option, DIM=-1):
    if DIM==-1:
        DIM = par.parval_from_str("DIM")
    IDX_OBJ_TMP = [[sp.sympify(objname + str(i) + str(j)) for j in range(DIM)] for i in range(DIM)]
    for i in range(DIM):
        for j in range(DIM):
            if symmetry_option == "sym01":
                if (j < i):
                    # j<i in g_{ij} would indicate, e.g., g_{21}.
                    #  By this convention, we must set
                    #  g_{21} = g_{12}:
                    IDX_OBJ_TMP[i][j] = IDX_OBJ_TMP[j][i]
            elif symmetry_option == "nosym":
                pass
            else:
                print("Error: symmetry option " + symmetry_option + " unsupported.")
                exit(1)
    return apply_symmetry_condition_to_derivatives(IDX_OBJ_TMP)

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
    gri.register_gridfunctions(gf_type,gf_list, is_indexed=True)

    # Step 3: Return array of SymPy variables
    return IDX_OBJ_TMP

def declarerank3(objname, symmetry_option, DIM=-1):
    if DIM==-1:
        DIM = par.parval_from_str("DIM")
    IDX_OBJ_TMP = [[[sp.sympify(objname + str(i) + str(j) + str(k)) for k in range(DIM)] for j in range(DIM)] for i in range(DIM)]
    for i in range(DIM):
        for j in range(DIM):
            for k in range(DIM):
                if symmetry_option == "sym01":
                    if j < i:
                        IDX_OBJ_TMP[i][j][k] = IDX_OBJ_TMP[j][i][k]
                if symmetry_option == "sym12":
                    if k < j:
                        IDX_OBJ_TMP[i][j][k] = IDX_OBJ_TMP[i][k][j]
                if not (symmetry_option == "sym01" or symmetry_option == "sym12" or symmetry_option == "nosym"):
                    print("Error: symmetry option " + symmetry_option + " unsupported.")
                    exit(1)
    return apply_symmetry_condition_to_derivatives(IDX_OBJ_TMP)

def declarerank4(objname, symmetry_option, DIM=-1):
    if DIM==-1:
        DIM = par.parval_from_str("DIM")
    IDX_OBJ_TMP = [[[[sp.sympify(objname + str(i) + str(j) + str(k) + str(l)) for l in range(DIM)] for k in range(DIM)] for j in range(DIM)] for i in range(DIM)]
    for i in range(DIM):
        for j in range(DIM):
            for k in range(DIM):
                for l in range(DIM):
                    IDX_OBJ_TMP[i][j][k][l] = sp.sympify(objname + str(i) + str(j) + str(k) + str(l))
                    if symmetry_option == "sym01" or symmetry_option == "sym01_sym23":
                        if(j < i):
                            IDX_OBJ_TMP[i][j][k][l] = IDX_OBJ_TMP[j][i][k][l]
                    if symmetry_option == "sym12":
                        if(k < j):
                            IDX_OBJ_TMP[i][j][k][l] = IDX_OBJ_TMP[i][k][j][l]
                    if symmetry_option == "sym23" or symmetry_option == "sym01_sym23":
                        if(l < k):
                            IDX_OBJ_TMP[i][j][k][l] = IDX_OBJ_TMP[i][j][l][k]
                    if not (symmetry_option=="sym01" or symmetry_option=="sym23" or symmetry_option=="sym01_sym23" or symmetry_option=="none"):
                        print("Error: symmetry option "+symmetry_option+" unsupported.")
                        exit(1)
    return apply_symmetry_condition_to_derivatives(IDX_OBJ_TMP)


# We use the following functions to evaluate 3-metric inverses
def symm_matrix_inverter3x3(in2Darray):
    # It is far more efficient to write out the matrix determinant and inverse by hand
    #   instead of using symmetry_optionPy's built-in functions, since the matrix is symmetric.
    outDET = -in2Darray[0][2]**2*in2Darray[1][1] + 2*in2Darray[0][1]*in2Darray[0][2]*in2Darray[1][2] - \
                in2Darray[0][0]*in2Darray[1][2]**2 - in2Darray[0][1]**2*in2Darray[2][2] + \
                in2Darray[0][0]*in2Darray[1][1]*in2Darray[2][2]

    outINV = [[sp.sympify(0) for i in range(3)] for j in range(3)]

    # First fill in the upper-triangle of the gPhysINV matrix...
    outINV[0][0] = (-in2Darray[1][2]**2              + in2Darray[1][1]*in2Darray[2][2])/outDET
    outINV[0][1] = (+in2Darray[0][2]*in2Darray[1][2] - in2Darray[0][1]*in2Darray[2][2])/outDET
    outINV[0][2] = (-in2Darray[0][2]*in2Darray[1][1] + in2Darray[0][1]*in2Darray[1][2])/outDET
    outINV[1][1] = (-in2Darray[0][2]**2              + in2Darray[0][0]*in2Darray[2][2])/outDET
    outINV[1][2] = (+in2Darray[0][1]*in2Darray[0][2] - in2Darray[0][0]*in2Darray[1][2])/outDET
    outINV[2][2] = (-in2Darray[0][1]**2              + in2Darray[0][0]*in2Darray[1][1])/outDET
    outINV[1][0] = outINV[0][1]
    outINV[2][0] = outINV[0][2]
    outINV[2][1] = outINV[1][2]
    return outINV, outDET

# symmetry_optionPy's generic matrix inverter is highly inefficient for 3x3 matrices, so here we have an optimized version.
def generic_matrix_inverter3x3(in2Darray):
    outDET = -in2Darray[0][2]*in2Darray[1][1]*in2Darray[2][0] + in2Darray[0][1]*in2Darray[1][2]*in2Darray[2][0] + \
              in2Darray[0][2]*in2Darray[1][0]*in2Darray[2][1] - in2Darray[0][0]*in2Darray[1][2]*in2Darray[2][1] - \
              in2Darray[0][1]*in2Darray[1][0]*in2Darray[2][2] + in2Darray[0][0]*in2Darray[1][1]*in2Darray[2][2]

    outINV = [[sp.sympify(0) for i in range(3)] for j in range(3)]

    outINV[0][0] = -in2Darray[1][2]*in2Darray[2][1] + in2Darray[1][1]*in2Darray[2][2]
    outINV[0][1] =  in2Darray[0][2]*in2Darray[2][1] - in2Darray[0][1]*in2Darray[2][2]
    outINV[0][2] = -in2Darray[0][2]*in2Darray[1][1] + in2Darray[0][1]*in2Darray[1][2]
    outINV[1][0] =  in2Darray[1][2]*in2Darray[2][0] - in2Darray[1][0]*in2Darray[2][2]
    outINV[1][1] = -in2Darray[0][2]*in2Darray[2][0] + in2Darray[0][0]*in2Darray[2][2]
    outINV[1][2] =  in2Darray[0][2]*in2Darray[1][0] - in2Darray[0][0]*in2Darray[1][2]
    outINV[2][0] = -in2Darray[1][1]*in2Darray[2][0] + in2Darray[1][0]*in2Darray[2][1]
    outINV[2][1] =  in2Darray[0][1]*in2Darray[2][0] - in2Darray[0][0]*in2Darray[2][1]
    outINV[2][2] = -in2Darray[0][1]*in2Darray[1][0] + in2Darray[0][0]*in2Darray[1][1]

    for i in range(3):
        for j in range(3):
            outINV[i][j] /= outDET

    return outINV, outDET
