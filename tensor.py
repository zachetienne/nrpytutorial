# tensor.py: functions related to tensors or, more generally, indexed objects:
# functions: registration of tensors
# Step 1: Load needed modules
import NRPy_param_funcs as par
import grid as gri
import sympy as sp

def zerorank1(DIM=-1):
    if DIM == -1:
        DIM = par.parval_from_str("DIM")
    return [sp.sympify(0) for i in range(DIM)]

def zerorank2(DIM=-1):
    if DIM == -1:
        DIM = par.parval_from_str("DIM")
    return [[sp.sympify(0) for i in range(DIM)] for j in range(DIM)]

def declarerank1(objname, DIM=-1):
    if DIM==-1:
        DIM = par.parval_from_str("DIM")
    IDX_OBJ_TMP = [sp.sympify(objname + str(i)) for i in range(DIM)]
    return IDX_OBJ_TMP

def register_gridfunctions_for_single_rank1(gf_type,gf_basename):
    IDX_OBJ_TMP = declarerank1(gf_basename)
    DIM = par.parval_from_str("DIM")
    stringIDX_OBJ_TMP = []
    for i in range(DIM):
        stringIDX_OBJ_TMP[i] = str(IDX_OBJ_TMP[i])
    gri.register_gridfunctions(gf_type,stringIDX_OBJ_TMP)
    return IDX_OBJ_TMP

def declarerank2(objname, symmetry_option, DIM=-1):
    if DIM==-1:
        DIM = par.parval_from_str("DIM")
    IDX_OBJ_TMP = [[sp.sympify(objname + str(i) + str(j)) for j in range(DIM)] for i in range(DIM)]
    for i in range(DIM):
        for j in range(DIM):
            if symmetry_option == "sym12":
                if (j < i):
                    # j<i in g_{ij} would indicate, e.g., g_{21}.
                    #  By this convention, we must set
                    #  g_{21} = g_{12}:
                    IDX_OBJ_TMP[i][j] = IDX_OBJ_TMP[j][i]
            elif symmetry_option == "none":
                pass
            else:
                print("Error: symmetry option " + symmetry_option + " unsupported.")
                exit(1)
    return IDX_OBJ_TMP

def register_gridfunctions_for_single_rank2(gf_type,gf_basename, symmetry_option):
    # First convert gf_names to a list if it's not already a list
    DIM = par.parval_from_str("DIM")
    IDX_OBJ_TMP = declarerank2(gf_basename,symmetry_option)
    # Now drop to a flat list, so we can register each gridfunction
    gf_list = []
    for i in range(DIM):
        for j in range(DIM):
            save = True
            for l in range(len(gf_list)):
                if gf_list[l] == str(IDX_OBJ_TMP[i][j]):
                    save = False
            if save == True:
                gf_list.append(str(IDX_OBJ_TMP[i][j]))
    gri.register_gridfunctions(gf_type,gf_list)
    return IDX_OBJ_TMP

def declarerank3(objname, symmetry_option, DIM=-1):
    if DIM==-1:
        DIM = par.parval_from_str("DIM")
    IDX_OBJ_TMP = [[[sp.sympify(objname + str(i) + str(j) + str(k)) for k in range(DIM)] for j in range(DIM)] for i in range(DIM)]
    for i in range(DIM):
        for j in range(DIM):
            for k in range(DIM):
                if symmetry_option == "sym12":
                    if j < i:
                        IDX_OBJ_TMP[i][j][k] = IDX_OBJ_TMP[j][i][k]
                if symmetry_option == "sym23":
                    if k < j:
                        IDX_OBJ_TMP[i][j][k] = IDX_OBJ_TMP[i][k][j]
                if not (symmetry_option == "sym12" or symmetry_option == "sym23" or symmetry_option == "none"):
                    print("Error: symmetry option " + symmetry_option + " unsupported.")
                    exit(1)
    return IDX_OBJ_TMP

def declarerank4(objname, symmetry_option, DIM=-1):
    if DIM==-1:
        DIM = par.parval_from_str("DIM")
    IDX_OBJ_TMP = [[[[sp.sympify(objname + str(i) + str(j) + str(k) + str(l)) for l in range(DIM)] for k in range(DIM)] for j in range(DIM)] for i in range(DIM)]
    for i in range(DIM):
        for j in range(DIM):
            for k in range(DIM):
                for l in range(DIM):
                    IDX_OBJ_TMP[i][j][k][l] = sp.sympify(objname + str(i) + str(j) + str(k) + str(l))
                    if symmetry_option == "sym12" or symmetry_option == "sym12_sym34":
                        if(j < i):
                            IDX_OBJ_TMP[i][j][k][l] = IDX_OBJ_TMP[j][i][k][l]
                    if symmetry_option == "sym23":
                        if(k < j):
                            IDX_OBJ_TMP[i][j][k][l] = IDX_OBJ_TMP[i][k][j][l]
                    if symmetry_option == "sym34" or symmetry_option == "sym12_sym34":
                        if(l < k):
                            IDX_OBJ_TMP[i][j][k][l] = IDX_OBJ_TMP[i][j][l][k]
                    if not (symmetry_option=="sym12" or symmetry_option=="sym34" or symmetry_option=="sym12_sym34" or symmetry_option=="none"):
                        print("Error: symmetry option "+symmetry_option+" unsupported.")
                        exit(1)
    return IDX_OBJ_TMP


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
def generic_matrix_inverter3x3(in2Darray, outINV):
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