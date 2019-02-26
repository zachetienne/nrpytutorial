# First we import needed core NRPy+ modules
from outputC import *
import NRPy_param_funcs as par
import grid as gri
import loop as lp
import indexedexp as ixp
import finite_difference as fin
import reference_metric as rfm

def Set_up_CurviBoundaryConditions():
    # Step 2.A.2. Output gridfunction #define aliases to file:
    evolved_variables_list, auxiliary_variables_list = \
        gri.output__gridfunction_defines_h__return_gf_lists("CurviBoundaryConditions")

    # Step 2.B.1: set the parity conditions on each gridfunction in gf_list,
    #       based on the digits at the end of its name
    # For example, if the gridfunction name ends with "01", then (based on the table
    #       in the Curvi. BCs tutorial module) the "set_parity_types()" function below
    #       will set the parity_type of that gridfunction to 5. We can be assured this
    #       is a robust algorithm, because gri.register_gridfunctions() in grid.py will
    #       throw an error if a gridfunction base name ends in an integer.
    #
    # After each parity type is found, we store the parity type of each gridfunction to
    #       const int8_t arrays evol_gf_parity and aux_gf_parity, appended to the end of
    #       "CurviBoundaryConditions/gridfunction_defines.h".
    def set_parity_types(gf_list):
        parity_type = []
        for i in range(len(gf_list)):
            varname = gf_list[i]
            parity_type__orig_len = len(parity_type)
            if len(varname) > 2:
                if varname[-2] == "0" and varname[-1] == "0":  # In Python, a[-1] points to the last
                    # element of a list; a[-2] the
                    # second-to-last element, etc.
                    parity_type.append(4)
                elif varname[-2] == "0" and varname[-1] == "1":
                    parity_type.append(5)
                elif varname[-2] == "0" and varname[-1] == "2":
                    parity_type.append(6)
                elif varname[-2] == "1" and varname[-1] == "1":
                    parity_type.append(7)
                elif varname[-2] == "1" and varname[-1] == "2":
                    parity_type.append(8)
                elif varname[-2] == "2" and varname[-1] == "2":
                    parity_type.append(9)
            if len(varname) > 1 and len(parity_type) == parity_type__orig_len:
                if varname[-1] == "0":
                    parity_type.append(1)
                elif varname[-1] == "1":
                    parity_type.append(2)
                elif varname[-1] == "2":
                    parity_type.append(3)
            if varname[len(varname) - 1].isdigit() == False:
                parity_type.append(0)

            if len(parity_type) == parity_type__orig_len:
                print("Error: Could not figure out parity type for evolved variable: " + varname)
                exit(1)
        return parity_type

    evol_parity_type = set_parity_types(evolved_variables_list)
    aux_parity_type = set_parity_types(auxiliary_variables_list)

    with open("CurviBoundaryConditions/gridfunction_defines.h", "a") as file:
        file.write("\n\n/* PARITY TYPES FOR ALL GRIDFUNCTIONS.\n")
        file.write("   SEE \"Tutorial-Start_to_Finish-BSSNCurvilinear-Two_BHs_Collide.ipynb\" FOR DEFINITIONS. */\n")
        if len(evolved_variables_list) > 0:
            file.write("const int8_t evol_gf_parity[" + str(len(evolved_variables_list)) + "] = { ")
            for i in range(len(evolved_variables_list) - 1):
                file.write(str(evol_parity_type[i]) + ", ")
            file.write(str(evol_parity_type[len(evolved_variables_list) - 1]) + " };\n")

        if len(auxiliary_variables_list) > 0:
            file.write("const int8_t aux_gf_parity[" + str(len(auxiliary_variables_list)) + "] = { ")
            for i in range(len(auxiliary_variables_list) - 1):
                file.write(str(aux_parity_type[i]) + ", ")
            file.write(str(aux_parity_type[len(auxiliary_variables_list) - 1]) + " };\n")
    print("Wrote to file \"CurviBoundaryConditions/gridfunction_defines.h\"")

    # Step 2.B.2: Set up unit-vector dot products (=parity) for each of the 10 parity condition types
    parity = ixp.zerorank1(DIM=10)
    UnitVectors_inner = ixp.zerorank2()
    xx0_inbounds, xx1_inbounds, xx2_inbounds = sp.symbols("xx0_inbounds xx1_inbounds xx2_inbounds", real=True)
    for i in range(3):
        for j in range(3):
            UnitVectors_inner[i][j] = rfm.UnitVectors[i][j].subs(rfm.xx[0], xx0_inbounds).subs(rfm.xx[1],
                                                                                               xx1_inbounds).subs(
                rfm.xx[2], xx2_inbounds)
    # Type 0: scalar
    parity[0] = sp.sympify(1)
    # Type 1: i0-direction vector or one-form
    # Type 2: i1-direction vector or one-form
    # Type 3: i2-direction vector or one-form
    for i in range(3):
        for Type in range(1, 4):
            parity[Type] += rfm.UnitVectors[Type - 1][i] * UnitVectors_inner[Type - 1][i]
    # Type 4: i0i0-direction rank-2 tensor
    # parity[4] = parity[1]*parity[1]
    # Type 5: i0i1-direction rank-2 tensor
    # Type 6: i0i2-direction rank-2 tensor
    # Type 7: i1i1-direction rank-2 tensor
    # Type 8: i1i2-direction rank-2 tensor
    # Type 9: i2i2-direction rank-2 tensor
    count = 4
    for i in range(3):
        for j in range(i, 3):
            parity[count] = parity[i + 1] * parity[j + 1]
            count = count + 1

    lhs_strings = []
    for i in range(10):
        lhs_strings.append("parity[" + str(i) + "]")
    outputC(parity, lhs_strings, "CurviBoundaryConditions/set_parity_conditions.h")

    # Step 2.C: Modified version of the scalar wave in curvilinear coordinates boundary
    #           condition ghost zone mapping routine, so that it also defines the parity
    #           conditions.

    # First output code needed for mapping from any given curvilinear coordinate gridpoint
    #  to the Cartesian coordinate in the grid interior (xxCart), and then find the
    #  corresponding gridpoint index in the grid interior (Cart_to_xx; xxminmax).
    # Generic coordinate NRPy+ file output, Part 1: output the conversion from (x0,x1,x2) to Cartesian (x,y,z)
    outputC([rfm.xxCart[0], rfm.xxCart[1], rfm.xxCart[2]], ["xCart[0]", "xCart[1]", "xCart[2]"],
            "CurviBoundaryConditions/xxCart.h")
    # Generic coordinate NRPy+ file output, Part 2: output the coordinate bounds xxmin[] and xxmax[]:
    with open("CurviBoundaryConditions/xxminmax.h", "w") as file:
        file.write(
            "const REAL xxmin[3] = {" + str(rfm.xxmin[0]) + "," + str(rfm.xxmin[1]) + "," + str(rfm.xxmin[2]) + "};\n")
        file.write(
            "const REAL xxmax[3] = {" + str(rfm.xxmax[0]) + "," + str(rfm.xxmax[1]) + "," + str(rfm.xxmax[2]) + "};\n")
    print("Wrote to file \"CurviBoundaryConditions/xxminmax.h\"")

    # Generic coordinate NRPy+ file output, Part 3: output the conversion from Cartesian (x,y,z) to interior/OB (x0,x1,x2)
    outputC([rfm.Cart_to_xx[0], rfm.Cart_to_xx[1], rfm.Cart_to_xx[2]],
            ["Cart_to_xx0_inbounds", "Cart_to_xx1_inbounds", "Cart_to_xx2_inbounds"],
            "CurviBoundaryConditions/Cart_to_xx.h")