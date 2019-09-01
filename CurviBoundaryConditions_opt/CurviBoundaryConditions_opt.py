# This module provides functions for setting up Curvilinear boundary conditions,
#     as documented in Tutorial-Start_to_Finish-Curvilinear_BCs.ipynb

# Author: Zachariah B. Etienne
#         zachetie **at** gmail **dot* com

# First we import needed core NRPy+ modules
from outputC import *
import NRPy_param_funcs as par
import grid as gri
import loop as lp
import indexedexp as ixp
import finite_difference as fin
import reference_metric as rfm
import sys

def Set_up_CurviBoundaryConditions(outdir="CurviBoundaryConditions/",verbose=True):
    # Step 0: Set up reference metric in case it hasn't already been set up.
    #         (Doing it twice hurts nothing).
    rfm.reference_metric()

    # Step 1: Set unit-vector dot products (=parity) for each of the 10 parity condition types
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
    outputC(parity, lhs_strings, outdir + "parity_conditions_symbolic_dot_products.h")


    # Step 2.a: Generate outdir+gridfunction_defines.h file,
    #       containing human-readable gridfunction aliases
    evolved_variables_list, auxiliary_variables_list, auxevol_variables_list = gri.output__gridfunction_defines_h__return_gf_lists(outdir)

    # Step 2.b: set the parity conditions on all gridfunctions in gf_list,
    #       based on how many digits are at the end of their names
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
                sys.exit(1)
        return parity_type

    evol_parity_type = set_parity_types(evolved_variables_list)
    aux_parity_type = set_parity_types(auxiliary_variables_list)
    auxevol_parity_type = set_parity_types(auxevol_variables_list)

    # Step 2.c: Output all gridfunctions to outdir+"/gridfunction_defines.h"
    # ... then append to the file the parity type for each gridfunction.
    with open(outdir + "/gridfunction_defines.h", "a") as file:
        file.write("\n\n/* PARITY TYPES FOR ALL GRIDFUNCTIONS.\n")
        file.write(
            "   SEE \"Tutorial-Start_to_Finish-Curvilinear_BCs.ipynb\" FOR DEFINITIONS. */\n")
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

        if len(auxevol_variables_list) > 0:
            file.write("const int8_t auxevol_gf_parity[" + str(len(auxevol_variables_list)) + "] = { ")
            for i in range(len(auxevol_variables_list) - 1):
                file.write(str(auxevol_parity_type[i]) + ", ")
            file.write(str(auxevol_parity_type[len(auxevol_variables_list) - 1]) + " };\n")

    if verbose == True:
        for i in range(len(evolved_variables_list)):
            print("Evolved gridfunction \"" + evolved_variables_list[i] + "\" has parity type " + str(
                evol_parity_type[i]) + ".")
        for i in range(len(auxiliary_variables_list)):
            print("Auxiliary gridfunction \"" + auxiliary_variables_list[i] + "\" has parity type " + str(
                aux_parity_type[i]) + ".")
        for i in range(len(auxevol_variables_list)):
            print("AuxEvol gridfunction \"" + auxevol_variables_list[i] + "\" has parity type " + str(
                auxevol_parity_type[i]) + ".")

    # Step 3: Find the Eigen-Coordinate and set up the Eigen-Coordinate's reference metric:
    CoordSystem_orig = par.parval_from_str("reference_metric::CoordSystem")
    par.set_parval_from_str("reference_metric::CoordSystem", rfm.get_EigenCoord())
    rfm.reference_metric()

    # Step 4: Output C code for the Eigen-Coordinate mapping from xx->Cartesian:
    rfm.xxCart_h("EigenCoord_xxCart", "../set_Cparameters.h", outdir + "EigenCoord_xxCart.h")

    # Step 5: Output the Eigen-Coordinate mapping from Cartesian->xx:
    # Step 5.a: Sanity check: First make sure that rfm.Cart_to_xx has been set. Error out if not!
    if rfm.Cart_to_xx[0] == 0 or rfm.Cart_to_xx[1] == 0 or rfm.Cart_to_xx[2] == 0:
        print("ERROR: rfm.Cart_to_xx[], which maps Cartesian -> xx, has not been set for")
        print("       reference_metric::CoordSystem = " + par.parval_from_str("reference_metric::CoordSystem"))
        print("       Boundary conditions in curvilinear coordinates REQUIRE this be set.")
        sys.exit(1)
    # Step 5.b: Output C code for the Eigen-Coordinate mapping from Cartesian->xx:
    outputC([rfm.Cart_to_xx[0], rfm.Cart_to_xx[1], rfm.Cart_to_xx[2]],
            ["Cart_to_xx0_inbounds", "Cart_to_xx1_inbounds", "Cart_to_xx2_inbounds"],
            outdir + "EigenCoord_Cart_to_xx.h")

    # Step 6: Restore reference_metric::CoordSystem back to the original CoordSystem
    par.set_parval_from_str("reference_metric::CoordSystem", CoordSystem_orig)
    rfm.reference_metric()