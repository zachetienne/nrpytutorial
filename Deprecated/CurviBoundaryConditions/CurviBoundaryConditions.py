# This module provides functions for setting up Curvilinear boundary conditions,
#     as documented in Tutorial-Start_to_Finish-Curvilinear_BCs.ipynb

# Authors: Zachariah B. Etienne
#         zachetie **at** gmail **dot* com
#          Terrence Pierre Jacques

# First we import needed core NRPy+ modules
from outputC import outputC      # NRPy+: Core C code output module
import NRPy_param_funcs as par   # NRPy+: Parameter interface
import sympy as sp               # SymPy: The Python computer algebra package upon which NRPy+ depends
import grid as gri               # NRPy+: Functions having to do with numerical grids
import deprecated_grid as evil_gri               # NRPy+: Functions having to do with numerical grids
import indexedexp as ixp         # NRPy+: Symbolic indexed expression (e.g., tensors, vectors, etc.) support
import reference_metric as rfm   # NRPy+: Reference metric support
import cmdline_helper as cmd     # NRPy+: Multi-platform Python command-line interface
import shutil, os, sys           # Standard Python modules for multiplatform OS-level functions
import deprecated_reference_metric as evil_rfm  # NRPy+: PLEASE DON'T USE.

def Set_up_CurviBoundaryConditions(Ccodesdir,verbose=True,Cparamspath=os.path.join("../"),
                                   enable_copy_of_static_Ccodes=True, BoundaryCondition="QuadraticExtrapolation",
                                   path_prefix=""):
    # Step P0: Check that Ccodesdir is not the same as CurviBoundaryConditions/boundary_conditions,
    #          to prevent trusted versions of these C codes from becoming contaminated.
    if os.path.join(Ccodesdir) == os.path.join("CurviBoundaryConditions", "boundary_conditions"):
        print("Error: Tried to output boundary conditions C code into CurviBoundaryConditions/boundary_conditions,"
              "       which is not allowed, to prevent trusted versions of these C codes from becoming contaminated.")
        sys.exit(1)

    # Step P1: Create the C codes output directory & copy static CurviBC files
    #          from CurviBoundaryConditions/boundary_conditions to Ccodesdir/
    if enable_copy_of_static_Ccodes:
        cmd.mkdir(os.path.join(Ccodesdir))

        # Choosing boundary condition drivers with in NRPy+
        #  - current options are Quadratic Polynomial Extrapolation for any coordinate system,
        #    and the Sommerfeld boundary condition for only cartesian coordinates
        if   str(BoundaryCondition) == "QuadraticExtrapolation":
            for file in ["apply_bcs_curvilinear.h", "BCs_data_structs.h", "bcstruct_freemem.h", "CurviBC_include_Cfunctions.h",
                         "driver_bcstruct.h", "set_bcstruct.h", "set_up__bc_gz_map_and_parity_condns.h"]:
                shutil.copy(os.path.join(path_prefix,"CurviBoundaryConditions", "boundary_conditions", file),
                            os.path.join(Ccodesdir))

            with open(os.path.join(Ccodesdir,"CurviBC_include_Cfunctions.h"),"a") as file:
                file.write("\n#include \"apply_bcs_curvilinear.h\"")


        elif str(BoundaryCondition) == "Sommerfeld":
            for file in ["BCs_data_structs.h", "bcstruct_freemem.h", "CurviBC_include_Cfunctions.h",
                         "driver_bcstruct.h", "set_bcstruct.h", "set_up__bc_gz_map_and_parity_condns.h"]:
                shutil.copy(os.path.join(path_prefix,"CurviBoundaryConditions", "boundary_conditions", file),
                            os.path.join(Ccodesdir))

            with open(os.path.join(Ccodesdir,"CurviBC_include_Cfunctions.h"),"a") as file:
                file.write("\n#include \"apply_bcs_sommerfeld.h\"")


        elif str(BoundaryCondition) == "QuadraticExtrapolation&Sommerfeld":
            for file in ["apply_bcs_curvilinear.h", "BCs_data_structs.h", "bcstruct_freemem.h", "CurviBC_include_Cfunctions.h",
                         "driver_bcstruct.h", "set_bcstruct.h", "set_up__bc_gz_map_and_parity_condns.h"]:
                shutil.copy(os.path.join(path_prefix,"CurviBoundaryConditions", "boundary_conditions", file),
                            os.path.join(Ccodesdir))

            with open(os.path.join(Ccodesdir,"CurviBC_include_Cfunctions.h"),"a") as file:
                file.write("\n#include \"apply_bcs_sommerfeld.h\"" +
                           "\n#include \"apply_bcs_curvilinear.h\"")


        else:
            print("ERROR: Only Quadratic Polynomial Extrapolation (QuadraticExtrapolation) and Sommerfeld boundary conditions are currently supported\n")
            sys.exit(1)


    # Step P2: Output correct #include for set_Cparameters.h to
    #          Ccodesdir/boundary_conditions/RELATIVE_PATH__set_Cparameters.h
    with open(os.path.join(Ccodesdir, "RELATIVE_PATH__set_Cparameters.h"), "w") as file:
        file.write("#include \"" + Cparamspath + "/set_Cparameters.h\"\n") # #include's may include forward slashes for paths, even in Windows.

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
    outputC(parity, lhs_strings, os.path.join(Ccodesdir, "parity_conditions_symbolic_dot_products.h"))


    # Step 2.a: Generate Ccodesdir/gridfunction_defines.h file,
    #       containing human-readable gridfunction aliases
    evolved_variables_list, auxiliary_variables_list, auxevol_variables_list = evil_gri.output__gridfunction_defines_h__return_gf_lists(Ccodesdir)[0:3]

    # Step 2.b: set the parity conditions on all gridfunctions in gf_list,
    #       based on how many digits are at the end of their names
    def set_parity_types(list_of_gf_names):
        parity_type = []
        for name in list_of_gf_names:
            for gf in gri.glb_gridfcs_list:
                if gf.name == name:
                    parity_type__orig_len = len(parity_type)
                    if gf.DIM < 3 or gf.DIM > 4:
                        print("Error: Cannot currently specify parity conditions on gridfunctions with DIM<3 or >4.")
                        sys.exit(1)
                    if gf.rank == 0:
                        parity_type.append(0)
                    elif gf.rank == 1:
                        if gf.DIM == 3:
                            parity_type.append(int(gf.name[-1]) + 1)  # = 1 for e.g., beta^0; = 2 for e.g., beta^1, etc.
                        elif gf.DIM == 4:
                            parity_type.append(int(gf.name[-1]))  # = 0 for e.g., b4^0; = 1 for e.g., beta^1, etc.
                    elif gf.rank == 2:
                        if gf.DIM == 3:
                            # element of a list; a[-2] the
                            # second-to-last element, etc.
                            idx0 = gf.name[-2]
                            idx1 = gf.name[-1]
                            if idx0 == "0" and idx1 == "0":
                                parity_type.append(4)
                            elif (idx0 == "0" and idx1 == "1") or (idx0 == "1" and idx1 == "0"):
                                parity_type.append(5)
                            elif (idx0 == "0" and idx1 == "2") or (idx0 == "2" and idx1 == "0"):
                                parity_type.append(6)
                            elif idx0 == "1" and idx1 == "1":
                                parity_type.append(7)
                            elif (idx0 == "1" and idx1 == "2") or (idx0 == "2" and idx1 == "1"):
                                parity_type.append(8)
                            elif idx0 == "2" and idx1 == "2":
                                parity_type.append(9)
                        elif gf.DIM == 4:
                            idx0 = gf.name[-2]
                            idx1 = gf.name[-1]
                            # g4DD00 = g_{tt} : parity type = 0
                            # g4DD01 = g_{tx} : parity type = 1
                            # g4DD02 = g_{ty} : parity type = 2
                            # g4DD0a = g_{ta} : parity type = a
                            if idx0 == "0":
                                parity_type.append(int(idx1))
                            elif idx1 == "0":
                                parity_type.append(int(idx0))
                            if idx0 == "1" and idx1 == "1":
                                parity_type.append(4)
                            elif (idx0 == "1" and idx1 == "2") or (idx0 == "2" and idx1 == "1"):
                                parity_type.append(5)
                            elif (idx0 == "1" and idx1 == "3") or (idx0 == "3" and idx1 == "1"):
                                parity_type.append(6)
                            elif idx0 == "2" and idx1 == "2":
                                parity_type.append(7)
                            elif (idx0 == "2" and idx1 == "3") or (idx0 == "3" and idx1 == "2"):
                                parity_type.append(8)
                            elif idx0 == "3" and idx1 == "3":
                                parity_type.append(9)
                    if len(parity_type) == parity_type__orig_len:
                        print("Error: Could not figure out parity type for "+gf.gftype+" gridfunction: " + gf.name,gf.DIM,gf.name[-2],gf.name[-1],gf.rank)
                        sys.exit(1)
        if len(parity_type) != len(list_of_gf_names):
            print("Error: For some reason the length of the parity types list did not match the length of the gf list.")
            sys.exit(1)
        return parity_type

    evol_parity_type = set_parity_types(evolved_variables_list)
    aux_parity_type = set_parity_types(auxiliary_variables_list)
    auxevol_parity_type = set_parity_types(auxevol_variables_list)

    # Step 2.c: Output all gridfunctions to Ccodesdir+"/gridfunction_defines.h"
    # ... then append to the file the parity type for each gridfunction.
    with open(os.path.join(Ccodesdir, "gridfunction_defines.h"), "a") as file:
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
        import textwrap
        wrapper = textwrap.TextWrapper(initial_indent="",subsequent_indent="    ",width=75)
        def print_parity_list(gf_type, variable_names,parity_types):
            outstr = ""
            if len(variable_names) != 0:
                outstr += gf_type+" parity: ( "
                for i in range(len(variable_names)):
                    outstr += variable_names[i] + ":" + str(parity_types[i])
                    if i != len(variable_names)-1:
                        outstr += ", "
                outstr += " )"
            print(wrapper.fill(outstr))
        print_parity_list("Evolved"  ,evolved_variables_list  ,evol_parity_type)
        print_parity_list("Auxiliary",auxiliary_variables_list,aux_parity_type)
        print_parity_list("AuxEvol"  ,auxevol_variables_list  ,auxevol_parity_type)

    # Step 3: Find the Eigen-Coordinate and set up the Eigen-Coordinate's reference metric:
    CoordSystem_orig = par.parval_from_str("reference_metric::CoordSystem")
    par.set_parval_from_str("reference_metric::CoordSystem", rfm.get_EigenCoord())
    rfm.reference_metric()

    # Step 4: Output C code for the Eigen-Coordinate mapping from xx->Cartesian:
    evil_rfm.xx_to_Cart_h("EigenCoord_xx_to_Cart", os.path.join(Cparamspath,"set_Cparameters.h"), os.path.join(Ccodesdir, "EigenCoord_xx_to_Cart.h"))

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
            os.path.join(Ccodesdir, "EigenCoord_Cart_to_xx.h"))

    # Step 6: Restore reference_metric::CoordSystem back to the original CoordSystem
    par.set_parval_from_str("reference_metric::CoordSystem", CoordSystem_orig)
    rfm.reference_metric()

# Sommerfeld boundary condition class; generates Sommerfeld parameters and functions to be used in subsequent
# C code
# Author: Terrence Pierre Jacques
class sommerfeld_boundary_condition_class():
    """
    Class for generating C code to apply Sommerfeld boundary conditions
    """
    # class variables should be the resulting dicts
    # Set class variable default values
    # radial falloff power n = 3 has been found to yield the best results
    #  - see Tutorial-SommerfeldBoundaryCondition.ipynb Step 2 for details
    def __init__(self, fd_order=2, vars_at_inf_default = 0., vars_radial_falloff_power_default = 3., vars_speed_default = 1.):
        evolved_variables_list, _, _ = gri.gridfunction_lists()

        # set class finite differencing order
        self.fd_order = fd_order

        NRPy_FD_order = par.parval_from_str("finite_difference::FD_CENTDERIVS_ORDER")

        if NRPy_FD_order < fd_order:
            print("ERROR: The global central finite differencing order within NRPy+ must be greater than or equal to the Sommerfeld boundary condition's finite differencing order")
            sys.exit(1)

        # Define class dictionaries to store sommerfeld parameters for each EVOL gridfunction

        # EVOL gridfunction asymptotic value at infinity
        self.vars_at_infinity = {}

        # EVOL gridfunction wave speed at outer boundaries
        self.vars_speed = {}

        # EVOL gridfunction radial falloff power
        self.vars_radial_falloff_power = {}

        # Set default values for each specific EVOL gridfunction
        for gf in evolved_variables_list:
            self.vars_at_infinity[gf.upper() + 'GF'] = vars_at_inf_default
            self.vars_radial_falloff_power[gf.upper() + 'GF'] = vars_radial_falloff_power_default
            self.vars_speed[gf.upper() + 'GF'] = vars_speed_default

    def sommerfeld_params(self):
        # Write parameters to C file

        # Creating array for EVOL gridfunction values at infinity
        var_at_inf_string = "{"
        for _gf,val in self.vars_at_infinity.items():
            var_at_inf_string += str(val) + ", "
        var_at_inf_string = var_at_inf_string[:-2] + "};"

        # Creating array for EVOL gridfunction values of radial falloff power
        vars_radial_falloff_power_string = "{"
        for _gf,val in self.vars_radial_falloff_power.items():
            vars_radial_falloff_power_string += str(val) + ", "
        vars_radial_falloff_power_string = vars_radial_falloff_power_string[:-2] + "};"

        # Creating array for EVOL gridfunction values of wave speed at outer boundaries
        var_speed_string = "{"
        for _gf,val in self.vars_speed.items():
            var_speed_string += str(val) + ", "
        var_speed_string = var_speed_string[:-2] + "};"

        # Writing to values to sommerfeld_params.h file
        out_str = """
// Sommerfeld EVOL grid function parameters
const REAL evolgf_at_inf[NUM_EVOL_GFS] = """+var_at_inf_string+"""
const REAL evolgf_radial_falloff_power[NUM_EVOL_GFS] = """+vars_radial_falloff_power_string+"""
const REAL evolgf_speed[NUM_EVOL_GFS] = """+var_speed_string+"""
"""
        return out_str

    @staticmethod
    def dfdr_function(fd_order):
        # function to write c code to calculate dfdr term in Sommerfeld boundary condition

        # Read what # of dimensions being used
        DIM = par.parval_from_str("grid::DIM")

        # Set up the chosen reference metric from chosen coordinate system, set within NRPy+
        CoordSystem = par.parval_from_str("reference_metric::CoordSystem")
        rfm.reference_metric()

        # Simplifying the results make them easier to interpret.
        do_simplify = True
        if "Sinh" in CoordSystem:
            # Simplification takes too long on Sinh* coordinate systems
            do_simplify = False

        # Construct Jacobian matrix, output Jac_dUSph_dDrfmUD[i][j] = \partial x_{Sph}^i / \partial x^j:
        Jac_dUSph_dDrfmUD = ixp.zerorank2()
        for i in range(3):
            for j in range(3):
                Jac_dUSph_dDrfmUD[i][j] = sp.diff(rfm.xxSph[i],rfm.xx[j])

        # Invert Jacobian matrix, output to Jac_dUrfm_dDSphUD.
        Jac_dUrfm_dDSphUD, dummyDET = ixp.generic_matrix_inverter3x3(Jac_dUSph_dDrfmUD)

        # Jac_dUrfm_dDSphUD[i][0] stores \partial x^i / \partial r
        if do_simplify:
            for i in range(3):
                Jac_dUrfm_dDSphUD[i][0] = sp.simplify(Jac_dUrfm_dDSphUD[i][0])

        # Declare \partial_i f, which is actually computed later on
        fdD = ixp.declarerank1("fdD") # = [fdD0, fdD1, fdD2]
        contraction = sp.sympify(0)
        for i in range(3):
            contraction += fdD[i]*Jac_dUrfm_dDSphUD[i][0]

        if do_simplify:
            contraction = sp.simplify(contraction)

        r_str_and_contraction_str = outputC([rfm.xxSph[0],contraction],
                 ["*_r","*_partial_i_f"],filename="returnstring",params="includebraces=False")


        def gen_central_2oFD_stencil_str(intdirn):
            if intdirn == 0:
                return "(gfs[IDX4S(which_gf,i0+1,i1,i2)]-gfs[IDX4S(which_gf,i0-1,i1,i2)])*0.5"  # Does not include the 1/dx multiplication
            if intdirn == 1:
                return "(gfs[IDX4S(which_gf,i0,i1+1,i2)]-gfs[IDX4S(which_gf,i0,i1-1,i2)])*0.5"  # Does not include the 1/dy multiplication
            return "(gfs[IDX4S(which_gf,i0,i1,i2+1)]-gfs[IDX4S(which_gf,i0,i1,i2-1)])*0.5"  # Does not include the 1/dz multiplication

        def gen_central_4oFD_stencil_str(intdirn):
            if intdirn == 0:
                return """(-c2*gfs[IDX4S(which_gf,i0+2,i1,i2)]
                         +c1*gfs[IDX4S(which_gf,i0+1,i1,i2)]
                         -c1*gfs[IDX4S(which_gf,i0-1,i1,i2)]
                         +c2*gfs[IDX4S(which_gf,i0-2,i1,i2)])"""  # Does not include the 1/dx multiplication
            if intdirn == 1:
                return """(-c2*gfs[IDX4S(which_gf,i0,i1+2,i2)]
                         +c1*gfs[IDX4S(which_gf,i0,i1+1,i2)]
                         -c1*gfs[IDX4S(which_gf,i0,i1-1,i2)]
                         +c2*gfs[IDX4S(which_gf,i0,i1-2,i2)])"""  # Does not include the 1/dy multiplication
            return """(-c2*gfs[IDX4S(which_gf,i0,i1,i2+2)]
                         +c1*gfs[IDX4S(which_gf,i0,i1,i2+1)]
                         -c1*gfs[IDX4S(which_gf,i0,i1,i2-1)]
                         +c2*gfs[IDX4S(which_gf,i0,i1,i2-2)])"""  # Does not include the 1/dz multiplication

        def gen_central_6oFD_stencil_str(intdirn):
            if intdirn == 0:
                return """( c3*gfs[IDX4S(which_gf,i0+3,i1,i2)]
                         -c2*gfs[IDX4S(which_gf,i0+2,i1,i2)]
                         +c1*gfs[IDX4S(which_gf,i0+1,i1,i2)]
                         -c1*gfs[IDX4S(which_gf,i0-1,i1,i2)]
                         +c2*gfs[IDX4S(which_gf,i0-2,i1,i2)]
                         -c3*gfs[IDX4S(which_gf,i0-3,i1,i2)])"""  # Does not include the 1/dx multiplication

            if intdirn == 1:
                return """( c3*gfs[IDX4S(which_gf,i0,i1+3,i2)]
                         -c2*gfs[IDX4S(which_gf,i0,i1+2,i2)]
                         +c1*gfs[IDX4S(which_gf,i0,i1+1,i2)]
                         -c1*gfs[IDX4S(which_gf,i0,i1-1,i2)]
                         +c2*gfs[IDX4S(which_gf,i0,i1-2,i2)]
                         -c3*gfs[IDX4S(which_gf,i0,i1-3,i2)])"""  # Does not include the 1/dy multiplication

            return """( c3*gfs[IDX4S(which_gf,i0,i1,i2+3)]
                         -c2*gfs[IDX4S(which_gf,i0,i1,i2+2)]
                         +c1*gfs[IDX4S(which_gf,i0,i1,i2+1)]
                         -c1*gfs[IDX4S(which_gf,i0,i1,i2-1)]
                         +c2*gfs[IDX4S(which_gf,i0,i1,i2-2)]
                         -c3*gfs[IDX4S(which_gf,i0,i1,i2-3)])"""  # Does not include the 1/dz multiplication

        def gen_central_fd_stencil_str(intdirn, fd_order):
            if fd_order==2:
                return gen_central_2oFD_stencil_str(intdirn)
            if fd_order==4:
                return gen_central_4oFD_stencil_str(intdirn)
            return gen_central_6oFD_stencil_str(intdirn)

        def output_dfdx(intdirn, fd_order):
            dirn = str(intdirn)
            dirnp1 = str((intdirn+1)%3)  # if dirn='0', then we want this to be '1'; '1' then '2'; and '2' then '0'
            dirnp2 = str((intdirn+2)%3)  # if dirn='0', then we want this to be '2'; '1' then '0'; and '2' then '1'

            preface = """
// On a +x"""+dirn+""" or -x"""+dirn+""" face, do up/down winding as appropriate:
if(abs(FACEXi["""+dirn+"""])==1 || i"""+dirn+"""+NGHOSTS >= Nxx_plus_2NGHOSTS"""+dirn+""" || i"""+dirn+"""-NGHOSTS <= 0) {
    int8_t SHIFTSTENCIL"""+dirn+""" = FACEXi["""+dirn+"""];
    if(i"""+dirn+"""+NGHOSTS >= Nxx_plus_2NGHOSTS"""+dirn+""") SHIFTSTENCIL"""+dirn+""" = -1;
    if(i"""+dirn+"""-NGHOSTS <= 0)                  SHIFTSTENCIL"""+dirn+""" = +1;
    SHIFTSTENCIL"""+dirnp1+""" = 0;
    SHIFTSTENCIL"""+dirnp2+""" = 0;
"""
            if fd_order == 2:
                return preface + """

    fdD"""+dirn+"""
        = SHIFTSTENCIL"""+dirn+"""*(-1.5*gfs[IDX4S(which_gf,i0+0*SHIFTSTENCIL0,i1+0*SHIFTSTENCIL1,i2+0*SHIFTSTENCIL2)]
                         +2.*gfs[IDX4S(which_gf,i0+1*SHIFTSTENCIL0,i1+1*SHIFTSTENCIL1,i2+1*SHIFTSTENCIL2)]
                         -0.5*gfs[IDX4S(which_gf,i0+2*SHIFTSTENCIL0,i1+2*SHIFTSTENCIL1,i2+2*SHIFTSTENCIL2)]
                        )*invdx"""+dirn+""";

// Not on a +x"""+dirn+""" or -x"""+dirn+""" face, using centered difference:
} else {
    fdD"""+dirn+""" = """+gen_central_fd_stencil_str(intdirn, 2)+"""*invdx"""+dirn+""";
}
"""
            if fd_order == 4:
                return preface + """

    fdD"""+dirn+"""
        = SHIFTSTENCIL"""+dirn+"""*(u0*gfs[IDX4S(which_gf,i0+0*SHIFTSTENCIL0,i1+0*SHIFTSTENCIL1,i2+0*SHIFTSTENCIL2)]
                         +u1*gfs[IDX4S(which_gf,i0+1*SHIFTSTENCIL0,i1+1*SHIFTSTENCIL1,i2+1*SHIFTSTENCIL2)]
                         +u2*gfs[IDX4S(which_gf,i0+2*SHIFTSTENCIL0,i1+2*SHIFTSTENCIL1,i2+2*SHIFTSTENCIL2)]
                         +u3*gfs[IDX4S(which_gf,i0+3*SHIFTSTENCIL0,i1+3*SHIFTSTENCIL1,i2+3*SHIFTSTENCIL2)]
                         +u4*gfs[IDX4S(which_gf,i0+4*SHIFTSTENCIL0,i1+4*SHIFTSTENCIL1,i2+4*SHIFTSTENCIL2)]
                        )*invdx"""+dirn+""";

// Not on a +x"""+dirn+""" or -x"""+dirn+""" face, using centered difference:
} else {
    fdD"""+dirn+""" = """+gen_central_fd_stencil_str(intdirn, 4)+"""*invdx"""+dirn+""";
}
"""
            if fd_order == 6:
                return preface + """

    fdD"""+dirn+"""
        = SHIFTSTENCIL"""+dirn+"""*(u0*gfs[IDX4S(which_gf,i0+0*SHIFTSTENCIL0,i1+0*SHIFTSTENCIL1,i2+0*SHIFTSTENCIL2)]
                         +u1*gfs[IDX4S(which_gf,i0+1*SHIFTSTENCIL0,i1+1*SHIFTSTENCIL1,i2+1*SHIFTSTENCIL2)]
                         +u2*gfs[IDX4S(which_gf,i0+2*SHIFTSTENCIL0,i1+2*SHIFTSTENCIL1,i2+2*SHIFTSTENCIL2)]
                         +u3*gfs[IDX4S(which_gf,i0+3*SHIFTSTENCIL0,i1+3*SHIFTSTENCIL1,i2+3*SHIFTSTENCIL2)]
                         +u4*gfs[IDX4S(which_gf,i0+4*SHIFTSTENCIL0,i1+4*SHIFTSTENCIL1,i2+4*SHIFTSTENCIL2)]
                         +u5*gfs[IDX4S(which_gf,i0+5*SHIFTSTENCIL0,i1+5*SHIFTSTENCIL1,i2+5*SHIFTSTENCIL2)]
                         +u6*gfs[IDX4S(which_gf,i0+6*SHIFTSTENCIL0,i1+6*SHIFTSTENCIL1,i2+6*SHIFTSTENCIL2)]
                        )*invdx"""+dirn+""";

// Not on a +x"""+dirn+""" or -x"""+dirn+""" face, using centered difference:
} else {
    fdD"""+dirn+""" = """+gen_central_fd_stencil_str(intdirn, 6)+"""*invdx"""+dirn+""";
}
"""

            print("Error: fd_order = "+str(fd_order)+" currently unsupported.")
            sys.exit(1)

        contraction_term_func = """

// Function to calculate the radial derivative of a grid function
void contraction_term(const paramstruct *restrict params, const int which_gf, const REAL *restrict gfs, REAL *restrict xx[3],
           const int8_t FACEXi[3], const int i0, const int i1, const int i2, REAL *restrict _r, REAL *restrict _partial_i_f) {

#include "RELATIVE_PATH__set_Cparameters.h" /* Header file containing correct #include for set_Cparameters.h;
                                             * accounting for the relative path */

// Initialize derivatives to crazy values, to ensure that
//   we will notice in case they aren't set properly.
REAL fdD0=1e100;
REAL fdD1=1e100;
REAL fdD2=1e100;

REAL xx0 = xx[0][i0];
REAL xx1 = xx[1][i1];
REAL xx2 = xx[2][i2];

int8_t SHIFTSTENCIL0;
int8_t SHIFTSTENCIL1;
int8_t SHIFTSTENCIL2;

"""
        if fd_order == 4:
            contraction_term_func +="""
// forward/backward finite difference coefficients
const REAL u0 =-25./12.;
const REAL u1 = 4.;
const REAL u2 = -3.;
const REAL u3 = 4./3.;
const REAL u4 = -1./4.;

// central finite difference coefficients
const REAL c1 = 2./3.;
const REAL c2 = 1./12.;

"""
        if fd_order == 6:
            contraction_term_func +="""
// forward/backward finite difference coefficients
const REAL u0 = -49./20.;
const REAL u1 =  6.;
const REAL u2 = -15./2.;
const REAL u3 =  20./3.;
const REAL u4 = -15./4.;
const REAL u5 =  6./5.;
const REAL u6 = -1./6.;

// central finite difference coefficients
const REAL c1 = 3./4.;
const REAL c2 = 3./20.;
const REAL c3 = 1./60;

"""
        for i in range(DIM):
            if "fdD"+str(i) in r_str_and_contraction_str:
                contraction_term_func += output_dfdx(i, fd_order)

        contraction_term_func += "\n" + r_str_and_contraction_str

        contraction_term_func +="""
} // END contraction_term function
"""
        return contraction_term_func

    def write_sommerfeld_main_Ccode(self, Ccodesdir):
        main_Ccode = """
// Boundary condition driver routine: Apply BCs to all
// boundary faces of the 3D numerical domain, filling in the
// outer boundary ghost zone layers, starting with the innermost
// layer and working outward.
"""
        main_Ccode += self.sommerfeld_params()
        main_Ccode += self.dfdr_function(self.fd_order)

        main_Ccode += """
void apply_bcs_sommerfeld(const paramstruct *restrict params, REAL *restrict xx[3],
                          const bc_struct *restrict bcstruct, const int NUM_GFS,
                          const int8_t *restrict gfs_parity, REAL *restrict gfs,
                          REAL *restrict rhs_gfs) {

    #pragma omp parallel for
        for(int which_gf=0;which_gf<NUM_GFS;which_gf++) {
          const REAL char_speed             = evolgf_speed[which_gf];
          const REAL var_at_infinity        = evolgf_at_inf[which_gf];
          const REAL radial_falloff_power = evolgf_radial_falloff_power[which_gf];


          #include "RELATIVE_PATH__set_Cparameters.h" /* Header file containing correct #include for set_Cparameters.h;
                                                       * accounting for the relative path */


            for(int which_gz = 0; which_gz < NGHOSTS; which_gz++) {
                for(int pt=0;pt<bcstruct->num_ob_gz_pts[which_gz];pt++) {
                    const int i0 = bcstruct->outer[which_gz][pt].outer_bc_dest_pt.i0;
                    const int i1 = bcstruct->outer[which_gz][pt].outer_bc_dest_pt.i1;
                    const int i2 = bcstruct->outer[which_gz][pt].outer_bc_dest_pt.i2;
                    const int8_t FACEX0 = bcstruct->outer[which_gz][pt].FACEi0;
                    const int8_t FACEX1 = bcstruct->outer[which_gz][pt].FACEi1;
                    const int8_t FACEX2 = bcstruct->outer[which_gz][pt].FACEi2;

                    const int8_t FACEXi[3] = {FACEX0, FACEX1, FACEX2};

                    // Initialize derivatives to crazy values, to ensure that
                    //   we will notice in case they aren't set properly.
                    REAL r = 1e100;
                    REAL partial_i_f = 1e100;

                    contraction_term(params, which_gf, gfs, xx, FACEXi, i0, i1, i2, &r, &partial_i_f);

                    const REAL invr = 1./r;

                    const REAL source_rhs = -char_speed*(partial_i_f + invr*(gfs[IDX4S(which_gf,i0,i1,i2)] - var_at_infinity));
                    rhs_gfs[IDX4S(which_gf,i0,i1,i2)] = source_rhs;

                    /************* For radial falloff and the extrapolated k term *************/
                    if (radial_falloff_power > 0) {

                      // Move one point away from gz point to compare pure advection to df/dt|interior

                      const int i0_offset = i0+FACEX0;
                      const int i1_offset = i1+FACEX1;
                      const int i2_offset = i2+FACEX2;

                      // Initialize derivatives to crazy values, to ensure that
                      //   we will notice in case they aren't set properly.
                      REAL r_offset = 1e100;
                      REAL partial_i_f_offset = 1e100;

                      contraction_term(params, which_gf, gfs, xx, FACEXi, i0_offset, i1_offset, i2_offset, &r_offset, &partial_i_f_offset);

                      const REAL invr_offset = 1./r_offset;

                      // Pure advection: [FIXME: Add equation (appearing in Jupyter notebook documentation)]
                      const REAL extrap_rhs = char_speed*(partial_i_f_offset + invr_offset*(gfs[IDX4S(which_gf,i0_offset,i1_offset,i2_offset)] - var_at_infinity));

                      // Take difference between pure advection and df/dt|interior
                      const REAL diff_between_advection_and_f_rhs =
                          rhs_gfs[IDX4S(which_gf,i0_offset,i1_offset,i2_offset)] + extrap_rhs;

                      // Solve for k/(r_gz)^n+1 term
                      rhs_gfs[IDX4S(which_gf,i0,i1,i2)] += diff_between_advection_and_f_rhs*pow(r_offset*invr,radial_falloff_power);

                  }
                } // END for(int pt=0;pt<num_ob_gz_pts[which_gz];pt++)

            // Apply INNER (parity) boundary conditions:
                for(int pt=0;pt<bcstruct->num_ib_gz_pts[which_gz];pt++) {
                    const int i0dest = bcstruct->inner[which_gz][pt].inner_bc_dest_pt.i0;
                    const int i1dest = bcstruct->inner[which_gz][pt].inner_bc_dest_pt.i1;
                    const int i2dest = bcstruct->inner[which_gz][pt].inner_bc_dest_pt.i2;
                    const int i0src  = bcstruct->inner[which_gz][pt].inner_bc_src_pt.i0;
                    const int i1src  = bcstruct->inner[which_gz][pt].inner_bc_src_pt.i1;
                    const int i2src  = bcstruct->inner[which_gz][pt].inner_bc_src_pt.i2;

                    rhs_gfs[IDX4S(which_gf,i0dest,i1dest,i2dest)] =
                            bcstruct->inner[which_gz][pt].parity[gfs_parity[which_gf]] * rhs_gfs[IDX4S(which_gf, i0src,i1src,i2src)];
            } // END for(int pt=0;pt<num_ib_gz_pts[which_gz];pt++)
        } // END for(int which_gz = 0; which_gz < NGHOSTS; which_gz++)
    } // END for(int which_gf=0;which_gf<NUM_GFS;which_gf++)
} // END function
"""
        with open(os.path.join(Ccodesdir,"boundary_conditions/apply_bcs_sommerfeld.h"),"w") as file:
            file.write(main_Ccode)

    def write_sommerfeld_file(self, Ccodesdir):
        self.write_sommerfeld_main_Ccode(Ccodesdir)
        print("""\nSuccessfully generated Sommerfeld boundary condition C code""")
