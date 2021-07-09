# grid.py: functions & parameters related to numerical grids
# functions: Automatic loop output, output C code needed for gridfunction memory I/O, gridfunction registration

# Author: Zachariah B. Etienne
#         zachetie **at** gmail **dot* com

import NRPy_param_funcs as par      # NRPy+: Parameter interface
import sympy as sp                  # Import SymPy, a computer algebra system written entirely in Python
from collections import namedtuple  # Standard Python `collections` module: defines named tuples data structure
import os                           # Standard Python module for multiplatform OS-level functions

# Initialize globals related to the grid
glb_gridfcs_list = []
glb_gridfc = namedtuple('gridfunction', 'gftype name rank DIM')

thismodule = __name__
par.initialize_param(par.glb_param("char", thismodule, "GridFuncMemAccess", "SENRlike"))
par.initialize_param(par.glb_param("char", thismodule, "MemAllocStyle", "210"))
par.initialize_param(par.glb_param("int",  thismodule, "DIM", 3))

Nxx = par.Cparameters("int", thismodule,["Nxx0","Nxx1","Nxx2"],[64,32,64]) # Default to 64x32x64 grid
Nxx_plus_2NGHOSTS = par.Cparameters("int", thismodule,
                      ["Nxx_plus_2NGHOSTS0","Nxx_plus_2NGHOSTS1","Nxx_plus_2NGHOSTS2"],
                      [                  70,                  38,                  70]) # Default to 64x32x64 grid w/ NGHOSTS=3
xx = par.Cparameters("REAL", thismodule, ["xx0", "xx1", "xx2"],1e300) # These are C variables, not parameters, and
                                                                      # will be overwritten; best to initialize to crazy
                                                                      # number to ensure they are overwritten!
dxx = par.Cparameters("REAL", thismodule, ["dxx0", "dxx1", "dxx2"], 0.1)
invdx = par.Cparameters("REAL", thismodule, ["invdx0", "invdx1", "invdx2"], 1.0)

def variable_type(var):
    var_is_gf = False
    for gf in range(len(glb_gridfcs_list)):
        if str(var) == glb_gridfcs_list[gf].name:
            var_is_gf = True
    var_is_parameter = False
#     for paramname in range(len(par.glb_params_list)):
#         if str(var) == par.glb_params_list[paramname].parname:
#             var_is_parameter = True
    for paramname in range(len(par.glb_Cparams_list)):
        if str(var) == par.glb_Cparams_list[paramname].parname:
            var_is_parameter = True
    if var_is_parameter and var_is_gf:
        print("Error: variable "+str(var)+" is registered both as a gridfunction and as a Cparameter.")
        sys.exit(1)
    if not (var_is_parameter or var_is_gf):
        return "other"
    if var_is_parameter:
        return "Cparameter"
    if var_is_gf:
        return "gridfunction"
    print("grid.py: Could not find variable_type.")
    sys.exit(1)

def find_gftype(varname):
    for gf in glb_gridfcs_list:
        if gf.name == varname:
            return gf.gftype
    print("grid.py: Could not find gftype.")
    sys.exit(1)

def gfaccess(gfarrayname = "", varname = "", ijklstring = ""):
    found_registered_gf = False
    for gf in glb_gridfcs_list:
        if gf.name == varname:
            if found_registered_gf:
                print("Error: found duplicate gridfunction name: "+gf.name)
                sys.exit(1)
            found_registered_gf = True

    if not found_registered_gf:
        print("Error: gridfunction \""+varname+"\" is not registered!")
        sys.exit(1)

    gftype = find_gftype(varname)

    DIM = par.parval_from_str("DIM")
    retstring = ""
    if par.parval_from_str("GridFuncMemAccess") == "SENRlike":
        if gfarrayname == "":
            print("Error: GridFuncMemAccess = SENRlike requires gfarrayname be passed to gfaccess()")
            sys.exit(1)
        # FIXME: if gftype == "AUX" then override gfarrayname to aux_gfs[].
        #        This enables expressions containing a mixture of AUX and EVOL
        #        gridfunctions, though in a slightly hacky way.
        if gftype == "AUX":
            gfarrayname = "aux_gfs"
        elif gftype == "AUXEVOL":
            gfarrayname = "auxevol_gfs"
        # Return gfarrayname[IDX3(varname,i0)] for DIM=1, gfarrayname[IDX3(varname,i0,i1)] for DIM=2, etc.
        retstring += gfarrayname + "[IDX" + str(DIM+1) + "(" + varname.upper()+"GF" + ", "
    elif par.parval_from_str("GridFuncMemAccess") == "ETK":
        # Return varname[CCTK_GFINDEX3D(i0,i1,i2)] for DIM=3. Error otherwise
        if DIM != 3:
            print("Error: GridFuncMemAccess = ETK currently requires that gridfunctions be 3D. Can be easily extended.")
            sys.exit(1)
        if gfarrayname == "rhs_gfs":
            retstring += varname + "_rhsGF" + "[CCTK_GFINDEX"+str(DIM)+"D(cctkGH, "
        else:
            retstring += varname + "GF" + "[CCTK_GFINDEX"+str(DIM)+"D(cctkGH, "
    else:
        print("grid::GridFuncMemAccess = "+par.parval_from_str("GridFuncMemAccess")+" not supported")
        sys.exit(1)
    if ijklstring == "":
        for i in range(DIM):
            retstring += "i"+str(i)
            if i != DIM-1:
                retstring += ', '
    else:
        retstring += ijklstring
    return retstring + ")]"

# Gridfunction basenames cannot end in integers.
#    For example, rank-1 gridfunction vecU1 has
#    basename "vecU". Thus this gridfunction has
#    a valid basename. If we instead defined a
#    scalar gridfunction "u2", this would have a
#    basename of "u2" -- not a valid basename.
#    Being so strict enables us to determine
#    quickly what a gridfunction is by its name
#    alone, which is useful, e.g., when setting
#    up parity boundary conditions.
def verify_gridfunction_basename_is_valid(gf_basename):
    # First check for zero-length basenames:
    if len(gf_basename) == 0:
        print("Error: tried to register gridfunction without a name!")
        sys.exit(1)

    # https://stackoverflow.com/questions/1303243/how-to-find-out-if-a-python-object-is-a-string
    if sys.version_info[0] < 3:
        if not isinstance(gf_basename,basestring):
            print("ERROR: gf_names must be strings")
            sys.exit(1)
    else:
        if not isinstance(gf_basename, str):
            print("ERROR: gf_names must be strings")
            sys.exit(1)

    if len(gf_basename) > 0 and gf_basename[-1].isdigit():
        print("Error: tried to register gridfunction with base name: "+gf_basename)
        print(" Gridfunctions with base names ending in an integer is forbidden; pick a new name.")
        sys.exit(1)

import sys
def register_gridfunctions(gf_type,gf_names,rank=0,is_indexed=False,DIM=3):
    # Step 0: Sanity check
    if (rank > 0 and is_indexed == False) or (rank == 0 and is_indexed == True):
        print("Error: Attempted to register *indexed* gridfunction(s) with rank 0, or *scalar* gridfunctions with rank>0.")
        print("       Gridfunctions = ",gf_names)
        sys.exit(1)

    # Step 1: convert gf_names to a list if it's not already a list
    if type(gf_names) is not list:
        gf_namestmp = [gf_names]
        gf_names = gf_namestmp

    # Step 2: if the gridfunction is not indexed, then
    #         gf_names == base names. Check that the
    #         gridfunction basenames are valid:
    if is_indexed==False:
        for i in range(len(gf_names)):
            verify_gridfunction_basename_is_valid(gf_names[i])

    # Step 3: Verify that gridfunction type is valid.
    if gf_type not in ('EVOL', 'AUX', 'AUXEVOL'):
        print("Error in registering gridfunction(s) with unsupported type "+gf_type+".")
        print("Supported gridfunction types include:")
        print("    \"EVOL\": for evolved quantities (i.e., quantities stepped forward in time),")
        print("    \"AUXEVOL\": for auxiliary quantities needed at all points by evolved quantities,")
        print("    \"AUX\": for all other quantities needed at all gridpoints.")
        sys.exit(1)

    # Step 4: Check for duplicate grid function registrations. If:
    #         a) A duplicate is found, error out. Otherwise
    #         b) Append to list of gridfunctions, stored in glb_gridfcs_list[].
    for i in range(len(gf_names)):
        for j in range(len(glb_gridfcs_list)):
            if gf_names[i] == glb_gridfcs_list[j].name:
                print("Error: Tried to register the gridfunction \""+gf_names[i]+"\" twice (ignored type)\n\n")
                sys.exit(1)
        # If no duplicate found, append to "gridfunctions" list:
        glb_gridfcs_list.append(glb_gridfc(gf_type,gf_names[i],rank,DIM))

    # Step 5: Return SymPy object corresponding to symbol or
    #         list of symbols representing gridfunction in
    #         SymPy expression
    OBJ_TMPS = []
    for i in range(len(gf_names)):
        OBJ_TMPS.append(sp.symbols(gf_names[i], real=True))
#        OBJ_TMPS.append(sp.sympify(gf_names[i]))
    if len(gf_names)==1:
        return OBJ_TMPS[0]
    return OBJ_TMPS

# Given output directory "outdir" as input, the
#   following function outputs a file called
#   "outdir/gridfunction_defines.h", which
#   #define's all the gridfunction aliases, and
#   returns two lists, corresponding to the
#   names (strings) of the evolved and auxiliary
#   gridfunction names respectively.
#
# For example, if we define only two gridfunctions uu and vv,
#   which are evolved quantities (i.e., represent
#   the solution of the PDEs we are solving and are
#   registered with gftype == "EVOL"), then
#   this function will create a file with the following
#   content:
#
# | /* This file is automatically generated by NRPy+. Do not edit. */
# | /* EVOLVED VARIABLES: */
# | #define NUM_EVOL_GFS 2
# | #define UUGF 0
# | #define VVGF 1
# |
# | /* AUXILIARY VARIABLES: */
# | #define NUM_AUX_GFS 0
#
# The function will return two lists: the lists of
#    EVOL and AUX gridfunction names, respectively.
#    For this example, the first list (all gridfunctions
#    registered as EVOL) will be ["uu","vv"], and the
#    second (all gridfunctions registered as AUX) will
#    be the empty list: []

def gridfunction_lists():
    evolved_variables_list   = []
    auxiliary_variables_list = []
    auxevol_variables_list = []
    for i in range(len(glb_gridfcs_list)):
        if glb_gridfcs_list[i].gftype == "EVOL":
            evolved_variables_list.append(glb_gridfcs_list[i].name)
        if glb_gridfcs_list[i].gftype == "AUX":
            auxiliary_variables_list.append(glb_gridfcs_list[i].name)
        if glb_gridfcs_list[i].gftype == "AUXEVOL":
            auxevol_variables_list.append(glb_gridfcs_list[i].name)

    # Next we alphabetize the lists
    evolved_variables_list.sort()
    auxiliary_variables_list.sort()
    auxevol_variables_list.sort()

    return evolved_variables_list,auxiliary_variables_list,auxevol_variables_list

def output__gridfunction_defines_h__return_gf_lists(outdir):

    evolved_variables_list, auxiliary_variables_list, auxevol_variables_list = gridfunction_lists()

    # Finally we set up the #define statements:
    with open(os.path.join(outdir,"gridfunction_defines.h"), "w") as file:
        file.write("/* This file is automatically generated by NRPy+. Do not edit. */\n\n")

        file.write("/* EVOLVED VARIABLES: */\n")
        file.write("#define NUM_EVOL_GFS "+str(len(evolved_variables_list))+"\n")
        for i in range(len(evolved_variables_list)):
            file.write("#define "+evolved_variables_list[i].upper()+"GF\t"+str(i)+"\n")

        file.write("\n\n /* AUXILIARY VARIABLES: */\n")
        file.write("#define NUM_AUX_GFS " + str(len(auxiliary_variables_list)) + "\n")
        for i in range(len(auxiliary_variables_list)):
            file.write("#define " + auxiliary_variables_list[i].upper() + "GF\t" + str(i) + "\n")

        file.write("\n\n /* AUXEVOL VARIABLES: */\n")
        file.write("#define NUM_AUXEVOL_GFS "+str(len(auxevol_variables_list))+"\n")
        for i in range(len(auxevol_variables_list)):
            file.write("#define "+auxevol_variables_list[i].upper()+"GF\t"+str(i)+"\n")

    return evolved_variables_list,auxiliary_variables_list,auxevol_variables_list
