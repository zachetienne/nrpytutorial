# grid.py: functions & parameters related to numerical grids
# functions: Automatic loop output, output C code needed for gridfunction memory I/O, gridfunction registration

# Author: Zachariah B. Etienne
#         zachetie **at** gmail **dot* com

import NRPy_param_funcs as par      # NRPy+: Parameter interface
import sympy as sp                  # Import SymPy, a computer algebra system written entirely in Python
from collections import namedtuple  # Standard Python `collections` module: defines named tuples data structure
import os                           # Standard Python module for multiplatform OS-level functions
from suffixes import setsuffix
import re
from fstr import f

# Initialize globals related to the grid
ET_driver = ""

# These next two are needed by Tutorial-Numerical_Grids.py
tile_tmp_variables_list = []
scalar_tmp_variables_list = []

glb_gridfc = namedtuple('gridfunction', 'gftype name rank DIM f_infinity wavespeed centering external_module')

glb_gridfcs_list = []
def glb_gridfcs_map():
    m = {}
    for gf in glb_gridfcs_list:
        m[gf.name] = gf
    assert len(glb_gridfcs_list) == len(m)
    return m

# Grids may have centerings. These will be C=cell-centered or V=vertex centered, with either one C or V per dimension
gf_centering = {}

# griddata_struct contains data needed by each grid
glb_griddata = namedtuple('griddata', 'module string')
glb_griddata_struct_list = []

thismodule = __name__
par.initialize_param(par.glb_param("char", thismodule, "GridFuncMemAccess", "SENRlike"))
par.initialize_param(par.glb_param("char", thismodule, "MemAllocStyle", "210"))
par.initialize_param(par.glb_param("int",  thismodule, "DIM", 3))

Nxx = par.Cparameters("int", thismodule, ["Nxx0", "Nxx1", "Nxx2"], [64, 32, 64])  # Default to 64x32x64 grid
Nxx_plus_2NGHOSTS = par.Cparameters("int", thismodule,
                      ["Nxx_plus_2NGHOSTS0","Nxx_plus_2NGHOSTS1","Nxx_plus_2NGHOSTS2"],
                      [                  70,                  38,                  70]) # Default to 64x32x64 grid w/ NGHOSTS=3
xx = par.Cparameters("REAL", thismodule, ["xx0", "xx1", "xx2"], 1e300) # These are C variables, not parameters, and
                                                                       # will be overwritten; best to initialize to crazy
                                                                       # number to ensure they are overwritten!
# TODO: dt = par.Cparameters("REAL", thismodule, ["dt"], 0.1)
dxx = par.Cparameters("REAL", thismodule, ["dxx0", "dxx1", "dxx2"], 0.1)
xxmin0, xxmin1, xxmin2 = par.Cparameters("REAL", thismodule, ["xxmin0", "xxmin1", "xxmin2"], 1.0)
xxmax0, xxmax1, xxmax2 = par.Cparameters("REAL", thismodule, ["xxmax0", "xxmax1", "xxmax2"], 2.0)
invdx = par.Cparameters("REAL", thismodule, ["invdx0", "invdx1", "invdx2"], 1.0)

# Origin of grid in Cartesian coordinates, relative to global grid
Cart_origin = par.Cparameters("REAL", thismodule, ["Cart_originx", "Cart_originy", "Cart_originz"], 0.0)
Cart_CoM_offset = par.Cparameters("REAL", thismodule, ["Cart_CoM_offsetx", "Cart_CoM_offsety", "Cart_CoM_offsetz"], 0.0)

def variable_type(var):
    var_data = glb_gridfcs_map().get(str(var),None)
    var_is_gf = var_data is not None
    var_is_parameter = False
#     for paramname in range(len(par.glb_params_list)):
#         if str(var) == par.glb_params_list[paramname].parname:
#             var_is_parameter = True
    for paramname in range(len(par.glb_Cparams_list)):
        if str(var) == par.glb_Cparams_list[paramname].parname:
            var_is_parameter = True
    if var_is_parameter and var_is_gf:
        raise Exception("Error: variable "+str(var)+" is registered both as a gridfunction and as a Cparameter.")
    elif not (var_is_parameter or var_is_gf):
        return "other"
    elif var_is_parameter:
        return "Cparameter"
    elif var_is_gf:
        return "gridfunction"
    else:
        raise Exception("grid.py: Could not find variable_type for '"+var+"'.")
    if not found_registered_gf:
        print("Error: gridfunction \""+varname+"\" is not registered!")
        print("Here's the list of registered gridfunctions:", grid.glb_gridfcs_list)
        sys.exit(1)

def find_gfnames():
    return sorted(list(glb_gridfcs_map().keys()))

def find_gftype(varname,die=True):
    assert type(varname) == str or type(varname) == unicode
    var_data = glb_gridfcs_map().get(varname, None)
    if var_data is None:
        if die:
            raise Exception("grid.py: Could not find gftype for '"+varname+"'.")
        else:
            return None
    else:
        return var_data.gftype

def find_gfmodule(varname,die=True):
    var_data = glb_gridfcs_map().get(varname, None)
    if var_data is None:
        if die:
            raise Exception("grid.py: Could not find module for '"+varname+"'.")
        else:
            return None
    else:
        return var_data.external_module

from var_access import var_from_access, set_access

def gfaccess(gfarrayname = "", varname = "", ijklstring = "", context = "DECL"):
    ret = _gfaccess(gfarrayname, varname, ijklstring, context)
    assert ret is not None, f("gfarrayname={gfarrayname}, varname={varname}, ijklstring={ijklstring}, context={context}")
    #from_access[ret] = varname
    set_access(ret, varname)
    return ret

def _gfaccess(gfarrayname, varname, ijklstring, context):
    var_data = glb_gridfcs_map().get(varname, None)

    assert context in ["DECL","USE"], "The context must be either DECL or USE, not '"+context+"'"

    if var_data is None:
        raise Exception("Error: gridfunction '"+varname+"' is not registered!")

    gftype = var_data.gftype

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
        elif gftype == "EXTERNAL":
            gfarrayname = "ext_gfs"
        elif gftype == "CORE":
            gfarrayname = "core_gfs"
        elif gftype == "TILE_TMP":
            gfarrayname = "tile_tmp_gfs"
        elif gftype == "SCALAR_TMP":
            gfarrayname = "scalar_tmp_gfs"
        elif gftype == "EVOL":
            pass
        # Return gfarrayname[IDX3(varname,i0)] for DIM=1, gfarrayname[IDX3(varname,i0,i1)] for DIM=2, etc.
        retstring += gfarrayname + "[IDX" + str(DIM+1) + "S(" + varname.upper()+"GF" + ", "
    elif par.parval_from_str("GridFuncMemAccess") == "ETK":
        # Return varname[CCTK_GFINDEX3D(i0,i1,i2)] for DIM=3. Error otherwise
        if DIM != 3:
            print("Error: GridFuncMemAccess = ETK currently requires that gridfunctions be 3D. Can be easily extended.")
            sys.exit(1)
        if ET_driver == "Carpet":
            if gfarrayname == "rhs_gfs":
                retstring += varname + "_rhsGF" + "[CCTK_GFINDEX"+str(DIM)+"D(cctkGH, "
            elif gftype == "EXTERNAL":
                retstring += varname + "[CCTK_GFINDEX"+str(DIM)+"D(cctkGH, "
            #elif gftype == "CORE":
            #    print("Error: gftype = CORE should only be used with the CarpetX driver.")
            #    sys.exit(1)
            elif gftype == "TILE_TMP":
                print("Error: gftype = TILE_TMP should only be used with the CarpetX driver.")
                sys.exit(1)
            elif gftype == "SCALAR_TMP":
                print("Error: gftype = SCALAR_TMP should only be used with the CarpetX driver.")
                sys.exit(1)
            else:
                retstring += varname + "GF" + "[CCTK_GFINDEX"+str(DIM)+"D(cctkGH, "
        elif ET_driver == "CarpetX":
            vartype = par.parval_from_str("PRECISION")
            mask = "mask, " if vartype == "CCTK_REALVEC" else ""
            if gfarrayname == "rhs_gfs":
                return retstring + varname + "_rhsGF" + "("+mask+find_centering(varname)+"_index)"
            elif gftype == "EXTERNAL":
                return retstring + varname + "("+mask+find_centering(varname)+"_index)"
            elif gftype == "CORE":
                if varname in ["x", "y", "z"]:
                    return retstring + "p." + varname
                elif varname == "regrid_error":
                    return retstring + varname + "(CCC_index)"
                else:
                    raise Exception("Unknown CORE variable: "+varname)
            elif gftype == "TILE_TMP":
                if ijklstring == "":
                    return retstring + varname + "("+mask+find_centering(varname)+"_tmp_index)"
                else:
                    return retstring + varname + "("+mask+find_centering(varname)+"_tmp_layout, p.I{ijklstring})"
            elif gftype == "SCALAR_TMP":
                if context == "DECL":
                    return retstring + "const " + vartype + " " + varname + " CCTK_ATTRIBUTE_UNUSED /* decl */"
                else:
                    return varname
            else:
                if ijklstring == "":
                    return retstring + varname + "GF("+mask+find_centering(varname)+"_index)"
                else:
                    return retstring + varname + "GF("+mask+find_centering(varname)+"_layout, p.I"+ijklstring+")"
        else:
            retstring += varname + "GF[CCTK_GFINDEX"+str(DIM)+"D(cctkGH, "
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
        if not isinstance(gf_basename, basestring):
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

def find_centering(gf_name):
    gf = glb_gridfcs_map().get(gf_name, None)
    if gf is not None:
        return gf.centering

import sys
def register_gridfunctions(gf_type,gf_names,rank=0,is_indexed=False,DIM=3, f_infinity=None, wavespeed=None, centering=None, external_module=None):
    global gf_centering

    # Step 0: Sanity check
    if (rank > 0 and not is_indexed) or (rank == 0 and is_indexed):
        print("Error: Attempted to register *indexed* gridfunction(s) with rank 0, or *scalar* gridfunctions with rank>0.")
        print("       Gridfunctions = ",gf_names)
        sys.exit(1)

    assert centering is None or type(centering) == str, \
        f("Centering, if supplied, should be of type str. Type was '{type(centering)}'")
    if centering is not None:
        assert len(centering) == DIM, f("len(centering) ({len(centering)}) is not equal to DIM ({DIM})")
        for c in centering:
            assert c in "CV", f("Centering should contain only 'C' or 'V'. The letter '{c}' was found.")

    # Step 1: convert gf_names to a list if it's not already a list
    if not isinstance(gf_names, list):
        gf_namestmp = [gf_names]
        gf_names = gf_namestmp

    for gf_name in gf_names:
        gf_centering[gf_name] = centering

    f_inf = []
    if f_infinity is None:
        f_infinity = 0.0  # set to default
    if not isinstance(f_infinity, list):
        for gf in gf_names:
            f_inf.append(f_infinity)
        f_infinity = f_inf


    wavespd = []
    if wavespeed is None:
        wavespeed = 1.0  # set to default
    if not isinstance(wavespeed, list):
        for gf in gf_names:
            wavespd.append(wavespeed)
        wavespeed = wavespd

    cent = []
    if centering is None:
        centering = "VVV"
    if not isinstance(centering, list):
        for gf in gf_names:
            cent.append(centering)
        centering = cent

    if len(f_infinity) != len(gf_names) or len(wavespeed) != len(gf_names) or len(centering) != len(gf_names):
        print("ERROR: Tried to register a list of gridfunctions with length " + str(len(gf_names))+" but f_infinity (len="+str(len(f_infinity))+") or wavespeed (len="+str(len(wavespeed))+") or centering (len="+str(len(centering))+") lists not of the same length.")
        sys.exit(1)

    # Step 2: if the gridfunction is not indexed, then
    #         gf_names == base names. Check that the
    #         gridfunction basenames are valid:
    if not is_indexed:
        for i in range(len(gf_names)):
            verify_gridfunction_basename_is_valid(gf_names[i])

    # Step 3: Verify that gridfunction type is valid.
    if gf_type not in ('EVOL', 'AUX', 'AUXEVOL', 'EXTERNAL', 'CORE', 'TILE_TMP', 'SCALAR_TMP'):
        print("Error in registering gridfunction(s) with unsupported type "+gf_type+".")
        print("Supported gridfunction types include:")
        print("    \"EVOL\": for evolved quantities (i.e., quantities stepped forward in time),")
        print("    \"AUXEVOL\": for auxiliary quantities needed at all points by evolved quantities,")
        print("    \"AUX\": for all other quantities needed at all gridpoints.")
        print("    \"EXTERNAL\": for all quantities defined in other modules.")
        print("    \"CORE\": for all quantities defined inside the CarpetX driver.")
        print("    \"TILE_TMP\": for all temporary quantities defined for CarpetX tiles.")
        print("    \"SCALAR_TMP\": for all temporary quantities defined for doubles.")
        #sys.exit(1)
        raise Exception("unsupported type")

    if gf_type == "EXTERNAL":
        for gf_name in gf_names:
            setsuffix(gf_name, "_ext")

    if gf_type == "CORE":
        for gf_name in gf_names:
            setsuffix(gf_name, "_core")

    if gf_type == "TILE_TMP":
        for gf_name in gf_names:
            setsuffix(gf_name, "_tile_tmp")

    if gf_type == "SCALAR_TMP":
        for gf_name in gf_names:
            setsuffix(gf_name, "")

    # Step 4: Check for duplicate grid function registrations. If:
    #         a) A duplicate is found, error out. Otherwise
    #         b) Add to map of gridfunctions, stored in glb_gridfcs_list
    for i in range(len(gf_names)):
        gf_name = gf_names[i]
        if gf_name in glb_gridfcs_map():
            assert gf_type == glb_gridfcs_map()[gf_name].gftype, \
                'Error: Tried to register the gridfunction "'+gf_names[i]+'" twice with different types'
        # If no duplicate found, append to "gridfunctions" list:
        var_data = (glb_gridfc(gf_type, gf_name, rank, DIM, f_infinity[i], wavespeed[i], centering[i], external_module))
        glb_gridfcs_map()[gf_name] = var_data
        if gf_name not in glb_gridfcs_map():
            glb_gridfcs_list.append(var_data)

    # Step 5: Return SymPy object corresponding to symbol or
    #         list of symbols representing gridfunction in
    #         SymPy expression
    OBJ_TMPS = []
    for i in range(len(gf_names)):
        OBJ_TMPS.append(sp.symbols(gf_names[i], real=True))
    if len(gf_names) == 1:
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
    external_variables_list = []
    core_variables_list = []
    tmp_variables_list = []
    for gf in glb_gridfcs_list:
        gf_type = gf.gftype
        gf_name = gf.name
        if gf_type == "EVOL":
            evolved_variables_list.append(gf_name)
        elif gf_type == "AUX":
            auxiliary_variables_list.append(gf_name)
        elif gf_type == "AUXEVOL":
            auxevol_variables_list.append(gf_name)
        elif gf_type == "EXTERNAL":
            external_variables_list.append(gf_name)
        elif gf_type == "CORE":
           core_variables_list.append(gf_name)
        elif gf_type == "TILE_TMP":
           tile_tmp_variables_list.append(gf_name)
        elif gf_type == "SCALAR_TMP":
           scalar_tmp_variables_list.append(gf_name)
        else:
           raise Exception("Bad gftype '"+gf_type+"'")

    # Next we alphabetize the lists
    evolved_variables_list.sort()
    auxiliary_variables_list.sort()
    auxevol_variables_list.sort()
    external_variables_list.sort()
    core_variables_list.sort()
    tile_tmp_variables_list.sort()
    scalar_tmp_variables_list.sort()

    return evolved_variables_list, auxiliary_variables_list, auxevol_variables_list, external_variables_list, core_variables_list, tile_tmp_variables_list, scalar_tmp_variables_list

def gridfunction_defines():
    evolved_variables_list, auxiliary_variables_list, auxevol_variables_list, external_variables_list, core_variables_list, tile_tmp_variables_list, scalar_tmp_variables_list  = gridfunction_lists()

    outstr  = """// EVOLVED VARIABLES:
#define NUM_EVOL_GFS """ + str(len(evolved_variables_list)) + "\n"
    for i in range(len(evolved_variables_list)):
        outstr += "#define " + evolved_variables_list[i].upper() + "GF\t" + str(i) + "\n"

    outstr += """\n\n// AUXILIARY VARIABLES:
#define NUM_AUX_GFS """ + str(len(auxiliary_variables_list)) + "\n"
    for i in range(len(auxiliary_variables_list)):
        outstr += "#define " + auxiliary_variables_list[i].upper() + "GF\t" + str(i) + "\n"

    outstr += """\n\n// AUXEVOL VARIABLES:
#define NUM_AUXEVOL_GFS """ + str(len(auxevol_variables_list)) + "\n"
    for i in range(len(auxevol_variables_list)):
        outstr += "#define " + auxevol_variables_list[i].upper() + "GF\t" + str(i) + "\n"

    outstr += """\n\n// EXTERNAL VARIABLES:
#define NUM_EXTERNAL_GFS """ + str(len(external_variables_list)) + "\n"
    for i in range(len(external_variables_list)):
        outstr += "#define " + external_variables_list[i].upper() + "GF\t" + str(i) + "\n"
#If CarpetX: Do I need to do this for CORE vars? I don't think so.

    if(len(evolved_variables_list)) > 0:
        outstr += """\n\n// SET gridfunctions_f_infinity[i] = value of gridfunction i in the limit r->infinity:
static const REAL gridfunctions_f_infinity[NUM_EVOL_GFS] = { """
        for evol_var in evolved_variables_list:  # This list is sorted
            #                                      We need to preserve the order to ensure consistency with the #defines
            for gf in glb_gridfcs_list:
                if gf.name == evol_var and gf.gftype == "EVOL":
                    outstr += str(gf.f_infinity) + ", "
        outstr = outstr[:-2] + " };\n"

        outstr += """\n\n// SET gridfunctions_wavespeed[i] = gridfunction i's characteristic wave speed:
static const REAL gridfunctions_wavespeed[NUM_EVOL_GFS] = { """
        for evol_var in evolved_variables_list:  # This list is sorted
            #                                      We need to preserve the order to ensure consistency with the #defines
            for gf in glb_gridfcs_list:
                if gf.name == evol_var and gf.gftype == "EVOL":
                    outstr += str(gf.wavespeed) + ", "
        outstr = outstr[:-2] + " };\n"

        # This code could never have worked
        outstr += """\n\n// SET gridfunctions_centering[i] = gridfunction i's vertex/cell centering:
//static const REAL gridfunctions_centering[NUM_EVOL_GFS] = { """
        for evol_var in evolved_variables_list:  # This list is sorted
            #                                      We need to preserve the order to ensure consistency with the #defines
            for gf in glb_gridfcs_list:
                if gf.name == evol_var and gf.gftype == "EVOL":
                    outstr += str(gf.centering) + ", "
        outstr = outstr[:-2] + " };\n"

    return outstr


####################
# TO BE DEPRECATED
def output__gridfunction_defines_h__return_gf_lists(outdir):
    with open(os.path.join(outdir, "gridfunction_defines.h"), "w") as file:
        file.write("/* This file is automatically generated by NRPy+. Do not edit. */\n\n")
        file.write(gridfunction_defines())
    return gridfunction_lists()
####################

#from outputC import outC_NRPy_basic_defines_h_dict
import defines_dict
def register_C_functions_and_NRPy_basic_defines(enable_griddata_struct=True,
                                                list_of_extras_in_griddata_struct=None):
    # First register C functions needed by grid

    # Then set up the dictionary entry for grid in NRPy_basic_defines
    Nbd_str  = gridfunction_defines()
    Nbd_str += r"""
// Declare the IDX4S(gf,i,j,k) macro, which enables us to store 4-dimensions of
//   data in a 1D array. In this case, consecutive values of "i"
//   (all other indices held to a fixed value) are consecutive in memory, where
//   consecutive values of "j" (fixing all other indices) are separated by
//   Nxx_plus_2NGHOSTS0 elements in memory. Similarly, consecutive values of
//   "k" are separated by Nxx_plus_2NGHOSTS0*Nxx_plus_2NGHOSTS1 in memory, etc.
#define IDX4S(g,i,j,k)                                                  \
  ( (i) + Nxx_plus_2NGHOSTS0 * ( (j) + Nxx_plus_2NGHOSTS1 * ( (k) + Nxx_plus_2NGHOSTS2 * (g) ) ) )
#define IDX4ptS(g,idx) ( (idx) + (Nxx_plus_2NGHOSTS0*Nxx_plus_2NGHOSTS1*Nxx_plus_2NGHOSTS2) * (g) )
#define IDX3S(i,j,k) ( (i) + Nxx_plus_2NGHOSTS0 * ( (j) + Nxx_plus_2NGHOSTS1 * ( (k) ) ) )
#define LOOP_REGION(i0min,i0max, i1min,i1max, i2min,i2max)              \
  for(int i2=i2min;i2<i2max;i2++) for(int i1=i1min;i1<i1max;i1++) for(int i0=i0min;i0<i0max;i0++)
#define LOOP_OMP(__OMP_PRAGMA__, i0,i0min,i0max, i1,i1min,i1max, i2,i2min,i2max) _Pragma(__OMP_PRAGMA__) \
    for(int (i2)=(i2min);(i2)<(i2max);(i2)++) for(int (i1)=(i1min);(i1)<(i1max);(i1)++) for(int (i0)=(i0min);(i0)<(i0max);(i0)++)
#define LOOP_NOOMP(i0,i0min,i0max, i1,i1min,i1max, i2,i2min,i2max)      \
  for(int (i2)=(i2min);(i2)<(i2max);(i2)++) for(int (i1)=(i1min);(i1)<(i1max);(i1)++) for(int (i0)=(i0min);(i0)<(i0max);(i0)++)
#define LOOP_BREAKOUT(i0,i1,i2, i0max,i1max,i2max) i0=(i0max); i1=(i1max); i2=(i2max); break;
#define IS_IN_GRID_INTERIOR(i0i1i2, Nxx_plus_2NGHOSTS0,Nxx_plus_2NGHOSTS1,Nxx_plus_2NGHOSTS2, NG) \
  ( i0i1i2[0] >= (NG) && i0i1i2[0] < (Nxx_plus_2NGHOSTS0)-(NG) &&       \
    i0i1i2[1] >= (NG) && i0i1i2[1] < (Nxx_plus_2NGHOSTS1)-(NG) &&       \
    i0i1i2[2] >= (NG) && i0i1i2[2] < (Nxx_plus_2NGHOSTS2)-(NG) )
"""
    if enable_griddata_struct:
        Nbd_str += """
typedef struct __griddata__ {
  // griddata_struct stores data needed on each grid
"""
        griddata_struct_list = [glb_griddata(__name__, "REAL *restrict xx[3];")]
        for item in griddata_struct_list + glb_griddata_struct_list:
            Nbd_str += "  " + item.string + "  // <- registered by NRPy+ module \"" + item.module + "\"\n"
        if isinstance(list_of_extras_in_griddata_struct, list):
            for extra in list_of_extras_in_griddata_struct:
                Nbd_str += "  " + extra + ";\n"
        Nbd_str += "} griddata_struct;\n"

    defines_dict.outC_NRPy_basic_defines_h_dict["grid"] = Nbd_str
