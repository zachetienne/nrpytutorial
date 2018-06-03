import NRPy_param_funcs as par
from collections import namedtuple
import sympy as sp

# grid.py: functions & parameters related to numerical grids:
# functions: Automatic loop output, output C code needed for gridfunction memory I/O, gridfunction registration

# Initialize globals related to the grid
glb_gridfcs_list = []
glb_gridfc  = namedtuple('gridfunction', 'gftype name')

thismodule = __name__
par.initialize_param(par.glb_param("char", thismodule, "GridFuncMemAccess",   "SENRlike"))
par.initialize_param(par.glb_param("char", thismodule, "MemAllocStyle","kji"))
par.initialize_param(par.glb_param("int",  thismodule, "DIM", 3))
par.initialize_param(par.glb_param("int",  thismodule, "Nx[DIM]", "SetAtCRuntime"))

def variable_type(var):
    var_is_gf = False
    for gf in range(len(glb_gridfcs_list)):
        if str(var) == glb_gridfcs_list[gf].name:
            var_is_gf = True
    var_is_parameter = False
    for paramname in range(len(par.glb_params_list)):
        if str(var) == par.glb_params_list[paramname].parname:
            var_is_parameter = True
    if var_is_parameter and var_is_gf:
        print("Error: variable "+str(var)+" is registered both as a gridfunction and as a Cparameter.")
        exit(1)
    if not (var_is_parameter or var_is_gf):
        return "other"
    if var_is_parameter:
        return "Cparameter"
    if var_is_gf:
        return "gridfunction"

def gfaccess(gfarrayname = "",varname = "",ijklstring = ""):
    DIM = par.parval_from_str("DIM")
    if par.parval_from_str("GridFuncMemAccess") == "SENRlike":
        if gfarrayname == "":
            print("Error: GridFuncMemAccess = SENRlike requires gfarrayname be passed to gfaccess()")
            exit(1)
        # Return gfarrayname[IDX3(varname,i0)] for DIM=1, gfarrayname[IDX3(varname,i0,i1)] for DIM=2, etc.
        retstring = gfarrayname + "[IDX" + str(DIM+1) + "(" + varname + ", "
    elif par.parval_from_str("GridFuncMemAccess") == "ETK":
        # Return varname[CCTK_GFINDEX3D(i0,i1,i2)] for DIM=3. Error otherwise
        if DIM != 3:
            print("Error: GridFuncMemAccess = ETK currently requires that gridfunctions be 3D. Can be easily extended.")
            exit(1)
        retstring = varname + "[CCTK_GFINDEX"+str(DIM)+"D(cctkGH, "
    if ijklstring == "":
        for i in range(DIM):
            retstring += "i"+str(i)
            if i != DIM-1:
                retstring += ', '
    else:
        retstring += ijklstring
    return retstring + ")]"


import sys
def register_gridfunctions(gf_type,gf_names):
    # First convert gf_names to a list if it's not already a list
    if type(gf_names) is not list:
        gf_namestmp = [gf_names]
        gf_names = gf_namestmp
    # Next check for duplicates, and error out if any are found.
    for i in range(len(gf_names)):
        # https://stackoverflow.com/questions/1303243/how-to-find-out-if-a-python-object-is-a-string
        if sys.version_info[0] < 3:
            if not isinstance(gf_names[i],basestring):
                print("ERROR: gf_names must be strings")
                exit(1)
        else:
            if not isinstance(gf_names[i], str):
                print("ERROR: gf_names must be strings")
                exit(1)

        for j in range(len(glb_gridfcs_list)):
            if gf_names[i] == glb_gridfcs_list[j].name:
                print("Error: Tried to register the gridfunction "+gf_names[i]+" twice (ignored type)")
                exit(1)
        # If no duplicate found, append to "gridfunctions" list:
        glb_gridfcs_list.append(glb_gridfc(gf_type,gf_names[i]))
    if not (gf_type == "EVOL" or gf_type == "AUX"):
        print("Error in registering gridfunction(s) with unsupported type "+gf_type+".")
        print("Supported types include \"EVOL\" for gridfunctions related to evolved quantities or \"AUX\" for all others.")
        exit(1)
    OBJ_TMPS = []
    for i in range(len(gf_names)):
        OBJ_TMPS.append(sp.sympify(gf_names[i]))
    return OBJ_TMPS
