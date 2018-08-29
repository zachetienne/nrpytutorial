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
par.initialize_param(par.glb_param("char", thismodule, "MemAllocStyle","210"))
par.initialize_param(par.glb_param("INT",  thismodule, "DIM", 3))
par.initialize_param(par.glb_param("INT",  thismodule, "Nx[DIM]", "SetAtCRuntime"))

xx = par.Cparameters("REALARRAY",thismodule,["xx0","xx1","xx2","xx3"])

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
    found_registered_gf = False
    for gf in glb_gridfcs_list:
        if gf.name == varname:
            if found_registered_gf:
                print("Error: found duplicate gridfunction name: "+gf.name)
                exit(1)
            found_registered_gf = True

    if not found_registered_gf:
        print("Error: gridfunction \""+varname+"\" is not registered!")
        exit(1)
       
    DIM = par.parval_from_str("DIM")
    retstring = ""
    if par.parval_from_str("GridFuncMemAccess") == "SENRlike":
        if gfarrayname == "":
            print("Error: GridFuncMemAccess = SENRlike requires gfarrayname be passed to gfaccess()")
            exit(1)
        # Return gfarrayname[IDX3(varname,i0)] for DIM=1, gfarrayname[IDX3(varname,i0,i1)] for DIM=2, etc.
        retstring += gfarrayname + "[IDX" + str(DIM+1) + "(" + varname.upper()+"GF" + ", "
    elif par.parval_from_str("GridFuncMemAccess") == "ETK":
        # Return varname[CCTK_GFINDEX3D(i0,i1,i2)] for DIM=3. Error otherwise
        if DIM != 3:
            print("Error: GridFuncMemAccess = ETK currently requires that gridfunctions be 3D. Can be easily extended.")
            exit(1)
        retstring += varname + "GF" + "[CCTK_GFINDEX"+str(DIM)+"D(cctkGH, "
    else:
        print("grid::GridFuncMemAccess = "+par.parval_from_str("GridFuncMemAccess")+" not supported")
        exit(1)
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
    # Step 1: convert gf_names to a list if it's not already a list
    if type(gf_names) is not list:
        gf_namestmp = [gf_names]
        gf_names = gf_namestmp

    # Step 2: check that all gridfunction names are strings. Python3 & Python2 have different behavior here:
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

    # Step 3: Verify that gridfunction type is valid.
    if not (gf_type == "EVOL" or gf_type == "AUX" or gf_type == "COORDARRAY"):
        print("Error in registering gridfunction(s) with unsupported type "+gf_type+".")
        print("Supported types include \"EVOL\" for gridfunctions related to evolved quantities or \"AUX\" for all others.")
        exit(1)

    # Step 4: Check for duplicate grid function registrations. If:
    #         a) A duplicate is found, error out. Otherwise
    #         b) Append to list of gridfunctions, stored in glb_gridfcs_list[].
    for i in range(len(gf_names)):
        for j in range(len(glb_gridfcs_list)):
            if gf_names[i] == glb_gridfcs_list[j].name:
                print("Error: Tried to register the gridfunction \""+gf_names[i]+"\" twice (ignored type)\n\n")
                exit(1)
        # If no duplicate found, append to "gridfunctions" list:
        glb_gridfcs_list.append(glb_gridfc(gf_type,gf_names[i]))

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
