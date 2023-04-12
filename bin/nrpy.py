import sys

# Step 0: Check Python version. NRPy+ is untested below Python 2.7-ish, but is compatible with 2.7+ and 3+
PYTHONVERSION3 = False
if sys.version_info[0]==2:
    if sys.version_info[1]<4:
        print("Sorry, NRPy won't work with Python < 2.4; sorting functs won't work. See https://docs.python.org/3/howto/sorting.html for details.")
        sys.exit(1)
if sys.version_info[0]==3:
    PYTHONVERSION3 = True

# Step 1: Print logo.
#from NRPy_logo import print_logo
#print_logo()

# Step 2: Initialize core parameter NRPy::MainModule,
#         which defines the desired main module.
#         E.g., scalarwave, BSSN_RHSs, BSSN_InitialData, etc.
# Contains parameter initialization, manipulation, and read-in routines
import NRPy_param_funcs as par

# Used to import needed modules dynamically
import importlib

# Initialize the MainModule parameter.
# This is the ONLY parameter initialized outside of a module!
MainModule = "scalarwave" # Default. To be overwritten later.
par.initialize_param(par.glb_param("char","NRPy","MainModule",MainModule))

# Step 4: Initialize NRPy+ as desired.

# Step 4a: Enter Interactive Mode if NRPy+ is run via
#   `python nrpy.py`
if(len(sys.argv) == 1):
    print("/* Run `python nrpy.py --help` for other command-line options */")
    print("/* Entering interactive mode                                  */\n")

#   Print help message if NRPy+ is run via
#    `python nrpy.py --help`
elif(len(sys.argv) == 2 and sys.argv[1] == "--help"):
    print("\n     \033[1m............................................\033[0m ")
    print("     -={ \033[1mNRPy+ supports multiple usage modes.\033[0m }=-\n")
    print("\033[1mUsage Mode 0\033[0m: `python nrpy.py`\n\t\t  initializes interactive session\n")
    print("\033[1mUsage Mode 1\033[0m: `python nrpy.py --help`\n\t\t  outputs this message\n")
    print("\033[1mUsage Mode 2\033[0m: `python nrpy.py --gen-defparam-file`\n\t\t generates default parameter file\n")
    print("\033[1mUsage Mode 3\033[0m: `python nrpy.py [PARAMETER FILE]`\n\t\t override default parameters with parameter file.")
    print("\t\t Parameter file takes the form of a list;")
    print("\t\t each list item has syntax `mainmodule:paramname=value`\n")
    print("\033[1mUsage Mode 4\033[0m: `python nrpy.py [PARAMETER FILE] [PARAMETER OVERRIDES]`")
    print("\t\t read parameter file, then override parameter file")
    print("\t\t settings with command-line parameters, with same")
    print("\t\t syntax (`mainmodule:paramname=value`)")
    print("     \033[1m............................................\033[0m ")
    sys.exit(0)

#    Run with parameter file & optional list of parameter overrides if NRPy+ is run via
#   `python nrpy.py [PARAMETER FILE] [(optional) PARAMETER OVERRIDES]`:
# Note that
# 1) parameters set in parameter file override parameter defaults, and
# 2) parameters set at command line override both parameter defaults
#    *and* parameter settings in parameter file
elif(len(sys.argv) >= 2):
    # When not in an interactive mode, the NRPy::MainModule parameter must be set,
    #   either in the param file (preferred!) or as a command line argument.
    with open(sys.argv[1], "r") as file:
        for line in file:
            par.set_paramsvals_value(line, sys.argv[1], FindMainModuleMode=True)
    # Search for NRPy::MainModule in the command line arguments
    for i in range(2,len(sys.argv)):
        par.set_paramsvals_value(sys.argv[i], "", FindMainModuleMode=True)
    # The NRPy::MainModule parameter has already been set,
    #   and parse_param_string__set__params_and_paramsvars()
    #   will error out with a "Critical Error" if it has not been set.
    idx = par.get_params_idx(par.glb_param("ignoretype", "NRPy", "MainModule", "ignoredefval"))
    MainModule = par.glb_paramsvals_list[idx]
    if MainModule == "NODEFAULT":
        print("Error: Could not find NRPy::MainModule defined in the parameter file \""+sys.argv[1]+"\" or on the command line!")
        sys.exit(1)

    # Next initialize all of MainModule's parameters.
    # Note that MainModule must also initialize parameters for modules
    # on which it depends, except outputC (which NRPy+ loads by default).
    # https://stackoverflow.com/questions/10675054/how-to-import-a-module-in-python-with-importlib-import-module
    importlib.import_module(MainModule+"."+MainModule)

    # Next overwrite default parameters with values specified in the parameter file.
    with open(sys.argv[1], "r") as file:
        for line in file:
            par.set_paramsvals_value(line, sys.argv[1])
    # Next overwrite default parameters and values specified in the parameter file with command-line parameter assignments.
    for i in range(2,len(sys.argv)):
        par.set_paramsvals_value(sys.argv[i], "")

# Next load the MainModule, if it hasn't been loaded already.
#importlib.import_module(MainModule+"."+MainModule)
getattr(importlib.import_module(MainModule+"."+MainModule), MainModule)()
# Initialize parameters for the core outputC module,
# which is the default MainModule.
# Call function initparams() from outputC module (in outputC.py):
# getattr(importlib.import_module(MainModule), 'initparams')()

# Step 3b (temporary): Set a SymPy expression to test processing
# import sympy as sp
#from outputC import *
#getattr(importlib.import_module("outputC"),'outputC')([sympify("3*a*b**4+c*sin(a*b**4)"),sympify("6*a*b**4")],["output1","output2"])
#getattr(importlib.import_module("outputC"),'outputC')(sp.sympify("3*a*b**4+2*c*sin(a*b**4)"),"output1")
