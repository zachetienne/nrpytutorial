# BaikalETK: NRPy+-Based BSSN Solvers for the Einstein Toolkit

# As documented in the NRPy+ tutorial notebook
#   Tutorial-BaikalETK.ipynb
#   This module generates Baikal and BaikalVacuum,
#   Einstein Toolkit thorns for solving Einstein's
#   equations in the BSSN formalism, in Cartesian
#   coordinates. These thorns are highly optimized
#   for modern CPU architectures, featuring SIMD
#   intrinsics and OpenMP support.
#   In particular,
#   Baikal contains T^{mu nu} source terms
#   and
#   BaikalVacuum does not.

# Author: Zachariah B. Etienne
#         zachetie **at** gmail **dot* com

# INSTRUCTIONS FOR RUNNING THIS:
#  Run it from the parent directory via
# python3 BaikalETK/BaikalETK_main_codegen_driver.py

###############################
# Step 0: Import needed core NRPy+ and standard Python modules
import os, sys                # Standard Python modules for multiplatform OS-level functions
import pickle                   # Standard Python module for converting arbitrary data structures to a uniform format.

nrpy_dir_path = os.path.join(".")
if nrpy_dir_path not in sys.path:
    sys.path.append(nrpy_dir_path)

import grid as gri              # NRPy+: Functions having to do with numerical grids
import NRPy_param_funcs as par  # NRPy+: Parameter interface
###############################

###############################
# Step 1: Set compile-time and runtime parameters for both Baikal & BaikalVacuum:

# Step 1.a: Set compile-time (i.e., NRPy+-time) parameters
LapseCondition = "OnePlusLog"
ShiftCondition = "GammaDriving2ndOrder_NoCovariant"

# Step 1.b: Set runtime parameters for Baikal and BaikalVacuum
#         Current runtime choices:
#         Baikal:       FD_orders = [2,4]  ; Gamma-driving eta parameter; Kreiss-Oliger dissipation strength
#         BaikalVacuum: FD_orders = [4,6,8]; Gamma-driving eta parameter; Kreiss-Oliger dissipation strength
paramslist = []
FD_orders = [2,4,6,8]
WhichParamSet = 0
for WhichPart in ["BSSN_RHSs","Ricci","BSSN_constraints","detgammabar_constraint"]:
    for FD_order in FD_orders:
        enable_stress_energy_list = [True,False]
        if FD_order == 2:
            enable_stress_energy_list = [True]
        elif FD_order >= 6:
            enable_stress_energy_list = [False]
        for enable_stress_energy in enable_stress_energy_list:
            ThornName = "Baikal"
            if enable_stress_energy == False:
                ThornName = "BaikalVacuum"
            paramstr = "WhichPart="+WhichPart+","
            paramstr+= "WhichParamSet="+str(WhichParamSet)+","
            paramstr+= "ThornName="+ThornName+","
            paramstr+= "FD_order="+str(FD_order)+","
            paramstr+= "LapseCondition="+LapseCondition+","
            paramstr+= "ShiftCondition="+ShiftCondition+","
            paramstr+= "enable_stress_energy_source_terms="+str(enable_stress_energy)
            paramslist.append(paramstr)
            WhichParamSet = WhichParamSet + 1
###############################

###############################
# Step 2: Generate all C-code kernels for Baikal and BaikalVacuum,
#         in parallel if supported by this OS.
nrpy_dir_path = os.path.join(".")
if nrpy_dir_path not in sys.path:
    sys.path.append(nrpy_dir_path)

NRPyEnvVars = []

import time   # Standard Python module for benchmarking
start = time.time()
if __name__ == "__main__":
    try:
        if os.name == 'nt':
            # Windows & Jupyter multiprocessing do not mix, so we run in serial on Windows.
            #  Here's why: https://stackoverflow.com/questions/45719956/python-multiprocessing-in-jupyter-on-windows-attributeerror-cant-get-attribut
            raise Exception("Parallel codegen currently not available in Windows")
        # Step 3.d.ii: Import the multiprocessing module.
        import multiprocessing
        
        # Step 3.d.iii: Define master function for parallelization.
        #           Note that lambdifying this doesn't work in Python 3
        def master_func(i):
            import BaikalETK_validate.BaikalETK_C_kernels_codegen as BCk
            return BCk.BaikalETK_C_kernels_codegen_onepart(NRPyDir=os.path.join("."),params=paramslist[i])
        
        # Step 3.d.iv: Evaluate list of functions in parallel if possible;
        #           otherwise fallback to serial evaluation:
        pool = multiprocessing.Pool() #processes=len(paramslist))
        NRPyEnvVars.append(pool.map(master_func,range(len(paramslist))))
        pool.terminate()
        pool.join()
    except:
        # Steps 3.d.ii-iv, alternate: As fallback, evaluate functions in serial.
        #       This will happen on Android and Windows systems
        import BaikalETK_validate.BaikalETK_C_kernels_codegen as BCk
        for param in paramslist:
            NRPyEnvVars.append(BCk.BaikalETK_C_kernels_codegen_onepart(NRPyDir=os.path.join("."),params=param))

print("Finished C kernel codegen for Baikal and BaikalVacuum in "+str(time.time()-start)+" seconds.")
###############################

###############################
# Step 3: Generate Einstein Toolkit ccl files

# Step 3.a: Create an appropriate NRPy+ environment based on
#           environment variables set ("pickled") in the above codegen
#           and stored in the NRPyEnvVars. Here we un-pickle the data
#           to construct the NRPy+ environment.

# Long description:
# The Einstein Toolkit (ETK) ccl files contain runtime parameters (param.ccl),
# registered gridfunctions (interface.ccl), and function scheduling (schedule.ccl).
# As parameters and gridfunctions are registered with NRPy+ when the C-code kernels
# are generated, and this generation occurs on separate processes in parallel, we
# store the entire NRPy+ environment for each process. This results in a tremendous
# amount of duplication, which is sorted out next. Once all duplicated environment
# variables (e.g., registered gridfunctions) are removed, we replace the current NRPy+
# environment with the new one, by setting
# gri.glb_gridfcs_list[],par.glb_params_list[],par.glb_Cparams_list[].

# Store all NRPy+ environment variables to file so NRPy+ environment from within this subprocess can be easily restored

# https://www.pythonforthelab.com/blog/storing-binary-data-and-serializing/
grfcs_list = []
param_list = []
Cparm_list = []

for WhichParamSet in NRPyEnvVars[0]:
    # gridfunctions
    i=0
    num_elements = pickle.loads(WhichParamSet[i]); i+=1
    for lst in range(num_elements):
        grfcs_list.append(gri.glb_gridfc(gftype=pickle.loads(WhichParamSet[i+0]),
                                         name  =pickle.loads(WhichParamSet[i+1]),
                                         rank  =pickle.loads(WhichParamSet[i+2]),
                                         DIM   =pickle.loads(WhichParamSet[i+3]))) ; i+=4
    # parameters
    num_elements = pickle.loads(WhichParamSet[i]); i+=1
    for lst in range(num_elements):
        param_list.append(par.glb_param(type      =pickle.loads(WhichParamSet[i+0]),
                                        module    =pickle.loads(WhichParamSet[i+1]),
                                        parname   =pickle.loads(WhichParamSet[i+2]),
                                        defaultval=pickle.loads(WhichParamSet[i+3]))) ; i+=4
    # Cparameters
    num_elements = pickle.loads(WhichParamSet[i]); i+=1
    for lst in range(num_elements):
        Cparm_list.append(par.glb_Cparam(type      =pickle.loads(WhichParamSet[i+0]),
                                         module    =pickle.loads(WhichParamSet[i+1]),
                                         parname   =pickle.loads(WhichParamSet[i+2]),
                                         defaultval=pickle.loads(WhichParamSet[i+3]))) ; i+=4
grfcs_list_uniq = []
for gf_ntuple_stored in grfcs_list:
    found_gf = False
    for gf_ntuple_new in grfcs_list_uniq:
        if gf_ntuple_new == gf_ntuple_stored:
            found_gf = True
    if found_gf == False:
        grfcs_list_uniq.append(gf_ntuple_stored)

param_list_uniq = []
for pr_ntuple_stored in param_list:
    found_pr = False
    for pr_ntuple_new in param_list_uniq:
        if pr_ntuple_new == pr_ntuple_stored:
            found_pr = True
    if found_pr == False:
        param_list_uniq.append(pr_ntuple_stored)

# Set glb_paramsvals_list:
# Step 1: Reset all paramsvals to their defaults
par.glb_paramsvals_list = []
for parm in param_list_uniq:
    par.glb_paramsvals_list.append(parm.defaultval)

Cparm_list_uniq = []
for Cp_ntuple_stored in Cparm_list:
    found_Cp = False
    for Cp_ntuple_new in Cparm_list_uniq:
        if Cp_ntuple_new == Cp_ntuple_stored:
            found_Cp = True
    if found_Cp == False:
        Cparm_list_uniq.append(Cp_ntuple_stored)

gri.glb_gridfcs_list = []
par.glb_params_list  = []
par.glb_Cparams_list = []

gri.glb_gridfcs_list = grfcs_list_uniq
par.glb_params_list  = param_list_uniq
par.glb_Cparams_list = Cparm_list_uniq

# Step 3.b: Override defaults with values used here.
#           Note that almost no NRPy+ parameters
#           are used after this point (DIM is definitely
#           one exception), so most of these lines have no
#           effect.
par.set_parval_from_str("grid::DIM",3)
import reference_metric as rfm
par.set_parval_from_str("reference_metric::CoordSystem","Cartesian")
rfm.reference_metric() # Create ReU, ReDD needed for rescaling B-L initial data, generating BSSN RHSs, etc.
par.set_parval_from_str("grid::GridFuncMemAccess","ETK")
par.set_parval_from_str("BSSN.BSSN_gauge_RHSs::ShiftEvolutionOption", ShiftCondition)
par.set_parval_from_str("BSSN.BSSN_gauge_RHSs::LapseEvolutionOption", LapseCondition)

# Step 3.c: Now that NRPy+ parameters and gridfunctions have been
#           defined, we now have all the information we need loaded
#           into our NRPy+ environment to generate the Einstein
#           Toolkit ccl files
import BaikalETK_validate.BaikalETK_ETK_ccl_files_codegen as cclgen

for enable_stress_energy_source_terms in [True,False]:
    ThornName="Baikal"
    if enable_stress_energy_source_terms==False:
        ThornName="BaikalVacuum"
    cclgen.output_param_ccl(ThornName,enable_stress_energy_source_terms)
    cclgen.output_interface_ccl(ThornName,enable_stress_energy_source_terms)
    cclgen.output_schedule_ccl(ThornName,enable_stress_energy_source_terms)
###############################

###############################
# Step 4: Generate C driver functions for ETK registration & NRPy+-generated kernels

# Now that we have constructed the basic C code kernels and the
# needed Einstein Toolkit ccl files, we next write the driver
# functions for registering BaikalETK within the Toolkit and the
# C code kernels. Each of these driver functions will be called
# directly from the thorn's schedule.ccl in the ETK.

# Step 4.a: First we call the functions in the `BaikalETK.BaikalETK_C_drivers_codegen`
#           Python module) to store all needed driver C files to a Python dictionary,
#           then we simply output the dictionary to the appropriate files.
import BaikalETK_validate.BaikalETK_C_drivers_codegen as driver

# The following Python dictionaries consist of a key, which is the filename
#    in the thorn's src/ directory (e.g., "driver_BSSN_constraints.c"),
#    and a value, which is the corresponding source code, stored as a
#    Python string.
Vac_Csrcdict = {}
Reg_Csrcdict = {}

# We'll need lists of gridfunctions for these driver functions
evol_gfs_list    = cclgen.evol_gfs_list
aux_gfs_list     = cclgen.aux_gfs_list
auxevol_gfs_list = cclgen.auxevol_gfs_list

# Generate driver codes for Baikal thorn (i.e., populate the Reg_Csrcdict dictionary)
driver.driver_C_codes(Reg_Csrcdict, "Baikal",
                               cclgen.rhs_list,cclgen.evol_gfs_list,cclgen.aux_gfs_list,cclgen.auxevol_gfs_list,
                               LapseCondition = LapseCondition, enable_stress_energy_source_terms=True)

# Generate driver codes for BaikalVacuum thorn (i.e., populate the Vac_Csrcdict dictionary)
driver.driver_C_codes(Vac_Csrcdict, "BaikalVacuum",
                               cclgen.rhs_list,cclgen.evol_gfs_list,cclgen.aux_gfs_list,cclgen.auxevol_gfs_list,
                               LapseCondition = LapseCondition, enable_stress_energy_source_terms=False)

# Next we output the contents of the Reg_Csrcdict and
#   Vac_Csrcdict dictionaries to files in the respective
#   thorns' directories.
for key,val in Reg_Csrcdict.items():
    with open(os.path.join("Baikal","src",key),"w") as file:
        file.write(val)
for key,val in Vac_Csrcdict.items():
    with open(os.path.join("BaikalVacuum","src",key),"w") as file:
        file.write(val)
###############################

###############################
# Step 5: Generate the make.code.defn file, needed for the Einstein Toolkit build system

# Finally output the thorns' make.code.defn files, consisting of
#   a list of all C codes in the above dictionaries. This is
#   part of the ETK build system so that these files are output.

def output_make_code_defn(dictionary, ThornName):
    with open(os.path.join(ThornName,"src","make.code.defn"), "w") as file:
        file.write("""
# Main make.code.defn file for thorn """+ThornName+"""

# Source files in this directory
SRCS =""")
        filestring = ""
        list_of_C_driver_files = list(dictionary.keys())
        for i in range(len(list_of_C_driver_files)):
            filestring += "      "+list_of_C_driver_files[i]
            if i != len(list_of_C_driver_files)-1:
                filestring += " \\\n"
            else:
                filestring += "\n"
        file.write(filestring)

output_make_code_defn(Reg_Csrcdict,"Baikal")
output_make_code_defn(Vac_Csrcdict,"BaikalVacuum")
###############################