#!/usr/bin/env python
# coding: utf-8

# <script async src="https://www.googletagmanager.com/gtag/js?id=UA-59152712-8"></script>
# <script>
#   window.dataLayer = window.dataLayer || [];
#   function gtag(){dataLayer.push(arguments);}
#   gtag('js', new Date());
# 
#   gtag('config', 'UA-59152712-8');
# </script>
# 
# # Start-to-Finish Example: Validating `GiRaFFEfood_NRPy` against original, trusted `GiRaFFEfood`: 
# 
# ## Author: Patrick Nelson
# 
# **Notebook Status:** <font color='green'><b>Validated</b></font>
# 
# **Validation Notes:** This module validates all expressions used to set up initial data in 
# * [Tutorial-GiRaFFEfood_NRPy_Aligned_Rotator](Tutorial-GiRaFFEfood_NRPy_Aligned_Rotator.ipynb)
# * [Tutorial-GiRaFFEfood_NRPy_1D_tests](Tutorial-GiRaFFEfood_NRPy_1D_tests.ipynb), 
# * [Tutorial-GiRaFFEfood_NRPy_1D_tests-fast_wave](Tutorial-GiRaFFEfood_NRPy_1D_tests-fast_wave.ipynb), 
# * [Tutorial-GiRaFFEfood_NRPy_1D_tests-degen_Alfven_wave](Tutorial-GiRaFFEfood_NRPy_1D_tests-degen_Alfven_wave.ipynb), 
# * [Tutorial-GiRaFFEfood_NRPy_1D_tests-three_waves](Tutorial-GiRaFFEfood_NRPy_1D_tests-three_waves.ipynb), and
# * [Tutorial-GiRaFFEfood_NRPy_1D_tests-FFE_breakdown](Tutorial-GiRaFFEfood_NRPy_1D_tests-FFE_breakdown.ipynb), 
# 
# against the C-code implementation of these expressions found in the original (trusted) [`GiRaFFEfood` Einstein Toolkit thorn](link), and confirms roundoff-level agreement.
# 
# ### NRPy+ Source Code for this module: 
# * [GiRaFFEfood_NRPy/GiRaFFEfood_NRPy_Aligned_Rotator.py](../../edit/in_progress/GiRaFFEfood_NRPy/GiRaFFEfood_NRPy_Aligned_Rotator.py) [\[**tutorial**\]](Tutorial-GiRaFFEfood_NRPy_Aligned_Rotator.ipynb) Generates Aligned Rotator initial data
# * [GiRaFFEfood_NRPy/GiRaFFEfood_NRPy_1D_tests.py](../../edit/in_progress/GiRaFFEfood_NRPy/GiRaFFEfood_NRPy_1D_tests.py) [\[**tutorial**\]](Tutorial-GiRaFFEfood_NRPy_1D_tests.ipynb) Generates Alfv&eacute;n Wave initial data.
# * [GiRaFFEfood_NRPy/GiRaFFEfood_NRPy_1D_tests_fast_wave.py](../../edit/in_progress/GiRaFFEfood_NRPy/GiRaFFEfood_NRPy_1D_tests_fast_wave.py) [\[**tutorial**\]](Tutorial-GiRaFFEfood_NRPy_1D_tests-fast_wave.ipynb) Generates Alfv&eacute;n Wave initial data.
# * [GiRaFFEfood_NRPy/GiRaFFEfood_NRPy_1D_tests_degen_Alfven_wave.py](../../edit/in_progress/GiRaFFEfood_NRPy/GiRaFFEfood_NRPy_1D_tests_degen_Alfven_wave.py) [\[**tutorial**\]](Tutorial-GiRaFFEfood_NRPy_1D_tests-degen_Alfven_wave.ipynb) Generates Degenerate Alfv&eacute;n Wave initial data.
# * [GiRaFFEfood_NRPy/GiRaFFEfood_NRPy_1D_tests_three_waves.py](../../edit/in_progress/GiRaFFEfood_NRPy/GiRaFFEfood_NRPy_1D_tests_three_waves.py) [\[**tutorial**\]](Tutorial-GiRaFFEfood_NRPy_1D_tests-three_waves.ipynb) Generates Three Waves initial data.
# * [GiRaFFEfood_NRPy/GiRaFFEfood_NRPy_1D_tests_FFE_breakdown.py](../../edit/in_progress/GiRaFFEfood_NRPy/GiRaFFEfood_NRPy_1D_tests_FFE_breakdown.py) [\[**tutorial**\]](Tutorial-GiRaFFEfood_NRPy_1D_tests-FFE_breakdown.ipynb) Generates FFE Breakdown initial data.
# 
# ## Introduction:
# 
# This notebook validates the initial data routines that we will use for `GiRaFFE_NRPy`, collectively referred to as `GiRaFFEfood_NRPy`. To do so, we will generate the initial data with both our code and the original `GiRaFFEfood` code. Then, we will directly compare the velocities and show round-off level agreement between the two. 
# 
# When this notebook is run, the significant digits of agreement between the old `GiRaFFE` and new `GiRaFFE_NRPy` versions of the algorithm will be evaluated. If the agreement falls below a thresold, the point, quantity, and level of agreement are reported [here](#compile_run).
# 

# <a id='toc'></a>
# 
# # Table of Contents
# $$\label{toc}$$
# 
# This notebook is organized as follows
# 
# 1. [Step 1](#setup): Set up core functions and parameters for unit testing the initial data algorithms
#     1. [Step 1.a](#initial_data) Generate the initial data C function
#     1. [Step 1.b](#download) Download original `GiRaFFE` files
#     1. [Step 1.c](#free_params) Output C codes needed for declaring and setting Cparameters; also set `free_parameters.h`
#     1. [Step 1.d](#interface) Create dummy files for the CCTK version of the code
# 1. [Step 2](#mainc): `GiRaFFEfood_NRPy_unit_test.c`: The Main C Code
#     1. [Step 2.a](#compile_run): Compile and run the code to validate the output
# 1. [Step 3](#drift_notes): Output this notebook to $\LaTeX$-formatted PDF file
# 1. [Step 4](#latex_pdf_output): Output this notebook to $\LaTeX$-formatted PDF file

# <a id='setup'></a>
# 
# # Step 1: Set up core functions and parameters for unit testing the initial data algorithms" \[Back to [top](#toc)\]
# 
# $$\label{setup}$$
# 
# We'll start by appending the relevant paths to `sys.path` so that we can access sympy modules in other places. Then, we'll import NRPy+ core functionality and set up a directory in which to carry out our test. We will also declare the gridfunctions that are needed for this portion of the code.

# In[1]:


# There are several initial data routines we need to test. We'll control which one we use with a string option
initial_data = "AlfvenWave" # Valid options: "AllTests", "AlignedRotator", "AlfvenWave", "FastWave",
                           # "DegenAlfvenWave", "ThreeWaves", "FFE_Breakdown"


# In[2]:


import os, sys           # Standard Python modules for multiplatform OS-level functions
# First, we'll add the parent directory to the list of directories Python will check for modules.
nrpy_dir_path = os.path.join("..")
if nrpy_dir_path not in sys.path:
    sys.path.append(nrpy_dir_path)
nrpy_dir_path = os.path.join("..","..")
if nrpy_dir_path not in sys.path:
    sys.path.append(nrpy_dir_path)

from outputC import outCfunction, lhrh # NRPy+: Core C code output module
import sympy as sp               # SymPy: The Python computer algebra package upon which NRPy+ depends
import finite_difference as fin  # NRPy+: Finite difference C code generation module
import NRPy_param_funcs as par   # NRPy+: Parameter interface
import grid as gri               # NRPy+: Functions having to do with numerical grids
import indexedexp as ixp         # NRPy+: Symbolic indexed expression (e.g., tensors, vectors, etc.) support
import reference_metric as rfm   # NRPy+: Reference metric support
import cmdline_helper as cmd     # NRPy+: Multi-platform Python command-line interface

out_dir = "Validation/"
cmd.mkdir(out_dir)

thismodule = "Start_to_Finish_UnitTest-GiRaFFEfood_NRPy"

# Register the gridfunctions we need for this function
AD = ixp.register_gridfunctions_for_single_rank1("EVOL","AD")
ValenciavU = ixp.register_gridfunctions_for_single_rank1("AUXEVOL","ValenciavU")
# gammaDD = ixp.register_gridfunctions_for_single_rank2("AUXEVOL","gammaDD","sym01")
betaU = ixp.register_gridfunctions_for_single_rank1("AUXEVOL","betaU")
alpha = gri.register_gridfunctions("AUXEVOL","alpha")


# <a id='initial_data'></a>
# 
# ## Step 1.a: Generate the initial data C function \[Back to [top](#toc)\]
# $$\label{initial_data}$$
# 
# First, we'll use NRPy+ to build the C function that will generate the initial data. There are several different cases here, one for each type of initial test data.

# In[3]:


if initial_data=="AlfvenWave":
    import GiRaFFEfood_NRPy.GiRaFFEfood_NRPy_1D_tests as gid
    gid.GiRaFFEfood_NRPy_1D_tests()
    desc = "Generate Alfven wave 1D initial test data for GiRaFFEfood_NRPy."
elif initial_data=="FastWave":
    import GiRaFFEfood_NRPy.GiRaFFEfood_NRPy_1D_tests_fast_wave as gid
    gid.GiRaFFEfood_NRPy_1D_tests_fast_wave()
    desc = "Generate fast wave 1D initial test data for GiRaFFEfood_NRPy."
elif initial_data=="DegenAlfvenWave":
    import GiRaFFEfood_NRPy.GiRaFFEfood_NRPy_1D_tests_degen_Alfven_wave as gid
    gid.GiRaFFEfood_NRPy_1D_tests_degen_Alfven_wave()
    desc = "Generate degenerate Alfven wave 1D initial test data for GiRaFFEfood_NRPy."
elif initial_data=="ThreeWaves":
    import GiRaFFEfood_NRPy.GiRaFFEfood_NRPy_1D_tests_three_waves as gid
    gid.GiRaFFEfood_NRPy_1D_tests_three_waves()
    desc = "Generate three waves 1D initial test data for GiRaFFEfood_NRPy."
elif initial_data=="FFE_Breakdown":
    import GiRaFFEfood_NRPy.GiRaFFEfood_NRPy_1D_tests_FFE_breakdown as gid
    gid.GiRaFFEfood_NRPy_1D_tests_FFE_breakdown()
    desc = "Generate FFE breakdown 1D initial test data for GiRaFFEfood_NRPy."
elif initial_data=="AlignedRotator":
    import GiRaFFEfood_NRPy.GiRaFFEfood_NRPy_Aligned_Rotator as gid
    gid.GiRaFFEfood_NRPy_Aligned_Rotator()
    desc = "Generate aligned rotator initial test data for GiRaFFEfood_NRPy."
elif initial_data=="ExactWald":
    import GiRaFFEfood_NRPy.GiRaFFEfood_NRPy_Exact_Wald as gid
    gid.GiRaFFEfood_NRPy_Exact_Wald(gammaDD,sks.M,sks.r0,stagger=True)
    desc = "Generate exact Wald initial test data for GiRaFFEfood_NRPy."
else:
    print("Unsupported Initial Data string "+initial_data+"! Supported ID: AllTests, AlfvenWave, FastWave, DegenAlfvenWave, ThreeWaves, FFE_Breakdown, AlignedRotator, or ExactWald")

name = "GiRaFFE_NRPy_initial_data"

values_to_print = [                   lhrh(lhs=gri.gfaccess("out_gfs","AD0"),rhs=gid.AD[0]),                   lhrh(lhs=gri.gfaccess("out_gfs","AD1"),rhs=gid.AD[1]),                   lhrh(lhs=gri.gfaccess("out_gfs","AD2"),rhs=gid.AD[2]),                   lhrh(lhs=gri.gfaccess("auxevol_gfs","ValenciavU0"),rhs=gid.ValenciavU[0]),                   lhrh(lhs=gri.gfaccess("auxevol_gfs","ValenciavU1"),rhs=gid.ValenciavU[1]),                   lhrh(lhs=gri.gfaccess("auxevol_gfs","ValenciavU2"),rhs=gid.ValenciavU[2])                  ]

outCfunction(
    outfile  = os.path.join(out_dir,name+".h"), desc=desc, name=name,
    params   ="const paramstruct *params,REAL *xx[3],REAL *auxevol_gfs,REAL *out_gfs",
    body     = fin.FD_outputC("returnstring",values_to_print,params="outCverbose=False").replace("IDX4","IDX4S"),
    loopopts ="AllPoints,Read_xxs")


# <a id='download'></a>
# 
# ## Step 1.b: Download original `GiRaFFE` files \[Back to [top](#toc)\]
# 
# $$\label{download}$$
# 
# Here, we download the relevant portion of the original `GiRaFFE` code from Bitbucket. 

# In[4]:


# First download the original GiRaFFE source code
import urllib

original_file_url  = [
                      "https://bitbucket.org/zach_etienne/wvuthorns/raw/231af720ccf3f1af50f7cce4a86b410fc8ea2e51/GiRaFFEfood/src/AlfvenWave.cc",
                      "https://bitbucket.org/zach_etienne/wvuthorns/raw/231af720ccf3f1af50f7cce4a86b410fc8ea2e51/GiRaFFEfood/src/FastWave.cc",
                      "https://bitbucket.org/zach_etienne/wvuthorns/raw/b826f6578b0c2b5a43ad9171e65a1b0af88d8b77/GiRaFFEfood/src/DegenAlfvenWave.cc",
                      "https://bitbucket.org/zach_etienne/wvuthorns/raw/b826f6578b0c2b5a43ad9171e65a1b0af88d8b77/GiRaFFEfood/src/ThreeAlfvenWave.cc",
                      "https://bitbucket.org/zach_etienne/wvuthorns/raw/b826f6578b0c2b5a43ad9171e65a1b0af88d8b77/GiRaFFEfood/src/FFEBreakdown.cc",
                      "https://bitbucket.org/zach_etienne/wvuthorns/raw/231af720ccf3f1af50f7cce4a86b410fc8ea2e51/GiRaFFEfood/src/AlignedRotator.cc",
                     ]
original_file_name = [
                      "AlfvenWave.cc",
                      "FastWave.cc",
                      "DegenAlfvenWave.cc",
                      "ThreeAlfvenWave.cc",
                      "FFEBreakdown.cc",
                      "AlignedRotator.cc",
                     ]

for i in range(len(original_file_url)):
    original_file_path = os.path.join(out_dir,original_file_name[i])

    # Then download the original GiRaFFE source code
    # We try it here in a couple of ways in an attempt to keep
    # the code more portable
    try:
        original_file_code = urllib.request.urlopen(original_file_url[i]).read().decode('utf-8')
    except:
        original_file_code = urllib.urlopen(original_file_url[i]).read().decode('utf-8')

    # Write down the file the original GiRaFFE source code
    with open(original_file_path,"w") as file:
        file.write(original_file_code)


# <a id='free_params'></a>
# 
# ## Step 1.c: Output C codes needed for declaring and setting Cparameters; also set `free_parameters.h` \[Back to [top](#toc)\]
# 
# $$\label{free_params}$$
# 
# Based on declared NRPy+ Cparameters, first we generate `declare_Cparameters_struct.h`, `set_Cparameters_default.h`, and `set_Cparameters[-SIMD].h`.
# 
# Then we output `free_parameters.h`, which sets some basic grid parameters as well as the speed limit parameter we need for this function.

# In[5]:


# Step 3.d
# Step 3.d.ii: Set free_parameters.h
with open(os.path.join(out_dir,"free_parameters.h"),"w") as file:
    file.write("""
// Set free-parameter values.

const int NGHOSTS = 3;

// Set free-parameter values for the initial data.
// Override parameter defaults with values based on command line arguments and NGHOSTS.
const int Nx0x1x2 = 5;
params.Nxx0 = Nx0x1x2;
params.Nxx1 = Nx0x1x2;
params.Nxx2 = Nx0x1x2;
params.Nxx_plus_2NGHOSTS0 = params.Nxx0 + 2*NGHOSTS;
params.Nxx_plus_2NGHOSTS1 = params.Nxx1 + 2*NGHOSTS;
params.Nxx_plus_2NGHOSTS2 = params.Nxx2 + 2*NGHOSTS;
// Step 0d: Set up space and time coordinates
// Step 0d.i: Declare \Delta x^i=dxx{0,1,2} and invdxx{0,1,2}, as well as xxmin[3] and xxmax[3]:
const REAL xxmin[3] = {-1.0,-1.0,-1.0};
const REAL xxmax[3] = { 1.0, 1.0, 1.0};

params.dxx0 = (xxmax[0] - xxmin[0]) / ((REAL)params.Nxx0);
params.dxx1 = (xxmax[1] - xxmin[1]) / ((REAL)params.Nxx1);
params.dxx2 = (xxmax[2] - xxmin[2]) / ((REAL)params.Nxx2);
params.invdx0 = 1.0 / params.dxx0;
params.invdx1 = 1.0 / params.dxx1;
params.invdx2 = 1.0 / params.dxx2;
\n""")

# Generates declare_Cparameters_struct.h, set_Cparameters_default.h, and set_Cparameters[-SIMD].h
par.generate_Cparameters_Ccodes(os.path.join(out_dir))


# <a id='interface'></a>
# 
# ## Step 1.d: Create dummy files for the CCTK version of the code \[Back to [top](#toc)\]
# 
# $$\label{interface}$$
# 
# The original `GiRaFFE` code depends on some functionalities of the CCTK. Since we only care about this one small function, we can get around this by creating some nearly-empty, non-functional files that can be included to satisfy the pre-processor without changing functionality. We will later replace what little functionality we need with some basic global variables and macros.

# In[6]:


#incldue "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"
#include "Symmetry.h"
with open(os.path.join(out_dir,"cctk.h"),"w") as file:
    file.write("""//""")

with open(os.path.join(out_dir,"cctk_Arguments.h"),"w") as file:
    file.write("""#define DECLARE_CCTK_ARGUMENTS //
#define CCTK_ARGUMENTS void
""")

with open(os.path.join(out_dir,"cctk_Parameters.h"),"w") as file:
    file.write("""#define DECLARE_CCTK_PARAMETERS //
""")

with open(os.path.join(out_dir,"Symmetry.h"),"w") as file:
    file.write("""//""")


# <a id='mainc'></a>
# 
# # Step 2: `GiRaFFEfood_NRPy_unit_test.c`: The Main C Code \[Back to [top](#toc)\]
# 
# $$\label{mainc}$$
# 
# Now that we have our vector potential and analytic magnetic field to compare against, we will start writing our unit test. We'll also import common C functionality, define `REAL`, the number of ghost zones, and the faces, and set the standard macros for NRPy+ style memory access.

# In[7]:


get_ipython().run_cell_magic('writefile', '$out_dir/GiRaFFEfood_NRPy_unit_test.C', '\n// These are common packages that we are likely to need.\n#include "stdio.h"\n#include "stdlib.h"\n#include "math.h"\n#include <string> // Needed for strncmp, etc.\n#include "stdint.h" // Needed for Windows GCC 6.x compatibility\n#include <time.h>   // Needed to set a random seed.\n\n#define REAL double\n#include "declare_Cparameters_struct.h"\n\n// Standard NRPy+ memory access:\n#define IDX4S(g,i,j,k) \\\n( (i) + Nxx_plus_2NGHOSTS0 * ( (j) + Nxx_plus_2NGHOSTS1 * ( (k) + Nxx_plus_2NGHOSTS2 * (g) ) ) )\n\n// Standard formula to calculate significant digits of agreement:\n#define SDA(a,b) 1.0-log10(2.0*fabs(a-b)/(fabs(a)+fabs(b)))\n\n// Memory access definitions for NRPy+\n#define BU0GF 0\n#define BU1GF 1\n#define BU2GF 2\n#define VALENCIAVU0GF 3\n#define VALENCIAVU1GF 4\n#define VALENCIAVU2GF 5\n#define BETAU0GF 6\n#define BETAU1GF 7\n#define BETAU2GF 8\n#define ALPHAGF 9\n#define NUM_AUXEVOL_GFS 10\n\n#define AD0GF 0\n#define AD1GF 1\n#define AD2GF 2\n#define NUM_EVOL_GFS 3\n\n// Include the functions that we want to test:\n#include "GiRaFFE_NRPy_initial_data.h"\n\n// Define CCTK macros\n#define CCTK_REAL double\n#define CCTK_INT int\nstruct cGH{};\nconst cGH* cctkGH;\n\n// GiRaFFE parameters in ETK\nconst CCTK_REAL min_radius_inside_of_which_conserv_to_prims_FFE_and_FFE_evolution_is_DISABLED = -1;\nconst int current_sheet_null_v = 1;\n\n// More definitions to interface with ETK code:\nconst int cctk_lsh[3] = {11,11,11};\nconst int grid_size = cctk_lsh[0]*cctk_lsh[1]*cctk_lsh[2];\nCCTK_REAL Avec[3*grid_size];\nCCTK_REAL vel[3*grid_size];\nCCTK_REAL Ax[grid_size];\nCCTK_REAL Ay[grid_size];\nCCTK_REAL Az[grid_size];\nCCTK_REAL vx[grid_size];\nCCTK_REAL vy[grid_size];\nCCTK_REAL vz[grid_size];\nCCTK_REAL Bx[grid_size];\nCCTK_REAL By[grid_size];\nCCTK_REAL Bz[grid_size];\nCCTK_REAL x[grid_size];\nCCTK_REAL y[grid_size];\nCCTK_REAL z[grid_size];\nCCTK_REAL r[grid_size];\nCCTK_REAL *alp;\nCCTK_REAL *betax;\nCCTK_REAL *betay;\nCCTK_REAL *betaz;\n\n// We need to declare these to compile functions we won\'t call:\nint Compute_Exact_Every;\nint cctk_iteration;\nCCTK_REAL *delpsi6phi;\nCCTK_REAL *psi6phi;\nCCTK_REAL *delAx;\nCCTK_REAL *delAy;\nCCTK_REAL *delAz;\nCCTK_REAL *exactBx;\nCCTK_REAL *exactBy;\nCCTK_REAL *exactBz;\nCCTK_REAL *delBx;\nCCTK_REAL *delBy;\nCCTK_REAL *delBz;\nCCTK_REAL *exactVx;\nCCTK_REAL *exactVy;\nCCTK_REAL *exactVz;\nCCTK_REAL *delvx;\nCCTK_REAL *delvy;\nCCTK_REAL *delvz;\nCCTK_REAL cctk_time = 0.0;\nCCTK_REAL *exactBx_ThreeWaves;\nCCTK_REAL *exactBy_ThreeWaves;\nCCTK_REAL *exactBz_ThreeWaves;\nCCTK_REAL *delBx_ThreeWaves;\nCCTK_REAL *delBy_ThreeWaves;\nCCTK_REAL *delBz_ThreeWaves;\nCCTK_REAL *mhd_st_x;\nCCTK_REAL *mhd_st_y;\nCCTK_REAL *mhd_st_z;\nCCTK_REAL *Ex;\nCCTK_REAL *Ey;\nCCTK_REAL *Ez;\nCCTK_REAL *B2mE2;\n\n// Set constants to default for comparison\nCCTK_REAL wave_speed = -0.5;\nCCTK_REAL Omega_aligned_rotator = 1e3;\nCCTK_REAL R_NS_aligned_rotator = 1.0;\nCCTK_REAL B_p_aligned_rotator = 1e-5;\nCCTK_REAL Wald_B0 = 1.0;\n\n// Define dz in CCTK\nCCTK_REAL cactus_dxx[3];\n#define CCTK_DELTA_SPACE(i) cactus_dxx[i]\n\n// Dummy ETK function:\n#define CCTK_GFINDEX3D(cctkGH,i,j,k) (i) + cctk_lsh[0] * ( (j) + cctk_lsh[1] * (k) )\n#define CCTK_GFINDEX4D(cctkGH,i,j,k,g) \\\n( (i) + cctk_lsh[0] * ( (j) + cctk_lsh[1] * ( (k) + cctk_lsh[2] * (g) ) ) )\n#define CCTK_VInfo(...) //\n//#define CCTK_VWarn(...) //\n\n#include "AlfvenWave.cc"\n#include "FastWave.cc"\n#include "AlignedRotator.cc"\n#include "DegenAlfvenWave.cc"\n#include "ThreeAlfvenWave.cc"\n#include "FFEBreakdown.cc"\n\nint main(int argc, char** argv) {\n    paramstruct params;\n#include "set_Cparameters_default.h"\n\n    // Step 0c: Set free parameters, overwriting Cparameters defaults\n    //          by hand or with command-line input, as desired.\n#include "free_parameters.h"\n#include "set_Cparameters-nopointer.h"\n\n    // Now that we\'ve calculated dxx2,  we can define a cactus equivalent\n    cactus_dxx[0] = dxx0;\n    cactus_dxx[1] = dxx1;\n    cactus_dxx[2] = dxx2;\n\n    // Step 0d.ii: Set up uniform coordinate grids\n    REAL *xx[3];\n    xx[0] = (REAL *)malloc(sizeof(REAL)*Nxx_plus_2NGHOSTS0);\n    xx[1] = (REAL *)malloc(sizeof(REAL)*Nxx_plus_2NGHOSTS1);\n    xx[2] = (REAL *)malloc(sizeof(REAL)*Nxx_plus_2NGHOSTS2);\n    for(int j=0;j<Nxx_plus_2NGHOSTS0;j++) xx[0][j] = xxmin[0] + (j-NGHOSTS)*dxx0;\n    for(int j=0;j<Nxx_plus_2NGHOSTS1;j++) xx[1][j] = xxmin[1] + (j-NGHOSTS)*dxx1;\n    for(int j=0;j<Nxx_plus_2NGHOSTS2;j++) xx[2][j] = xxmin[2] + (j-NGHOSTS)*dxx2;\n\n    for(int k=0;k<Nxx_plus_2NGHOSTS2;k++)\n        for(int j=0;j<Nxx_plus_2NGHOSTS1;j++)\n            for(int i=0;i<Nxx_plus_2NGHOSTS0;i++) {\n                int index = CCTK_GFINDEX3D(cctkGH,i,j,k);\n                x[index] = xx[0][i];\n                y[index] = xx[1][j];\n                z[index] = xx[2][k];\n                r[index] = sqrt(x[index]*x[index] + y[index]*y[index] + z[index]*z[index]);\n    }\n\n    //for(int j=0;j<Nxx_plus_2NGHOSTS0;j++) printf("x[%d] = %.5e\\n",j,xx[0][j]);\n\n    // This is the array to which we\'ll write the NRPy+ variables.\n    REAL *auxevol_gfs  = (REAL *)malloc(sizeof(REAL) * NUM_AUXEVOL_GFS * Nxx_plus_2NGHOSTS2 * Nxx_plus_2NGHOSTS1 * Nxx_plus_2NGHOSTS0);\n    REAL *evol_gfs  = (REAL *)malloc(sizeof(REAL) * NUM_EVOL_GFS * Nxx_plus_2NGHOSTS2 * Nxx_plus_2NGHOSTS1 * Nxx_plus_2NGHOSTS0);\n\n    GiRaFFE_NRPy_initial_data(&params,xx,auxevol_gfs,evol_gfs);\n    if(atoi(argv[1])==0) GiRaFFEfood_AlfvenWave();\n    else if(atoi(argv[1])==1) GiRaFFEfood_AlignedRotator();\n    else if(atoi(argv[1])==2) GiRaFFEfood_FastWave();\n    else if(atoi(argv[1])==3) GiRaFFEfood_DegenAlfvenWave();\n    else if(atoi(argv[1])==4) GiRaFFEfood_ThreeAlfvenWave();\n    else if(atoi(argv[1])==5) GiRaFFEfood_FFEBreakdown();\n\n    int all_agree = 1;\n\n    for(int i0=0;i0<Nxx_plus_2NGHOSTS0;i0++){\n        for(int i1=0;i1<Nxx_plus_2NGHOSTS1;i1++){\n            for(int i2=0;i2<Nxx_plus_2NGHOSTS2;i2++){\n                if(SDA(auxevol_gfs[IDX4S(VALENCIAVU0GF, i0,i1,i2)],vel[CCTK_GFINDEX4D(cctkGH,i0,i1,i2,0)])<10.0){\n                    printf("Quantity ValenciavU0 only agrees with the original GiRaFFE to %.2f digits at i0,i1,i2=%d,%d,%d!\\n",\n                           SDA(auxevol_gfs[IDX4S(ALPHAGF, i0,i1,i2)]*auxevol_gfs[IDX4S(VALENCIAVU0GF, i0,i1,i2)]-auxevol_gfs[IDX4S(BETAU0GF, i0,i1,i2)],vel[CCTK_GFINDEX4D(cctkGH,i0,i1,i2,0)]),i0,i1,i2);\n                    all_agree=0;\n                }\n                if(SDA(auxevol_gfs[IDX4S(VALENCIAVU1GF, i0,i1,i2)],vel[CCTK_GFINDEX4D(cctkGH,i0,i1,i2,1)])<10.0){\n                    printf("Quantity ValenciavU1 only agrees with the original GiRaFFE to %.2f digits at i0,i1,i2=%d,%d,%d!\\n",\n                           SDA(auxevol_gfs[IDX4S(ALPHAGF, i0,i1,i2)]*auxevol_gfs[IDX4S(VALENCIAVU1GF, i0,i1,i2)]-auxevol_gfs[IDX4S(BETAU1GF, i0,i1,i2)],vel[CCTK_GFINDEX4D(cctkGH,i0,i1,i2,1)]),i0,i1,i2);\n                    all_agree=0;\n                }\n                if(SDA(auxevol_gfs[IDX4S(VALENCIAVU2GF, i0,i1,i2)],vel[CCTK_GFINDEX4D(cctkGH,i0,i1,i2,2)])<10.0){\n                    printf("Quantity ValenciavU2 only agrees with the original GiRaFFE to %.2f digits at i0,i1,i2=%d,%d,%d!\\n",\n                           SDA(auxevol_gfs[IDX4S(ALPHAGF, i0,i1,i2)]*auxevol_gfs[IDX4S(VALENCIAVU2GF, i0,i1,i2)]-auxevol_gfs[IDX4S(BETAU2GF, i0,i1,i2)],vel[CCTK_GFINDEX4D(cctkGH,i0,i1,i2,2)]),i0,i1,i2);\n                    all_agree=0;\n                }\n                //printf("NRPy: %.15e,%.15e,%.15e\\n",auxevol_gfs[IDX4S(VALENCIAVU0GF, i0,i1,i2)],auxevol_gfs[IDX4S(VALENCIAVU1GF, i0,i1,i2)],auxevol_gfs[IDX4S(VALENCIAVU2GF, i0,i1,i2)]);\n                //printf("CCTK: %.15e,%.15e,%.15e\\n",vel[CCTK_GFINDEX4D(cctkGH,i0,i1,i2,0)],vel[CCTK_GFINDEX4D(cctkGH,i0,i1,i2,1)],vel[CCTK_GFINDEX4D(cctkGH,i0,i1,i2,2)]);\n            }\n        }\n    }\n\n    // Shift the grid to compare A_x\n    for(int j=0;j<Nxx_plus_2NGHOSTS1;j++) xx[1][j] += 0.5*dxx1;\n    for(int j=0;j<Nxx_plus_2NGHOSTS2;j++) xx[2][j] += 0.5*dxx2;\n\n    GiRaFFE_NRPy_initial_data(&params,xx,auxevol_gfs,evol_gfs);\n    if(atoi(argv[1])==0) GiRaFFEfood_AlfvenWave();\n    else if(atoi(argv[1])==1) GiRaFFEfood_AlignedRotator();\n    else if(atoi(argv[1])==2) GiRaFFEfood_FastWave();\n    else if(atoi(argv[1])==3) GiRaFFEfood_DegenAlfvenWave();\n    else if(atoi(argv[1])==4) GiRaFFEfood_ThreeAlfvenWave();\n    else if(atoi(argv[1])==5) GiRaFFEfood_FFEBreakdown();\n    for(int i0=0;i0<Nxx_plus_2NGHOSTS0;i0++){\n        for(int i1=0;i1<Nxx_plus_2NGHOSTS1;i1++){\n            for(int i2=0;i2<Nxx_plus_2NGHOSTS2;i2++){\n                if(SDA(evol_gfs[IDX4S(AD0GF, i0,i1,i2)],Avec[CCTK_GFINDEX4D(cctkGH,i0,i1,i2,0)])<10.0){\n                    printf("Quantity AD0 only agrees with the original GiRaFFE to %.2f digits at i0,i1,i2=%d,%d,%d!\\n",\n                           SDA(evol_gfs[IDX4S(AD0GF, i0,i1,i2)],Avec[CCTK_GFINDEX4D(cctkGH,i0,i1,i2,0)]),i0,i1,i2);\n                    all_agree=0;\n                }\n            }\n        }\n    }\n\n    // Shift the grid to compare A_y\n    for(int j=0;j<Nxx_plus_2NGHOSTS0;j++) xx[0][j] += 0.5*dxx0;\n    for(int j=0;j<Nxx_plus_2NGHOSTS1;j++) xx[1][j] -= 0.5*dxx1;\n    GiRaFFE_NRPy_initial_data(&params,xx,auxevol_gfs,evol_gfs);\n    if(atoi(argv[1])==0) GiRaFFEfood_AlfvenWave();\n    else if(atoi(argv[1])==1) GiRaFFEfood_AlignedRotator();\n    else if(atoi(argv[1])==2) GiRaFFEfood_FastWave();\n    else if(atoi(argv[1])==3) GiRaFFEfood_DegenAlfvenWave();\n    else if(atoi(argv[1])==4) GiRaFFEfood_ThreeAlfvenWave();\n    else if(atoi(argv[1])==5) GiRaFFEfood_FFEBreakdown();\n    for(int i0=0;i0<Nxx_plus_2NGHOSTS0;i0++){\n        for(int i1=0;i1<Nxx_plus_2NGHOSTS1;i1++){\n            for(int i2=0;i2<Nxx_plus_2NGHOSTS2;i2++){\n                if(SDA(evol_gfs[IDX4S(AD1GF, i0,i1,i2)],Avec[CCTK_GFINDEX4D(cctkGH,i0,i1,i2,1)])<10.0){\n                    printf("Quantity AD1 only agrees with the original GiRaFFE to %.2f digits at i0,i1,i2=%d,%d,%d!\\n",\n                           SDA(evol_gfs[IDX4S(AD1GF, i0,i1,i2)],Avec[CCTK_GFINDEX4D(cctkGH,i0,i1,i2,1)]),i0,i1,i2);\n                    all_agree=0;\n                printf("NRPy: %.15e\\n",evol_gfs[IDX4S(AD1GF, i0,i1,i2)]);\n                printf("CCTK: %.15e\\n\\n",Avec[CCTK_GFINDEX4D(cctkGH,i0,i1,i2,1)]);\n                }\n            }\n        }\n    }\n\n    // Shift the grid to compare A_z\n    for(int j=0;j<Nxx_plus_2NGHOSTS1;j++) xx[1][j] += 0.5*dxx1;\n    for(int j=0;j<Nxx_plus_2NGHOSTS2;j++) xx[2][j] -= 0.5*dxx2;\n    GiRaFFE_NRPy_initial_data(&params,xx,auxevol_gfs,evol_gfs);\n    if(atoi(argv[1])==0) GiRaFFEfood_AlfvenWave();\n    else if(atoi(argv[1])==1) GiRaFFEfood_AlignedRotator();\n    else if(atoi(argv[1])==2) GiRaFFEfood_FastWave();\n    else if(atoi(argv[1])==3) GiRaFFEfood_DegenAlfvenWave();\n    else if(atoi(argv[1])==4) GiRaFFEfood_ThreeAlfvenWave();\n    else if(atoi(argv[1])==5) GiRaFFEfood_FFEBreakdown();\n    for(int i0=0;i0<Nxx_plus_2NGHOSTS0;i0++){\n        for(int i1=0;i1<Nxx_plus_2NGHOSTS1;i1++){\n            for(int i2=0;i2<Nxx_plus_2NGHOSTS2;i2++){\n                if(SDA(evol_gfs[IDX4S(AD2GF, i0,i1,i2)],Avec[CCTK_GFINDEX4D(cctkGH,i0,i1,i2,2)])<10.0){\n                    printf("Quantity AD2 only agrees with the original GiRaFFE to %.2f digits at i0,i1,i2=%d,%d,%d!\\n",\n                           SDA(evol_gfs[IDX4S(AD2GF, i0,i1,i2)],Avec[CCTK_GFINDEX4D(cctkGH,i0,i1,i2,2)]),i0,i1,i2);\n                    all_agree=0;\n                }\n            }\n        }\n    }\n\n    //printf("NRPy: %.15e,%.15e,%.15e\\n",evol_gfs[IDX4S(AD0GF, i0,i1,i2)],evol_gfs[IDX4S(AD1GF, i0,i1,i2)],evol_gfs[IDX4S(AD2GF, i0,i1,i2)]);\n    //printf("CCTK: %.15e,%.15e,%.15e\\n",Avec[CCTK_GFINDEX4D(cctkGH,i0,i1,i2,0)],Avec[CCTK_GFINDEX4D(cctkGH,i0,i1,i2,1)],Avec[CCTK_GFINDEX4D(cctkGH,i0,i1,i2,2)]);\n\n    if(all_agree) printf("All quantities agree at all points!\\n");\n}')


# <a id='compile_run'></a>
# 
# ## Step 2.a: Compile and run the code to validate the output \[Back to [top](#toc)\]
# 
# $$\label{compile_run}$$
# 
# Finally, we can compile and run the code we have written. Once run, this code will output the level of agreement between the two codes and some information to help interpret those numbers.

# In[8]:


import time

print("Now compiling, should take ~2 seconds...\n")
start = time.time()
# cmd.C_compile(os.path.join(out_dir,"GiRaFFEfood_NRPy_unit_test.C"), os.path.join(out_dir,"GiRaFFEfood_NRPy_unit_test"))
get_ipython().system('g++ -Ofast -fopenmp -march=native -funroll-loops Validation/GiRaFFEfood_NRPy_unit_test.C -o Validation/GiRaFFEfood_NRPy_unit_test -lstdc++')
end = time.time()
print("Finished in "+str(end-start)+" seconds.\n\n")

results_file = "out_GiRaFFEfood_NRPy_test.txt"

# os.chdir(out_dir)
os.chdir(out_dir)
# cmd.Execute(os.path.join("GiRaFFEfood_NRPy_unit_test"))
if initial_data=="AlfvenWave":
    cmd.Execute("GiRaFFEfood_NRPy_unit_test","0",results_file)
elif initial_data=="AlignedRotator":
    cmd.Execute("GiRaFFEfood_NRPy_unit_test","1",results_file)
elif initial_data=="FastWave":
    cmd.Execute("GiRaFFEfood_NRPy_unit_test","2",results_file)
elif initial_data=="DegenAlfvenWave":
    cmd.Execute("GiRaFFEfood_NRPy_unit_test","3",results_file)
elif initial_data=="ThreeWaves":
    cmd.Execute("GiRaFFEfood_NRPy_unit_test","4",results_file)
elif initial_data=="FFE_Breakdown":
    cmd.Execute("GiRaFFEfood_NRPy_unit_test","5",results_file)
os.chdir(os.path.join("../"))


# Here, we add some emergency brakes so that if the output from the test isn't good, we throw an error to stop the notebook dead in its tracks. This way, our automatic testing infrastructure can let us know if something goes wrong. We will also print the output from the test for convenience's sake.

# In[9]:


with open(os.path.join(out_dir,results_file),"r") as file:
    output = file.readline()
    print(output)
    if output!="All quantities agree at all points!\n": # If this isn't the first line of this file, something went wrong!
        sys.exit(1)


# <a id='latex_pdf_output'></a>
# 
# # Step 4: Output this notebook to $\LaTeX$-formatted PDF file \[Back to [top](#toc)\]
# $$\label{latex_pdf_output}$$
# 
# The following code cell converts this Jupyter notebook into a proper, clickable $\LaTeX$-formatted PDF file. After the cell is successfully run, the generated PDF may be found in the root NRPy+ tutorial directory, with filename
# [Tutorial-Start_to_Finish_UnitTest-GiRaFFEfood_NRPy.pdf](Tutorial-Start_to_Finish_UnitTest-GiRaFFEfood_NRPy.pdf) (Note that clicking on this link may not work; you may need to open the PDF file through another means.)

# In[10]:


import cmdline_helper as cmd    # NRPy+: Multi-platform Python command-line interface
cmd.output_Jupyter_notebook_to_LaTeXed_PDF("Tutorial-Start_to_Finish_UnitTest-GiRaFFEfood_NRPy")

