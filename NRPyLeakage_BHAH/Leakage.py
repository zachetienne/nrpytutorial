import shutil, os, sys                # Standard Python modules for multiplatform OS-level functions

# Step P1: Import needed NRPy+ core modules:
from outputC import outputC, lhrh, outCfunction, add_to_Cfunction_dict, outC_function_master_list, outC_function_dict
# NRPy+: Core C code output module
from pickling import unpickle_NRPy_env  # NRPy+: Pickle/unpickle NRPy+ environment, for parallel codegen
import finite_difference as fin       # NRPy+: Finite difference C code generation module
import NRPy_param_funcs as par        # NRPy+: Parameter interface
import grid as gri                    # NRPy+: Functions having to do with numerical grids
import reference_metric as rfm        # NRPy+: Reference metric support
import indexedexp as ixp         # NRPy+: Symbolic indexed expression (e.g., tensors, vectors, etc.) support
import cmdline_helper as cmd          # NRPy+: Multi-platform Python command-line interface
import sympy as sp

thismodule = __name__

# restricted:
# REAL tauChangeThreshold "Threshold to stop PathOfLeastResistance algorithm during initialization." STEERABLE = ALWAYS
# {
#   *:* :: "Anything goes"
# } 1e-8
tauChangeThreshold = par.Cparameters("REAL", thismodule, ["tauChangeThreshold"], 1e-8)

# restricted:
# REAL psi6_threshold "If psi^6 > psi6_threshold we assume to be insidee the horizon and do not leak." STEERABLE = ALWAYS
# {
#   *:* :: "Anything goes"
# } 1e100
psi6_threshold = par.Cparameters("REAL", thismodule, ["psi6_threshold"], 1e100)

# restricted:
# REAL rho_min_threshold "Only leak for rho > rho_min_threshold" STEERABLE = ALWAYS
# {
#   0:* :: "Any non-negative number"
# } 0.0
rho_min_threshold = par.Cparameters("REAL", thismodule, ["rho_min_threshold"], 0.0)

# restricted:
# REAL rho_max_threshold "Only leak for rho < rho_max_threshold" STEERABLE = ALWAYS
# {
#   0:* :: "Any non-negative number"
# } 1e100
rho_max_threshold = par.Cparameters("REAL", thismodule, ["rho_max_threshold"], 1e100)

# restricted:
# CCTK_INT max_iterations "Maximum number of iterations to initialize the optical depths"
# {
#   1:* :: "Should be greater than one"
# } 2048
max_iterations = par.Cparameters("int", thismodule, ["max_iterations"], 2048)


CoordSystem     = "Cartesian"
par.set_parval_from_str("reference_metric::CoordSystem", CoordSystem)
rfm.reference_metric()

import BSSN.BSSN_quantities as Bq
Bq.BSSN_basic_tensors()

Bq.phi_and_derivs()

psi6 = Bq.exp_m4phi**(-0.5*3)

psi6_str = outputC(psi6, "const REAL psi6L", filename='returnstring', params="includebraces=false")


def add_to_Cfunction_dict_NRPyLeakageET_Initialize():
        
  name = "NRPyLeakageET_Initialize"
  includes = ["NRPy_basic_defines.h", "NRPy_function_prototypes.h"]
  c_type = "void"
  desc = "initialize opacities and optical depths"
  params = "const paramstruct *restrict params, const ghl_eos_parameters *restrict ghl_eos, const REAL *restrict evol_gfs, REAL *restrict auxevol_gfs"

  body = r"""
// Step 0: set pointers to GFs not in use at the moment in initialization

const int Nxx_plus_2NGHOSTS_TOT = Nxx_plus_2NGHOSTS0*Nxx_plus_2NGHOSTS1*Nxx_plus_2NGHOSTS2;

REAL *tau_0_nue_aux  = auxevol_gfs + Nxx_plus_2NGHOSTS_TOT*T4UU00GF;
REAL *tau_1_nue_aux  = auxevol_gfs + Nxx_plus_2NGHOSTS_TOT*T4UU01GF;
REAL *tau_0_anue_aux = auxevol_gfs + Nxx_plus_2NGHOSTS_TOT*T4UU02GF;
REAL *tau_1_anue_aux = auxevol_gfs + Nxx_plus_2NGHOSTS_TOT*T4UU03GF;
REAL *tau_0_nux_aux  = auxevol_gfs + Nxx_plus_2NGHOSTS_TOT*T4UU11GF;
REAL *tau_1_nux_aux  = auxevol_gfs + Nxx_plus_2NGHOSTS_TOT*T4UU12GF;

REAL *tau_0_nue_p  = auxevol_gfs + Nxx_plus_2NGHOSTS_TOT*KDD00GF;
REAL *tau_1_nue_p  = auxevol_gfs + Nxx_plus_2NGHOSTS_TOT*KDD01GF;
REAL *tau_0_anue_p = auxevol_gfs + Nxx_plus_2NGHOSTS_TOT*KDD02GF;
REAL *tau_1_anue_p = auxevol_gfs + Nxx_plus_2NGHOSTS_TOT*KDD11GF;
REAL *tau_0_nux_p  = auxevol_gfs + Nxx_plus_2NGHOSTS_TOT*KDD12GF;
REAL *tau_1_nux_p  = auxevol_gfs + Nxx_plus_2NGHOSTS_TOT*KDD22GF;

REAL *kappa_0_nue_p  = auxevol_gfs + Nxx_plus_2NGHOSTS_TOT*RBARDD00GF;
REAL *kappa_1_nue_p  = auxevol_gfs + Nxx_plus_2NGHOSTS_TOT*RBARDD01GF;
REAL *kappa_0_anue_p = auxevol_gfs + Nxx_plus_2NGHOSTS_TOT*RBARDD02GF;
REAL *kappa_1_anue_p = auxevol_gfs + Nxx_plus_2NGHOSTS_TOT*RBARDD11GF;
REAL *kappa_0_nux_p  = auxevol_gfs + Nxx_plus_2NGHOSTS_TOT*RBARDD12GF;
REAL *kappa_1_nux_p  = auxevol_gfs + Nxx_plus_2NGHOSTS_TOT*RBARDD22GF;

// Step 1: Initialize all optical depth gridfunctions to zero
NRPyLeakageET_optical_depths_initialize_to_zero(params, auxevol_gfs,
                                                tau_0_nue_aux,  tau_1_nue_aux,
                                                tau_0_anue_aux, tau_1_anue_aux,
                                                tau_0_nux_aux,  tau_1_nux_aux,
                                                tau_0_nue_p,    tau_1_nue_p,
                                                tau_0_anue_p,   tau_1_anue_p,
                                                tau_0_nux_p,    tau_1_nux_p,
                                                kappa_0_nue_p,  kappa_1_nue_p,
                                                kappa_0_anue_p, kappa_1_anue_p,
                                                kappa_0_nux_p,  kappa_1_nux_p);

// Step 2: Now perform iterations of the path of least resistance algorithm
int counter = 0;
int RemainingIterations = max_iterations;

for(int i=1;i<=max_iterations;i++) {
  // Step 2.a: First perform the iteration on refinement level startRefLev
  printf("Beginning iteration %d...\n",counter+1);

  // Step 2.b: Now loop over remaining refinement levels, performing the
  //           path of least resistance (POLR) algorithm
  // Step 2.b.i: Compute opacities
  NRPyLeakageET_compute_neutrino_opacities(params, ghl_eos, 
                                           evol_gfs, auxevol_gfs, 
                                           tau_0_nue_p,  tau_1_nue_p, 
                                           tau_0_anue_p, tau_1_anue_p, 
                                           tau_0_nux_p,  tau_1_nux_p);
  // Step 2.b.ii: Perform POLR algorithm on every grid component at this level
  NRPyLeakageET_optical_depths_PathOfLeastResistance(params, evol_gfs, auxevol_gfs, 
                                                     tau_0_nue_p,  tau_1_nue_p, 
                                                     tau_0_anue_p, tau_1_anue_p, 
                                                     tau_0_nux_p,  tau_1_nux_p);

  // Step 2.c: Perform a reduction of the L2-norm among all processors
  REAL l2norm_global;
  NRPyLeakageET_getGlobalL2norm(params, evol_gfs, auxevol_gfs, 
                                tau_0_nue_aux,  tau_1_nue_aux, 
                                tau_0_anue_aux, tau_1_anue_aux, 
                                tau_0_nux_aux,  tau_1_nux_aux,
                                &l2norm_global);

  // Step 2.d: Increment counter; check convergence
  counter++;
  if( l2norm_global < tauChangeThreshold ) {
    // Step 2.d.i: Converged!
    RemainingIterations = 0;
    i = max_iterations*10;
    printf("Algorithm converged! l2norm_tau_change = %e (threshold = %e)\n",l2norm_global,tauChangeThreshold);
  }
  else {
    // Step 2.d.ii: Have not converged. Restrict current solution from the
    // finest to the coarsest initialized refinement levels; Decrement
    // number of iterations remaining; continue.
    if( RemainingIterations > 1 ) {
      RemainingIterations--;
      printf("Completed iteration %d. Global change: %e. Remaining iterations: %4d\n",counter,l2norm_global,RemainingIterations);
    }
    else {
      printf("Maximum iterations reached, but convergence threshold not met.\n");
    }
  }
}
  """

  add_to_Cfunction_dict(
                      includes=includes,
                      desc=desc,
                      c_type=c_type, name=name, params=params,
                      body=body,
                      rel_path_to_Cparams=os.path.join("."))





def add_to_Cfunction_dict_NRPyLeakageET_optical_depths_initialize_to_zero():
    
  includes = ["NRPy_basic_defines.h", "NRPy_function_prototypes.h"]
  c_type = "void"
  name = "NRPyLeakageET_optical_depths_initialize_to_zero"
  desc = "initialize GFs to 0.0"
  params = r"""const paramstruct *restrict params, REAL *restrict auxevol_gfs,
          REAL *restrict tau_0_nue_aux,
          REAL *restrict tau_1_nue_aux,
          REAL *restrict tau_0_anue_aux,
          REAL *restrict tau_1_anue_aux,
          REAL *restrict tau_0_nux_aux,
          REAL *restrict tau_1_nux_aux,
          REAL *restrict tau_0_nue_p,
          REAL *restrict tau_1_nue_p,
          REAL *restrict tau_0_anue_p,
          REAL *restrict tau_1_anue_p,
          REAL *restrict tau_0_nux_p,
          REAL *restrict tau_1_nux_p,
          REAL *restrict kappa_0_nue_p,
          REAL *restrict kappa_1_nue_p,
          REAL *restrict kappa_0_anue_p,
          REAL *restrict kappa_1_anue_p,
          REAL *restrict kappa_0_nux_p,
          REAL *restrict kappa_1_nux_p"""
  
  body = r"""  
// Step 1: Initialize all time levels of the optical
//         depth gridfunctions to zero.
const int index = IDX3S(i0, i1, i2);
auxevol_gfs[IDX4ptS(TAU_0_NUEGF,  index)] = auxevol_gfs[IDX4ptS(KAPPA_0_NUEGF,  index)] = 0.0;
auxevol_gfs[IDX4ptS(TAU_1_NUEGF,  index)] = auxevol_gfs[IDX4ptS(KAPPA_1_NUEGF,  index)] = 0.0;
auxevol_gfs[IDX4ptS(TAU_0_ANUEGF, index)] = auxevol_gfs[IDX4ptS(KAPPA_0_ANUEGF, index)] = 0.0;
auxevol_gfs[IDX4ptS(TAU_1_ANUEGF, index)] = auxevol_gfs[IDX4ptS(KAPPA_1_ANUEGF, index)] = 0.0;
auxevol_gfs[IDX4ptS(TAU_0_NUXGF,  index)] = auxevol_gfs[IDX4ptS(KAPPA_0_NUXGF,  index)] = 0.0;
auxevol_gfs[IDX4ptS(TAU_1_NUXGF,  index)] = auxevol_gfs[IDX4ptS(KAPPA_1_NUXGF,  index)] = 0.0;

tau_0_nue_p [index] = kappa_0_nue_p [index] = 0.0;
tau_1_nue_p [index] = kappa_1_nue_p [index] = 0.0;
tau_0_anue_p[index] = kappa_0_anue_p[index] = 0.0;
tau_1_anue_p[index] = kappa_1_anue_p[index] = 0.0;
tau_0_nux_p [index] = kappa_0_nux_p [index] = 0.0;
tau_1_nux_p [index] = kappa_1_nux_p [index] = 0.0;

tau_0_nue_aux [index] = 0.0;
tau_1_nue_aux [index] = 0.0;
tau_0_anue_aux[index] = 0.0;
tau_1_anue_aux[index] = 0.0;
tau_0_nux_aux [index] = 0.0;
tau_1_nux_aux [index] = 0.0;
"""
  
  add_to_Cfunction_dict(
                  includes=includes,
                  desc=desc,
                  c_type=c_type, name=name, params=params,
                  body=body,
                  rel_path_to_Cparams=os.path.join("."), 
                  loopopts="AllPoints")


def add_to_Cfunction_dict_NRPyLeakageET_compute_neutrino_opacities():
        
  includes = ["NRPy_basic_defines.h", "NRPy_function_prototypes.h"]
  c_type = "void"
  name = "NRPyLeakageET_compute_neutrino_opacities"
  desc = r"""
(c) Leo Werneck
Compute opacities
https://adsabs.harvard.edu/pdf/1996A%26A...311..532R
"""
  params = r"""const paramstruct *restrict params, const ghl_eos_parameters *restrict ghl_eos,
  const REAL *restrict evol_gfs, REAL *restrict auxevol_gfs, 
  REAL *restrict tau_0_nue_p,
  REAL *restrict tau_1_nue_p, 
  REAL *restrict tau_0_anue_p, 
  REAL *restrict tau_1_anue_p, 
  REAL *restrict tau_0_nux_p, 
  REAL *restrict tau_1_nux_p"""
  
  
  body = r"""
  
// Step 1: Set local gridfunction index
const int index = IDX3S(i0, i1, i2);

// Step 2: Declare variables
const REAL rhoL = auxevol_gfs[IDX4ptS(RHOBGF, index)];

ghl_neutrino_opacities kappa;
ghl_neutrino_optical_depths tau;

// Step 3: Check density threshold, compute opacities
if( rhoL < rho_min_threshold || rhoL > rho_max_threshold ) {
  // Step 3.a: Below density threshold; set opacities and optical
  // depths to zero
  kappa.nue[0] = kappa.anue[0] = kappa.nux[0] = 0.0;
  kappa.nue[1] = kappa.anue[1] = kappa.nux[1] = 0.0;
  tau.nue  [0] = tau.anue  [0] = tau.nux  [0] = 0.0;
  tau.nue  [1] = tau.anue  [1] = tau.nux  [1] = 0.0;
}
else {
const REAL cf = evol_gfs[IDX4ptS(CFGF, index)];
"""+psi6_str+r"""
  if( psi6L > psi6_threshold ) {
    kappa.nue[0] = kappa.anue[0] = kappa.nux[0] = 0.0;
    kappa.nue[1] = kappa.anue[1] = kappa.nux[1] = 0.0;
    tau.nue  [0] = tau.anue  [0] = tau.nux  [0] = 0.0;
    tau.nue  [1] = tau.anue  [1] = tau.nux  [1] = 0.0;
  }
  else {
    // Step 3.b: Above density threshold, compute opacities
    // Step 3.b.i: Read from main memory
    const REAL Y_eL         = auxevol_gfs[IDX4ptS(YEGF, index)];
    const REAL temperatureL = auxevol_gfs[IDX4ptS(TEMPERATUREGF, index)];
    tau.nue [0]             = auxevol_gfs[IDX4ptS(TAU_0_NUEGF,  index)];
    tau.nue [1]             = auxevol_gfs[IDX4ptS(TAU_1_NUEGF,  index)];
    tau.anue[0]             = auxevol_gfs[IDX4ptS(TAU_0_ANUEGF, index)];
    tau.anue[1]             = auxevol_gfs[IDX4ptS(TAU_1_ANUEGF, index)];
    tau.nux [0]             = auxevol_gfs[IDX4ptS(TAU_0_NUXGF,  index)];
    tau.nux [1]             = auxevol_gfs[IDX4ptS(TAU_1_NUXGF,  index)];

    // Step 3.b.ii: Compute opacities
    NRPyLeakage_compute_neutrino_opacities(ghl_eos, rhoL, Y_eL, temperatureL, &tau, &kappa);
  }
}

// Step 4: Write to main memory
auxevol_gfs[IDX4ptS(KAPPA_0_NUEGF,  index)] = kappa.nue [0];
auxevol_gfs[IDX4ptS(KAPPA_1_NUEGF,  index)] = kappa.nue [1];
auxevol_gfs[IDX4ptS(KAPPA_0_ANUEGF, index)] = kappa.anue[0];
auxevol_gfs[IDX4ptS(KAPPA_1_ANUEGF, index)] = kappa.anue[1];
auxevol_gfs[IDX4ptS(KAPPA_0_NUXGF,  index)] = kappa.nux [0];
auxevol_gfs[IDX4ptS(KAPPA_1_NUXGF,  index)] = kappa.nux [1];

// Step 5: Copy the optical depths to the previous time level as a placeholder
tau_0_nue_p [index] = tau.nue [0];
tau_1_nue_p [index] = tau.nue [1];
tau_0_anue_p[index] = tau.anue[0];
tau_1_anue_p[index] = tau.anue[1];
tau_0_nux_p [index] = tau.nux [0];
tau_1_nux_p [index] = tau.nux [1];
"""
  
  add_to_Cfunction_dict(
                      includes=includes,
                      desc=desc,
                      c_type=c_type, name=name, params=params,
                      body=body,
                      rel_path_to_Cparams=os.path.join("."), 
                      loopopts="AllPoints")

def add_to_Cfunction_dict_NRPyLeakageET_optical_depths_PathOfLeastResistance():

  import BSSN.ADM_in_terms_of_BSSN as AtoB
  AtoB.ADM_in_terms_of_BSSN()
  import BSSN.BSSN_quantities as Bq
  
  gammaDD = AtoB.gammaDD
  KDD = AtoB.KDD
  betaU = Bq.betaU
  
  pre_string = r"""  
const REAL hDD00 = evol_gfs[IDX4ptS(HDD00GF, index)];
const REAL hDD11 = evol_gfs[IDX4ptS(HDD11GF, index)];
const REAL hDD22 = evol_gfs[IDX4ptS(HDD22GF, index)];
const REAL cf    = evol_gfs[IDX4ptS(CFGF, index)];  
"""
  
  POLR_metric = outputC([gammaDD[0][0], gammaDD[1][1], gammaDD[2][2]], 
              ["*gammaDD00","*gammaDD11", "*gammaDD22"], "returnstring", prestring=pre_string, params="includebraces=False")
  
  POLR_metric_str_declare = r"""
static void compute_gDDii(const paramstruct *restrict params, const int index, const REAL *restrict evol_gfs, 
                          double *restrict gammaDD00, double *restrict gammaDD11, double *restrict gammaDD22)"""
  
  POLR_metric_str = POLR_metric_str_declare + "\n{\n" +r"""#include "./set_Cparameters.h" """+"\n"+ POLR_metric +"\n\n}"
        
  prefunc = r"""
static void set_opacity_struct_from_gfs(
      const int index,
      const REAL *restrict kappa_0_nue,
      const REAL *restrict kappa_1_nue,
      const REAL *restrict kappa_0_anue,
      const REAL *restrict kappa_1_anue,
      const REAL *restrict kappa_0_nux,
      const REAL *restrict kappa_1_nux,
      ghl_neutrino_opacities *restrict kappa) {
  kappa->nue [0] = kappa_0_nue [index];
  kappa->nue [1] = kappa_1_nue [index];
  kappa->anue[0] = kappa_0_anue[index];
  kappa->anue[1] = kappa_1_anue[index];
  kappa->nux [0] = kappa_0_nux [index];
  kappa->nux [1] = kappa_1_nux [index];
}

static void set_optical_depths_struct_from_gfs(
      const int index,
      const REAL *restrict tau_0_nue,
      const REAL *restrict tau_1_nue,
      const REAL *restrict tau_0_anue,
      const REAL *restrict tau_1_anue,
      const REAL *restrict tau_0_nux,
      const REAL *restrict tau_1_nux,
      ghl_neutrino_optical_depths *restrict tau) {
  tau->nue [0] = tau_0_nue [index];
  tau->nue [1] = tau_1_nue [index];
  tau->anue[0] = tau_0_anue[index];
  tau->anue[1] = tau_1_anue[index];
  tau->nux [0] = tau_0_nux [index];
  tau->nux [1] = tau_1_nux [index];
}
""" + POLR_metric_str_declare + ";"
  
  includes = ["NRPy_basic_defines.h", "NRPy_function_prototypes.h"]
  c_type = "void"
  name = "NRPyLeakageET_optical_depths_PathOfLeastResistance"
  desc = r"""
(c) Leo Werneck
Compute opacities
https://adsabs.harvard.edu/pdf/1996A%26A...311..532R
"""
  params = r"""const paramstruct *restrict params, const REAL *restrict evol_gfs, REAL *restrict auxevol_gfs, 
  REAL *restrict tau_0_nue_p,
  REAL *restrict tau_1_nue_p, 
  REAL *restrict tau_0_anue_p, 
  REAL *restrict tau_1_anue_p, 
  REAL *restrict tau_0_nux_p, 
  REAL *restrict tau_1_nux_p"""
  
  preloop = "const REAL dxx[3] = {dxx0, dxx1, dxx2};"
  
  body = r""" 

const int Nxx_plus_2NGHOSTS_TOT = Nxx_plus_2NGHOSTS0*Nxx_plus_2NGHOSTS1*Nxx_plus_2NGHOSTS2;

REAL *kappa_0_nue  = auxevol_gfs + Nxx_plus_2NGHOSTS_TOT*KAPPA_0_NUEGF;
REAL *kappa_1_nue  = auxevol_gfs + Nxx_plus_2NGHOSTS_TOT*KAPPA_1_NUEGF;
REAL *kappa_0_anue = auxevol_gfs + Nxx_plus_2NGHOSTS_TOT*KAPPA_0_ANUEGF;
REAL *kappa_1_anue = auxevol_gfs + Nxx_plus_2NGHOSTS_TOT*KAPPA_1_ANUEGF;
REAL *kappa_0_nux  = auxevol_gfs + Nxx_plus_2NGHOSTS_TOT*KAPPA_0_NUXGF;
REAL *kappa_1_nux  = auxevol_gfs + Nxx_plus_2NGHOSTS_TOT*KAPPA_1_NUXGF;

// Step 1: Set gridpoint indices
const int i_j_k   = IDX3S(i0  , i1,   i2  );
const int ip1_j_k = IDX3S(i0+1, i1,   i2  );
const int im1_j_k = IDX3S(i0-1, i1,   i2  );
const int i_jp1_k = IDX3S(i0,   i1+1, i2  );
const int i_jm1_k = IDX3S(i0,   i1-1, i2  );
const int i_j_kp1 = IDX3S(i0,   i1  , i2+1);
const int i_j_km1 = IDX3S(i0,   i1  , i2-1);

const REAL rhoL = auxevol_gfs[IDX4ptS(RHOBGF, i_j_k)];
// Only compute optical depths if density is above a threshold
if( rhoL < rho_min_threshold || rhoL > rho_max_threshold ) {
  auxevol_gfs[IDX4ptS(TAU_0_NUEGF,  i_j_k)] = 0.0;
  auxevol_gfs[IDX4ptS(TAU_1_NUEGF,  i_j_k)] = 0.0;
  auxevol_gfs[IDX4ptS(TAU_0_ANUEGF, i_j_k)] = 0.0;
  auxevol_gfs[IDX4ptS(TAU_1_ANUEGF, i_j_k)] = 0.0;
  auxevol_gfs[IDX4ptS(TAU_0_NUXGF,  i_j_k)] = 0.0;
  auxevol_gfs[IDX4ptS(TAU_1_NUXGF,  i_j_k)] = 0.0;
}
else {
  const REAL cf = evol_gfs[IDX4ptS(CFGF, i_j_k)];
    """+psi6_str+r"""
  if( psi6L > psi6_threshold ) {
    auxevol_gfs[IDX4ptS(TAU_0_NUEGF,  i_j_k)] = 0.0;
    auxevol_gfs[IDX4ptS(TAU_1_NUEGF,  i_j_k)] = 0.0;
    auxevol_gfs[IDX4ptS(TAU_0_ANUEGF, i_j_k)] = 0.0;
    auxevol_gfs[IDX4ptS(TAU_1_ANUEGF, i_j_k)] = 0.0;
    auxevol_gfs[IDX4ptS(TAU_0_NUXGF,  i_j_k)] = 0.0;
    auxevol_gfs[IDX4ptS(TAU_1_NUXGF,  i_j_k)] = 0.0;
  }
  else {
    // Step 2: Read in metric gfs from main memory

    REAL stencil_gxx[3];
    REAL stencil_gyy[3];
    REAL stencil_gzz[3];

    REAL dummy1, dummy2;

    compute_gDDii(params, i_j_k, evol_gfs, 
                             &stencil_gxx[1], &stencil_gyy[1], &stencil_gzz[1]);
            
            
    compute_gDDii(params, im1_j_k, evol_gfs, 
                 &stencil_gxx[0], &dummy1, &dummy2);

    compute_gDDii(params, ip1_j_k, evol_gfs, 
                 &stencil_gxx[2], &dummy1, &dummy2);


    compute_gDDii(params, i_jm1_k, evol_gfs, 
                 &dummy1, &stencil_gyy[0], &dummy2);

    compute_gDDii(params, i_jp1_k, evol_gfs, 
                 &dummy1, &stencil_gyy[2], &dummy2);


    compute_gDDii(params, i_j_km1, evol_gfs, 
                 &dummy1, &dummy2, &stencil_gzz[0]);

    compute_gDDii(params, i_j_kp1, evol_gfs, 
                 &dummy1, &dummy2, &stencil_gzz[2]);
                  
    //const REAL stencil_gxx[3] = {gxx[im1_j_k], gxx[i_j_k], gxx[ip1_j_k]};
    //const REAL stencil_gyy[3] = {gyy[i_jm1_k], gyy[i_j_k], gyy[i_jp1_k]};
    //const REAL stencil_gzz[3] = {gzz[i_j_km1], gzz[i_j_k], gzz[i_j_kp1]};

    // Step 3: Read in opacity gfs from main memory
    ghl_neutrino_opacities kappa_i_j_k;
    ghl_neutrino_opacities kappa_ip1_j_k, kappa_im1_j_k;
    ghl_neutrino_opacities kappa_i_jp1_k, kappa_i_jm1_k;
    ghl_neutrino_opacities kappa_i_j_kp1, kappa_i_j_km1;
    set_opacity_struct_from_gfs(i_j_k  , kappa_0_nue, kappa_1_nue, kappa_0_anue, kappa_1_anue, kappa_0_nux, kappa_1_nux, &kappa_i_j_k  );
    set_opacity_struct_from_gfs(ip1_j_k, kappa_0_nue, kappa_1_nue, kappa_0_anue, kappa_1_anue, kappa_0_nux, kappa_1_nux, &kappa_ip1_j_k);
    set_opacity_struct_from_gfs(im1_j_k, kappa_0_nue, kappa_1_nue, kappa_0_anue, kappa_1_anue, kappa_0_nux, kappa_1_nux, &kappa_im1_j_k);
    set_opacity_struct_from_gfs(i_jp1_k, kappa_0_nue, kappa_1_nue, kappa_0_anue, kappa_1_anue, kappa_0_nux, kappa_1_nux, &kappa_i_jp1_k);
    set_opacity_struct_from_gfs(i_jm1_k, kappa_0_nue, kappa_1_nue, kappa_0_anue, kappa_1_anue, kappa_0_nux, kappa_1_nux, &kappa_i_jm1_k);
    set_opacity_struct_from_gfs(i_j_kp1, kappa_0_nue, kappa_1_nue, kappa_0_anue, kappa_1_anue, kappa_0_nux, kappa_1_nux, &kappa_i_j_kp1);
    set_opacity_struct_from_gfs(i_j_km1, kappa_0_nue, kappa_1_nue, kappa_0_anue, kappa_1_anue, kappa_0_nux, kappa_1_nux, &kappa_i_j_km1);

    // Step 4: Read in optical depth gfs from main memory
    ghl_neutrino_optical_depths tau_ip1_j_k, tau_im1_j_k;
    ghl_neutrino_optical_depths tau_i_jp1_k, tau_i_jm1_k;
    ghl_neutrino_optical_depths tau_i_j_kp1, tau_i_j_km1;
    set_optical_depths_struct_from_gfs(ip1_j_k, tau_0_nue_p, tau_1_nue_p, tau_0_anue_p, tau_1_anue_p, tau_0_nux_p, tau_1_nux_p, &tau_ip1_j_k);
    set_optical_depths_struct_from_gfs(im1_j_k, tau_0_nue_p, tau_1_nue_p, tau_0_anue_p, tau_1_anue_p, tau_0_nux_p, tau_1_nux_p, &tau_im1_j_k);
    set_optical_depths_struct_from_gfs(i_jp1_k, tau_0_nue_p, tau_1_nue_p, tau_0_anue_p, tau_1_anue_p, tau_0_nux_p, tau_1_nux_p, &tau_i_jp1_k);
    set_optical_depths_struct_from_gfs(i_jm1_k, tau_0_nue_p, tau_1_nue_p, tau_0_anue_p, tau_1_anue_p, tau_0_nux_p, tau_1_nux_p, &tau_i_jm1_k);
    set_optical_depths_struct_from_gfs(i_j_kp1, tau_0_nue_p, tau_1_nue_p, tau_0_anue_p, tau_1_anue_p, tau_0_nux_p, tau_1_nux_p, &tau_i_j_kp1);
    set_optical_depths_struct_from_gfs(i_j_km1, tau_0_nue_p, tau_1_nue_p, tau_0_anue_p, tau_1_anue_p, tau_0_nux_p, tau_1_nux_p, &tau_i_j_km1);

    // Step 5: Compute the new optical depths
    ghl_neutrino_optical_depths tau_i_j_k;
    NRPyLeakage_optical_depths_PathOfLeastResistance(dxx, stencil_gxx, stencil_gyy, stencil_gzz,
                                                     &kappa_ip1_j_k, &kappa_im1_j_k,
                                                     &kappa_i_jp1_k, &kappa_i_jm1_k,
                                                     &kappa_i_j_kp1, &kappa_i_j_km1,
                                                     &tau_ip1_j_k  , &tau_im1_j_k,
                                                     &tau_i_jp1_k  , &tau_i_jm1_k,
                                                     &tau_i_j_kp1  , &tau_i_j_km1,
                                                     &kappa_i_j_k  , &tau_i_j_k);

    // Step 6: Write to main memory
    auxevol_gfs[IDX4ptS(TAU_0_NUEGF,  i_j_k)] = tau_i_j_k.nue [0];
    auxevol_gfs[IDX4ptS(TAU_1_NUEGF,  i_j_k)] = tau_i_j_k.nue [1];
    auxevol_gfs[IDX4ptS(TAU_0_ANUEGF, i_j_k)] = tau_i_j_k.anue[0];
    auxevol_gfs[IDX4ptS(TAU_1_ANUEGF, i_j_k)] = tau_i_j_k.anue[1];
    auxevol_gfs[IDX4ptS(TAU_0_NUXGF,  i_j_k)] = tau_i_j_k.nux [0];
    auxevol_gfs[IDX4ptS(TAU_1_NUXGF,  i_j_k)] = tau_i_j_k.nux [1];
  }
}
"""

  add_to_Cfunction_dict(
                        includes=includes,
                        desc=desc,
                        c_type=c_type, name=name, params=params,
                        prefunc=prefunc,
                        preloop=preloop,
                        body=body,
                        rel_path_to_Cparams=os.path.join("."), 
                        loopopts="InteriorPoints")

  outC_function_dict[name] = outC_function_dict[name]+POLR_metric_str



def add_to_Cfunction_dict_NRPyLeakageET_getGlobalL2norm():
        
  name = "NRPyLeakageET_getGlobalL2norm"
  includes = ["NRPy_basic_defines.h", "NRPy_function_prototypes.h"]
  c_type = "void"
  desc = "compute L2norm of change in optical depth"
  
  params = r"""const paramstruct *restrict params, const REAL *restrict evol_gfs, REAL *restrict auxevol_gfs, 
  REAL *restrict tau_0_nue_aux,
  REAL *restrict tau_1_nue_aux, 
  REAL *restrict tau_0_anue_aux, 
  REAL *restrict tau_1_anue_aux, 
  REAL *restrict tau_0_nux_aux, 
  REAL *restrict tau_1_nux_aux,
  REAL *restrict l2norm_global"""
  
  body = r"""
// Step 1: Compute the differences in the optical depths between two consecutive iterations
NRPyLeakageET_compute_optical_depth_change(params, evol_gfs, auxevol_gfs, 
                                           tau_0_nue_aux,  tau_1_nue_aux,
                                           tau_0_anue_aux, tau_1_anue_aux,
                                           tau_0_nux_aux,  tau_1_nux_aux);

const REAL *tau_auxs[6] = {tau_0_nue_aux,  tau_1_nue_aux,
                           tau_0_anue_aux, tau_1_anue_aux,
                           tau_0_nux_aux,  tau_1_nux_aux};
// Step 6: Compute the global L2-norm:
// l2norm = sqrt(dTau_0_nue + dTau_0_anue + dTau_0_nux
//             + dTau_1_nue + dTau_1_anue + dTau_1_nux)                         
  
REAL l2norm = 0.0;
#pragma omp parallel for reduction(+:l2norm)
LOOP_REGION(0, Nxx_plus_2NGHOSTS2,
            0, Nxx_plus_2NGHOSTS1,
            0, Nxx_plus_2NGHOSTS2){

    const int index = IDX3S(i0, i1, i2);
    for (int i = 0; i<6; i++) l2norm += tau_auxs[i][index];
}

*l2norm_global = sqrt(l2norm);

// Step 3: Update auxiliary variables
NRPyLeakage_CopyOpticalDepthsToAux(params, auxevol_gfs, 
                                   tau_0_nue_aux,  tau_1_nue_aux,
                                   tau_0_anue_aux, tau_1_anue_aux,
                                   tau_0_nux_aux,  tau_1_nux_aux);
"""
    
  add_to_Cfunction_dict(
                      includes=includes,
                      desc=desc,
                      c_type=c_type, name=name, params=params,
                      body=body,
                      rel_path_to_Cparams=os.path.join("."))




def add_to_Cfunction_dict_NRPyLeakage_CopyOpticalDepthsToAux():
        
  includes = ["NRPy_basic_defines.h", "NRPy_function_prototypes.h"]
  c_type = "void"
  name = "NRPyLeakage_CopyOpticalDepthsToAux"
  desc = "copy optical depth GFs to auxillary GFs"
  params = r"""const paramstruct *restrict params, REAL *restrict auxevol_gfs, 
  REAL *restrict tau_0_nue_aux,
  REAL *restrict tau_1_nue_aux, 
  REAL *restrict tau_0_anue_aux, 
  REAL *restrict tau_1_anue_aux, 
  REAL *restrict tau_0_nux_aux, 
  REAL *restrict tau_1_nux_aux"""
  
  body = r"""
const int index = IDX3S(i0, i1, i2);
tau_0_nue_aux [index] = auxevol_gfs[IDX4ptS(TAU_0_NUEGF,  index)];
tau_1_nue_aux [index] = auxevol_gfs[IDX4ptS(TAU_1_NUEGF,  index)];
tau_0_anue_aux[index] = auxevol_gfs[IDX4ptS(TAU_0_ANUEGF, index)];
tau_1_anue_aux[index] = auxevol_gfs[IDX4ptS(TAU_1_ANUEGF, index)];
tau_0_nux_aux [index] = auxevol_gfs[IDX4ptS(TAU_0_NUXGF,  index)];
tau_1_nux_aux [index] = auxevol_gfs[IDX4ptS(TAU_1_NUXGF,  index)];
"""    
  add_to_Cfunction_dict(
          includes=includes,
          desc=desc,
          c_type=c_type, name=name, params=params,
          body=body,
          rel_path_to_Cparams=os.path.join("."), 
          loopopts="AllPoints")



def add_to_Cfunction_dict_NRPyLeakageET_compute_optical_depth_change():
        
  prefunc = r"""
static inline REAL relative_difference(const REAL a, const REAL b) {
  if     (a!=0.0) return 1.0-b/a;
  else if(b!=0.0) return 1.0-a/b;
  else            return 0.0;
}
  """
    
  includes = ["NRPy_basic_defines.h", "NRPy_function_prototypes.h"]
  c_type = "void"
  name = "NRPyLeakageET_compute_optical_depth_change"
  desc = "compute change in optical depth at some iteration"
  params = r"""const paramstruct *restrict params, const REAL *restrict evol_gfs, REAL *restrict auxevol_gfs, 
  REAL *restrict tau_0_nue_aux,
  REAL *restrict tau_1_nue_aux, 
  REAL *restrict tau_0_anue_aux, 
  REAL *restrict tau_1_anue_aux, 
  REAL *restrict tau_0_nux_aux, 
  REAL *restrict tau_1_nux_aux"""
  
  body = r"""
// Step 1: Set local index
const int index = IDX3S(i0, i1, i2);
const REAL rhoL = auxevol_gfs[IDX4ptS(RHOBGF, index)];
if( rhoL < rho_min_threshold || rhoL > rho_max_threshold ) {
  tau_0_nue_aux [index] = 0.0;
  tau_1_nue_aux [index] = 0.0;
  tau_0_anue_aux[index] = 0.0;
  tau_1_anue_aux[index] = 0.0;
  tau_0_nux_aux [index] = 0.0;
  tau_1_nux_aux [index] = 0.0;
}
else {
  const REAL cf = evol_gfs[IDX4ptS(CFGF, index)];
    """+psi6_str+r"""
    
  if( psi6L > psi6_threshold ) {
    tau_0_nue_aux [index] = 0.0;
    tau_1_nue_aux [index] = 0.0;
    tau_0_anue_aux[index] = 0.0;
    tau_1_anue_aux[index] = 0.0;
    tau_0_nux_aux [index] = 0.0;
    tau_1_nux_aux [index] = 0.0;
  }
  else {
    // Step 2: Compute relative differences in tau
    const REAL rel_diff_tau_0_nue  = relative_difference(auxevol_gfs[IDX4ptS(TAU_0_NUEGF,  index)],
                                                                                  tau_0_nue_aux [index]);
    const REAL rel_diff_tau_1_nue  = relative_difference(auxevol_gfs[IDX4ptS(TAU_1_NUEGF,  index)],
                                                                                  tau_1_nue_aux [index]);
    const REAL rel_diff_tau_0_anue = relative_difference(auxevol_gfs[IDX4ptS(TAU_0_ANUEGF, index)],
                                                                                  tau_0_anue_aux[index]);
    const REAL rel_diff_tau_1_anue = relative_difference(auxevol_gfs[IDX4ptS(TAU_1_ANUEGF, index)],
                                                                                  tau_1_anue_aux[index]);
    const REAL rel_diff_tau_0_nux  = relative_difference(auxevol_gfs[IDX4ptS(TAU_0_NUXGF,  index)],
                                                                                  tau_0_nux_aux [index]);
    const REAL rel_diff_tau_1_nux  = relative_difference(auxevol_gfs[IDX4ptS(TAU_1_NUXGF,  index)],
                                                                                  tau_1_nux_aux [index]);

    // Step 3: Write to main memory
    tau_0_nue_aux [index] = rel_diff_tau_0_nue *rel_diff_tau_0_nue;
    tau_1_nue_aux [index] = rel_diff_tau_1_nue *rel_diff_tau_1_nue;
    tau_0_anue_aux[index] = rel_diff_tau_0_anue*rel_diff_tau_0_anue;
    tau_1_anue_aux[index] = rel_diff_tau_1_anue*rel_diff_tau_1_anue;
    tau_0_nux_aux [index] = rel_diff_tau_0_nux *rel_diff_tau_0_nux;
    tau_1_nux_aux [index] = rel_diff_tau_1_nux *rel_diff_tau_1_nux;
  }
}
"""
    
  add_to_Cfunction_dict(
              includes=includes,
              desc=desc,
              c_type=c_type, name=name, params=params,
              prefunc=prefunc,
              body=body,
              rel_path_to_Cparams=os.path.join("."), 
              loopopts="InteriorPoints")


# For evolutions
def add_to_Cfunction_dict_NRPyLeakageET_compute_neutrino_opacities_and_add_GRHD_source_terms():
        
  includes = ["NRPy_basic_defines.h", "NRPy_function_prototypes.h"]
  c_type = "void"
  name = "NRPyLeakageET_compute_neutrino_opacities_and_add_GRHD_source_terms"
  desc = r"""
(c) Leo Werneck
Compute opacities
https://adsabs.harvard.edu/pdf/1996A%26A...311..532R
"""
  params = r"""const paramstruct *restrict params, const ghl_eos_parameters *restrict ghl_eos,
  const REAL *restrict evol_gfs, REAL *restrict rhs_gfs, REAL *restrict auxevol_gfs"""

  R_sourceL, Q_sourceL = sp.symbols(["R_sourceL", "Q_sourceL"])

  vU = ixp.declarerank1("vU")
  u4Ut = sp.symbols("u4Ut")

  u4U = ixp.zerorank1(DIM=4)

  u4U[0] = u4Ut
  for i in range(3):
      u4U[i+1] = u4Ut * vU[i]

  import BSSN.ADM_in_terms_of_BSSN as AtoB
  AtoB.ADM_in_terms_of_BSSN()
  import BSSN.BSSN_quantities as Bq

  Bq.BSSN_basic_tensors()

  import BSSN.ADMBSSN_tofrom_4metric as AB4m
  AB4m.g4DD_ito_BSSN_or_ADM("BSSN")

  u4D = ixp.zerorank1(DIM=4)

  for mu in range(4):
      for nu in range(4):
          u4D[mu] += AB4m.g4DD[mu][nu] * u4U[nu]

  for mu in range(4):
      u4D[mu] = sp.simplify(u4D[mu])

  Bq.phi_and_derivs()
  psi6 = Bq.exp_m4phi**(-0.5*3)

  Y_e_star_sourceL =  Bq.alpha     * psi6 * R_sourceL
  tau_sourceL      = (Bq.alpha**2) * psi6 * Q_sourceL * u4Ut

  StildeD_sourceL  = ixp.zerorank1()
  for i in range(3):
      StildeD_sourceL[i] = sp.simplify(Bq.alpha * psi6 * Q_sourceL * u4D[i+1])


  sourceL_str = outputC([Y_e_star_sourceL, tau_sourceL, *StildeD_sourceL], 
                        ["Y_e_star_sourceL", "tau_sourceL", "StildeD_sourceL[0]", "StildeD_sourceL[1]", "StildeD_sourceL[2]"], 
                        filename='returnstring', params="includebraces=true")  
  
  body = r"""  
// Step 1: Set local gridfunction index
const int index = IDX3S(i0, i1, i2);

// Step 2: Declare variables
const REAL rho = auxevol_gfs[IDX4ptS(RHOBGF, index)];
const REAL vU0 = auxevol_gfs[IDX4ptS(VU0GF,  index)];
const REAL vU1 = auxevol_gfs[IDX4ptS(VU1GF,  index)];
const REAL vU2 = auxevol_gfs[IDX4ptS(VU2GF,  index)];
const REAL u4Ut = auxevol_gfs[IDX4ptS(U4UTGF,  index)];

const REAL alpha = evol_gfs[IDX4ptS(ALPHAGF, index)];
const REAL vetU0 = evol_gfs[IDX4ptS(VETU0GF, index)];
const REAL vetU1 = evol_gfs[IDX4ptS(VETU1GF, index)];
const REAL vetU2 = evol_gfs[IDX4ptS(VETU2GF, index)];

const REAL hDD00 = evol_gfs[IDX4ptS(HDD00GF, index)];
const REAL hDD01 = evol_gfs[IDX4ptS(HDD01GF, index)];
const REAL hDD02 = evol_gfs[IDX4ptS(HDD02GF, index)];
const REAL hDD11 = evol_gfs[IDX4ptS(HDD11GF, index)];
const REAL hDD12 = evol_gfs[IDX4ptS(HDD12GF, index)];
const REAL hDD22 = evol_gfs[IDX4ptS(HDD22GF, index)];
const REAL cf    = evol_gfs[IDX4ptS(CFGF,    index)];

ghl_neutrino_opacities kappa;
ghl_neutrino_optical_depths tau;

// Step 3: Check density threshold, compute opacities
if( rho < rho_min_threshold || rho > rho_max_threshold ) {
  // Step 3.a: Below density threshold; set opacities and optical
  // depths to zero
  kappa.nue[0] = kappa.anue[0] = kappa.nux[0] = 0.0;
  kappa.nue[1] = kappa.anue[1] = kappa.nux[1] = 0.0;
}
else {
"""+psi6_str+r"""
  if( psi6L > psi6_threshold ) {
    kappa.nue[0] = kappa.anue[0] = kappa.nux[0] = 0.0;
    kappa.nue[1] = kappa.anue[1] = kappa.nux[1] = 0.0;
  }
  else {
    // Step 3.b: Above density threshold, compute opacities
    // Step 3.b.i: Read from main memory
    const REAL Y_e         = auxevol_gfs[IDX4ptS(YEGF, index)];
    const REAL temperature = auxevol_gfs[IDX4ptS(TEMPERATUREGF, index)];
    tau.nue [0]            = auxevol_gfs[IDX4ptS(TAU_0_NUEGF,  index)];
    tau.nue [1]            = auxevol_gfs[IDX4ptS(TAU_1_NUEGF,  index)];
    tau.anue[0]            = auxevol_gfs[IDX4ptS(TAU_0_ANUEGF, index)];
    tau.anue[1]            = auxevol_gfs[IDX4ptS(TAU_1_ANUEGF, index)];
    tau.nux [0]            = auxevol_gfs[IDX4ptS(TAU_0_NUXGF,  index)];
    tau.nux [1]            = auxevol_gfs[IDX4ptS(TAU_1_NUXGF,  index)];

    // Compute R, Q, and the neutrino opacities
    REAL R_sourceL, Q_sourceL;
    NRPyLeakage_compute_neutrino_opacities_and_GRMHD_source_terms(ghl_eos,
                                                                  rho, Y_e, temperature,
                                                                  &tau, &kappa, &R_sourceL, &Q_sourceL);
    // Now add to RHSs
    REAL Y_e_star_sourceL, tau_sourceL;
    REAL StildeD_sourceL[3];
""" + sourceL_str + r"""
    rhs_gfs[IDX4ptS(YE_STARGF, index)]   += Y_e_star_sourceL;
    rhs_gfs[IDX4ptS(TAU_TILDEGF, index)] += tau_sourceL;
    rhs_gfs[IDX4ptS(STILDED0GF, index)]  += StildeD_sourceL[0];
    rhs_gfs[IDX4ptS(STILDED1GF, index)]  += StildeD_sourceL[1];
    rhs_gfs[IDX4ptS(STILDED2GF, index)]  += StildeD_sourceL[2];
  }
}

// Step 4: Write to main memory
auxevol_gfs[IDX4ptS(KAPPA_0_NUEGF,  index)] = kappa.nue [0];
auxevol_gfs[IDX4ptS(KAPPA_1_NUEGF,  index)] = kappa.nue [1];
auxevol_gfs[IDX4ptS(KAPPA_0_ANUEGF, index)] = kappa.anue[0];
auxevol_gfs[IDX4ptS(KAPPA_1_ANUEGF, index)] = kappa.anue[1];
auxevol_gfs[IDX4ptS(KAPPA_0_NUXGF,  index)] = kappa.nux [0];
auxevol_gfs[IDX4ptS(KAPPA_1_NUXGF,  index)] = kappa.nux [1];
"""  
  add_to_Cfunction_dict(
                      includes=includes,
                      desc=desc,
                      c_type=c_type, name=name, params=params,
                      body=body,
                      rel_path_to_Cparams=os.path.join("."),
                      loopopts="InteriorPoints")



def add_to_Cfunction_dict_NRPyLeakageET_update_opacities_and_optical_depths__add_to_RHSs():
          
  name = "NRPyLeakageET_update_opacities_and_optical_depths__add_to_RHSs"
  includes = ["NRPy_basic_defines.h", "NRPy_function_prototypes.h"]
  c_type = "void"
  desc = "update opacities and optical depths to current timestep, add to GRHD RHSs"
  params = "const paramstruct *restrict params, const ghl_eos_parameters *restrict ghl_eos, const REAL *restrict evol_gfs, REAL *restrict rhs_gfs, REAL *restrict auxevol_gfs"

  body = r"""
// Step 0: set pointers to GFs not in use at the moment

const int Nxx_plus_2NGHOSTS_TOT = Nxx_plus_2NGHOSTS0*Nxx_plus_2NGHOSTS1*Nxx_plus_2NGHOSTS2;

REAL *tau_0_nue_p  = auxevol_gfs + Nxx_plus_2NGHOSTS_TOT*RBARDD00GF;
REAL *tau_1_nue_p  = auxevol_gfs + Nxx_plus_2NGHOSTS_TOT*RBARDD01GF;
REAL *tau_0_anue_p = auxevol_gfs + Nxx_plus_2NGHOSTS_TOT*RBARDD02GF;
REAL *tau_1_anue_p = auxevol_gfs + Nxx_plus_2NGHOSTS_TOT*RBARDD11GF;
REAL *tau_0_nux_p  = auxevol_gfs + Nxx_plus_2NGHOSTS_TOT*RBARDD12GF;
REAL *tau_1_nux_p  = auxevol_gfs + Nxx_plus_2NGHOSTS_TOT*RBARDD22GF;

// Now perform update for opacities and source terms
NRPyLeakageET_compute_neutrino_opacities_and_add_GRHD_source_terms(params, ghl_eos, 
                                                                   evol_gfs, rhs_gfs,
                                                                   auxevol_gfs);

// Now copy valid data to these tmp storage
NRPyLeakage_CopyOpticalDepthsToAux(params, auxevol_gfs, 
                                   tau_0_nue_p,  tau_1_nue_p,
                                   tau_0_anue_p, tau_1_anue_p,
                                   tau_0_nux_p,  tau_1_nux_p);

// Finally, update optical depths
NRPyLeakageET_optical_depths_PathOfLeastResistance(params, evol_gfs, auxevol_gfs,
                                                   tau_0_nue_p,  tau_1_nue_p, 
                                                   tau_0_anue_p, tau_1_anue_p, 
                                                   tau_0_nux_p,  tau_1_nux_p);                                                  
"""
  add_to_Cfunction_dict(
                      includes=includes,
                      desc=desc,
                      c_type=c_type, name=name, params=params,
                      body=body,
                      rel_path_to_Cparams=os.path.join("."))


def add_to_all_Cfunctions():
  add_to_Cfunction_dict_NRPyLeakageET_Initialize()
  add_to_Cfunction_dict_NRPyLeakageET_compute_optical_depth_change()
  add_to_Cfunction_dict_NRPyLeakage_CopyOpticalDepthsToAux()
  add_to_Cfunction_dict_NRPyLeakageET_getGlobalL2norm()
  add_to_Cfunction_dict_NRPyLeakageET_optical_depths_initialize_to_zero()
  add_to_Cfunction_dict_NRPyLeakageET_compute_neutrino_opacities()
  add_to_Cfunction_dict_NRPyLeakageET_optical_depths_PathOfLeastResistance()
  add_to_Cfunction_dict_NRPyLeakageET_compute_neutrino_opacities_and_add_GRHD_source_terms()
  add_to_Cfunction_dict_NRPyLeakageET_update_opacities_and_optical_depths__add_to_RHSs()