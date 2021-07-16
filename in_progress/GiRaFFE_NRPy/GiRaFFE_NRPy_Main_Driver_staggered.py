# Step 0: Add NRPy's directory to the path
# https://stackoverflow.com/questions/16780014/import-file-from-parent-directory
import os,sys
nrpy_dir_path = os.path.join("..")
if nrpy_dir_path not in sys.path:
    sys.path.append(nrpy_dir_path)

from outputC import outCfunction, lhrh # NRPy+: Core C code output module
import finite_difference as fin  # NRPy+: Finite difference C code generation module
import NRPy_param_funcs as par   # NRPy+: Parameter interface
import grid as gri               # NRPy+: Functions having to do with numerical grids
import indexedexp as ixp         # NRPy+: Symbolic indexed expression (e.g., tensors, vectors, etc.) support
import reference_metric as rfm   # NRPy+: Reference metric support
import cmdline_helper as cmd     # NRPy+: Multi-platform Python command-line interface
import GiRaFFE_NRPy.GiRaFFE_NRPy_Metric_Face_Values as FCVAL
import GiRaFFE_NRPy.GiRaFFE_NRPy_PPM as PPM
import GiRaFFE_NRPy.GiRaFFE_NRPy_staggered_Afield_flux as Af
import GiRaFFE_NRPy.Stilde_flux as Sf
import GiRaFFE_NRPy.GiRaFFE_NRPy_BCs as BC
import GiRaFFE_NRPy.GiRaFFE_NRPy_staggered_A2B as A2B
import GiRaFFE_NRPy.GiRaFFE_NRPy_C2P_P2C as C2P_P2C
import GiRaFFE_NRPy.GiRaFFE_NRPy_Source_Terms as source
import GiRaFFE_NRPy.GiRaFFE_NRPy_staggered_Source_Terms as stgsrc

thismodule = "GiRaFFE_NRPy_Main_Driver"

CoordSystem = "Cartesian"

# Generating the CSE for GRFFE expressions
# operation in this notebook, and much of the CSE
# time is spent sorting CSE expressions. Disabling
# this sorting makes the C codegen 3-4x faster,
# but the tradeoff is that every time this is
# run, the CSE patterns will be different
# (though they should result in mathematically
# *identical* expressions). You can expect
# roundoff-level differences as a result.

# Once everything is debugged, we'll want to
# remove the CSE_sorting=none option so that
# the C code output is completely identical
# every time the below Python function is run.
outCparams = "outCverbose=False,CSE_sorting=none"

par.set_parval_from_str("reference_metric::CoordSystem",CoordSystem)
rfm.reference_metric() # Create ReU, ReDD needed for rescaling B-L initial data, generating BSSN RHSs, etc.

xi_damping = par.Cparameters("REAL",thismodule,"xi_damping",0.1)

# Default Kreiss-Oliger dissipation strength
default_KO_strength = 0.1
diss_strength = par.Cparameters("REAL", thismodule, "diss_strength", default_KO_strength)

def GiRaFFE_NRPy_Main_Driver_generate_all(out_dir):
    cmd.mkdir(out_dir)

    gammaDD = ixp.register_gridfunctions_for_single_rank2("AUXEVOL","gammaDD","sym01",DIM=3)
    betaU = ixp.register_gridfunctions_for_single_rank1("AUXEVOL","betaU",DIM=3)
    alpha = gri.register_gridfunctions("AUXEVOL","alpha")
    ixp.register_gridfunctions_for_single_rank1("EVOL","AD")
    BU = ixp.register_gridfunctions_for_single_rank1("AUXEVOL","BU")
    ixp.register_gridfunctions_for_single_rank1("AUXEVOL","BstaggerU")
    ValenciavU = ixp.register_gridfunctions_for_single_rank1("AUXEVOL","ValenciavU")
    gri.register_gridfunctions("EVOL","psi6Phi")
    StildeD = ixp.register_gridfunctions_for_single_rank1("EVOL","StildeD")

    gri.register_gridfunctions("AUXEVOL","psi6_temp")
    gri.register_gridfunctions("AUXEVOL","psi6center")
    gri.register_gridfunctions("AUXEVOL","cmax_x")
    gri.register_gridfunctions("AUXEVOL","cmin_x")
    gri.register_gridfunctions("AUXEVOL","cmax_y")
    gri.register_gridfunctions("AUXEVOL","cmin_y")
    gri.register_gridfunctions("AUXEVOL","cmax_z")
    gri.register_gridfunctions("AUXEVOL","cmin_z")

    subdir = "RHSs"
    stgsrc.GiRaFFE_NRPy_Source_Terms(os.path.join(out_dir,subdir))

    # Declare this symbol:
    sqrt4pi = par.Cparameters("REAL",thismodule,"sqrt4pi","sqrt(4.0*M_PI)")

    cmd.mkdir(os.path.join(out_dir,subdir))
    source.write_out_functions_for_StildeD_source_term(os.path.join(out_dir,subdir),outCparams,gammaDD,betaU,alpha,
                                                       ValenciavU,BU,sqrt4pi)

    subdir = "FCVAL"
    cmd.mkdir(os.path.join(out_dir, subdir))
    FCVAL.GiRaFFE_NRPy_FCVAL(os.path.join(out_dir,subdir))

    subdir = "PPM"
    cmd.mkdir(os.path.join(out_dir, subdir))
    PPM.GiRaFFE_NRPy_PPM(os.path.join(out_dir,subdir))


    # We will pass values of the gridfunction on the cell faces into the function. This requires us
    # to declare them as C parameters in NRPy+. We will denote this with the _face infix/suffix.
    alpha_face = gri.register_gridfunctions("AUXEVOL","alpha_face")
    gamma_faceDD = ixp.register_gridfunctions_for_single_rank2("AUXEVOL","gamma_faceDD","sym01")
    beta_faceU = ixp.register_gridfunctions_for_single_rank1("AUXEVOL","beta_faceU")

    # We'll need some more gridfunctions, now, to represent the reconstructions of BU and ValenciavU
    # on the right and left faces
    Valenciav_rU = ixp.register_gridfunctions_for_single_rank1("AUXEVOL","Valenciav_rU",DIM=3)
    B_rU = ixp.register_gridfunctions_for_single_rank1("AUXEVOL","B_rU",DIM=3)
    Valenciav_lU = ixp.register_gridfunctions_for_single_rank1("AUXEVOL","Valenciav_lU",DIM=3)
    B_lU = ixp.register_gridfunctions_for_single_rank1("AUXEVOL","B_lU",DIM=3)

    ixp.register_gridfunctions_for_single_rank1("AUXEVOL","Stilde_flux_HLLED")

    ixp.register_gridfunctions_for_single_rank1("AUXEVOL","Valenciav_rrU",DIM=3)
    ixp.register_gridfunctions_for_single_rank1("AUXEVOL","Valenciav_rlU",DIM=3)
    ixp.register_gridfunctions_for_single_rank1("AUXEVOL","Valenciav_lrU",DIM=3)
    ixp.register_gridfunctions_for_single_rank1("AUXEVOL","Valenciav_llU",DIM=3)
    ixp.register_gridfunctions_for_single_rank1("AUXEVOL","Bstagger_rU",DIM=3)
    ixp.register_gridfunctions_for_single_rank1("AUXEVOL","Bstagger_lU",DIM=3)

    subdir = "RHSs"
    Af.GiRaFFE_NRPy_Afield_flux(os.path.join(out_dir, subdir))

    Sf.generate_C_code_for_Stilde_flux(os.path.join(out_dir,subdir), True, alpha_face,gamma_faceDD,beta_faceU,
                                       Valenciav_rU,B_rU,Valenciav_lU,B_lU, sqrt4pi, write_cmax_cmin=True)

    subdir = "boundary_conditions"
    cmd.mkdir(os.path.join(out_dir,subdir))
    BC.GiRaFFE_NRPy_BCs(os.path.join(out_dir,subdir))

    subdir = "A2B"
    cmd.mkdir(os.path.join(out_dir,subdir))
    A2B.GiRaFFE_NRPy_A2B(os.path.join(out_dir,subdir))

    C2P_P2C.GiRaFFE_NRPy_C2P(StildeD,BU,gammaDD,betaU,alpha)

    values_to_print = [
                       lhrh(lhs=gri.gfaccess("in_gfs","StildeD0"),rhs=C2P_P2C.outStildeD[0]),
                       lhrh(lhs=gri.gfaccess("in_gfs","StildeD1"),rhs=C2P_P2C.outStildeD[1]),
                       lhrh(lhs=gri.gfaccess("in_gfs","StildeD2"),rhs=C2P_P2C.outStildeD[2]),
                       lhrh(lhs=gri.gfaccess("auxevol_gfs","ValenciavU0"),rhs=C2P_P2C.ValenciavU[0]),
                       lhrh(lhs=gri.gfaccess("auxevol_gfs","ValenciavU1"),rhs=C2P_P2C.ValenciavU[1]),
                       lhrh(lhs=gri.gfaccess("auxevol_gfs","ValenciavU2"),rhs=C2P_P2C.ValenciavU[2])
                      ]

    subdir = "C2P"
    cmd.mkdir(os.path.join(out_dir,subdir))
    desc = "Apply fixes to \tilde{S}_i and recompute the velocity to match with current sheet prescription."
    name = "GiRaFFE_NRPy_cons_to_prims"
    outCfunction(
        outfile  = os.path.join(out_dir,subdir,name+".h"), desc=desc, name=name,
        params   ="const paramstruct *params,REAL *xx[3],REAL *auxevol_gfs,REAL *in_gfs",
        body     = fin.FD_outputC("returnstring",values_to_print,params=outCparams).replace("IDX4","IDX4S"),
        loopopts ="AllPoints,Read_xxs",
        rel_path_to_Cparams=os.path.join("../"))

    C2P_P2C.GiRaFFE_NRPy_P2C(gammaDD,betaU,alpha,  ValenciavU,BU, sqrt4pi)

    values_to_print = [
                       lhrh(lhs=gri.gfaccess("in_gfs","StildeD0"),rhs=C2P_P2C.StildeD[0]),
                       lhrh(lhs=gri.gfaccess("in_gfs","StildeD1"),rhs=C2P_P2C.StildeD[1]),
                       lhrh(lhs=gri.gfaccess("in_gfs","StildeD2"),rhs=C2P_P2C.StildeD[2]),
                      ]

    desc = "Recompute StildeD after current sheet fix to Valencia 3-velocity to ensure consistency between conservative & primitive variables."
    name = "GiRaFFE_NRPy_prims_to_cons"
    outCfunction(
        outfile  = os.path.join(out_dir,subdir,name+".h"), desc=desc, name=name,
        params   ="const paramstruct *params,REAL *auxevol_gfs,REAL *in_gfs",
        body     = fin.FD_outputC("returnstring",values_to_print,params=outCparams).replace("IDX4","IDX4S"),
        loopopts ="AllPoints",
        rel_path_to_Cparams=os.path.join("../"))

    # Write out the main driver itself:
    with open(os.path.join(out_dir,"GiRaFFE_NRPy_Main_Driver.h"),"w") as file:
        file.write(r"""// Structure to track ghostzones for PPM:
typedef struct __gf_and_gz_struct__ {
  REAL *gf;
  int gz_lo[4],gz_hi[4];
} gf_and_gz_struct;
// Some additional constants needed for PPM:
static const int VX=0,VY=1,VZ=2,
  BX_CENTER=3,BY_CENTER=4,BZ_CENTER=5,BX_STAGGER=6,BY_STAGGER=7,BZ_STAGGER=8,
  VXR=9,VYR=10,VZR=11,VXL=12,VYL=13,VZL=14;  //<-- Be _sure_ to define MAXNUMVARS appropriately!
const int NUM_RECONSTRUCT_GFS = 15;

#define WORKAROUND_ENABLED

// Include ALL functions needed for evolution
#include "PPM/reconstruct_set_of_prims_PPM_GRFFE_NRPy.c"
#include "FCVAL/interpolate_metric_gfs_to_cell_faces.h"
#include "RHSs/calculate_StildeD0_source_term.h"
#include "RHSs/calculate_StildeD1_source_term.h"
#include "RHSs/calculate_StildeD2_source_term.h"
#include "RHSs/Lorenz_psi6phi_rhs__add_gauge_terms_to_A_i_rhs.h"
#include "RHSs/A_i_rhs_no_gauge_terms.h"
#include "A2B/compute_B_and_Bstagger_from_A.h"
#include "RHSs/calculate_Stilde_flux_D0.h"
#include "RHSs/calculate_Stilde_flux_D1.h"
#include "RHSs/calculate_Stilde_flux_D2.h"
#include "RHSs/calculate_Stilde_rhsD.h"
#include "boundary_conditions/GiRaFFE_boundary_conditions.h"
#include "C2P/GiRaFFE_NRPy_cons_to_prims.h"
#include "C2P/GiRaFFE_NRPy_prims_to_cons.h"

void workaround_Valencia_to_Drift_velocity(const paramstruct *params, REAL *vU0, const REAL *alpha, const REAL *betaU0, const REAL flux_dirn) {
#include "set_Cparameters.h"
    // Converts Valencia 3-velocities to Drift 3-velocities for testing. The variable argument
    // vu0 is any Valencia 3-velocity component or reconstruction thereof.
#pragma omp parallel for
    for (int i2 = 2*(flux_dirn==3);i2 < Nxx_plus_2NGHOSTS2-1*(flux_dirn==3);i2++) for (int i1 = 2*(flux_dirn==2);i1 < Nxx_plus_2NGHOSTS1-1*(flux_dirn==2);i1++) for (int i0 = 2*(flux_dirn==1);i0 < Nxx_plus_2NGHOSTS0-1*(flux_dirn==1);i0++) {
        int ii = IDX3S(i0,i1,i2);
        // Here, we modify the velocity in place.
        vU0[ii] = alpha[ii]*vU0[ii]-betaU0[ii];
    }
}

void workaround_Drift_to_Valencia_velocity(const paramstruct *params, REAL *vU0, const REAL *alpha, const REAL *betaU0, const REAL flux_dirn) {
#include "set_Cparameters.h"
    // Converts Drift 3-velocities to Valencia 3-velocities for testing. The variable argument
    // vu0 is any drift (i.e. IllinoisGRMHD's definition) 3-velocity component or reconstruction thereof.
#pragma omp parallel for
    for (int i2 = 2*(flux_dirn==3);i2 < Nxx_plus_2NGHOSTS2-1*(flux_dirn==3);i2++) for (int i1 = 2*(flux_dirn==2);i1 < Nxx_plus_2NGHOSTS1-1*(flux_dirn==2);i1++) for (int i0 = 2*(flux_dirn==1);i0 < Nxx_plus_2NGHOSTS0-1*(flux_dirn==1);i0++) {
        int ii = IDX3S(i0,i1,i2);
        // Here, we modify the velocity in place.
        vU0[ii] = (vU0[ii]+betaU0[ii])/alpha[ii];
    }
}

void GiRaFFE_NRPy_RHSs(const paramstruct *restrict params,REAL *restrict auxevol_gfs,REAL *restrict in_gfs,REAL *restrict rhs_gfs) {
#include "set_Cparameters.h"
    // First thing's first: initialize the RHSs to zero!
#pragma omp parallel for
    for(int ii=0;ii<Nxx_plus_2NGHOSTS0*Nxx_plus_2NGHOSTS1*Nxx_plus_2NGHOSTS2*NUM_EVOL_GFS;ii++) {
        rhs_gfs[ii] = 0.0;
    }

    // Now, we set up a bunch of structs of pointers to properly guide the PPM algorithm.
    // They also count the number of ghostzones available.
    gf_and_gz_struct in_prims[NUM_RECONSTRUCT_GFS], out_prims_r[NUM_RECONSTRUCT_GFS], out_prims_l[NUM_RECONSTRUCT_GFS];
    int which_prims_to_reconstruct[NUM_RECONSTRUCT_GFS],num_prims_to_reconstruct;
    const int Nxxp2NG012 = Nxx_plus_2NGHOSTS0*Nxx_plus_2NGHOSTS1*Nxx_plus_2NGHOSTS2;

    REAL *temporary = auxevol_gfs + Nxxp2NG012*PSI6_TEMPGF; // Using dedicated temporary variables for the staggered grid
    REAL *psi6center = auxevol_gfs + Nxxp2NG012*PSI6CENTERGF; // Because the prescription requires more acrobatics.
    // This sets pointers to the portion of auxevol_gfs containing the relevant gridfunction.
    int ww=0;
    in_prims[ww].gf      = auxevol_gfs + Nxxp2NG012*VALENCIAVU0GF;
      out_prims_r[ww].gf = auxevol_gfs + Nxxp2NG012*VALENCIAV_RU0GF;
      out_prims_l[ww].gf = auxevol_gfs + Nxxp2NG012*VALENCIAV_LU0GF;
    ww++;
    in_prims[ww].gf      = auxevol_gfs + Nxxp2NG012*VALENCIAVU1GF;
      out_prims_r[ww].gf = auxevol_gfs + Nxxp2NG012*VALENCIAV_RU1GF;
      out_prims_l[ww].gf = auxevol_gfs + Nxxp2NG012*VALENCIAV_LU1GF;
    ww++;
    in_prims[ww].gf      = auxevol_gfs + Nxxp2NG012*VALENCIAVU2GF;
      out_prims_r[ww].gf = auxevol_gfs + Nxxp2NG012*VALENCIAV_RU2GF;
      out_prims_l[ww].gf = auxevol_gfs + Nxxp2NG012*VALENCIAV_LU2GF;
    ww++;
    in_prims[ww].gf      = auxevol_gfs + Nxxp2NG012*BU0GF;
      out_prims_r[ww].gf = auxevol_gfs + Nxxp2NG012*B_RU0GF;
      out_prims_l[ww].gf = auxevol_gfs + Nxxp2NG012*B_LU0GF;
    ww++;
    in_prims[ww].gf      = auxevol_gfs + Nxxp2NG012*BU1GF;
      out_prims_r[ww].gf = auxevol_gfs + Nxxp2NG012*B_RU1GF;
      out_prims_l[ww].gf = auxevol_gfs + Nxxp2NG012*B_LU1GF;
    ww++;
    in_prims[ww].gf      = auxevol_gfs + Nxxp2NG012*BU2GF;
      out_prims_r[ww].gf = auxevol_gfs + Nxxp2NG012*B_RU2GF;
      out_prims_l[ww].gf = auxevol_gfs + Nxxp2NG012*B_LU2GF;
    ww++;
    in_prims[ww].gf      = auxevol_gfs + Nxxp2NG012*BSTAGGERU0GF;
      out_prims_r[ww].gf = auxevol_gfs + Nxxp2NG012*BSTAGGER_RU0GF;
      out_prims_l[ww].gf = auxevol_gfs + Nxxp2NG012*BSTAGGER_LU0GF;
    ww++;
    in_prims[ww].gf      = auxevol_gfs + Nxxp2NG012*BSTAGGERU1GF;
      out_prims_r[ww].gf = auxevol_gfs + Nxxp2NG012*BSTAGGER_RU1GF;
      out_prims_l[ww].gf = auxevol_gfs + Nxxp2NG012*BSTAGGER_LU1GF;
    ww++;
    in_prims[ww].gf      = auxevol_gfs + Nxxp2NG012*BSTAGGERU2GF;
      out_prims_r[ww].gf = auxevol_gfs + Nxxp2NG012*BSTAGGER_RU2GF;
      out_prims_l[ww].gf = auxevol_gfs + Nxxp2NG012*BSTAGGER_LU2GF;
    ww++;
    in_prims[ww].gf      = auxevol_gfs + Nxxp2NG012*VALENCIAV_RU0GF;
      out_prims_r[ww].gf = auxevol_gfs + Nxxp2NG012*VALENCIAV_RRU0GF;
      out_prims_l[ww].gf = auxevol_gfs + Nxxp2NG012*VALENCIAV_RLU0GF;
    ww++;
    in_prims[ww].gf      = auxevol_gfs + Nxxp2NG012*VALENCIAV_RU1GF;
      out_prims_r[ww].gf = auxevol_gfs + Nxxp2NG012*VALENCIAV_RRU1GF;
      out_prims_l[ww].gf = auxevol_gfs + Nxxp2NG012*VALENCIAV_RLU1GF;
    ww++;
    in_prims[ww].gf      = auxevol_gfs + Nxxp2NG012*VALENCIAV_RU2GF;
      out_prims_r[ww].gf = auxevol_gfs + Nxxp2NG012*VALENCIAV_RRU2GF;
      out_prims_l[ww].gf = auxevol_gfs + Nxxp2NG012*VALENCIAV_RLU2GF;
    ww++;
    in_prims[ww].gf      = auxevol_gfs + Nxxp2NG012*VALENCIAV_LU0GF;
      out_prims_r[ww].gf = auxevol_gfs + Nxxp2NG012*VALENCIAV_LRU0GF;
      out_prims_l[ww].gf = auxevol_gfs + Nxxp2NG012*VALENCIAV_LLU0GF;
    ww++;
    in_prims[ww].gf      = auxevol_gfs + Nxxp2NG012*VALENCIAV_LU1GF;
      out_prims_r[ww].gf = auxevol_gfs + Nxxp2NG012*VALENCIAV_LRU1GF;
      out_prims_l[ww].gf = auxevol_gfs + Nxxp2NG012*VALENCIAV_LLU1GF;
    ww++;
    in_prims[ww].gf      = auxevol_gfs + Nxxp2NG012*VALENCIAV_LU2GF;
      out_prims_r[ww].gf = auxevol_gfs + Nxxp2NG012*VALENCIAV_LRU2GF;
      out_prims_l[ww].gf = auxevol_gfs + Nxxp2NG012*VALENCIAV_LLU2GF;
    ww++;

    // Prims are defined AT ALL GRIDPOINTS, so we set the # of ghostzones to zero:
    for(int i=0;i<NUM_RECONSTRUCT_GFS;i++) for(int j=1;j<=3;j++) { in_prims[i].gz_lo[j]=0; in_prims[i].gz_hi[j]=0; }
    // Left/right variables are not yet defined, yet we set the # of gz's to zero by default:
    for(int i=0;i<NUM_RECONSTRUCT_GFS;i++) for(int j=1;j<=3;j++) { out_prims_r[i].gz_lo[j]=0; out_prims_r[i].gz_hi[j]=0; }
    for(int i=0;i<NUM_RECONSTRUCT_GFS;i++) for(int j=1;j<=3;j++) { out_prims_l[i].gz_lo[j]=0; out_prims_l[i].gz_hi[j]=0; }

  int flux_dirn;
  flux_dirn=0;
  interpolate_metric_gfs_to_cell_faces(params,auxevol_gfs,flux_dirn+1);
  // ftilde = 0 in GRFFE, since P=rho=0.

  /* There are two stories going on here:
   * 1) Computation of \partial_x on RHS of \partial_t {mhd_st_{x,y,z}},
   *    via PPM reconstruction onto (i-1/2,j,k), so that
   *    \partial_y F = [ F(i+1/2,j,k) - F(i-1/2,j,k) ] / dx
   * 2) Computation of \partial_t A_i, where A_i are *staggered* gridfunctions,
   *    where A_x is defined at (i,j+1/2,k+1/2), A_y at (i+1/2,j,k+1/2), etc.
   *    Ai_rhs = \partial_t A_i = \epsilon_{ijk} \psi^{6} (v^j B^k - v^j B^k),
   *    where \epsilon_{ijk} is the flat-space antisymmetric operator.
   * 2A) Az_rhs is defined at (i+1/2,j+1/2,k), and it depends on {Bx,By,vx,vy},
   *     so the trick is to reconstruct {Bx,By,vx,vy} cleverly to get to these
   *     staggered points. For example:
   * 2Aa) vx and vy are at (i,j,k), and we reconstruct them to (i-1/2,j,k) below. After
   *      this, we'll reconstruct again in the y-dir'n to get {vx,vy} at (i-1/2,j-1/2,k)
   * 2Ab) By_stagger is at (i,j+1/2,k), and we reconstruct below to (i-1/2,j+1/2,k). */
  ww=0;
  which_prims_to_reconstruct[ww]=VX;        ww++;
  which_prims_to_reconstruct[ww]=VY;        ww++;
  which_prims_to_reconstruct[ww]=VZ;        ww++;
  //which_prims_to_reconstruct[ww]=BX_CENTER; ww++;
  which_prims_to_reconstruct[ww]=BY_CENTER; ww++;
  which_prims_to_reconstruct[ww]=BZ_CENTER; ww++;
  which_prims_to_reconstruct[ww]=BY_STAGGER;ww++;
  num_prims_to_reconstruct=ww;
  // This function is housed in the file: "reconstruct_set_of_prims_PPM_GRFFE.C"
  reconstruct_set_of_prims_PPM_GRFFE_NRPy(params, auxevol_gfs, flux_dirn+1, num_prims_to_reconstruct,
                                          which_prims_to_reconstruct, in_prims, out_prims_r, out_prims_l, temporary);
  //Right and left face values of BI_CENTER are used in GRFFE__S_i__flux computation (first to compute b^a).
  //   Instead of reconstructing, we simply set B^x face values to be consistent with BX_STAGGER.
#pragma omp parallel for
  for(int k=0;k<Nxx_plus_2NGHOSTS2;k++) for(int j=0;j<Nxx_plus_2NGHOSTS1;j++) for(int i=0;i<Nxx_plus_2NGHOSTS0;i++) {
        const int index=IDX3S(i,j,k), indexim1=IDX3S(i-1+(i==0),j,k); /* indexim1=0 when i=0 */
        out_prims_r[BX_CENTER].gf[index]=out_prims_l[BX_CENTER].gf[index]=in_prims[BX_STAGGER].gf[indexim1]; }
  // Then add fluxes to RHS for hydro variables {vx,vy,vz}:
  // This function is housed in the file: "add_fluxes_and_source_terms_to_hydro_rhss.C"
  calculate_StildeD0_source_term(params,auxevol_gfs,rhs_gfs);
  calculate_Stilde_flux_D0(params,auxevol_gfs,rhs_gfs);
  calculate_Stilde_rhsD(flux_dirn+1,params,auxevol_gfs,rhs_gfs);

  // Note that we have already reconstructed vx and vy along the x-direction,
  //   at (i-1/2,j,k). That result is stored in v{x,y}{r,l}.  Bx_stagger data
  //   are defined at (i+1/2,j,k).
  // Next goal: reconstruct Bx, vx and vy at (i+1/2,j+1/2,k).
  flux_dirn=1;
  // ftilde = 0 in GRFFE, since P=rho=0.

  // in_prims[{VXR,VXL,VYR,VYL}].gz_{lo,hi} ghostzones are set to all zeros, which
  //    is incorrect. We fix this below.
  // [Note that this is a cheap operation, copying only 8 integers and a pointer.]
  in_prims[VXR]=out_prims_r[VX];
  in_prims[VXL]=out_prims_l[VX];
  in_prims[VYR]=out_prims_r[VY];
  in_prims[VYL]=out_prims_l[VY];

  /* There are two stories going on here:
   * 1) Computation of \partial_y on RHS of \partial_t {mhd_st_{x,y,z}},
   *    via PPM reconstruction onto (i,j-1/2,k), so that
   *    \partial_y F = [ F(i,j+1/2,k) - F(i,j-1/2,k) ] / dy
   * 2) Computation of \partial_t A_i, where A_i are *staggered* gridfunctions,
   *    where A_x is defined at (i,j+1/2,k+1/2), A_y at (i+1/2,j,k+1/2), etc.
   *    Ai_rhs = \partial_t A_i = \epsilon_{ijk} \psi^{6} (v^j B^k - v^j B^k),
   *    where \epsilon_{ijk} is the flat-space antisymmetric operator.
   * 2A) Az_rhs is defined at (i+1/2,j+1/2,k), and it depends on {Bx,By,vx,vy},
   *     so the trick is to reconstruct {Bx,By,vx,vy} cleverly to get to these
   *     staggered points. For example:
   * 2Aa) VXR = [right-face of vx reconstructed along x-direction above] is at (i-1/2,j,k),
   *      and we reconstruct it to (i-1/2,j-1/2,k) below. Similarly for {VXL,VYR,VYL}
   * 2Ab) Bx_stagger is at (i+1/2,j,k), and we reconstruct to (i+1/2,j-1/2,k) below
   * 2Ac) By_stagger is at (i-1/2,j+1/2,k) already for Az_rhs, from the previous step.
   * 2B) Ax_rhs is defined at (i,j+1/2,k+1/2), and it depends on {By,Bz,vy,vz}.
   *     Again the trick is to reconstruct these onto these staggered points.
   * 2Ba) Bz_stagger is at (i,j,k+1/2), and we reconstruct to (i,j-1/2,k+1/2) below */
  ww=0;
  // NOTE! The order of variable reconstruction is important here,
  //   as we don't want to overwrite {vxr,vxl,vyr,vyl}!
  which_prims_to_reconstruct[ww]=VXR;       ww++;
  which_prims_to_reconstruct[ww]=VYR;       ww++;
  which_prims_to_reconstruct[ww]=VXL;       ww++;
  which_prims_to_reconstruct[ww]=VYL;       ww++;
  num_prims_to_reconstruct=ww;
#ifdef WORKAROUND_ENABLED
  workaround_Valencia_to_Drift_velocity(params,auxevol_gfs+Nxxp2NG012*VALENCIAV_RU0GF,auxevol_gfs+Nxxp2NG012*ALPHA_FACEGF,auxevol_gfs+Nxxp2NG012*BETA_FACEU0GF,flux_dirn+1);
  workaround_Valencia_to_Drift_velocity(params,auxevol_gfs+Nxxp2NG012*VALENCIAV_RU1GF,auxevol_gfs+Nxxp2NG012*ALPHA_FACEGF,auxevol_gfs+Nxxp2NG012*BETA_FACEU1GF,flux_dirn+1);
  workaround_Valencia_to_Drift_velocity(params,auxevol_gfs+Nxxp2NG012*VALENCIAV_LU0GF,auxevol_gfs+Nxxp2NG012*ALPHA_FACEGF,auxevol_gfs+Nxxp2NG012*BETA_FACEU0GF,flux_dirn+1);
  workaround_Valencia_to_Drift_velocity(params,auxevol_gfs+Nxxp2NG012*VALENCIAV_LU1GF,auxevol_gfs+Nxxp2NG012*ALPHA_FACEGF,auxevol_gfs+Nxxp2NG012*BETA_FACEU1GF,flux_dirn+1);
#endif /*WORKAROUND_ENABLED*/
  // This function is housed in the file: "reconstruct_set_of_prims_PPM_GRFFE.C"
  reconstruct_set_of_prims_PPM_GRFFE_NRPy(params, auxevol_gfs, flux_dirn+1, num_prims_to_reconstruct,
                                          which_prims_to_reconstruct, in_prims, out_prims_r, out_prims_l, temporary);
#ifdef WORKAROUND_ENABLED
  workaround_Drift_to_Valencia_velocity(params,auxevol_gfs+Nxxp2NG012*VALENCIAV_RU0GF,auxevol_gfs+Nxxp2NG012*ALPHA_FACEGF,auxevol_gfs+Nxxp2NG012*BETA_FACEU0GF,flux_dirn+1);
  workaround_Drift_to_Valencia_velocity(params,auxevol_gfs+Nxxp2NG012*VALENCIAV_RU1GF,auxevol_gfs+Nxxp2NG012*ALPHA_FACEGF,auxevol_gfs+Nxxp2NG012*BETA_FACEU1GF,flux_dirn+1);
  workaround_Drift_to_Valencia_velocity(params,auxevol_gfs+Nxxp2NG012*VALENCIAV_LU0GF,auxevol_gfs+Nxxp2NG012*ALPHA_FACEGF,auxevol_gfs+Nxxp2NG012*BETA_FACEU0GF,flux_dirn+1);
  workaround_Drift_to_Valencia_velocity(params,auxevol_gfs+Nxxp2NG012*VALENCIAV_LU1GF,auxevol_gfs+Nxxp2NG012*ALPHA_FACEGF,auxevol_gfs+Nxxp2NG012*BETA_FACEU1GF,flux_dirn+1);
#endif /*WORKAROUND_ENABLED*/
  interpolate_metric_gfs_to_cell_faces(params,auxevol_gfs,flux_dirn+1);
  ww=0;
  // Reconstruct other primitives last!
  which_prims_to_reconstruct[ww]=VX;        ww++;
  which_prims_to_reconstruct[ww]=VY;        ww++;
  which_prims_to_reconstruct[ww]=VZ;        ww++;
  which_prims_to_reconstruct[ww]=BX_CENTER; ww++;
  //which_prims_to_reconstruct[ww]=BY_CENTER; ww++;
  which_prims_to_reconstruct[ww]=BZ_CENTER; ww++;
  which_prims_to_reconstruct[ww]=BX_STAGGER;ww++;
  which_prims_to_reconstruct[ww]=BZ_STAGGER;ww++;
  num_prims_to_reconstruct=ww;
  // This function is housed in the file: "reconstruct_set_of_prims_PPM_GRFFE.C"
  reconstruct_set_of_prims_PPM_GRFFE_NRPy(params, auxevol_gfs, flux_dirn+1, num_prims_to_reconstruct,
                                          which_prims_to_reconstruct, in_prims, out_prims_r, out_prims_l, temporary);
  //Right and left face values of BI_CENTER are used in GRFFE__S_i__flux computation (first to compute b^a).
  //   Instead of reconstructing, we simply set B^y face values to be consistent with BY_STAGGER.
#pragma omp parallel for
  for(int k=0;k<Nxx_plus_2NGHOSTS2;k++) for(int j=0;j<Nxx_plus_2NGHOSTS1;j++) for(int i=0;i<Nxx_plus_2NGHOSTS0;i++) {
        const int index=IDX3S(i,j,k), indexjm1=IDX3S(i,j-1+(j==0),k); /* indexjm1=0 when j=0 */
        out_prims_r[BY_CENTER].gf[index]=out_prims_l[BY_CENTER].gf[index]=in_prims[BY_STAGGER].gf[indexjm1]; }
  // Then add fluxes to RHS for hydro variables {vx,vy,vz}:
  // This function is housed in the file: "add_fluxes_and_source_terms_to_hydro_rhss.C"
  calculate_StildeD1_source_term(params,auxevol_gfs,rhs_gfs);
  calculate_Stilde_flux_D1(params,auxevol_gfs,rhs_gfs);
  calculate_Stilde_rhsD(flux_dirn+1,params,auxevol_gfs,rhs_gfs);

  /*****************************************
   * COMPUTING RHS OF A_z, BOOKKEEPING NOTE:
   * We want to compute
   * \partial_t A_z - [gauge terms] = \psi^{6} (v^x B^y - v^y B^x).
   * A_z is defined at (i+1/2,j+1/2,k).
   * ==========================
   * Where defined  | Variables
   * (i-1/2,j-1/2,k)| {vxrr,vxrl,vxlr,vxll,vyrr,vyrl,vylr,vyll}
   * (i+1/2,j-1/2,k)| {Bx_stagger_r,Bx_stagger_l} (see Table 1 in arXiv:1007.2848)
   * (i-1/2,j+1/2,k)| {By_stagger_r,By_stagger_l} (see Table 1 in arXiv:1007.2848)
   * (i,j,k)        | {phi}
   * ==========================
   ******************************************/
  // Interpolates to i+1/2
#define IPH(METRICm1,METRICp0,METRICp1,METRICp2) (-0.0625*((METRICm1) + (METRICp2)) + 0.5625*((METRICp0) + (METRICp1)))
  // Next compute sqrt(gamma)=psi^6 at (i+1/2,j+1/2,k):
  // To do so, we first compute the sqrt of the metric determinant at all points:
#pragma omp parallel for
  for(int k=0;k<Nxx_plus_2NGHOSTS2;k++) for(int j=0;j<Nxx_plus_2NGHOSTS1;j++) for(int i=0;i<Nxx_plus_2NGHOSTS0;i++) {
        const int index=IDX3S(i,j,k);
        const REAL gxx = auxevol_gfs[IDX4ptS(GAMMADD00GF,index)];
        const REAL gxy = auxevol_gfs[IDX4ptS(GAMMADD01GF,index)];
        const REAL gxz = auxevol_gfs[IDX4ptS(GAMMADD02GF,index)];
        const REAL gyy = auxevol_gfs[IDX4ptS(GAMMADD11GF,index)];
        const REAL gyz = auxevol_gfs[IDX4ptS(GAMMADD12GF,index)];
        const REAL gzz = auxevol_gfs[IDX4ptS(GAMMADD22GF,index)];
        psi6center[index] = sqrt( gxx*gyy*gzz
                               -  gxx*gyz*gyz
                               +2*gxy*gxz*gyz
                               -  gyy*gxz*gxz
                               -  gzz*gxy*gxy );
      }
#pragma omp parallel for
  for(int k=0;k<Nxx_plus_2NGHOSTS2;k++) for(int j=1;j<Nxx_plus_2NGHOSTS1-2;j++) for(int i=1;i<Nxx_plus_2NGHOSTS0-2;i++) {
        temporary[IDX3S(i,j,k)]=
          IPH(IPH(psi6center[IDX3S(i-1,j-1,k)],psi6center[IDX3S(i,j-1,k)],psi6center[IDX3S(i+1,j-1,k)],psi6center[IDX3S(i+2,j-1,k)]),
              IPH(psi6center[IDX3S(i-1,j  ,k)],psi6center[IDX3S(i,j  ,k)],psi6center[IDX3S(i+1,j  ,k)],psi6center[IDX3S(i+2,j  ,k)]),
              IPH(psi6center[IDX3S(i-1,j+1,k)],psi6center[IDX3S(i,j+1,k)],psi6center[IDX3S(i+1,j+1,k)],psi6center[IDX3S(i+2,j+1,k)]),
              IPH(psi6center[IDX3S(i-1,j+2,k)],psi6center[IDX3S(i,j+2,k)],psi6center[IDX3S(i+1,j+2,k)],psi6center[IDX3S(i+2,j+2,k)]));
      }

  int A_directionz=3;
  A_i_rhs_no_gauge_terms(A_directionz,params,out_prims_r,out_prims_l,temporary,
                         auxevol_gfs+Nxxp2NG012*CMAX_XGF,
                         auxevol_gfs+Nxxp2NG012*CMIN_XGF,
                         auxevol_gfs+Nxxp2NG012*CMAX_YGF,
                         auxevol_gfs+Nxxp2NG012*CMIN_YGF,
                         rhs_gfs+Nxxp2NG012*AD2GF);

  // in_prims[{VYR,VYL,VZR,VZL}].gz_{lo,hi} ghostzones are not correct, so we fix
  //    this below.
  // [Note that this is a cheap operation, copying only 8 integers and a pointer.]
  in_prims[VYR]=out_prims_r[VY];
  in_prims[VYL]=out_prims_l[VY];
  in_prims[VZR]=out_prims_r[VZ];
  in_prims[VZL]=out_prims_l[VZ];

  flux_dirn=2;
  // ftilde = 0 in GRFFE, since P=rho=0.

  /* There are two stories going on here:
   * 1) Single reconstruction to (i,j,k-1/2) for {vx,vy,vz,Bx,By,Bz} to compute
   *    z-dir'n advection terms in \partial_t {mhd_st_{x,y,z}} at (i,j,k)
   * 2) Multiple reconstructions for *staggered* gridfunctions A_i:
   *    \partial_t A_i = \epsilon_{ijk} \psi^{6} (v^j B^k - v^j B^k),
   *    where \epsilon_{ijk} is the flat-space antisymmetric operator.
   * 2A) Ax_rhs is defined at (i,j+1/2,k+1/2), depends on v{y,z} and B{y,z}
   * 2Aa) v{y,z}{r,l} are at (i,j-1/2,k), so we reconstruct here to (i,j-1/2,k-1/2)
   * 2Ab) Bz_stagger{r,l} are at (i,j-1/2,k+1/2) already.
   * 2Ac) By_stagger is at (i,j+1/2,k), and below we reconstruct its value at (i,j+1/2,k-1/2)
   * 2B) Ay_rhs is defined at (i+1/2,j,k+1/2), depends on v{z,x} and B{z,x}.
   * 2Ba) v{x,z} are reconstructed to (i,j,k-1/2). Later we'll reconstruct again to (i-1/2,j,k-1/2).
   * 2Bb) Bz_stagger is at (i,j,k+1/2). Later we will reconstruct to (i-1/2,j,k+1/2).
   * 2Bc) Bx_stagger is at (i+1/2,j,k), and below we reconstruct its value at (i+1/2,j,k-1/2)
   */
  ww=0;
  // NOTE! The order of variable reconstruction is important here,
  //   as we don't want to overwrite {vxr,vxl,vyr,vyl}!
  which_prims_to_reconstruct[ww]=VYR;       ww++;
  which_prims_to_reconstruct[ww]=VZR;       ww++;
  which_prims_to_reconstruct[ww]=VYL;       ww++;
  which_prims_to_reconstruct[ww]=VZL;       ww++;
  num_prims_to_reconstruct=ww;
#ifdef WORKAROUND_ENABLED
  workaround_Valencia_to_Drift_velocity(params,auxevol_gfs+Nxxp2NG012*VALENCIAV_RU1GF,auxevol_gfs+Nxxp2NG012*ALPHA_FACEGF,auxevol_gfs+Nxxp2NG012*BETA_FACEU1GF,flux_dirn+1);
  workaround_Valencia_to_Drift_velocity(params,auxevol_gfs+Nxxp2NG012*VALENCIAV_RU2GF,auxevol_gfs+Nxxp2NG012*ALPHA_FACEGF,auxevol_gfs+Nxxp2NG012*BETA_FACEU2GF,flux_dirn+1);
  workaround_Valencia_to_Drift_velocity(params,auxevol_gfs+Nxxp2NG012*VALENCIAV_LU1GF,auxevol_gfs+Nxxp2NG012*ALPHA_FACEGF,auxevol_gfs+Nxxp2NG012*BETA_FACEU1GF,flux_dirn+1);
  workaround_Valencia_to_Drift_velocity(params,auxevol_gfs+Nxxp2NG012*VALENCIAV_LU2GF,auxevol_gfs+Nxxp2NG012*ALPHA_FACEGF,auxevol_gfs+Nxxp2NG012*BETA_FACEU2GF,flux_dirn+1);
#endif /*WORKAROUND_ENABLED*/
  // This function is housed in the file: "reconstruct_set_of_prims_PPM_GRFFE.C"
  reconstruct_set_of_prims_PPM_GRFFE_NRPy(params, auxevol_gfs, flux_dirn+1, num_prims_to_reconstruct,
                                          which_prims_to_reconstruct, in_prims, out_prims_r, out_prims_l, temporary);
#ifdef WORKAROUND_ENABLED
  workaround_Drift_to_Valencia_velocity(params,auxevol_gfs+Nxxp2NG012*VALENCIAV_RU1GF,auxevol_gfs+Nxxp2NG012*ALPHA_FACEGF,auxevol_gfs+Nxxp2NG012*BETA_FACEU1GF,flux_dirn+1);
  workaround_Drift_to_Valencia_velocity(params,auxevol_gfs+Nxxp2NG012*VALENCIAV_RU2GF,auxevol_gfs+Nxxp2NG012*ALPHA_FACEGF,auxevol_gfs+Nxxp2NG012*BETA_FACEU2GF,flux_dirn+1);
  workaround_Drift_to_Valencia_velocity(params,auxevol_gfs+Nxxp2NG012*VALENCIAV_LU1GF,auxevol_gfs+Nxxp2NG012*ALPHA_FACEGF,auxevol_gfs+Nxxp2NG012*BETA_FACEU1GF,flux_dirn+1);
  workaround_Drift_to_Valencia_velocity(params,auxevol_gfs+Nxxp2NG012*VALENCIAV_LU2GF,auxevol_gfs+Nxxp2NG012*ALPHA_FACEGF,auxevol_gfs+Nxxp2NG012*BETA_FACEU2GF,flux_dirn+1);
#endif /*WORKAROUND_ENABLED*/
  interpolate_metric_gfs_to_cell_faces(params,auxevol_gfs,flux_dirn+1);
  // Reconstruct other primitives last!
  ww=0;
  which_prims_to_reconstruct[ww]=VX;        ww++;
  which_prims_to_reconstruct[ww]=VY;        ww++;
  which_prims_to_reconstruct[ww]=VZ;        ww++;
  which_prims_to_reconstruct[ww]=BX_CENTER; ww++;
  which_prims_to_reconstruct[ww]=BY_CENTER; ww++;
  //which_prims_to_reconstruct[ww]=BZ_CENTER; ww++;
  which_prims_to_reconstruct[ww]=BX_STAGGER; ww++;
  which_prims_to_reconstruct[ww]=BY_STAGGER; ww++;
  num_prims_to_reconstruct=ww;
  // This function is housed in the file: "reconstruct_set_of_prims_PPM_GRFFE.C"
  reconstruct_set_of_prims_PPM_GRFFE_NRPy(params, auxevol_gfs, flux_dirn+1, num_prims_to_reconstruct,
                                          which_prims_to_reconstruct, in_prims, out_prims_r, out_prims_l, temporary);
  //Right and left face values of BI_CENTER are used in GRFFE__S_i__flux computation (first to compute b^a).
  //   Instead of reconstructing, we simply set B^z face values to be consistent with BZ_STAGGER.
#pragma omp parallel for
  for(int k=0;k<Nxx_plus_2NGHOSTS2;k++) for(int j=0;j<Nxx_plus_2NGHOSTS1;j++) for(int i=0;i<Nxx_plus_2NGHOSTS0;i++) {
        const int index=IDX3S(i,j,k), indexkm1=IDX3S(i,j,k-1+(k==0)); /* indexkm1=0 when k=0 */
        out_prims_r[BZ_CENTER].gf[index]=out_prims_l[BZ_CENTER].gf[index]=in_prims[BZ_STAGGER].gf[indexkm1]; }
  // Then add fluxes to RHS for hydro variables {vx,vy,vz}:
  // This function is housed in the file: "add_fluxes_and_source_terms_to_hydro_rhss.C"
  calculate_StildeD2_source_term(params,auxevol_gfs,rhs_gfs);
  calculate_Stilde_flux_D2(params,auxevol_gfs,rhs_gfs);
  calculate_Stilde_rhsD(flux_dirn+1,params,auxevol_gfs,rhs_gfs);

  // in_prims[{VYR,VYL,VZR,VZL}].gz_{lo,hi} ghostzones are not set correcty.
  //    We fix this below.
  // [Note that this is a cheap operation, copying only 8 integers and a pointer.]
  in_prims[VXR]=out_prims_r[VX];
  in_prims[VZR]=out_prims_r[VZ];
  in_prims[VXL]=out_prims_l[VX];
  in_prims[VZL]=out_prims_l[VZ];
  // FIXME: lines above seem to be inconsistent with lines below.... Possible bug, not major enough to affect evolutions though.
  in_prims[VZR].gz_lo[1]=in_prims[VZR].gz_hi[1]=0;
  in_prims[VXR].gz_lo[1]=in_prims[VXR].gz_hi[1]=0;
  in_prims[VZL].gz_lo[1]=in_prims[VZL].gz_hi[1]=0;
  in_prims[VXL].gz_lo[1]=in_prims[VXL].gz_hi[1]=0;
  /*****************************************
   * COMPUTING RHS OF A_x, BOOKKEEPING NOTE:
   * We want to compute
   * \partial_t A_x - [gauge terms] = \psi^{6} (v^y B^z - v^z B^y).
   * A_x is defined at (i,j+1/2,k+1/2).
   * ==========================
   * Where defined  | Variables
   * (i,j-1/2,k-1/2)| {vyrr,vyrl,vylr,vyll,vzrr,vzrl,vzlr,vzll}
   * (i,j+1/2,k-1/2)| {By_stagger_r,By_stagger_l} (see Table 1 in arXiv:1007.2848)
   * (i,j-1/2,k+1/2)| {Bz_stagger_r,Bz_stagger_l} (see Table 1 in arXiv:1007.2848)
   * (i,j,k)        | {phi}
   * ==========================
   ******************************************/
  // Next compute phi at (i,j+1/2,k+1/2):
#pragma omp parallel for
  for(int k=1;k<Nxx_plus_2NGHOSTS2-2;k++) for(int j=1;j<Nxx_plus_2NGHOSTS1-2;j++) for(int i=0;i<Nxx_plus_2NGHOSTS0;i++) {
        temporary[IDX3S(i,j,k)]=
          IPH(IPH(psi6center[IDX3S(i,j-1,k-1)],psi6center[IDX3S(i,j,k-1)],psi6center[IDX3S(i,j+1,k-1)],psi6center[IDX3S(i,j+2,k-1)]),
              IPH(psi6center[IDX3S(i,j-1,k  )],psi6center[IDX3S(i,j,k  )],psi6center[IDX3S(i,j+1,k  )],psi6center[IDX3S(i,j+2,k  )]),
              IPH(psi6center[IDX3S(i,j-1,k+1)],psi6center[IDX3S(i,j,k+1)],psi6center[IDX3S(i,j+1,k+1)],psi6center[IDX3S(i,j+2,k+1)]),
              IPH(psi6center[IDX3S(i,j-1,k+2)],psi6center[IDX3S(i,j,k+2)],psi6center[IDX3S(i,j+1,k+2)],psi6center[IDX3S(i,j+2,k+2)]));
      }

  int A_directionx=1;
  A_i_rhs_no_gauge_terms(A_directionx,params,out_prims_r,out_prims_l,temporary,
                         auxevol_gfs+Nxxp2NG012*CMAX_YGF,
                         auxevol_gfs+Nxxp2NG012*CMIN_YGF,
                         auxevol_gfs+Nxxp2NG012*CMAX_ZGF,
                         auxevol_gfs+Nxxp2NG012*CMIN_ZGF,
                         rhs_gfs+Nxxp2NG012*AD0GF);

  // We reprise flux_dirn=0 to finish up computations of Ai_rhs's!
  flux_dirn=0;
  // ftilde = 0 in GRFFE, since P=rho=0.

  ww=0;
  // NOTE! The order of variable reconstruction is important here,
  //   as we don't want to overwrite {vxr,vxl,vyr,vyl}!
  which_prims_to_reconstruct[ww]=VXR;       ww++;
  which_prims_to_reconstruct[ww]=VZR;       ww++;
  which_prims_to_reconstruct[ww]=VXL;       ww++;
  which_prims_to_reconstruct[ww]=VZL;       ww++;
  which_prims_to_reconstruct[ww]=BZ_STAGGER;ww++;
  num_prims_to_reconstruct=ww;
#ifdef WORKAROUND_ENABLED
  workaround_Valencia_to_Drift_velocity(params,auxevol_gfs+Nxxp2NG012*VALENCIAV_RU0GF,auxevol_gfs+Nxxp2NG012*ALPHA_FACEGF,auxevol_gfs+Nxxp2NG012*BETA_FACEU0GF,flux_dirn+1);
  workaround_Valencia_to_Drift_velocity(params,auxevol_gfs+Nxxp2NG012*VALENCIAV_RU2GF,auxevol_gfs+Nxxp2NG012*ALPHA_FACEGF,auxevol_gfs+Nxxp2NG012*BETA_FACEU2GF,flux_dirn+1);
  workaround_Valencia_to_Drift_velocity(params,auxevol_gfs+Nxxp2NG012*VALENCIAV_LU0GF,auxevol_gfs+Nxxp2NG012*ALPHA_FACEGF,auxevol_gfs+Nxxp2NG012*BETA_FACEU0GF,flux_dirn+1);
  workaround_Valencia_to_Drift_velocity(params,auxevol_gfs+Nxxp2NG012*VALENCIAV_LU2GF,auxevol_gfs+Nxxp2NG012*ALPHA_FACEGF,auxevol_gfs+Nxxp2NG012*BETA_FACEU2GF,flux_dirn+1);
#endif /*WORKAROUND_ENABLED*/
  // This function is housed in the file: "reconstruct_set_of_prims_PPM_GRFFE.C"
  reconstruct_set_of_prims_PPM_GRFFE_NRPy(params, auxevol_gfs, flux_dirn+1, num_prims_to_reconstruct,
                                          which_prims_to_reconstruct, in_prims, out_prims_r, out_prims_l, temporary);
#ifdef WORKAROUND_ENABLED
  workaround_Drift_to_Valencia_velocity(params,auxevol_gfs+Nxxp2NG012*VALENCIAV_RU0GF,auxevol_gfs+Nxxp2NG012*ALPHA_FACEGF,auxevol_gfs+Nxxp2NG012*BETA_FACEU0GF,flux_dirn+1);
  workaround_Drift_to_Valencia_velocity(params,auxevol_gfs+Nxxp2NG012*VALENCIAV_RU2GF,auxevol_gfs+Nxxp2NG012*ALPHA_FACEGF,auxevol_gfs+Nxxp2NG012*BETA_FACEU2GF,flux_dirn+1);
  workaround_Drift_to_Valencia_velocity(params,auxevol_gfs+Nxxp2NG012*VALENCIAV_LU0GF,auxevol_gfs+Nxxp2NG012*ALPHA_FACEGF,auxevol_gfs+Nxxp2NG012*BETA_FACEU0GF,flux_dirn+1);
  workaround_Drift_to_Valencia_velocity(params,auxevol_gfs+Nxxp2NG012*VALENCIAV_LU2GF,auxevol_gfs+Nxxp2NG012*ALPHA_FACEGF,auxevol_gfs+Nxxp2NG012*BETA_FACEU2GF,flux_dirn+1);
#endif /*WORKAROUND_ENABLED*/

  /*****************************************
   * COMPUTING RHS OF A_y, BOOKKEEPING NOTE:
   * We want to compute
   * \partial_t A_y - [gauge terms] = \psi^{6} (v^z B^x - v^x B^z).
   * A_y is defined at (i+1/2,j,k+1/2).
   * ==========================
   * Where defined  | Variables
   * (i-1/2,j,k-1/2)| {vyrr,vyrl,vylr,vyll,vzrr,vzrl,vzlr,vzll}
   * (i+1/2,j,k-1/2)| {By_stagger_r,By_stagger_l} (see Table 1 in arXiv:1007.2848)
   * (i-1/2,j,k+1/2)| {Bz_stagger_r,Bz_stagger_l} (see Table 1 in arXiv:1007.2848)
   * (i,j,k)        | {phi}
   * ==========================
   ******************************************/
  // Next compute phi at (i+1/2,j,k+1/2):
#pragma omp parallel for
  for(int k=1;k<Nxx_plus_2NGHOSTS2-2;k++) for(int j=0;j<Nxx_plus_2NGHOSTS1;j++) for(int i=1;i<Nxx_plus_2NGHOSTS0-2;i++) {
        temporary[IDX3S(i,j,k)]=
          IPH(IPH(psi6center[IDX3S(i-1,j,k-1)],psi6center[IDX3S(i,j,k-1)],psi6center[IDX3S(i+1,j,k-1)],psi6center[IDX3S(i+2,j,k-1)]),
              IPH(psi6center[IDX3S(i-1,j,k  )],psi6center[IDX3S(i,j,k  )],psi6center[IDX3S(i+1,j,k  )],psi6center[IDX3S(i+2,j,k  )]),
              IPH(psi6center[IDX3S(i-1,j,k+1)],psi6center[IDX3S(i,j,k+1)],psi6center[IDX3S(i+1,j,k+1)],psi6center[IDX3S(i+2,j,k+1)]),
              IPH(psi6center[IDX3S(i-1,j,k+2)],psi6center[IDX3S(i,j,k+2)],psi6center[IDX3S(i+1,j,k+2)],psi6center[IDX3S(i+2,j,k+2)]));
      }

  int A_directiony=2;
  A_i_rhs_no_gauge_terms(A_directiony,params,out_prims_r,out_prims_l,temporary,
                         auxevol_gfs+Nxxp2NG012*CMAX_ZGF,
                         auxevol_gfs+Nxxp2NG012*CMIN_ZGF,
                         auxevol_gfs+Nxxp2NG012*CMAX_XGF,
                         auxevol_gfs+Nxxp2NG012*CMIN_XGF,
                         rhs_gfs+Nxxp2NG012*AD1GF);

  // Next compute psi6phi_rhs, and add gauge terms to A_i_rhs terms!
  //   Note that in the following function, we don't bother with reconstruction, instead interpolating.
  // We need A^i, but only have A_i. So we add gtupij to the list of input variables.
  REAL *interp_vars[MAXNUMINTERP];
  ww=0;
  interp_vars[ww]=auxevol_gfs+Nxxp2NG012*BETAU0GF;   ww++;
  interp_vars[ww]=auxevol_gfs+Nxxp2NG012*BETAU1GF;   ww++;
  interp_vars[ww]=auxevol_gfs+Nxxp2NG012*BETAU2GF;   ww++;
  interp_vars[ww]=auxevol_gfs+Nxxp2NG012*GAMMADD00GF;  ww++;
  interp_vars[ww]=auxevol_gfs+Nxxp2NG012*GAMMADD01GF;  ww++;
  interp_vars[ww]=auxevol_gfs+Nxxp2NG012*GAMMADD02GF;  ww++;
  interp_vars[ww]=auxevol_gfs+Nxxp2NG012*GAMMADD11GF;  ww++;
  interp_vars[ww]=auxevol_gfs+Nxxp2NG012*GAMMADD12GF;  ww++;
  interp_vars[ww]=auxevol_gfs+Nxxp2NG012*GAMMADD22GF;  ww++;
  interp_vars[ww]=temporary;ww++;
  interp_vars[ww]=auxevol_gfs+Nxxp2NG012*ALPHAGF;   ww++;
  interp_vars[ww]=in_gfs+Nxxp2NG012*AD0GF;      ww++;
  interp_vars[ww]=in_gfs+Nxxp2NG012*AD1GF;      ww++;
  interp_vars[ww]=in_gfs+Nxxp2NG012*AD2GF;      ww++;
  const int max_num_interp_variables=ww;
//   if(max_num_interp_variables>MAXNUMINTERP) {CCTK_VError(VERR_DEF_PARAMS,"Error: Didn't allocate enough space for interp_vars[]."); }
  // We are FINISHED with v{x,y,z}{r,l} and P{r,l} so we use these 8 gridfunctions' worth of space as temp storage.
  Lorenz_psi6phi_rhs__add_gauge_terms_to_A_i_rhs(params,interp_vars,
                                                 in_gfs+Nxxp2NG012*PSI6PHIGF,
                                                 auxevol_gfs+Nxxp2NG012*VALENCIAV_RU0GF,  // WARNING:
                                                 auxevol_gfs+Nxxp2NG012*VALENCIAV_RU1GF,  // ALL VARIABLES
                                                 auxevol_gfs+Nxxp2NG012*VALENCIAV_RU2GF,  // ON THESE LINES
                                                 auxevol_gfs+Nxxp2NG012*VALENCIAV_LU0GF,  // ARE OVERWRITTEN
                                                 auxevol_gfs+Nxxp2NG012*VALENCIAV_LU1GF,  // FOR TEMP STORAGE
                                                 auxevol_gfs+Nxxp2NG012*VALENCIAV_LU2GF,  // .
                                                 auxevol_gfs+Nxxp2NG012*VALENCIAV_RRU0GF, // .
                                                 auxevol_gfs+Nxxp2NG012*VALENCIAV_RLU0GF, // .
                                                 rhs_gfs+Nxxp2NG012*PSI6PHIGF,
                                                 rhs_gfs+Nxxp2NG012*AD0GF,
                                                 rhs_gfs+Nxxp2NG012*AD1GF,
                                                 rhs_gfs+Nxxp2NG012*AD2GF);
/*#pragma omp parallel for
  for(int k=0;k<Nxx_plus_2NGHOSTS2;k++) for(int j=0;j<Nxx_plus_2NGHOSTS1;j++) for(int i=0;i<Nxx_plus_2NGHOSTS0;i++) {
        REAL x = xx[0][i];
        REAL y = xx[1][j];
        REAL z = xx[2][k];
        if(sqrt(x*x+y*y+z*z)<min_radius_inside_of_which_conserv_to_prims_FFE_and_FFE_evolution_is_DISABLED) {
          rhs_gfs[IDX4S(STILDED0GF,i,j,k)] = 0.0;
          rhs_gfs[IDX4S(STILDED1GF,i,j,k)] = 0.0;
          rhs_gfs[IDX4S(STILDED2GF,i,j,k)] = 0.0;

          rhs_gfs[IDX4S(PSI6PHIGF,i,j,k)] = 0.0;
          rhs_gfs[IDX4S(AD0GF,i,j,k)] = 0.0;
          rhs_gfs[IDX4S(AD1GF,i,j,k)] = 0.0;
          rhs_gfs[IDX4S(AD2GF,i,j,k)] = 0.0;
        }
      }*/
}

void GiRaFFE_NRPy_post_step(const paramstruct *restrict params,REAL *xx[3],REAL *restrict auxevol_gfs,REAL *restrict evol_gfs,const int n) {
#include "set_Cparameters.h"
    // First, apply BCs to AD and psi6Phi. Then calculate BU from AD
    apply_bcs_potential(params,evol_gfs);
    const int Nxxp2NG012 = Nxx_plus_2NGHOSTS0*Nxx_plus_2NGHOSTS1*Nxx_plus_2NGHOSTS2;
    GiRaFFE_compute_B_and_Bstagger_from_A(params,
                                          auxevol_gfs+Nxxp2NG012*GAMMADD00GF,
                                          auxevol_gfs+Nxxp2NG012*GAMMADD01GF,
                                          auxevol_gfs+Nxxp2NG012*GAMMADD02GF,
                                          auxevol_gfs+Nxxp2NG012*GAMMADD11GF,
                                          auxevol_gfs+Nxxp2NG012*GAMMADD12GF,
                                          auxevol_gfs+Nxxp2NG012*GAMMADD22GF,
                                          auxevol_gfs + Nxxp2NG012*PSI6_TEMPGF, /* Temporary storage */
                                          evol_gfs+Nxxp2NG012*AD0GF,
                                          evol_gfs+Nxxp2NG012*AD1GF,
                                          evol_gfs+Nxxp2NG012*AD2GF,
                                          auxevol_gfs+Nxxp2NG012*BU0GF,
                                          auxevol_gfs+Nxxp2NG012*BU1GF,
                                          auxevol_gfs+Nxxp2NG012*BU2GF,
                                          auxevol_gfs+Nxxp2NG012*BSTAGGERU0GF,
                                          auxevol_gfs+Nxxp2NG012*BSTAGGERU1GF,
                                          auxevol_gfs+Nxxp2NG012*BSTAGGERU2GF);
    //override_BU_with_old_GiRaFFE(params,auxevol_gfs,n);
    // Apply fixes to StildeD, then recompute the velocity at the new timestep.
    // Apply the current sheet prescription to the velocities
    GiRaFFE_NRPy_cons_to_prims(params,xx,auxevol_gfs,evol_gfs);
    // Then, recompute StildeD to be consistent with the new velocities
    //GiRaFFE_NRPy_prims_to_cons(params,auxevol_gfs,evol_gfs);
    // Finally, apply outflow boundary conditions to the velocities.
    apply_bcs_velocity(params,auxevol_gfs);
}
""")
