# The A-to-B driver
import shutil, os, sys           # Standard Python modules for multiplatform OS-level functions
# First, we'll add the parent directory to the list of directories Python will check for modules.
nrpy_dir_path = os.path.join("..")
if nrpy_dir_path not in sys.path:
    sys.path.append(nrpy_dir_path)

import cmdline_helper as cmd     # NRPy+: Multi-platform Python command-line interface

# Step 1: The A-to-B driver
from outputC import *            # NRPy+: Core C code output module
import finite_difference as fin  # NRPy+: Finite difference C code generation module
import NRPy_param_funcs as par   # NRPy+: Parameter interface
import grid as gri               # NRPy+: Functions having to do with numerical grids
import loop as lp                # NRPy+: Generate C code loops
import indexedexp as ixp         # NRPy+: Symbolic indexed expression (e.g., tensors, vectors, etc.) support
import reference_metric as rfm   # NRPy+: Reference metric support

global GiRaFFE_NRPy_A2B
def GiRaFFE_NRPy_A2B(outdir,gammaDD,AD,BU):
    cmd.mkdir(outdir)
    # Set spatial dimension (must be 3 for BSSN)
    DIM = 3
    par.set_parval_from_str("grid::DIM",DIM)
    # Compute the sqrt of the three metric determinant.
    import GRHD.equations as gh
    gh.compute_sqrtgammaDET(gammaDD)

    # Import the Levi-Civita symbol and build the corresponding tensor.
    # We already have a handy function to define the Levi-Civita symbol in WeylScalars
    import WeylScal4NRPy.WeylScalars_Cartesian as weyl
    LeviCivitaDDD = weyl.define_LeviCivitaSymbol_rank3()
    LeviCivitaUUU = ixp.zerorank3()
    for i in range(DIM):
        for j in range(DIM):
            for k in range(DIM):
                LCijk = LeviCivitaDDD[i][j][k]
                #LeviCivitaDDD[i][j][k] = LCijk * sp.sqrt(gho.gammadet)
                LeviCivitaUUU[i][j][k] = LCijk / gh.sqrtgammaDET

    AD_dD = ixp.declarerank2("AD_dD","nosym")
    BU = ixp.zerorank1() 
    for i in range(DIM):
        for j in range(DIM):
            for k in range(DIM):
                BU[i] += LeviCivitaUUU[i][j][k] * AD_dD[k][j]

    # Write the code to compute derivatives with shifted stencils as needed.
    with open(os.path.join(outdir,"driver_AtoB.h"),"w") as file:
        file.write("""void compute_A2B_in_ghostzones(const paramstruct *restrict params,REAL *restrict in_gfs,REAL *restrict auxevol_gfs,
                                      const int i0min,const int i0max, 
                                      const int i1min,const int i1max, 
                                      const int i2min,const int i2max) {
#include "../set_Cparameters.h"
    for(int i2=i2min;i2<i2max;i2++) for(int i1=i1min;i1<i1max;i1++) for(int i0=i0min;i0<i0max;i0++) {
        REAL dx_Ay,dx_Az,dy_Ax,dy_Az,dz_Ax,dz_Ay;
        // Check to see if we're on the +x or -x face. If so, use a downwinded- or upwinded-stencil, respectively.
        // Otherwise, use a centered stencil.
        if (i0 > 0 && i0 < Nxx_plus_2NGHOSTS0-1) {
            dx_Ay = invdx0*(-1.0/2.0*in_gfs[IDX4S(AD1GF, i0-1,i1,i2)] + (1.0/2.0)*in_gfs[IDX4S(AD1GF, i0+1,i1,i2)]);
            dx_Az = invdx0*(-1.0/2.0*in_gfs[IDX4S(AD2GF, i0-1,i1,i2)] + (1.0/2.0)*in_gfs[IDX4S(AD2GF, i0+1,i1,i2)]);
        }
        else if (i0==0) {
            dx_Ay = invdx0*(-3.0/2.0*in_gfs[IDX4S(AD1GF, i0,i1,i2)] + 2*in_gfs[IDX4S(AD1GF, i0+1,i1,i2)] - 1.0/2.0*in_gfs[IDX4S(AD1GF, i0+2,i1,i2)]);
            dx_Az = invdx0*(-3.0/2.0*in_gfs[IDX4S(AD2GF, i0,i1,i2)] + 2*in_gfs[IDX4S(AD2GF, i0+1,i1,i2)] - 1.0/2.0*in_gfs[IDX4S(AD2GF, i0+2,i1,i2)]);
        }
        else {
            dx_Ay = invdx0*((3.0/2.0)*in_gfs[IDX4S(AD1GF, i0,i1,i2)] - 2*in_gfs[IDX4S(AD1GF, i0-1,i1,i2)] + (1.0/2.0)*in_gfs[IDX4S(AD1GF, i0-2,i1,i2)]);
            dx_Az = invdx0*((3.0/2.0)*in_gfs[IDX4S(AD2GF, i0,i1,i2)] - 2*in_gfs[IDX4S(AD2GF, i0-1,i1,i2)] + (1.0/2.0)*in_gfs[IDX4S(AD2GF, i0-2,i1,i2)]);
        }
        // As above, but in the y direction.
        if (i1 > 0 && i1 < Nxx_plus_2NGHOSTS1-1) {
            dy_Ax = invdx1*(-1.0/2.0*in_gfs[IDX4S(AD0GF, i0,i1-1,i2)] + (1.0/2.0)*in_gfs[IDX4S(AD0GF, i0,i1+1,i2)]);
            dy_Az = invdx1*(-1.0/2.0*in_gfs[IDX4S(AD2GF, i0,i1-1,i2)] + (1.0/2.0)*in_gfs[IDX4S(AD2GF, i0,i1+1,i2)]);
        }
        else if (i1==0) {
            dy_Ax = invdx1*(-3.0/2.0*in_gfs[IDX4S(AD0GF, i0,i1,i2)] + 2*in_gfs[IDX4S(AD0GF, i0,i1+1,i2)] - 1.0/2.0*in_gfs[IDX4S(AD0GF, i0,i1+2,i2)]);
            dy_Az = invdx1*(-3.0/2.0*in_gfs[IDX4S(AD2GF, i0,i1,i2)] + 2*in_gfs[IDX4S(AD2GF, i0,i1+1,i2)] - 1.0/2.0*in_gfs[IDX4S(AD2GF, i0,i1+2,i2)]);
        }
        else {
            dy_Ax = invdx1*((3.0/2.0)*in_gfs[IDX4S(AD0GF, i0,i1,i2)] - 2*in_gfs[IDX4S(AD0GF, i0,i1-1,i2)] + (1.0/2.0)*in_gfs[IDX4S(AD0GF, i0,i1-2,i2)]);
            dy_Az = invdx1*((3.0/2.0)*in_gfs[IDX4S(AD2GF, i0,i1,i2)] - 2*in_gfs[IDX4S(AD2GF, i0,i1-1,i2)] + (1.0/2.0)*in_gfs[IDX4S(AD2GF, i0,i1-2,i2)]);
        }
        // As above, but in the z direction.
        if (i2 > 0 && i2 < Nxx_plus_2NGHOSTS2-1) {
            dz_Ax = invdx2*(-1.0/2.0*in_gfs[IDX4S(AD0GF, i0,i1,i2-1)] + (1.0/2.0)*in_gfs[IDX4S(AD0GF, i0,i1,i2+1)]);
            dz_Ay = invdx2*(-1.0/2.0*in_gfs[IDX4S(AD1GF, i0,i1,i2-1)] + (1.0/2.0)*in_gfs[IDX4S(AD1GF, i0,i1,i2+1)]);
        }
        else if (i2==0) {
            dz_Ax = invdx2*(-3.0/2.0*in_gfs[IDX4S(AD0GF, i0,i1,i2)] + 2*in_gfs[IDX4S(AD0GF, i0,i1,i2+1)] - 1.0/2.0*in_gfs[IDX4S(AD0GF, i0,i1,i2+2)]);
            dz_Ay = invdx2*(-3.0/2.0*in_gfs[IDX4S(AD1GF, i0,i1,i2)] + 2*in_gfs[IDX4S(AD1GF, i0,i1,i2+1)] - 1.0/2.0*in_gfs[IDX4S(AD1GF, i0,i1,i2+2)]);
        }
        else {
            dz_Ax = invdx2*((3.0/2.0)*in_gfs[IDX4S(AD0GF, i0,i1,i2)] - 2*in_gfs[IDX4S(AD0GF, i0,i1,i2-1)] + (1.0/2.0)*in_gfs[IDX4S(AD0GF, i0,i1,i2-2)]);
            dz_Ay = invdx2*((3.0/2.0)*in_gfs[IDX4S(AD1GF, i0,i1,i2)] - 2*in_gfs[IDX4S(AD1GF, i0,i1,i2-1)] + (1.0/2.0)*in_gfs[IDX4S(AD1GF, i0,i1,i2-2)]);
        }
        // Compute the magnetic field in the normal way, using the previously calculated derivatives.
        const double gammaDD00 = auxevol_gfs[IDX4S(GAMMADD00GF, i0,i1,i2)];
        const double gammaDD01 = auxevol_gfs[IDX4S(GAMMADD01GF, i0,i1,i2)];
        const double gammaDD02 = auxevol_gfs[IDX4S(GAMMADD02GF, i0,i1,i2)];
        const double gammaDD11 = auxevol_gfs[IDX4S(GAMMADD11GF, i0,i1,i2)];
        const double gammaDD12 = auxevol_gfs[IDX4S(GAMMADD12GF, i0,i1,i2)];
        const double gammaDD22 = auxevol_gfs[IDX4S(GAMMADD22GF, i0,i1,i2)];
        /* 
        * NRPy+ Finite Difference Code Generation, Step 2 of 1: Evaluate SymPy expressions and write to main memory:
        */
        const double invsqrtg = 1.0/sqrt(gammaDD00*gammaDD11*gammaDD22 
                                       - gammaDD00*gammaDD12*gammaDD12
                                       + 2*gammaDD01*gammaDD02*gammaDD12
                                       - gammaDD11*gammaDD02*gammaDD02
                                       - gammaDD22*gammaDD01*gammaDD01);
        auxevol_gfs[IDX4S(BU0GF, i0,i1,i2)] = (dy_Az-dz_Ay)*invsqrtg;
        auxevol_gfs[IDX4S(BU1GF, i0,i1,i2)] = (dz_Ax-dx_Az)*invsqrtg;
        auxevol_gfs[IDX4S(BU2GF, i0,i1,i2)] = (dx_Ay-dy_Ax)*invsqrtg;
    }
}
""")
        
    # Now, we'll also write some more auxiliary functions to handle the order-lowering method for A2B
    with open(os.path.join(outdir,"driver_AtoB.h"),"a") as file:
        file.write("""REAL relative_error(REAL a, REAL b) {
    if((a+b)!=0.0) {
        return 2.0*fabs(a-b)/fabs(a+b);
    }
    else {
        return 0.0;
    }
}

#define M2 0
#define M1 1
#define P0 2
#define P1 3
#define P2 4
#define FOURTH 0
#define SECOND 1
#define FORWARD 2
#define BACKWARD 3
void compute_Bx_pointwise(REAL *Bx, const REAL invdy, const REAL *Ay, const REAL invdz, const REAL *Az) {
    REAL dz_Ay,dy_Az;
    dz_Ay = invdz*((Ay[P1]-Ay[M1])*2.0/3.0 - (Ay[P2]-Ay[M2])/12.0);
    dy_Az = invdy*((Az[P1]-Az[M1])*2.0/3.0 - (Az[P2]-Az[M2])/12.0);
    Bx[FOURTH] = dy_Az - dz_Ay;
    
    dz_Ay = invdz*(Ay[P1]-Ay[M1])/2.0;
    dy_Az = invdy*(Az[P1]-Az[M1])/2.0;
    Bx[SECOND] = dy_Az - dz_Ay;
    
    dz_Ay = invdz*(Ay[P1]-Ay[P0]);
    dy_Az = invdy*(Az[P1]-Az[P0]);
    Bx[FORWARD] = dy_Az - dz_Ay;
    
    dz_Ay = invdz*(Ay[P0]-Ay[M1]);
    dy_Az = invdy*(Az[P0]-Az[M1]);
    Bx[BACKWARD] = dy_Az - dz_Ay;
}

#define TOLERANCE_A2B 1.0e-4
REAL find_accepted_Bx_order(REAL *Bx) {
    REAL accepted_val = Bx[FOURTH];
    REAL Rel_error_o1F_vs_o2 = relative_error(Bx[FORWARD], Bx[SECOND]);
    REAL Rel_error_o1B_vs_o2 = relative_error(Bx[BACKWARD],Bx[SECOND]);
    REAL Rel_error_o1_vs_o2 = 0.5*(Rel_error_o1F_vs_o2 + Rel_error_o1B_vs_o2 );
    if(relative_error(Bx[SECOND],Bx[FOURTH]) > 100*TOLERANCE_A2B) {
        accepted_val = Bx[SECOND];
        if(Rel_error_o1_vs_o2 > TOLERANCE_A2B) {
            if(Rel_error_o1_vs_o2 > Rel_error_o1B_vs_o2) {
                accepted_val = Bx[BACKWARD];
            }
            else {
                accepted_val = Bx[FORWARD];
            }
        }
    }
    return accepted_val;
}
""")

    order_lowering_body = """REAL AD0_1[5],AD0_2[5],AD1_2[5],AD1_0[5],AD2_0[5],AD2_1[5];
const double gammaDD00 = auxevol_gfs[IDX4S(GAMMADD00GF, i0,i1,i2)];
const double gammaDD01 = auxevol_gfs[IDX4S(GAMMADD01GF, i0,i1,i2)];
const double gammaDD02 = auxevol_gfs[IDX4S(GAMMADD02GF, i0,i1,i2)];
const double gammaDD11 = auxevol_gfs[IDX4S(GAMMADD11GF, i0,i1,i2)];
const double gammaDD12 = auxevol_gfs[IDX4S(GAMMADD12GF, i0,i1,i2)];
const double gammaDD22 = auxevol_gfs[IDX4S(GAMMADD22GF, i0,i1,i2)];
AD0_2[M2] = in_gfs[IDX4S(AD0GF, i0,i1,i2-2)];
AD0_2[M1] = in_gfs[IDX4S(AD0GF, i0,i1,i2-1)];
AD0_1[M2] = in_gfs[IDX4S(AD0GF, i0,i1-2,i2)];
AD0_1[M1] = in_gfs[IDX4S(AD0GF, i0,i1-1,i2)];
AD0_1[P0] = AD0_2[P0] = in_gfs[IDX4S(AD0GF, i0,i1,i2)];
AD0_1[P1] = in_gfs[IDX4S(AD0GF, i0,i1+1,i2)];
AD0_1[P2] = in_gfs[IDX4S(AD0GF, i0,i1+2,i2)];
AD0_2[P1] = in_gfs[IDX4S(AD0GF, i0,i1,i2+1)];
AD0_2[P2] = in_gfs[IDX4S(AD0GF, i0,i1,i2+2)];
AD1_2[M2] = in_gfs[IDX4S(AD1GF, i0,i1,i2-2)];
AD1_2[M1] = in_gfs[IDX4S(AD1GF, i0,i1,i2-1)];
AD1_0[M2] = in_gfs[IDX4S(AD1GF, i0-2,i1,i2)];
AD1_0[M1] = in_gfs[IDX4S(AD1GF, i0-1,i1,i2)];
AD1_2[P0] = AD1_0[P0] = in_gfs[IDX4S(AD1GF, i0,i1,i2)];
AD1_0[P1] = in_gfs[IDX4S(AD1GF, i0+1,i1,i2)];
AD1_0[P2] = in_gfs[IDX4S(AD1GF, i0+2,i1,i2)];
AD1_2[P1] = in_gfs[IDX4S(AD1GF, i0,i1,i2+1)];
AD1_2[P2] = in_gfs[IDX4S(AD1GF, i0,i1,i2+2)];
AD2_1[M2] = in_gfs[IDX4S(AD2GF, i0,i1-2,i2)];
AD2_1[M1] = in_gfs[IDX4S(AD2GF, i0,i1-1,i2)];
AD2_0[M2] = in_gfs[IDX4S(AD2GF, i0-2,i1,i2)];
AD2_0[M1] = in_gfs[IDX4S(AD2GF, i0-1,i1,i2)];
AD2_0[P0] = AD2_1[P0] = in_gfs[IDX4S(AD2GF, i0,i1,i2)];
AD2_0[P1] = in_gfs[IDX4S(AD2GF, i0+1,i1,i2)];
AD2_0[P2] = in_gfs[IDX4S(AD2GF, i0+2,i1,i2)];
AD2_1[P1] = in_gfs[IDX4S(AD2GF, i0,i1+1,i2)];
AD2_1[P2] = in_gfs[IDX4S(AD2GF, i0,i1+2,i2)];
const double invsqrtg = 1.0/sqrt(gammaDD00*gammaDD11*gammaDD22 
                               - gammaDD00*gammaDD12*gammaDD12
                               + 2*gammaDD01*gammaDD02*gammaDD12
                               - gammaDD11*gammaDD02*gammaDD02
                               - gammaDD22*gammaDD01*gammaDD01);

REAL BU0[4],BU1[4],BU2[4];
compute_Bx_pointwise(BU0,invdx2,AD1_2,invdx1,AD2_1);
compute_Bx_pointwise(BU1,invdx0,AD2_0,invdx2,AD0_2);
compute_Bx_pointwise(BU2,invdx1,AD0_1,invdx0,AD1_0);

auxevol_gfs[IDX4S(BU0GF, i0,i1,i2)] = find_accepted_Bx_order(BU0)*invsqrtg;
auxevol_gfs[IDX4S(BU1GF, i0,i1,i2)] = find_accepted_Bx_order(BU1)*invsqrtg;
auxevol_gfs[IDX4S(BU2GF, i0,i1,i2)] = find_accepted_Bx_order(BU2)*invsqrtg;
"""
    
    # Here, we'll use the outCfunction() function to output a function that will compute the magnetic field
    # on the interior. Then, we'll add postloop code to handle the ghostzones.    
    desc="Compute the magnetic field from the vector potential everywhere, including ghostzones"
    name="driver_A_to_B"
    driver_Ccode = outCfunction(
        outfile  = "returnstring", desc=desc, name=name,
        params   = "const paramstruct *restrict params,REAL *restrict in_gfs,REAL *restrict auxevol_gfs",
        body     = order_lowering_body,
        postloop = """
    int imin[3] = { NGHOSTS_A2B, NGHOSTS_A2B, NGHOSTS_A2B };
    int imax[3] = { NGHOSTS+Nxx0, NGHOSTS+Nxx1, NGHOSTS+Nxx2 };
    // Now, we loop over the ghostzones to calculate the magnetic field there. 
    for(int which_gz = 0; which_gz < NGHOSTS_A2B; which_gz++) {
        // After updating each face, adjust imin[] and imax[] 
        //   to reflect the newly-updated face extents.
        compute_A2B_in_ghostzones(params,in_gfs,auxevol_gfs,imin[0]-1,imin[0], imin[1],imax[1], imin[2],imax[2]); imin[0]--;
        compute_A2B_in_ghostzones(params,in_gfs,auxevol_gfs,imax[0],imax[0]+1, imin[1],imax[1], imin[2],imax[2]); imax[0]++;

        compute_A2B_in_ghostzones(params,in_gfs,auxevol_gfs,imin[0],imax[0], imin[1]-1,imin[1], imin[2],imax[2]); imin[1]--;
        compute_A2B_in_ghostzones(params,in_gfs,auxevol_gfs,imin[0],imax[0], imax[1],imax[1]+1, imin[2],imax[2]); imax[1]++;

        compute_A2B_in_ghostzones(params,in_gfs,auxevol_gfs,imin[0],imax[0], imin[1],imax[1], imin[2]-1,imin[2]); imin[2]--;
        compute_A2B_in_ghostzones(params,in_gfs,auxevol_gfs,imin[0],imax[0], imin[1],imax[1], imax[2],imax[2]+1); imax[2]++;
    }
""",
        loopopts="InteriorPoints",
        rel_path_for_Cparams=os.path.join("../")).replace("=NGHOSTS","=NGHOSTS_A2B").replace("NGHOSTS+Nxx0","Nxx_plus_2NGHOSTS0-NGHOSTS_A2B").replace("NGHOSTS+Nxx1","Nxx_plus_2NGHOSTS1-NGHOSTS_A2B").replace("NGHOSTS+Nxx2","Nxx_plus_2NGHOSTS2-NGHOSTS_A2B")

    with open(os.path.join(outdir,"driver_AtoB.h"),"a") as file:
        file.write(driver_Ccode)
