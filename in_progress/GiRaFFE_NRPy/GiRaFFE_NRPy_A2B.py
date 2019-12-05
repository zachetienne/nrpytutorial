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
def GiRaFFE_NRPy_A2B(outdir):
    cmd.mkdir(outdir)
    # Set spatial dimension (must be 3 for BSSN)
    DIM = 3
    par.set_parval_from_str("grid::DIM",DIM)
    # Register the gridfunction gammadet. This determinant will be calculated separately
    gammadet = gri.register_gridfunctions("AUXEVOL","gammadet")

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
                LeviCivitaUUU[i][j][k] = LCijk / sp.sqrt(gammadet)
    # Step 1.a: A useful function
    AD = ixp.register_gridfunctions_for_single_rank1("EVOL","AD")
    BU = ixp.register_gridfunctions_for_single_rank1("AUXEVOL","BU")
    AD_dD = ixp.declarerank2("AD_dD","nosym")
    BU = ixp.zerorank1() # BU is already registered as a gridfunction, but we need to zero its values and declare it in this scope.

    for i in range(DIM):
        for j in range(DIM):
            for k in range(DIM):
                BU[i] += LeviCivitaUUU[i][j][k] * AD_dD[k][j]

    # Write the basic #includes and #defines the C code needs
#     with open(os.path.join(outdir,"driver_AtoB.h"),"w") as file:
#         file.write("""#include <math.h>
# #include <stdio.h>
# #include <stdlib.h>

# #ifndef REAL
# #define REAL double
# #define IDX4S(g,i,j,k) \\
# ( (i) + Nxx_plus_2NGHOSTS0 * ( (j) + Nxx_plus_2NGHOSTS1 * ( (k) + Nxx_plus_2NGHOSTS2 * (g) ) ) )
# #endif""")
        
    # Write the code to compute derivatives with shifted stencils as needed.
    with open(os.path.join(outdir,"driver_AtoB.h"),"a") as file:
        file.write("""
inline void compute_A2B_in_ghostzones(const paramstruct *restrict params,REAL *restrict in_gfs,REAL *restrict out_gfs,i0min,i0max, i1min,i1max, i2min,i2max) {
#include "set_Cparameters.h"
    for(int i2=i2min;i2<i2max;i2++) for(int i1=i1min;i1<i1max;i1++) for(int i0=i0min;i0<i0max;i0++) {
        REAL dx_Ay,dx_Az,dy_Ax,dy_Az,dz_Ax,dz_Ay;
        // Check to see if we're on the +x or -x face. If so, use a downwinded- or upwinded-stencil, respectively.
        // Otherwise, use a centered stencil.
        switch(i0) {
            case 0:
                dx_Ay = invdx0*(-3.0/2.0*in_gfs[IDX4S(AD1GF, i0,i1,i2)] + 2*in_gfs[IDX4S(AD1GF, i0+1,i1,i2)] - 1.0/2.0*in_gfs[IDX4S(AD1GF, i0+2,i1,i2)]);
                dx_Az = invdx0*(-3.0/2.0*in_gfs[IDX4S(AD2GF, i0,i1,i2)] + 2*in_gfs[IDX4S(AD2GF, i0+1,i1,i2)] - 1.0/2.0*in_gfs[IDX4S(AD2GF, i0+2,i1,i2)]);
                break;
            case Nxx_plus_2NGHOSTS0-1:
                dx_Ay = invdx0*((3.0/2.0)*in_gfs[IDX4S(AD1GF, i0,i1,i2)] - 2*in_gfs[IDX4S(AD1GF, i0-1,i1,i2)] + (1.0/2.0)*in_gfs[IDX4S(AD1GF, i0-2,i1,i2)]);
                dx_Az = invdx0*((3.0/2.0)*in_gfs[IDX4S(AD2GF, i0,i1,i2)] - 2*in_gfs[IDX4S(AD2GF, i0-1,i1,i2)] + (1.0/2.0)*in_gfs[IDX4S(AD2GF, i0-2,i1,i2)]);
                break;
            default:
                dx_Ay = invdx0*(-1.0/2.0*in_gfs[IDX4S(AD1GF, i0-1,i1,i2)] + (1.0/2.0)*in_gfs[IDX4S(AD1GF, i0+1,i1,i2)]);
                dx_Az = invdx0*(-1.0/2.0*in_gfs[IDX4S(AD2GF, i0-1,i1,i2)] + (1.0/2.0)*in_gfs[IDX4S(AD2GF, i0+1,i1,i2)]);
        }
        // As above, but in the y direction.
        switch(i1) {
            case 0:
                dy_Ax = invdx1*(-3.0/2.0*in_gfs[IDX4S(AD0GF, i0,i1,i2)] + 2*in_gfs[IDX4S(AD0GF, i0,i1+1,i2)] - 1.0/2.0*in_gfs[IDX4S(AD0GF, i0,i1+2,i2)]);
                dy_Az = invdx1*(-3.0/2.0*in_gfs[IDX4S(AD2GF, i0,i1,i2)] + 2*in_gfs[IDX4S(AD2GF, i0,i1+1,i2)] - 1.0/2.0*in_gfs[IDX4S(AD2GF, i0,i1+2,i2)]);
                break;
            case Nxx_plus_2NGHOSTS1-1:
                dy_Ax = invdx1*((3.0/2.0)*in_gfs[IDX4S(AD0GF, i0,i1,i2)] - 2*in_gfs[IDX4S(AD0GF, i0,i1-1,i2)] + (1.0/2.0)*in_gfs[IDX4S(AD0GF, i0,i1-2,i2)]);
                dy_Az = invdx1*((3.0/2.0)*in_gfs[IDX4S(AD2GF, i0,i1,i2)] - 2*in_gfs[IDX4S(AD2GF, i0,i1-1,i2)] + (1.0/2.0)*in_gfs[IDX4S(AD2GF, i0,i1-2,i2)]);
                break;
            default:
                dy_Ax = invdx1*(-1.0/2.0*in_gfs[IDX4S(AD0GF, i0,i1-1,i2)] + (1.0/2.0)*in_gfs[IDX4S(AD1GF, i0,i1+1,i2)]);
                dy_Az = invdx1*(-1.0/2.0*in_gfs[IDX4S(AD2GF, i0,i1-1,i2)] + (1.0/2.0)*in_gfs[IDX4S(AD2GF, i0,i1+1,i2)]);
        }
        // As above, but in the z direction.
        switch(i2) {
            case 0:
                dz_Ax = invdx2*(-3.0/2.0*in_gfs[IDX4S(AD0GF, i0,i1,i2)] + 2*in_gfs[IDX4S(AD0GF, i0,i1,i2+1)] - 1.0/2.0*in_gfs[IDX4S(AD0GF, i0,i1,i2+2)]);
                dz_Ay = invdx2*(-3.0/2.0*in_gfs[IDX4S(AD1GF, i0,i1,i2)] + 2*in_gfs[IDX4S(AD1GF, i0,i1,i2+1)] - 1.0/2.0*in_gfs[IDX4S(AD1GF, i0,i1,i2+2)]);
                break;
            case Nxx_plus_2NGHOSTS2-1:
                dz_Ax = invdx2*((3.0/2.0)*in_gfs[IDX4S(AD0GF, i0,i1,i2)] - 2*in_gfs[IDX4S(AD0GF, i0,i1,i2-1)] + (1.0/2.0)*in_gfs[IDX4S(AD0GF, i0,i1,i2-2)]);
                dz_Ay = invdx2*((3.0/2.0)*in_gfs[IDX4S(AD1GF, i0,i1,i2)] - 2*in_gfs[IDX4S(AD1GF, i0,i1,i2-1)] + (1.0/2.0)*in_gfs[IDX4S(AD1GF, i0,i1,i2-2)]);
                break;
            default:
                dz_Ax = invdx2*(-1.0/2.0*in_gfs[IDX4S(AD0GF, i0,i1,i2-1)] + (1.0/2.0)*in_gfs[IDX4S(AD0GF, i0,i1,i2+1)]);
                dz_Ay = invdx2*(-1.0/2.0*in_gfs[IDX4S(AD1GF, i0,i1,i2-1)] + (1.0/2.0)*in_gfs[IDX4S(AD1GF, i0,i1,i2+1)]);
        }
        // Compute the magnetic field in the normal way, using the previously calculated derivatives.
        const REAL gammadet = out_gfs[IDX4S(GAMMADETGF, i0,i1,i2)];
        out_gfs[IDX4S(BU0GF, i0,i1,i2)] = (dy_Az-dz_Ay)/gammadet;
        out_gfs[IDX4S(BU1GF, i0,i1,i2)] = (dz_Ax-dx_Az)/gammadet;
        out_gfs[IDX4S(BU2GF, i0,i1,i2)] = (dx_Ay-dy_Ax)/gammadet;
    }
}
""")
    
    # Here, we'll use the outCfunction() function to output a function that will compute the magnetic field
    # on the interior. Then, we'll add postloop code to handle the ghostzones.    
    desc="Compute the magnetic field from the vector potential everywhere, including ghostzones"
    name="driver_A_to_B"
    driver_Ccode = outCfunction(
        outfile  = "returnstring", desc=desc, name=name,
        params   = "const paramstruct *restrict params,REAL *restrict in_gfs,REAL *restrict out_gfs",
        preloop  = """
    int imin[3] = { NGHOSTS, NGHOSTS, NGHOSTS };
    int imax[3] = { Nxx_plus_2NGHOSTS0-NGHOSTS, Nxx_plus_2NGHOSTS1-NGHOSTS, Nxx_plus_2NGHOSTS2-NGHOSTS };
""",
        body     = fin.FD_outputC("returnstring",[lhrh(lhs=gri.gfaccess("out_gfs","BU0"),rhs=BU[0]),\
                                                  lhrh(lhs=gri.gfaccess("out_gfs","BU1"),rhs=BU[1]),\
                                                  lhrh(lhs=gri.gfaccess("out_gfs","BU2"),rhs=BU[2])],
                                  params="outCverbose=False").replace("IDX4","IDX4S"),
        postloop = """
        // Now, we loop over the ghostzones to calculate the magnetic field there. 
        for(int which_gz = 0; which_gz < NGHOSTS; which_gz++) {
            // After updating each face, adjust imin[] and imax[] 
            //   to reflect the newly-updated face extents.
            compute_A2B_in_ghostzones(params,in_gfs,out_gfs,imin[0]-1,imin[0], imin[1],imax[1], imin[2],imax[2]); imin[0]--;
            compute_A2B_in_ghostzones(params,in_gfs,out_gfs,imax[0],imax[0]+1, imin[1],imax[1], imin[2],imax[2]); imax[0]++;

            compute_A2B_in_ghostzones(params,in_gfs,out_gfs,imin[0],imax[0], imin[1]-1,imin[1], imin[2],imax[2]); imin[1]--;
            compute_A2B_in_ghostzones(params,in_gfs,out_gfs,imin[0],imax[0], imax[1],imax[1]+1, imin[2],imax[2]); imax[1]++;

            compute_A2B_in_ghostzones(params,in_gfs,out_gfs,imin[0],imax[0], imin[1],imax[1], imin[2]-1,imin[2]); imin[2]--;
            compute_A2B_in_ghostzones(params,in_gfs,out_gfs,imin[0],imax[0], imin[1],imax[1], imax[2],imax[2]+1); imax[2]++;
        }
""",
        loopopts="InteriorPoints")

    with open(os.path.join(outdir,"driver_AtoB.h"),"a") as file:
        file.write(driver_Ccode)
