# Step 0: Add NRPy's directory to the path
# https://stackoverflow.com/questions/16780014/import-file-from-parent-directory
import os,sys
nrpy_dir_path = os.path.join("..")
if nrpy_dir_path not in sys.path:
    sys.path.append(nrpy_dir_path)

import cmdline_helper as cmd     # NRPy+: Multi-platform Python command-line interface
global GiRaFFE_NRPy_FCVAL

def GiRaFFE_NRPy_FCVAL(Ccodesdir):
    cmd.mkdir(Ccodesdir)
    # Write out the code to a file.
    with open(os.path.join(Ccodesdir,"interpolate_metric_gfs_to_cell_faces.h"),"w") as file:
        file.write("""// Side note: the following values could be used for cell averaged gfs: 
//     am2=-1.0/12.0, am1=7.0/12.0, a0=7.0/12.0, a1=-1.0/12.0
// However, since the metric gfs store the grid point values instead of the cell average, 
//     the following coefficients should be used: 
//     am2 = -1/16, am1 = 9/16, a0 = 9/16, a1 = -1/16
// This will yield the third-order-accurate face values at m-1/2, 
//      using values specified at {m-2,m-1,m,m+1}
#define AM2 -0.0625
#define AM1  0.5625
#define A0   0.5625
#define A1  -0.0625
#define COMPUTE_FCVAL(METRICm2,METRICm1,METRIC,METRICp1) (AM2*(METRICm2) + AM1*(METRICm1) + A0*(METRIC) + A1*(METRICp1))

const int metric_gfs_list[10] = {GAMMADD00GF,
                                 GAMMADD01GF,
                                 GAMMADD02GF,
                                 GAMMADD11GF,
                                 GAMMADD12GF,
                                 GAMMADD22GF,
                                 BETAU0GF,
                                 BETAU1GF,
                                 BETAU2GF,
                                 ALPHAGF}

const int metric_gfs_face_list[10] = {GAMMA_FACEDD00GF,
                                      GAMMA_FACEDD01GF,
                                      GAMMA_FACEDD02GF,
                                      GAMMA_FACEDD11GF,
                                      GAMMA_FACEDD12GF,
                                      GAMMA_FACEDD22GF,
                                      BETA_FACEU0GF,
                                      BETA_FACEU1GF,
                                      BETA_FACEU2GF,
                                      ALPHA_FACEGF}

const int num_metric_gfs = 10;
// Kronecker delta: all zero-offset
const int kronecker_delta[3][3] = { { 1,0,0 },
                                    { 0,1,0 },
                                    { 0,0,1 } };

void interpolate_metric_gfs_to_cell_faces(const paramstruct *params,REAL *auxevol_gfs,const int flux_dirn) {
#include "set_Cparameters.h"
 
#pragma omp parallel for
    for(int gf = 0;gf < num_metric_gfs;gf++) {
        int in_gf  = metric_gfs_list[gf];
        int out_gf = metric_gfs_face_list[gf];
        for (int i2 = NGHOSTS;i2 < Nxx2+NGHOSTS+1;i2++) {
            for (int i1 = NGHOSTS;i1 < Nxx2+NGHOSTS+1;i1++) {
                for (int i0 = NGHOSTS;i0 < Nxx2+NGHOSTS+1;i0++) {
                    REAL Qm2 = auxevol_gfs[IDX4S(in_gf,i0-2*kronecker_delta[0,flux_dirn],i1-2*kronecker_delta[1,flux_dirn],i2-2*kronecker_delta[2,flux_dirn])];
                    REAL Qm1 = auxevol_gfs[IDX4S(in_gf,i0-kronecker_delta[0,flux_dirn],i1-kronecker_delta[1,flux_dirn],i2-kronecker_delta[2,flux_dirn])];
                    REAL Qp0 = auxevol_gfs[IDX4S(in_gf,i0,i1,i2)];
                    REAL Qp1 = auxevol_gfs[IDX4S(in_gf,i0+kronecker_delta[0,flux_dirn],i1+kronecker_delta[1,flux_dirn],i2+kronecker_delta[2,flux_dirn])];
                    auxevol_gfs[IDX4S(out_gf,i0,i1,i2)] = COMPUTE_FCVAL(Qm2,Qm1,Qp0,Qp1);
                }
            }
        }
    }
}
""")