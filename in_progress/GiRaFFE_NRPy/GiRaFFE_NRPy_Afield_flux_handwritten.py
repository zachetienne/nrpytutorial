# Step 0: Add NRPy's directory to the path
# https://stackoverflow.com/questions/16780014/import-file-from-parent-directory
import os,sys
nrpy_dir_path = os.path.join("..")
if nrpy_dir_path not in sys.path:
    sys.path.append(nrpy_dir_path)

import cmdline_helper as cmd     # NRPy+: Multi-platform Python command-line interface
Ccodesdir = "GiRaFFE_standalone_Ccodes/RHSs"
cmd.mkdir(os.path.join(Ccodesdir))

from outputC import outputC # NRPy+: Core C code output module
import sympy as sp               # SymPy: The Python computer algebra package upon which NRPy+ depends
import indexedexp as ixp         # NRPy+: Symbolic indexed expression (e.g., tensors, vectors, etc.) support
import GiRaFFE_NRPy.GiRaFFE_NRPy_Characteristic_Speeds as chsp

def GiRaFFE_NRPy_Afield_flux(Ccodesdir):
    cmd.mkdir(Ccodesdir)
    # Write out the code to a file.

    gammaDD = ixp.declarerank2("gammaDD","sym01",DIM=3)
    betaU = ixp.declarerank1("betaU",DIM=3)
    alpha = sp.sympify("alpha")

    for flux_dirn in range(3):
        chsp.find_cmax_cmin(flux_dirn,gammaDD,betaU,alpha)
        Ccode_kernel = outputC([chsp.cmax,chsp.cmin],["cmax","cmin"],"returnstring",params="outCverbose=False,CSE_sorting=none")
        Ccode_kernel = Ccode_kernel.replace("cmax","*cmax").replace("cmin","*cmin")
        Ccode_kernel = Ccode_kernel.replace("betaU0","betaUi").replace("betaU1","betaUi").replace("betaU2","betaUi")

        with open(os.path.join(Ccodesdir,"compute_cmax_cmin_dirn"+str(flux_dirn)+".h"),"w") as file:
            file.write(Ccode_kernel)

    with open(os.path.join(Ccodesdir,"calculate_E_field_flat_all_in_one.h"),"w") as file:
        file.write("""void find_cmax_cmin(const REAL gammaDD00, const REAL gammaDD01, const REAL gammaDD02,
                    const REAL gammaDD11, const REAL gammaDD12, const REAL gammaDD22,
                    const REAL betaUi, const REAL alpha, const int flux_dirn,
                    REAL *cmax, REAL *cmin) {
    switch(flux_dirn) {
        case 0:
#include "compute_cmax_cmin_dirn0.h"
            break;
        case 1:
#include "compute_cmax_cmin_dirn1.h"
            break;
        case 2:
#include "compute_cmax_cmin_dirn2.h"
            break;
        default:
            printf("Invalid parameter flux_dirn!"); *cmax = 1.0/0.0; *cmin = 1.0/0.0;
            break;
    }
}

REAL HLLE_solve(REAL F0B1_r, REAL F0B1_l, REAL U_r, REAL U_l, REAL cmin, REAL cmax) {
  // Eq. 3.15 of https://epubs.siam.org/doi/abs/10.1137/1025002?journalCode=siread
  // F_HLLE = (c_min F_R + c_max F_L - c_min c_max (U_R-U_L)) / (c_min + c_max)
  return (cmin*F0B1_r + cmax*F0B1_l - cmin*cmax*(U_r-U_l)) / (cmin+cmax);
}

/*
Calculate the electric flux on both faces in the input direction.
The input count is an integer that is either 0 or 1. If it is 0, this implies
that the components are input in order of a backwards permutation  and the final
results will need to be multiplied by -1.0. If it is 1, then the permutation is forwards.
 */
void calculate_E_field_flat_all_in_one(const paramstruct *params,
                                       const REAL *Vr0,const REAL *Vr1,
                                       const REAL *Vl0,const REAL *Vl1,
                                       const REAL *Br0,const REAL *Br1,
                                       const REAL *Bl0,const REAL *Bl1,
                                       const REAL *Brflux_dirn,
                                       const REAL *Blflux_dirn,
                                       const REAL *gamma_faceDD00, const REAL *gamma_faceDD01, const REAL *gamma_faceDD02,
                                       const REAL *gamma_faceDD11, const REAL *gamma_faceDD12, const REAL *gamma_faceDD22,
                                       const REAL *beta_faceU0, const REAL *beta_faceU1, const REAL *alpha_face,
                                       REAL *A2_rhs,const REAL SIGN,const int flux_dirn) {
    // This function is written to be generic and compute the contribution for all three AD RHSs.
    // However, for convenience, the notation used in the function itself is for the contribution
    // to AD2, specifically the [F_HLL^x(B^y)]_z term, with reconstructions in the x direction. This
    // corresponds to flux_dirn=0 and count=1 (which corresponds to SIGN=+1.0).
    // Thus, Az(i,j,k) += 0.25 ( [F_HLL^x(B^y)]_z(i+1/2,j,k)+[F_HLL^x(B^y)]_z(i-1/2,j,k)) are solved here.
    // The other terms are computed by cyclically permuting the indices when calling this function.
#include "../set_Cparameters.h"

#pragma omp parallel for
    for(int i2=NGHOSTS; i2<NGHOSTS+Nxx2; i2++) {
        for(int i1=NGHOSTS; i1<NGHOSTS+Nxx1; i1++) {
            for(int i0=NGHOSTS; i0<NGHOSTS+Nxx0; i0++) {
                // First, we set the index from which we will read memory. indexp1 is incremented by
                // one point in the direction of reconstruction. These correspond to the faces at at
                // i-1/2 and i+1/2, respectively.

                // Now, we read in memory. We need the x and y components of velocity and magnetic field on both
                // the left and right sides of the interface at *both* faces.
                // Here, the point (i0,i1,i2) corresponds to the point (i-1/2,j,k)
                const int index           = IDX3S(i0,i1,i2);
                const double alpha        = alpha_face[index];
                const double betaU0       = beta_faceU0[index];
                const double betaU1       = beta_faceU1[index];
                const double v_rU0        = alpha*Vr0[index]-betaU0;
                const double v_rU1        = alpha*Vr1[index]-betaU1;
                const double B_rU0        = Br0[index];
                const double B_rU1        = Br1[index];
                const double B_rflux_dirn = Brflux_dirn[index];
                const double v_lU0        = alpha*Vl0[index]-betaU0;
                const double v_lU1        = alpha*Vl1[index]-betaU1;
                const double B_lU0        = Bl0[index];
                const double B_lU1        = Bl1[index];
                const double B_lflux_dirn = Blflux_dirn[index];
                // We will also need need the square root of the metric determinant here at this point:
                const REAL gxx = gamma_faceDD00[index];
                const REAL gxy = gamma_faceDD01[index];
                const REAL gxz = gamma_faceDD02[index];
                const REAL gyy = gamma_faceDD11[index];
                const REAL gyz = gamma_faceDD12[index];
                const REAL gzz = gamma_faceDD22[index];
                const REAL sqrtgammaDET = sqrt( gxx*gyy*gzz
                                             -  gxx*gyz*gyz
                                             +2*gxy*gxz*gyz
                                             -  gyy*gxz*gxz
                                             -  gzz*gxy*gxy );

                // *******************************
                // REPEAT ABOVE, but at i+1, which corresponds to point (i+1/2,j,k)
                //     Recall that the documentation here assumes flux_dirn==0, but the
                //     algorithm is generalized so that any flux_dirn or velocity/magnetic
                //     field component can be computed via permuting the inputs into this
                //     function.
                const int indexp1            = IDX3S(i0+(flux_dirn==0),i1+(flux_dirn==1),i2+(flux_dirn==2));
                const double alpha_p1        = alpha_face[indexp1];
                const double betaU0_p1       = beta_faceU0[indexp1];
                const double betaU1_p1       = beta_faceU1[indexp1];
                const double v_rU0_p1        = alpha_p1*Vr0[indexp1]-betaU0_p1;
                const double v_rU1_p1        = alpha_p1*Vr1[indexp1]-betaU1_p1;
                const double B_rU0_p1        = Br0[indexp1];
                const double B_rU1_p1        = Br1[indexp1];
                const double B_rflux_dirn_p1 = Brflux_dirn[indexp1];
                const double v_lU0_p1        = alpha_p1*Vl0[indexp1]-betaU0_p1;
                const double v_lU1_p1        = alpha_p1*Vl1[indexp1]-betaU1_p1;
                const double B_lU0_p1        = Bl0[indexp1];
                const double B_lU1_p1        = Bl1[indexp1];
                const double B_lflux_dirn_p1 = Blflux_dirn[indexp1];
                // We will also need need the square root of the metric determinant here at this point:
                const REAL gxx_p1 = gamma_faceDD00[indexp1];
                const REAL gxy_p1 = gamma_faceDD01[indexp1];
                const REAL gxz_p1 = gamma_faceDD02[indexp1];
                const REAL gyy_p1 = gamma_faceDD11[indexp1];
                const REAL gyz_p1 = gamma_faceDD12[indexp1];
                const REAL gzz_p1 = gamma_faceDD22[indexp1];
                const REAL sqrtgammaDET_p1 = sqrt( gxx_p1*gyy_p1*gzz_p1
                                                -  gxx_p1*gyz_p1*gyz_p1
                                                +2*gxy_p1*gxz_p1*gyz_p1
                                                -  gyy_p1*gxz_p1*gxz_p1
                                                -  gzz_p1*gxy_p1*gxy_p1 );

                // *******************************

                // DEBUGGING:
//                 if(flux_dirn==0 && SIGN>0 && i1==Nxx_plus_2NGHOSTS1/2 && i2==Nxx_plus_2NGHOSTS2/2) {
//                     printf("index=%d & indexp1=%d\\n",index,indexp1);
//                 }

                // Since we are computing A_z, the relevant equation here is:
                // -E_z(x_i,y_j,z_k) = 0.25 ( [F_HLL^x(B^y)]_z(i+1/2,j,k)+[F_HLL^x(B^y)]_z(i-1/2,j,k)
                //                           -[F_HLL^y(B^x)]_z(i,j+1/2,k)-[F_HLL^y(B^x)]_z(i,j-1/2,k) )
                // We will construct the above sum one half at a time, first with SIGN=+1, which
                // corresponds to flux_dirn = 0, count=1, and
                //  takes care of the terms:
                //  [F_HLL^x(B^y)]_z(i+1/2,j,k)+[F_HLL^x(B^y)]_z(i-1/2,j,k)

                // ( Note that we will repeat the above with flux_dirn = 1, count = 0, with SIGN=-1
                //   AND with the input components switched (x->y,y->x) so that we get the term
                // -[F_HLL^y(B^x)]_z(i,j+1/2,k)-[F_HLL^y(B^x)]_z(i,j-1/2,k)
                // thus completing the above sum. )

                // Here, [F_HLL^i(B^j)]_k = (v^i B^j - v^j B^i) in general.

                // Calculate the flux vector on each face for each component of the E-field:
                // The F(B) terms are as Eq. 6 in Giacomazzo: https://arxiv.org/pdf/1009.2468.pdf
                // [F^i(B^j)]_k = \sqrt{\gamma} (v^i B^j - v^j B^i)
                // Therefore since we want [F_HLL^x(B^y)]_z,
                // we will code     (v^x           B^y   - v^y           B^x) on both left and right faces.
                const REAL F0B1_r = sqrtgammaDET*(v_rU0*B_rU1 - v_rU1*B_rU0);
                const REAL F0B1_l = sqrtgammaDET*(v_lU0*B_lU1 - v_lU1*B_lU0);

                // Compute the state vector for these terms:
                const REAL U_r = B_rflux_dirn;
                const REAL U_l = B_lflux_dirn;

                REAL cmin,cmax;
                // Basic HLLE solver:
                find_cmax_cmin(gxx,gxy,gxz,
                               gyy,gyz,gzz,
                               betaU0,alpha,flux_dirn,
                               &cmax, &cmin);
                const REAL FHLL_0B1 = HLLE_solve(F0B1_r, F0B1_l, U_r, U_l, cmin, cmax);

                // ************************************
                // ************************************
                // REPEAT ABOVE, but at point i+1
                // Calculate the flux vector on each face for each component of the E-field:
                const REAL F0B1_r_p1 = sqrtgammaDET_p1*(v_rU0_p1*B_rU1_p1 - v_rU1_p1*B_rU0_p1);
                const REAL F0B1_l_p1 = sqrtgammaDET_p1*(v_lU0_p1*B_lU1_p1 - v_lU1_p1*B_lU0_p1);

                // Compute the state vector for this flux direction
                const REAL U_r_p1 = B_rflux_dirn_p1;
                const REAL U_l_p1 = B_lflux_dirn_p1;
                //const REAL U_r_p1 = B_rU1_p1;
                //const REAL U_l_p1 = B_lU1_p1;
                // Basic HLLE solver, but at the next point:
                find_cmax_cmin(gxx_p1,gxy_p1,gxz_p1,
                               gyy_p1,gyz_p1,gzz_p1,
                               betaU0_p1,alpha_p1,flux_dirn,
                               &cmax, &cmin);
                const REAL FHLL_0B1p1 = HLLE_solve(F0B1_r_p1, F0B1_l_p1, U_r_p1, U_l_p1, cmin, cmax);
                // ************************************
                // ************************************


                // With the Riemann problem solved, we add the contributions to the RHSs:
                // -E_z(x_i,y_j,z_k) &= 0.25 ( [F_HLL^x(B^y)]_z(i+1/2,j,k)+[F_HLL^x(B^y)]_z(i-1/2,j,k)
                //                            -[F_HLL^y(B^x)]_z(i,j+1/2,k)-[F_HLL^y(B^x)]_z(i,j-1/2,k) )
                // (Eq. 11 in https://arxiv.org/pdf/1009.2468.pdf)
                // This code, as written, solves the first two terms for flux_dirn=0. Calling this function for count=0
                // and flux_dirn=1 flips x for y to solve the latter two, switching to SIGN=-1 as well.

                // Here, we finally add together the output of the HLLE solver at i-1/2 and i+1/2
                // We also multiply by the SIGN dictated by the order of the input vectors and divide by 4.
                A2_rhs[index] += SIGN*0.25*(FHLL_0B1 + FHLL_0B1p1);
                // flux dirn = 0 ===================>   i-1/2       i+1/2
                //               Eq 11 in Giacomazzo:
                //               -FxBy(avg over i-1/2 and i+1/2) + FyBx(avg over j-1/2 and j+1/2)
                //               Eq 6 in Giacomazzo:
                //               FxBy = vxBy - vyBx
                //             ->
                //               FHLL_0B1 = vyBx - vxBy

            } // END LOOP: for(int i0=NGHOSTS; i0<NGHOSTS+Nxx0; i0++)
        } // END LOOP: for(int i1=NGHOSTS; i1<NGHOSTS+Nxx1; i1++)
    } // END LOOP: for(int i2=NGHOSTS; i2<NGHOSTS+Nxx2; i2++)
}
""")
