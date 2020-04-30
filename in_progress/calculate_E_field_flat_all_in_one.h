REAL HLLE_solve(REAL F0B1_r, REAL F0B1_l, REAL U_r, REAL U_l) {
  // Eq. 15 of https://epubs.siam.org/doi/abs/10.1137/1025002?journalCode=siread
  // F_HLLE = (c_min F_R + c_max F_L - c_min c_max (U_R-U_L)) / (c_min + c_max)
  return 0.5*(F0B1_r+F0B1_l-(U_r-U_l));
  // FIXME: Curved space implementation!
}

/*
Calculate the electric flux on both faces in the input direction.
The input count is an integer that is either 0 or 1. If it is 0, this implies
that the components are input in order of a forward permutation. If it is 1, 
then the permutation is backwards and the final results will need to be 
multiplied by -1.0
 */
void calculate_E_field_flat_all_in_one(const paramstruct *params,
                                       const REAL *Vr0,const REAL *Vr1,
                                       const REAL *Vl0,const REAL *Vl1,
                                       const REAL *Br0,const REAL *Br1,
                                       const REAL *Bl0,const REAL *Bl1,
                                       REAL *A2_rhs,const REAL SIGN,const int flux_dirn) {
    // FIXME: include metric functions!
    // This function is written to be generic and compute the contribution for all three AD RHSs.
    // However, for convenience, the notation used in the function itself is for the contribution
    // to AD2, specifically the [F_HLL^x(B^y)]_z term, with reconstructions in the x direction. This
    // corresponds to flux_dirn=0 and count=1 (which corresponds to SIGN=+1.0).
    // Thus, Az(i,j,k) += 0.25 ( [F_HLL^x(B^y)]_z(i+1/2,j,k)+[F_HLL^x(B^y)]_z(i-1/2,j,k)) are solved here.
    // The other terms are computed by cyclically permuting the indices when calling this function.
#include "GiRaFFE_standalone_Ccodes/set_Cparameters.h"

#pragma omp parallel for
    for(int i2=NGHOSTS; i2<NGHOSTS+Nxx2; i2++) {
        for(int i1=NGHOSTS; i1<NGHOSTS+Nxx1; i1++) {
            for(int i0=NGHOSTS; i0<NGHOSTS+Nxx0; i0++) {
                // First, we set the index from which we will read memory. indexp1 is incremented by
                // one point in the direction of reconstruction. These correspond to the faces at at 
                // i-1/2 and i+1/2
                int index   = IDX3S(i0,i1,i2);
                int indexp1 = IDX3S(i0+(flux_dirn==0),i1+(flux_dirn==1),i2+(flux_dirn==2));
                if(flux_dirn==0 && SIGN>0 && i1==Nxx_plus_2NGHOSTS1/2 && i2==Nxx_plus_2NGHOSTS2/2) {
                    printf("index=%d & indexp1=%d\n",index,indexp1);
                }
                
                // Now, we read in memory. We need all components of velocity and magnetic field on both 
                // the left and right sides of the interface at *both* faces.
                const double Valenciav_rU0 = Vr0[index];
                const double Valenciav_rU1 = Vr1[index];
                const double B_rU0 = Br0[index];
                const double B_rU1 = Br1[index];
                const double Valenciav_lU0 = Vl0[index];
                const double Valenciav_lU1 = Vl1[index];
                const double B_lU0 = Bl0[index];
                const double B_lU1 = Bl1[index];
                
                const double Valenciav_rU0_p1 = Vr0[indexp1];
                const double Valenciav_rU1_p1 = Vr1[indexp1];
                const double B_rU0_p1 = Br0[indexp1];
                const double B_rU1_p1 = Br1[indexp1];
                const double Valenciav_lU0_p1 = Vl0[indexp1];
                const double Valenciav_lU1_p1 = Vl1[indexp1];
                const double B_lU0_p1 = Bl0[indexp1];
                const double B_lU1_p1 = Bl1[indexp1];

                // -E_z(x_i,y_j,z_k) &= 0.25 ( [F_HLL^x(B^y)]_z(i+1/2,j,k)+[F_HLL^x(B^y)]_z(i-1/2,j,k)
                //                            -[F_HLL^y(B^x)]_z(i,j+1/2,k)-[F_HLL^y(B^x)]_z(i,j-1/2,k) )
                // We will construct the above sum one half at a time, first with SIGN=+1, which corresponds to flux_dirn = ???, and
                //  takes care of the term:
                //  [...............]

                // ( Note that we will repeat the above with flux_dirn = ???, with SIGN = ????
                //   so that we get the term
                //  [...............]
                // thus completing the above sum. )

                // Calculate the flux vector on each face for each component of the E-field:
                // The F(B) terms are as Eq. 6 in Giacomazzo: https://arxiv.org/pdf/1009.2468.pdf
                // [F^i(B^j)]_k = \sqrt{\gamma} (v^i B^j - v^j B^i)
                // Therefore for blah
                const REAL F0B1_r = (Valenciav_rU0*B_rU1 - Valenciav_rU1*B_rU0);
                const REAL F0B1_l = (Valenciav_lU0*B_lU1 - Valenciav_lU1*B_lU0);

                // Compute the state vector for this flux direction
                const REAL U_r = B_rU1;
                const REAL U_l = B_lU1;

                // Repeat at i+1
                // Calculate the flux vector on each face for each component of the E-field:
                const REAL F0B1_r_p1 = (Valenciav_rU0_p1*B_rU1_p1 - Valenciav_rU1_p1*B_rU0_p1);
                const REAL F0B1_l_p1 = (Valenciav_lU0_p1*B_lU1_p1 - Valenciav_lU1_p1*B_lU0_p1);
                
                // Compute the state vector for this flux direction
                const REAL U_r_p1 = B_rU1_p1;
                const REAL U_l_p1 = B_lU1_p1;

                // Basic HLLE solver: 
                const REAL FHLL_0B1 = HLLE_solve(F0B1_r, F0B1_l, U_r, U_l);

                // Basic HLLE solver, but at the next point: 
                const REAL FHLL_0B1p1 = HLLE_solve(F0B1_r_p1, F0B1_l_p1, U_r_p1, U_l_p1);

                // With the Riemann problem solved, we add the contributions to the RHSs:
                // -E_z(x_i,y_j,z_k) &= 0.25 ( [F_HLL^x(B^y)]_z(i+1/2,j,k)+[F_HLL^x(B^y)]_z(i-1/2,j,k)
                //                            -[F_HLL^y(B^x)]_z(i,j+1/2,k)-[F_HLL^y(B^x)]_z(i,j-1/2,k) )
                // (Eq. 11 in https://arxiv.org/pdf/1009.2468.pdf)
                // This code, as written, solves the first two terms for flux_dirn=0. Calling this function for count=1
                // flips x for y to solve the latter two, switching to SIGN=-1 as well.

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
