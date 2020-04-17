REAL HLLE_solve(F0B1_r, F0B1_l, U_r, U_l) {
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
				       REAL *Az_rhs,const int count,const int flux_dirn) {
  // FIXME: include metric functions!
#include "GiRaFFE_standalone_Ccodes/set_Cparameters.h"

  REAL SIGN = 1.0-2.0*((REAl) count); // 1.0 if count=0, -1.0 if count=1
  
#pragma omp parallel for
    for(int i2=NGHOSTS; i2<NGHOSTS+Nxx2; i2++) {
        for(int i1=NGHOSTS; i1<NGHOSTS+Nxx1; i1++) {
            for(int i0=NGHOSTS; i0<NGHOSTS+Nxx0; i0++) {
                // First, we set the index from which we will read memory. indexp1 is incremented by
                // one point in the direction of reconstruction. These correspond to the faces at at 
                // i-1/2 and i+1/2
                int index   = IDX3S(i0,i1,i2);
                int indexp1 = IDX3S(i0+k_delta[flux_dirn][0],i1+k_delta[flux_dirn][1],i2+k_delta[flux_dirn][2]);
                
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
                
                const double Valenciav_rU0_p1 = Vr0[index];
                const double Valenciav_rU1_p1 = Vr1[indexp1];
                const double B_rU0_p1 = Br0[indexp1];
                const double B_rU1_p1 = Br1[indexp1];
                const double Valenciav_lU0_p1 = Vl0[indexp1];
                const double Valenciav_lU1_p1 = Vl1[indexp1];
                const double B_lU0_p1 = Bl0[indexp1];
                const double B_lU1_p1 = Bl1[indexp1];

                // Calculate the flux vector on each face for each component of the E-field:
                const REAL F0B1_r = (Valenciav_rU0*B_rU1 - Valenciav_rU1*B_rU0);
                const REAL F0B1_l = (Valenciav_lU0*B_lU1 - Valenciav_lU1*B_lU0);

                // Compute the state vector for this flux direction
                const REAL U_r = B_rU0;
                const REAL U_l = B_lU0;

                // Repeat at i+1
                // Calculate the flux vector on each face for each component of the E-field:
                const REAL F0B1_r_p1 = (Valenciav_rU0_p1*B_rU1_p1 - Valenciav_rU1_p1*B_rU0_p1);
                const REAL F0B1_l_p1 = (Valenciav_lU0_p1*B_lU1_p1 - Valenciav_lU1_p1*B_lU0_p1);
                
                // Compute the state vector for this flux direction
                const REAL U_r_p1 = B_rU0_p1;
                const REAL U_l_p1 = B_lU0_p1;

                // Basic HLLE solver: 
                const REAL FHLL_0B1 = HLLE_solve(F0B1_r, F0B1_l, U_r, U_l);

                // Basic HLLE solver, but at the next point: 
                const REAL FHLL_0B1p1 = HLLE_solve(F0B1_r_p1, F0B1_l_p1, U_r_p1, U_l_p1);

                Az_rhs[index]+= SIGN*0.25*(FHLL_0B1 + FHLL_0B1p1)*(flux_dirn!=2);
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
