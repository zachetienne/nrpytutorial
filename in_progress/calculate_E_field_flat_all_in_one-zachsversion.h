/*
int k_delta[3][3] = {{1,0,0},
                     {0,1,0},
                     {0,0,1}};
*/
/*
Calculate the electric flux on both faces in the input direction.
 */

/*
  TO CALL THIS FROM MAIN DRIVER:

        for(int count = 1;count <= 2;count++) {
          int Ai = (flux_dirn+count)%3; // flux_dirn=1, count=1; Ai = 2: A2_rhs += - 0.25*[+FyBx(jp)+FyBx(jm)] : SIGN=-1
          //                               flux_dirn=1, count=2; Ai = 0
          //                               flux_dirn=0, count=1; Ai = 1
          //                               flux_dirn=0, count=2; Ai = 2: A2_rhs += - 0.25*[-FxBy(ip)-FxBy(im)] : SIGN=+1
          //                               flux_dirn=2, count=1; Ai = 0
          //                               flux_dirn=2, count=2; Ai = 1

          REAL SIGN=1.0;
          if(Ai == 0 && flux_dirn == 2) SIGN=-1.0;
          if(Ai == 1 && flux_dirn == 0) SIGN=-1.0;
          if(Ai == 2 && flux_dirn == 1) SIGN=-1.0;

          //printf("hey flux_dirn=%d ; Ai=%d ; SIGN=%e\n",flux_dirn,Ai,SIGN);
          if(SIGN==1.0) {
            calculate_E_field_flat_all_in_one_zachsversion(params,
                                                           &auxevol_gfs[IDX4ptS(VALENCIAV_RU0GF+(Ai+1)%3, 0)],&auxevol_gfs[IDX4ptS(VALENCIAV_RU0GF+(Ai+2)%3, 0)],
                                                           &auxevol_gfs[IDX4ptS(VALENCIAV_LU0GF+(Ai+1)%3, 0)],&auxevol_gfs[IDX4ptS(VALENCIAV_LU0GF+(Ai+2)%3, 0)],
                                                           &auxevol_gfs[IDX4ptS(B_RU0GF        +(Ai+1)%3, 0)],&auxevol_gfs[IDX4ptS(B_RU0GF        +(Ai+2)%3,0)],&auxevol_gfs[IDX4ptS(B_RU0GF+Ai, 0)],
                                                           &auxevol_gfs[IDX4ptS(B_LU0GF        +(Ai+1)%3, 0)],&auxevol_gfs[IDX4ptS(B_LU0GF        +(Ai+2)%3,0)],&auxevol_gfs[IDX4ptS(B_LU0GF+Ai, 0)],
                                                           &rhs_gfs[IDX4ptS(AD0GF+Ai,0)], SIGN, flux_dirn);
          } else {
            calculate_E_field_flat_all_in_one_zachsversion(params,
                                                           &auxevol_gfs[IDX4ptS(VALENCIAV_RU0GF+(Ai+2)%3, 0)],&auxevol_gfs[IDX4ptS(VALENCIAV_RU0GF+(Ai+1)%3, 0)],
                                                           &auxevol_gfs[IDX4ptS(VALENCIAV_LU0GF+(Ai+2)%3, 0)],&auxevol_gfs[IDX4ptS(VALENCIAV_LU0GF+(Ai+1)%3, 0)],
                                                           &auxevol_gfs[IDX4ptS(B_RU0GF        +(Ai+2)%3, 0)],&auxevol_gfs[IDX4ptS(B_RU0GF        +(Ai+1)%3, 0)],&auxevol_gfs[IDX4ptS(B_RU0GF+Ai, 0)],
                                                           &auxevol_gfs[IDX4ptS(B_LU0GF        +(Ai+2)%3, 0)],&auxevol_gfs[IDX4ptS(B_LU0GF        +(Ai+1)%3, 0)],&auxevol_gfs[IDX4ptS(B_LU0GF+Ai, 0)],
                                                           &rhs_gfs[IDX4ptS(AD0GF+Ai,0)], SIGN, flux_dirn);
          }
        }
 */

void calculate_E_field_flat_all_in_one_zachsversion(const paramstruct *params,
                                                    const REAL *Vr0,const REAL *Vr1,
                                                    const REAL *Vl0,const REAL *Vl1,

                                                    const REAL *Br0,const REAL *Br1,const REAL *Br2,
                                                    const REAL *Bl0,const REAL *Bl1,const REAL *Bl2,

                                                    REAL *Az_rhs,const REAL SIGN,const int flux_dirn) {
#include "GiRaFFE_standalone_Ccodes/set_Cparameters.h"
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
                const double B_rU0         = Br0[index];
                const double B_rU1         = Br1[index];
                const double B_rU2         = Br2[index];

                const double Valenciav_lU0 = Vl0[index];
                const double Valenciav_lU1 = Vl1[index];
                const double B_lU0         = Bl0[index];
                const double B_lU1         = Bl1[index];
                const double B_lU2         = Bl2[index];

                // Calculate the flux vector on each face for each component of the E-field:
                const REAL F0B1_r = (Valenciav_rU0*B_rU1 - Valenciav_rU1*B_rU0);
                const REAL F0B1_l = (Valenciav_lU0*B_lU1 - Valenciav_lU1*B_lU0);

                // Compute the state vector for this flux direction
                const REAL U_r = B_rU2;
                const REAL U_l = B_lU2;

                // Basic HLLE solver:
                const REAL FHLL_0B1 = 0.5*(F0B1_r + F0B1_l - (U_r-U_l));

                // Repeat at i+1
                // Now, we read in memory. We need all components of velocity and magnetic field on both
                // the left and right sides of the interface at *both* faces.
                const double Valenciav_rU0_p1 = Vr0[indexp1];
                const double Valenciav_rU1_p1 = Vr1[indexp1];
                const double B_rU0_p1         = Br0[indexp1];
                const double B_rU1_p1         = Br1[indexp1];
                //const double B_rU2_p1         = Br2[indexp1];

                const double Valenciav_lU0_p1 = Vl0[indexp1];
                const double Valenciav_lU1_p1 = Vl1[indexp1];
                const double B_lU0_p1         = Bl0[indexp1];
                const double B_lU1_p1         = Bl1[indexp1];
                //const double B_lU2_p1         = Bl2[indexp1];

                // Calculate the flux vector on each face for each component of the E-field:
                const REAL F0B1_r_p1 = (Valenciav_rU0_p1*B_rU1_p1 - Valenciav_rU1_p1*B_rU0_p1);
                const REAL F0B1_l_p1 = (Valenciav_lU0_p1*B_lU1_p1 - Valenciav_lU1_p1*B_lU0_p1);

                // WRONG (fixme):  Compute the state vector for this flux direction
                const REAL U_r_p1 = B_rU2_p1;
                printf("DO NOT USE calculate_E_field_flat_all_in_one-zachsversion.h!\n");
                exit(1);
                const REAL U_l_p1 = B_lU2_p1;

                // Basic HLLE solver
                const REAL FHLL_0B1_p1 = 0.5*(F0B1_r_p1 + F0B1_l_p1 - (U_r_p1-U_l_p1));

                Az_rhs[index] += SIGN*0.25*(FHLL_0B1 + FHLL_0B1_p1);
            } // END LOOP: for(int i0=NGHOSTS; i0<NGHOSTS+Nxx0; i0++)
        } // END LOOP: for(int i1=NGHOSTS; i1<NGHOSTS+Nxx1; i1++)
    } // END LOOP: for(int i2=NGHOSTS; i2<NGHOSTS+Nxx2; i2++)
}
