int k_delta[3][3] = {{1,0,0},
                     {0,1,0},
                     {0,0,1}};

/*
Calculate the electric flux on both faces in the input direction.
 */
void calculate_E_field_flat_all_in_one(const paramstruct *params,const REAL *auxevol_gfs,REAL *rhs_gfs,const int flux_dirn) {
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
                const double Valenciav_rU0 = auxevol_gfs[IDX4ptS(VALENCIAV_RU0GF, index)];
                const double Valenciav_rU1 = auxevol_gfs[IDX4ptS(VALENCIAV_RU1GF, index)];
                const double Valenciav_rU2 = auxevol_gfs[IDX4ptS(VALENCIAV_RU2GF, index)];
                const double B_rU0 = auxevol_gfs[IDX4ptS(B_RU0GF, index)];
                const double B_rU1 = auxevol_gfs[IDX4ptS(B_RU1GF, index)];
                const double B_rU2 = auxevol_gfs[IDX4ptS(B_RU2GF, index)];
                const double Valenciav_lU0 = auxevol_gfs[IDX4ptS(VALENCIAV_LU0GF, index)];
                const double Valenciav_lU1 = auxevol_gfs[IDX4ptS(VALENCIAV_LU1GF, index)];
                const double Valenciav_lU2 = auxevol_gfs[IDX4ptS(VALENCIAV_LU2GF, index)];
                const double B_lU0 = auxevol_gfs[IDX4ptS(B_LU0GF, index)];
                const double B_lU1 = auxevol_gfs[IDX4ptS(B_LU1GF, index)];
                const double B_lU2 = auxevol_gfs[IDX4ptS(B_LU2GF, index)];
                
                const double Valenciav_rU0_p1 = auxevol_gfs[IDX4ptS(VALENCIAV_RU0GF, indexp1)];
                const double Valenciav_rU1_p1 = auxevol_gfs[IDX4ptS(VALENCIAV_RU1GF, indexp1)];
                const double Valenciav_rU2_p1 = auxevol_gfs[IDX4ptS(VALENCIAV_RU2GF, indexp1)];
                const double B_rU0_p1 = auxevol_gfs[IDX4ptS(B_RU0GF, indexp1)];
                const double B_rU1_p1 = auxevol_gfs[IDX4ptS(B_RU1GF, indexp1)];
                const double B_rU2_p1 = auxevol_gfs[IDX4ptS(B_RU2GF, indexp1)];
                const double Valenciav_lU0_p1 = auxevol_gfs[IDX4ptS(VALENCIAV_LU0GF, indexp1)];
                const double Valenciav_lU1_p1 = auxevol_gfs[IDX4ptS(VALENCIAV_LU1GF, indexp1)];
                const double Valenciav_lU2_p1 = auxevol_gfs[IDX4ptS(VALENCIAV_LU2GF, indexp1)];
                const double B_lU0_p1 = auxevol_gfs[IDX4ptS(B_LU0GF, indexp1)];
                const double B_lU1_p1 = auxevol_gfs[IDX4ptS(B_LU1GF, indexp1)];
                const double B_lU2_p1 = auxevol_gfs[IDX4ptS(B_LU2GF, indexp1)];
                
                // Calculate the flux vector on each face for each component of the E-field:
                const REAL F1B2_r = (Valenciav_rU1*B_rU2 - Valenciav_rU2*B_rU1);
                const REAL F1B2_l = (Valenciav_lU1*B_lU2 - Valenciav_lU2*B_lU1);
                
                const REAL F2B0_r = (Valenciav_rU2*B_rU0 - Valenciav_rU0*B_rU2);
                const REAL F2B0_l = (Valenciav_lU2*B_lU0 - Valenciav_lU0*B_lU2);
                
                const REAL F0B1_r = (Valenciav_rU0*B_rU1 - Valenciav_rU1*B_rU0);
                const REAL F0B1_l = (Valenciav_lU0*B_lU1 - Valenciav_lU1*B_lU0);
                
                // Compute the state vector for this flux direction
                const REAL U_r = B_rU0*k_delta[flux_dirn][0] + B_rU1*k_delta[flux_dirn][1] + B_rU2*k_delta[flux_dirn][2];
                const REAL U_l = B_lU0*k_delta[flux_dirn][0] + B_lU1*k_delta[flux_dirn][1] + B_lU2*k_delta[flux_dirn][2];
                
                // Repeat at i+1
                // Calculate the flux vector on each face for each component of the E-field:
                const REAL F1B2_r_p1 = (Valenciav_rU1_p1*B_rU2_p1 - Valenciav_rU2_p1*B_rU1_p1);
                const REAL F1B2_l_p1 = (Valenciav_lU1_p1*B_lU2_p1 - Valenciav_lU2_p1*B_lU1_p1);
                
                const REAL F2B0_r_p1 = (Valenciav_rU2_p1*B_rU0_p1 - Valenciav_rU0_p1*B_rU2_p1);
                const REAL F2B0_l_p1 = (Valenciav_lU2_p1*B_lU0_p1 - Valenciav_lU0_p1*B_lU2_p1);
                
                const REAL F0B1_r_p1 = (Valenciav_rU0_p1*B_rU1_p1 - Valenciav_rU1_p1*B_rU0_p1);
                const REAL F0B1_l_p1 = (Valenciav_lU0_p1*B_lU1_p1 - Valenciav_lU1_p1*B_lU0_p1);
                
                // Compute the state vector for this flux direction
                const REAL U_r_p1 = B_rU0_p1*k_delta[flux_dirn][0] + B_rU1_p1*k_delta[flux_dirn][1] + B_rU2_p1*k_delta[flux_dirn][2];
                const REAL U_l_p1 = B_lU0_p1*k_delta[flux_dirn][0] + B_lU1_p1*k_delta[flux_dirn][1] + B_lU2_p1*k_delta[flux_dirn][2];

                // Basic HLLE solver: 
                const REAL FHLL_1B2 = 0.5*(F1B2_r + F1B2_l - (U_r-U_l));
                const REAL FHLL_2B0 = 0.5*(F2B0_r + F2B0_l - (U_r-U_l));
                const REAL FHLL_0B1 = 0.5*(F0B1_r + F0B1_l - (U_r-U_l));
                
                // Basic HLLE solver, but at the next point: 
                const REAL FHLL_1B2p1 = 0.5*(F1B2_r_p1 + F1B2_l_p1 - (U_r_p1-U_l_p1));
                const REAL FHLL_2B0p1 = 0.5*(F2B0_r_p1 + F2B0_l_p1 - (U_r_p1-U_l_p1));
                const REAL FHLL_0B1p1 = 0.5*(F0B1_r_p1 + F0B1_l_p1 - (U_r_p1-U_l_p1));
                
                
                rhs_gfs[IDX4ptS(AD0GF,index)] += 0.25*(FHLL_1B2 + FHLL_1B2p1)*(flux_dirn!=0); // Set to zero for the component in flux_dirn. Is it more efficient to do this sooner? An array-based implementation might be better, too.
                rhs_gfs[IDX4ptS(AD1GF,index)] += 0.25*(FHLL_2B0 + FHLL_2B0p1)*(flux_dirn!=1); 
                rhs_gfs[IDX4ptS(AD2GF,index)] += 0.25*(FHLL_0B1 + FHLL_0B1p1)*(flux_dirn!=2);
                
            } // END LOOP: for(int i0=NGHOSTS; i0<NGHOSTS+Nxx0; i0++)
        } // END LOOP: for(int i1=NGHOSTS; i1<NGHOSTS+Nxx1; i1++)
    } // END LOOP: for(int i2=NGHOSTS; i2<NGHOSTS+Nxx2; i2++)
}
