/* We evolve forward in time a set of functions called the 
 * "conservative variables" (magnetic field and Poynting vector), 
 * and any time the conserv's are updated, we must recover the 
 * primitive variables (velocities), before reconstructing & evaluating 
 * the RHSs of the MHD equations again. 
 *
 * This file contains the routine for this algebraic calculation. 
 * The velocity is calculated with formula (85), arXiv:1310.3274v2
 * $v^i = 4 \pi \alpha \gamma^{ij} {\tilde S}_j \gamma{-1/2} B^{-2} - \beta^i$ 
 * The force-free condition: $B^2>E^2$ is checked before computing the velocity.
 * and after imposing the constraint ${\tilde B}^i {\tilde S}_i = 0$
 
 * The procedure is as described in arXiv:1310.3274v2: 
 * 1. ${\tilde S}_i ->{\tilde S}_i - ({\tilde S}_j {\tilde B}^j) {\tilde B}_i/{\tilde B}^2$
 * 2. $f = \sqrt{(1-\gamma_{max}^{-2}){\tilde B}^4/(16 \pi^2 \gamma {\tilde S}^2)}$ 
 * 3. ${\tilde S}_i -> {\tilde S}_i min(1,f)
 * 4. $v^i = 4 \pi \alpha \gamma^{ij} {\tilde S}_j \gamma{-1/2} B^{-2} - \beta^i$
 * 5. ${\tilde n}_i v^i = 0$
 *
 * All equations are from: http://arxiv.org/pdf/1310.3274.pdf (v2)
 * */

//#include <iostream>
//#include <iomanip>
//#include <fstream>
#include <sys/time.h>
//#include <cmath>
//#include <ctime>
//#include <cstdlib>
//#include "Symmetry.h"

#ifndef M_PI
#define M_PI 3.141592653589793238463
#endif
#define MIN(a,b) ( ((a) < (b)) ? (a) : (b) )
void GiRaFFE_NRPy_compute_conservatives(const REAL gxxL,const REAL gxyL,const REAL gxzL,const REAL gyyL,const REAL gyzL,const REAL gzzL,
                                      const REAL BxL, const REAL ByL, const REAL BzL, const REAL vxL, const REAL vyL, const REAL vzL,
                                      //const REAL betaxL, const REAL betayL, const REAL betazL, const REAL alpL,
                                      const REAL sqrtg,REAL *StildeD0L, REAL *StildeD1L, REAL *StildeD2L);
#include "compute_conservatives_FFE.C"

#define REAL double

#define IDX4(g,i,j,k) \
( (i) + Nxx_plus_2NGHOSTS[0] * ( (j) + Nxx_plus_2NGHOSTS[1] * ( (k) + Nxx_plus_2NGHOSTS[2] * (g) ) ) )
#define IDX3(i,j,k) ( (i) + Nxx_plus_2NGHOSTS[0] * ( (j) + Nxx_plus_2NGHOSTS[1] * (k) ) )
// Assuming idx = IDX3(i,j,k). Much faster if idx can be reused over and over:
#define IDX4pt(g,idx)   ( (idx) + (Nxx_plus_2NGHOSTS[0]*Nxx_plus_2NGHOSTS[1]*Nxx_plus_2NGHOSTS[2]) * (g) )

#include "metric_quantities.h"
void GiRaFFE_NRPy_conserv_to_prims_FFE(const int Nxx[3],const int Nxx_plus_2NGHOSTS[3],const REAL dxx[3],
                                     REAL *xx[3], REAL *in_gfs, REAL *aux_gfs) {
  //printf("Starting conservative-to-primitive solver...\n");

  /*// We use proper C++ here, for file I/O later.
  using namespace std;*/

  const int imin=NGHOSTS,jmin=NGHOSTS,kmin=NGHOSTS;
  const int imax=Nxx_plus_2NGHOSTS[0]-NGHOSTS,jmax=Nxx_plus_2NGHOSTS[1]-NGHOSTS,kmax=Nxx_plus_2NGHOSTS[2]-NGHOSTS;
  
  const REAL dz = dxx[2];

  REAL error_int_numer=0,error_int_denom=0;

  int num_vel_limits=0,num_vel_nulls_current_sheet=0;

  GiRaFFE_NRPy_update_metric_det_inverse(Nxx_plus_2NGHOSTS, dxx, xx,aux_gfs);

#pragma omp parallel for reduction(+:error_int_numer,error_int_denom,num_vel_limits,num_vel_nulls_current_sheet) schedule(static)
  for(int k=kmin;k<kmax;k++)
    for(int j=jmin;j<jmax;j++)
      for(int i=imin;i<imax;i++) {
        const int index = IDX3(i,j,k);
        const REAL xx0 = xx[0][i];
        const REAL xx1 = xx[1][j];
        const REAL xx2 = xx[2][k];
        const REAL rL = sqrt(xx0*xx0+xx1*xx1+xx2*xx2);
        if(rL>min_radius_inside_of_which_conserv_to_prims_FFE_and_FFE_evolution_is_DISABLED) {
          const REAL sqrtg = sqrt(aux_gfs[IDX4pt(GAMMADETGF, index)]); // Determinant of 3-metric

          // \gamma_{ij}, computed from \tilde{\gamma}_{ij}
          const REAL gxxL = aux_gfs[IDX4pt(GAMMADD00GF, index)];
          const REAL gxyL = aux_gfs[IDX4pt(GAMMADD01GF, index)];
          const REAL gxzL = aux_gfs[IDX4pt(GAMMADD02GF, index)];
          const REAL gyyL = aux_gfs[IDX4pt(GAMMADD11GF, index)];
          const REAL gyzL = aux_gfs[IDX4pt(GAMMADD12GF, index)];
          const REAL gzzL = aux_gfs[IDX4pt(GAMMADD22GF, index)];

          // \gamma^{ij} = psim4 * \tilde{\gamma}^{ij}
          const REAL gupxxL = aux_gfs[IDX4pt(GAMMAUU00GF, index)];
          const REAL gupxyL = aux_gfs[IDX4pt(GAMMAUU01GF, index)];
          const REAL gupxzL = aux_gfs[IDX4pt(GAMMAUU02GF, index)];
          const REAL gupyyL = aux_gfs[IDX4pt(GAMMAUU11GF, index)];
          const REAL gupyzL = aux_gfs[IDX4pt(GAMMAUU12GF, index)];
          const REAL gupzzL = aux_gfs[IDX4pt(GAMMAUU22GF, index)];

          // Read in magnetic field and momentum variables once from memory, since memory access is expensive:
          const REAL BU0L = aux_gfs[IDX4pt(BU0GF, index)];
          const REAL BU1L = aux_gfs[IDX4pt(BU1GF, index)];
          const REAL BU2L = aux_gfs[IDX4pt(BU2GF, index)];

          // End of page 7 on http://arxiv.org/pdf/1310.3274.pdf
          const REAL BtildexL = BU0L*sqrtg;
          const REAL BtildeyL = BU1L*sqrtg;
          const REAL BtildezL = BU2L*sqrtg;

          const REAL Btilde_xL = gxxL*BtildexL + gxyL*BtildeyL + gxzL*BtildezL;
          const REAL Btilde_yL = gxyL*BtildexL + gyyL*BtildeyL + gyzL*BtildezL;
          const REAL Btilde_zL = gxzL*BtildexL + gyzL*BtildeyL + gzzL*BtildezL;

          REAL StildeD0L = in_gfs[IDX4pt(STILDED0GF, index)];
          REAL StildeD1L = in_gfs[IDX4pt(STILDED1GF, index)];
          REAL StildeD2L = in_gfs[IDX4pt(STILDED2GF, index)];

          const REAL StildeD0_orig = StildeD0L;
          const REAL StildeD1_orig = StildeD1L;
          const REAL StildeD2_orig = StildeD2L;

          const REAL ValenciavU0_orig = aux_gfs[IDX4pt(VALENCIAVU0GF, index)];
          const REAL ValenciavU1_orig = aux_gfs[IDX4pt(VALENCIAVU1GF, index)];
          const REAL ValenciavU2_orig = aux_gfs[IDX4pt(VALENCIAVU2GF, index)];

          //const REAL alpL = alp[index];
          //const REAL fourpialpha = 4.0*M_PI*alpL;
          const REAL fourpi = 4.0*M_PI;

          //const REAL betaxL = betax[index];
          //const REAL betayL = betay[index];
          //const REAL betazL = betaz[index];
          //* 1. Just below Eq 90: Enforce orthogonality of B^i & S^i, so that B^i S_i = 0
          //*    Correction ${\tilde S}_i ->{\tilde S}_i - ({\tilde S}_j {\tilde B}^j) {\tilde B}_i/{\tilde B}^2$
          //*    NOTICE THAT THE {\tilde B}_i IS LOWERED, AS IT SHOULD BE. THIS IS A TYPO IN PASCHALIDIS ET AL.

          // First compute Btilde^i Stilde_i:
          const REAL BtildeiSt_i = StildeD0L*BtildexL + StildeD1L*BtildeyL + StildeD2L*BtildezL;
          //printf("xterm = %f ; yterm = %f ; zterm = %f\n",StildeD0L*BtildexL,StildeD1L*BtildeyL,StildeD2L*BtildezL);

          // Then compute (Btilde)^2
          const REAL Btilde2 = gxxL*BtildexL*BtildexL + gyyL*BtildeyL*BtildeyL + gzzL*BtildezL*BtildezL
            + 2.0*(gxyL*BtildexL*BtildeyL + gxzL*BtildexL*BtildezL + gyzL*BtildeyL*BtildezL);

#define APPLY_GRFFE_FIXES

          // Now apply constraint: Stilde_i = Stilde_i - (Btilde^i Stilde_i) / (Btilde)^2
#ifdef APPLY_GRFFE_FIXES
          StildeD0L -= BtildeiSt_i*Btilde_xL/Btilde2;
          StildeD1L -= BtildeiSt_i*Btilde_yL/Btilde2;
          StildeD2L -= BtildeiSt_i*Btilde_zL/Btilde2;
          //printf("BtildeiSt_i = %f ; Btilde2 = %f\n",BtildeiSt_i,Btilde2);
#endif
          // Now that tildeS_i has been fixed, let's compute tildeS^i:
          REAL mhd_st_upx = gupxxL*StildeD0L + gupxyL*StildeD1L + gupxzL*StildeD2L;
          REAL mhd_st_upy = gupxyL*StildeD0L + gupyyL*StildeD1L + gupyzL*StildeD2L;
          REAL mhd_st_upz = gupxzL*StildeD0L + gupyzL*StildeD1L + gupzzL*StildeD2L;

          // Just below Eq. 86 in http://arxiv.org/pdf/1310.3274.pdf:
          REAL St2 = StildeD0L*mhd_st_upx + StildeD1L*mhd_st_upy + StildeD2L*mhd_st_upz;
          //* 2. Eq. 92: Factor $f = \sqrt{(1-\gamma_{max}^{-2}){\tilde B}^4/(16 \pi^2 \gamma {\tilde S}^2)}$ 

#ifdef APPLY_GRFFE_FIXES
          const REAL gmax = GAMMA_SPEED_LIMIT;
          if(St2 > (1.0 - 1.0/(gmax*gmax))*Btilde2*Btilde2/ (16.0*M_PI*M_PI*sqrtg*sqrtg)) {
            const REAL fact = sqrt((1.0 - 1.0/(gmax*gmax))/St2)*Btilde2/(4.0*M_PI*sqrtg);

            //* 3. ${\tilde S}_i -> {\tilde S}_i min(1,f)
            StildeD0L *= MIN(1.0,fact);
            StildeD1L *= MIN(1.0,fact);
            StildeD2L *= MIN(1.0,fact);

            // Recompute S^i
            mhd_st_upx = gupxxL*StildeD0L + gupxyL*StildeD1L + gupxzL*StildeD2L;
            mhd_st_upy = gupxyL*StildeD0L + gupyyL*StildeD1L + gupyzL*StildeD2L;
            mhd_st_upz = gupxzL*StildeD0L + gupyzL*StildeD1L + gupzzL*StildeD2L;
            /*
            printf("%e %e %e | %e %e %e | %e %e %e | oldgamma: %e %e should be > %e vfix\n",x[index],y[index],z[index],
                   BU0L,BU1L,BU2L,
                   St2,(1.0 - 1.0/(gmax*gmax))*Btilde2*Btilde2/ (16.0*M_PI*M_PI*sqrtg*sqrtg),gmax,
                   sqrt(Btilde2 / (Btilde2 - 16*M_PI*M_PI*sqrtg*sqrtg * St2 / Btilde2) ) , Btilde2,16*M_PI*M_PI*sqrtg*sqrtg * St2 / Btilde2  );
            //exit(1);
            */
            // Recompute Stilde^2:
            St2 = StildeD0L*mhd_st_upx + StildeD1L*mhd_st_upy + StildeD2L*mhd_st_upz;

            if( St2 >= Btilde2*Btilde2/ (16.0*M_PI*M_PI*sqrtg*sqrtg) ) {
              printf("ERROR: Velocity cap fix wasn't effective; still have B^2 > E^2\n"); exit(1);
            }
            num_vel_limits++;
          }
#endif
          //* 4. Eq. 85: $v^i = 4 pi \alpha \gamma^{ij} {\tilde S}_j \gamma{-1/2} B^{-2} - \beta^i$: 

          // See, e.g., Eq 71 in http://arxiv.org/pdf/1310.3274.pdf
          // ... or end of page 7 on http://arxiv.org/pdf/1310.3274.pdf:
          const REAL B2 = Btilde2/(sqrtg*sqrtg);
          /* 
             Eq. 75: 
             v^i = \alpha \gamma^{ij} S_j / \mathcal{B}^2 - \beta^i
             Eq. 7: \mathcal{B}^{\mu} = B^{\mu}/\sqrt{4 \pi}
             -> v^i = 4 \pi \alpha \gamma^{ij} S_j / B^2 - \beta^i
             Eq. 79: \tilde{S_i} = \sqrt{\gamma} S_i
             -> v^i = 4 \pi \alpha \gamma^{ij} \tilde{S}_j / (\sqrt{\gamma} B^2) - \beta^i
          */
          // Modified from the original GiRaFFE to use Valencia, not drift velocity
          const REAL ValenciavU0L = fourpi*mhd_st_upx/(sqrtg*B2);
          const REAL ValenciavU1L = fourpi*mhd_st_upy/(sqrtg*B2);
          /* ValenciavU2L not necessarily const! See below. */
          REAL ValenciavU2L = fourpi*mhd_st_upz/(sqrtg*B2);
          //* 5. Eq. 94: ${\tilde n}_i v^i = 0$ in the current sheet region
          //     n^i is defined as the normal from the current sheet, which lies in the 
          //     xy-plane (z=0). So n = (0,0,1) 
#ifdef APPLY_GRFFE_FIXES
          if(current_sheet_null_v) {
            if (fabs(xx2) <= (4.0 + 1.0e-2)*dz ) {
              //ValenciavU2L = 0.0;
              ValenciavU2L = - (ValenciavU0L*gxzL + ValenciavU1L*gyzL) / gzzL;
                // FIXME: This is probably not right, but also definitely not the problem. 
            
              // ValenciavU2L reset: TYPICALLY WOULD RESET CONSERVATIVES TO BE CONSISTENT. LET'S NOT DO THAT, TO AVOID MESSING UP B-FIELDS

              if(1==1) {
                GiRaFFE_NRPy_compute_conservatives(gxxL, gxyL, gxzL, gyyL, gyzL, gzzL,
                                                 BU0L, BU1L, BU2L, ValenciavU0L, ValenciavU1L, ValenciavU2L,
                                                 /*const REAL betaxL, const REAL betayL, const REAL betazL, const REAL alpL,*/
                                                 sqrtg, &StildeD0L, &StildeD1L, &StildeD2L);
              }
              num_vel_nulls_current_sheet++;
            }
          }
#endif
          aux_gfs[IDX4pt(VALENCIAVU0GF, index)] = ValenciavU0L;
          aux_gfs[IDX4pt(VALENCIAVU1GF, index)] = ValenciavU1L;      
          aux_gfs[IDX4pt(VALENCIAVU2GF, index)] = ValenciavU2L;      
          //Now we compute the difference between original & new conservatives, for diagnostic purposes:
          //error_int_numer += fabs(StildeD0L - StildeD0_orig) + fabs(StildeD1L - StildeD1_orig) + fabs(StildeD2L - StildeD2_orig);
          //error_int_denom += fabs(StildeD0_orig) + fabs(StildeD1_orig) + fabs(StildeD2_orig);
          /*
            if(fabs(ValenciavU0_orig) > 1e-13 && fabs(ValenciavU0L-ValenciavU0_orig)/ValenciavU0_orig > 1e-2) printf("BAD ValenciavU0: %e %e | %e %e %e\n",ValenciavU0L,ValenciavU0_orig,x[index],y[index],z[index]);
            if(fabs(ValenciavU1_orig) > 1e-13 && fabs(ValenciavU1L-ValenciavU1_orig)/ValenciavU1_orig > 1e-2) printf("BAD ValenciavU1: %e %e | %e %e %e\n",ValenciavU1L,ValenciavU1_orig,x[index],y[index],z[index]);
            if(fabs(ValenciavU2_orig) > 1e-13 && fabs(ValenciavU2L-ValenciavU2_orig)/ValenciavU2_orig > 1e-2) printf("BAD ValenciavU2: %e %e | %e %e %e\n",ValenciavU2L,ValenciavU2_orig,x[index],y[index],z[index]);
          */
            error_int_numer += fabs(ValenciavU0L - ValenciavU0_orig) + fabs(ValenciavU1L - ValenciavU1_orig) + fabs(ValenciavU2L - ValenciavU2_orig);
            error_int_denom += fabs(ValenciavU0_orig) + fabs(ValenciavU1_orig) + fabs(ValenciavU2_orig);
          


          in_gfs[IDX4pt(STILDED0GF, index)] = StildeD0L;
          in_gfs[IDX4pt(STILDED1GF, index)] = StildeD1L;
          in_gfs[IDX4pt(STILDED2GF, index)] = StildeD2L;
        }
      }
}
