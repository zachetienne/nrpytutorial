/* Compute the part of A_i_rhs that excludes the gauge terms. I.e., we set
 *   A_i_rhs = \partial_t A_i = \psi^{6} (v^z B^x - v^x B^z)   here.
 */
static void A_i_rhs_no_gauge_terms(const int A_dirn,const paramstruct *params,gf_and_gz_struct *out_prims_r,gf_and_gz_struct *out_prims_l,
                                   REAL *psi6_pointer,REAL *cmax_1,REAL *cmin_1,REAL *cmax_2,REAL *cmin_2, REAL *A3_rhs) {
  #include "../set_Cparameters.h"

  // If A_dirn=1, then v1_offset=1 (v1=VY) and v2_offset=2 (v2=VZ)
  // If A_dirn=2, then v1_offset=2 (v1=VZ) and v2_offset=0 (v2=VX)
  // If A_dirn=3, then v1_offset=0 (v1=VX) and v2_offset=1 (v2=VY)
  const int v1_offset  = ((A_dirn-1)+1)%3,        v2_offset = ((A_dirn-1)+2)%3;

  const REAL *v1rr=out_prims_r[VXR+v1_offset].gf, *v2rr=out_prims_r[VXR+v2_offset].gf;
  const REAL *v1rl=out_prims_l[VXR+v1_offset].gf, *v2rl=out_prims_l[VXR+v2_offset].gf;
  const REAL *v1lr=out_prims_r[VXL+v1_offset].gf, *v2lr=out_prims_r[VXL+v2_offset].gf;
  const REAL *v1ll=out_prims_l[VXL+v1_offset].gf, *v2ll=out_prims_l[VXL+v2_offset].gf;

  const REAL *B1r=out_prims_r[BX_STAGGER+v1_offset].gf, *B1l=out_prims_l[BX_STAGGER+v1_offset].gf;
  const REAL *B2r=out_prims_r[BX_STAGGER+v2_offset].gf, *B2l=out_prims_l[BX_STAGGER+v2_offset].gf;

  /**** V DEPENDENCIES ****/
  /* In the case of Ax_rhs, we need v{y,z}{r,l} at (i,j+1/2,k+1/2).
   *    However, v{y,z}{r,l}{r,l} are defined at (i,j-1/2,k-1/2), so
   *    v{y,z}{r,l} at (i,j+1/2,k+1/2) is stored at v{y,z}{r,l}{r,l}(i,j+1,k+1).
   * In the case of Ay_rhs, we need v{x,z}{r,l} at (i+1/2,j,k+1/2).
   *    However, v{x,z}{r,l}{r,l} are defined at (i-1/2,j,k-1/2), so
   *    v{x,z}{r,l} at (i+1/2,j,k+1/2) is stored at v{x,z}{r,l}{r,l}(i+1,j,k+1).
   * In the case of Az_rhs, we need v{x,y}{r,l} at (i+1/2,j+1/2,k).
   *    However, v{x,y}{r,l}{r,l} are defined at (i-1/2,j-1/2,k), so
   *    v{x,y}{r,l} at (i+1/2,j+1/2,k) is stored at v{x,y}{r,l}{r,l}(i+1,j+1,k). */
  static const int vs_ijk_offset[4][3] = { {0,0,0} , {0,1,1} , {1,0,1} , {1,1,0} }; // Note that vs_ijk_offset[0] is UNUSED; we choose a 1-offset for convenience.

  /**** B DEPENDENCIES ****/
  /* In the case of Ax_rhs, we need B{y,z}{r,l} at (i,j+1/2,k+1/2).
   *    However, By_stagger{r,l} is defined at (i,j+1/2,k-1/2), and
   *             Bz_stagger{r,l} is defined at (i,j-1/2,k+1/2), so
   *             By_stagger{r,l} at (i,j+1/2,k+1/2) is stored at By_stagger{r,l}(i,j,k+1), and
   *             Bz_stagger{r,l} at (i,j+1/2,k+1/2) is stored at Bz_stagger{r,l}(i,j+1,k).
   * In the case of Ay_rhs, we need B{z,x}_stagger{r,l} at (i+1/2,j,k+1/2).
   *    However, Bz_stagger{r,l} is defined at (i-1/2,j,k+1/2), and
   *             Bx_stagger{r,l} is defined at (i+1/2,j,k-1/2), so
   *             Bz_stagger{r,l} at (i+1/2,j,k+1/2) is stored at Bz_stagger{r,l}(i+1,j,k), and
   *             Bx_stagger{r,l} at (i+1/2,j,k+1/2) is stored at Bx_stagger{r,l}(i,j,k+1).
   * In the case of Az_rhs, we need B{x,y}_stagger{r,l} at (i+1/2,j+1/2,k).
   *    However, Bx_stagger{r,l} is defined at (i+1/2,j-1/2,k), and
   *             By_stagger{r,l} is defined at (i-1/2,j+1/2,k), so
   *             Bx_stagger{r,l} at (i+1/2,j+1/2,k) is stored at Bx_stagger{r,l}(i,j+1,k), and
   *             By_stagger{r,l} at (i+1/2,j+1/2,k) is stored at By_stagger{r,l}(i+1,j,k).
   */
  static const int B1_ijk_offset[4][3] = { {0,0,0} , {0,0,1} , {1,0,0} , {0,1,0} }; // Note that B1_ijk_offset[0] is UNUSED; we choose a 1-offset for convenience.
  static const int B2_ijk_offset[4][3] = { {0,0,0} , {0,1,0} , {0,0,1} , {1,0,0} }; // Note that B2_ijk_offset[0] is UNUSED; we choose a 1-offset for convenience.

#pragma omp parallel for
  for(int k=NGHOSTS;k<Nxx_plus_2NGHOSTS2-NGHOSTS;k++) for(int j=NGHOSTS;j<Nxx_plus_2NGHOSTS1-NGHOSTS;j++) for(int i=NGHOSTS;i<Nxx_plus_2NGHOSTS0-NGHOSTS;i++) {
        const int index=IDX3S(i,j,k);
        // The following lines set the indices appropriately. See justification in exorbitant comments above.
        const int index_v =IDX3S(i+vs_ijk_offset[A_dirn][0],j+vs_ijk_offset[A_dirn][1],k+vs_ijk_offset[A_dirn][2]);
        const int index_B1=IDX3S(i+B1_ijk_offset[A_dirn][0],j+B1_ijk_offset[A_dirn][1],k+B1_ijk_offset[A_dirn][2]);
        const int index_B2=IDX3S(i+B2_ijk_offset[A_dirn][0],j+B2_ijk_offset[A_dirn][1],k+B2_ijk_offset[A_dirn][2]);

        // Stores 1/sqrt(gamma)==exp(6 phi) at (i+1/2,j+1/2,k) for Az, (i+1/2,j,k+1/2) for Ay, and (i,j+1/2,k+1/2) for Az.
        const REAL psi6_interped=psi6_pointer[index];

        const REAL B1lL = B1l[index_B1];
        const REAL B1rL = B1r[index_B1];
        const REAL B2lL = B2l[index_B2];
        const REAL B2rL = B2r[index_B2];

        const REAL A3_rhs_rr = psi6_interped*(v1rr[index_v]*B2rL - v2rr[index_v]*B1rL);
        const REAL A3_rhs_rl = psi6_interped*(v1rl[index_v]*B2rL - v2rl[index_v]*B1lL);
        const REAL A3_rhs_lr = psi6_interped*(v1lr[index_v]*B2lL - v2lr[index_v]*B1rL);
        const REAL A3_rhs_ll = psi6_interped*(v1ll[index_v]*B2lL - v2ll[index_v]*B1lL);


        // All variables for the A_i_rhs computation are now at the appropriate staggered point,
        //   so it's time to compute the HLL flux!

        // Note that with PPM, cmin and cmax are defined between ijk=3 and ijk<cctk_lsh[]-2 for all directions.
        const REAL cmax_1L = cmax_1[index_B2];
        const REAL cmin_1L = cmin_1[index_B2];
        const REAL cmax_2L = cmax_2[index_B1];
        const REAL cmin_2L = cmin_2[index_B1];

        const REAL B1tilder_minus_B1tildel = psi6_interped*( B1rL - B1lL );
        const REAL B2tilder_minus_B2tildel = psi6_interped*( B2rL - B2lL );

        /*---------------------------
         * Implement 2D HLL flux
         * [see Del Zanna, Bucciantini & Londrillo A&A 400, 397 (2003), Eq. (44)]
         *
         * Note that cmax/cmin (\alpha^{\pm}  as defined in that paper) is at a slightly DIFFERENT
         * point than that described in the Del Zanna et al paper (e.g., (i+1/2,j,k) instead of
         * (i+1/2,j+1/2,k) for F3).  Yuk Tung Liu discussed this point with M. Shibata,
         * who found that the effect is negligible.
         ---------------------------*/
        A3_rhs[index] = (cmax_1L*cmax_2L*A3_rhs_ll + cmax_1L*cmin_2L*A3_rhs_lr +
                         cmin_1L*cmax_2L*A3_rhs_rl + cmin_1L*cmin_2L*A3_rhs_rr)
          /( (cmax_1L+cmin_1L)*(cmax_2L+cmin_2L) )
          - cmax_1L*cmin_1L*(B2tilder_minus_B2tildel)/(cmax_1L+cmin_1L)
          + cmax_2L*cmin_2L*(B1tilder_minus_B1tildel)/(cmax_2L+cmin_2L);
      }
}
