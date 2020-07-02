//          Boundary condtion driver routine: Apply BCs to all
//          boundary faces of the 3D numerical domain, filling in the
//          outer boundary ghost zone layers, starting with the innermost
//          layer and working outward.

#include "sommerfeld_params.h"
#include "radial_derivative.h"
#include <string.h>


void apply_bcs_sommerfeld(const paramstruct *restrict params, REAL *restrict xx[3],
                          const bc_struct *restrict bcstruct, const int NUM_GFS,
                          const int8_t *restrict gfs_parity, REAL *restrict gfs,
                          REAL *restrict rhs_gfs) {

    #pragma omp parallel for
        for(int which_gf=0;which_gf<NUM_GFS;which_gf++) {
          REAL var_at_infinity = evolgf_at_inf[which_gf];
          REAL radpower        = evolgf_radpower[which_gf];
          REAL char_speed      = evolgf_speed[which_gf];

          #include "RELATIVE_PATH__set_Cparameters.h" /* Header file containing correct #include for set_Cparameters.h;
                                                       * accounting for the relative path */


            for(int which_gz = 0; which_gz < NGHOSTS; which_gz++) {
                for(int pt=0;pt<bcstruct->num_ob_gz_pts[which_gz];pt++) {
                    int i0 = bcstruct->outer[which_gz][pt].outer_bc_dest_pt.i0;
                    int i1 = bcstruct->outer[which_gz][pt].outer_bc_dest_pt.i1;
                    int i2 = bcstruct->outer[which_gz][pt].outer_bc_dest_pt.i2;
                    int8_t FACEX0 = bcstruct->outer[which_gz][pt].FACEi0;
                    int8_t FACEX1 = bcstruct->outer[which_gz][pt].FACEi1;
                    int8_t FACEX2 = bcstruct->outer[which_gz][pt].FACEi2;

                    int8_t FACEXi[3] = {FACEX0, FACEX1, FACEX2};

                    // Initialize derivatives to crazy values, to ensure that
                    //   we will notice in case they aren't set properly.
                    REAL r = 1e100;
                    REAL dxU_fdD = 1e100;

                    contraction_term(params, which_gf, gfs, xx, FACEXi, i0, i1, i2, &r, &dxU_fdD);

                    REAL invr = 1./r;

                    REAL source_rhs = -char_speed*(dxU_fdD + invr*(gfs[IDX4S(which_gf,i0,i1,i2)] - var_at_infinity));
                    rhs_gfs[IDX4S(which_gf,i0,i1,i2)] = source_rhs;

                    /************* For radial falloff and the extrapolated k term *************/
                    if (radpower > 0) {

                      // Move one point away from gz point to compare pure advection to df/dt|interior

                      int i0_offset = i0+FACEX0;
                      int i1_offset = i1+FACEX1;
                      int i2_offset = i2+FACEX2;

                      // Initialize derivatives to crazy values, to ensure that
                      //   we will notice in case they aren't set properly.
                      REAL r_offset = 1e100;
                      REAL dxU_fdD_offset = 1e100;

                      contraction_term(params, which_gf, gfs, xx, FACEXi, i0_offset, i1_offset, i2_offset, &r_offset, &dxU_fdD_offset);

                      REAL invr_offset = 1./r_offset;

                      // Pure advection
                      REAL extrap_rhs = char_speed*(dxU_fdD_offset + invr_offset*(gfs[IDX4S(which_gf,i0_offset,i1_offset,i2_offset)] - var_at_infinity));

                      // Take difference between pure advection and df/dt|interior
                      REAL aux = rhs_gfs[IDX4S(which_gf,i0_offset,i1_offset,i2_offset)] + extrap_rhs;

                      // Solve for k/(r_gz)^n+1 term
                      rhs_gfs[IDX4S(which_gf,i0,i1,i2)] += aux*pow(r_offset*invr,radpower);
                    }
                }// END for(int pt=0;pt<num_ob_gz_pts[which_gz];pt++)

            // Then apply INNER (parity) boundary conditions:
            for(int pt=0;pt<bcstruct->num_ib_gz_pts[which_gz];pt++) {
                const int i0dest = bcstruct->inner[which_gz][pt].inner_bc_dest_pt.i0;
                const int i1dest = bcstruct->inner[which_gz][pt].inner_bc_dest_pt.i1;
                const int i2dest = bcstruct->inner[which_gz][pt].inner_bc_dest_pt.i2;
                const int i0src  = bcstruct->inner[which_gz][pt].inner_bc_src_pt.i0;
                const int i1src  = bcstruct->inner[which_gz][pt].inner_bc_src_pt.i1;
                const int i2src  = bcstruct->inner[which_gz][pt].inner_bc_src_pt.i2;
                const int8_t *prty= bcstruct->inner[which_gz][pt].parity;
//                printf("%d\n",bcstruct->inner_bc_parity[which_gz][pt].parity[gfs_parity[which_gf]]);
                gfs[IDX4S(which_gf,i0dest,i1dest,i2dest)] =
                        bcstruct->inner[which_gz][pt].parity[gfs_parity[which_gf]] * gfs[IDX4S(which_gf, i0src,i1src,i2src)];
            }// END for(int pt=0;pt<num_ib_gz_pts[which_gz];pt++)
        } // END for(int which_gz = 0; which_gz < NGHOSTS; which_gz++)
    } // END for(int which_gf=0;which_gf<NUM_GFS;which_gf++)
} // END function
