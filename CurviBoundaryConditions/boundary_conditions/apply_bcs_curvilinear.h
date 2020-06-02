
// Declare boundary condition BC_UPDATE_OUTER macro,
//          which updates a single outer boundary face
//          of the 3D grid cube using quadratic polynomial
//          extrapolation.

#define BC_UPDATE_OUTER(which_gf, i0,i1,i2, FACEX0,FACEX1,FACEX2) {     \
    const int idx3 = IDX3S(i0,i1,i2);                                   \
    gfs[IDX4S(which_gf,i0,i1,i2)] =                                     \
        +3.0*gfs[IDX4S(which_gf,i0+1*FACEX0,i1+1*FACEX1,i2+1*FACEX2)]   \
        -3.0*gfs[IDX4S(which_gf,i0+2*FACEX0,i1+2*FACEX1,i2+2*FACEX2)]   \
        +1.0*gfs[IDX4S(which_gf,i0+3*FACEX0,i1+3*FACEX1,i2+3*FACEX2)];  \
    }

// Curvilinear boundary condition driver routine: Apply BCs to all six
//          boundary faces of the 3D numerical domain, filling in the
//          innermost ghost zone layer first, and moving outward.

void apply_bcs_curvilinear_single_gf(const paramstruct *restrict params, const bc_struct *restrict bcstruct, const int8_t *restrict gfs_parity, const int which_gf, REAL *restrict gfs) {

#include "RELATIVE_PATH__set_Cparameters.h" /* Header file containing correct #include for set_Cparameters.h;
                                             * accounting for the relative path */

  for(int which_gz = 0; which_gz < NGHOSTS; which_gz++) {
    // First apply OUTER boundary conditions,
    //   in case an INNER (parity) boundary point
    //   needs data at the outer boundary:
    // After updating each face, adjust imin[] and imax[]
    //   to reflect the newly-updated face extents.
    for(int pt=0;pt<bcstruct->num_ob_gz_pts[which_gz];pt++) {
      BC_UPDATE_OUTER(which_gf,
                      bcstruct->outer[which_gz][pt].outer_bc_dest_pt.i0,
                      bcstruct->outer[which_gz][pt].outer_bc_dest_pt.i1,
                      bcstruct->outer[which_gz][pt].outer_bc_dest_pt.i2,
                      bcstruct->outer[which_gz][pt].FACEi0,
                      bcstruct->outer[which_gz][pt].FACEi1,
                      bcstruct->outer[which_gz][pt].FACEi2);
    }

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
    } // END for(int pt=0;pt<num_ib_gz_pts[which_gz];pt++)
  } // END for(int which_gz = 0; which_gz < NGHOSTS; which_gz++)
} // END function

void apply_bcs_curvilinear(const paramstruct *restrict params, const bc_struct *restrict bcstruct,
                           const int NUM_GFS, const int8_t *restrict gfs_parity, REAL *restrict gfs) {
#pragma omp parallel for
    for(int which_gf=0;which_gf<NUM_GFS;which_gf++) {

        apply_bcs_curvilinear_single_gf(params, bcstruct, gfs_parity, which_gf, gfs);

    } // END for(int which_gf=0;which_gf<NUM_GFS;which_gf++)
} // END function
