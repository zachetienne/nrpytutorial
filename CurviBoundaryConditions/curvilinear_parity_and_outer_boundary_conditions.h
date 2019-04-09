
// First we define the struct that will be used to store the 10 parity conditions at all gridpoints:
// We store the 10 parity conditions in a struct consisting of 10 integers, one for each condition.
// Note that these conditions can only take one of two values: +1 or -1.
typedef struct parity_conditions {
  int8_t parity[10];
} parity_condition;

typedef struct ghostzone_map {
  short i0,i1,i2;
} gz_map;

void set_bc_parity_conditions(REAL parity[10], const REAL xx0,const REAL xx1,const REAL xx2, 
                              const REAL xx0_inbounds,const REAL xx1_inbounds,const REAL xx2_inbounds) {
    #include "set_parity_conditions.h"
}

void set_up_bc_gz_map_and_parity_conditions(const int Nxx_plus_2NGHOSTS[3], REAL *xx[3], 
                                            const REAL dxx[3], const REAL xxmin[3], const REAL xxmax[3], 
                                            gz_map *bc_gz_map, parity_condition *bc_parity_conditions) {
  LOOP_REGION(0,Nxx_plus_2NGHOSTS[0],0,Nxx_plus_2NGHOSTS[1],0,Nxx_plus_2NGHOSTS[2]) {
    // First find Cartesian coordinate corresponding to (x_0,x_1,x_2)=(xx[0][i0],xx[1][i1],xx[2][i2]):
    REAL xCart[3];
    xxCart(xx, i0,i1,i2, xCart);
    REAL Cartx = xCart[0];
    REAL Carty = xCart[1];
    REAL Cartz = xCart[2];
    
    // Next find the (i0_inbounds,i1_inbounds,i2_inbounds) corresponding to the above Cartesian coordinate.
    //   If (i0_inbounds,i1_inbounds,i2_inbounds) is in a ghost zone, then it must equal (i0,i1,i2), and
    //      the point is an outer boundary point.
    //   Otherwise (i0_inbounds,i1_inbounds,i2_inbounds) is in the grid interior, and data at (i0,i1,i2)
    //      must be replaced with data at (i0_inbounds,i1_inbounds,i2_inbounds), but multiplied by the
    //      appropriate parity condition (+/- 1).
    REAL Cart_to_xx0_inbounds,Cart_to_xx1_inbounds,Cart_to_xx2_inbounds;
#include "Cart_to_xx.h"
    int i0_inbounds = (int)( (Cart_to_xx0_inbounds - xxmin[0] - (1.0/2.0)*dxx[0] + ((REAL)NGHOSTS)*dxx[0])/dxx[0] + 0.5 ); 
    int i1_inbounds = (int)( (Cart_to_xx1_inbounds - xxmin[1] - (1.0/2.0)*dxx[1] + ((REAL)NGHOSTS)*dxx[1])/dxx[1] + 0.5 );
    int i2_inbounds = (int)( (Cart_to_xx2_inbounds - xxmin[2] - (1.0/2.0)*dxx[2] + ((REAL)NGHOSTS)*dxx[2])/dxx[2] + 0.5 );

    REAL xCart_orig[3]; for(int ii=0;ii<3;ii++) xCart_orig[ii] = xCart[ii];
    xxCart(xx, i0_inbounds,i1_inbounds,i2_inbounds, xCart);

#define EPS_ABS 1e-8
    if(fabs( (double)(xCart_orig[0] - xCart[0]) ) > EPS_ABS ||
       fabs( (double)(xCart_orig[1] - xCart[1]) ) > EPS_ABS ||
       fabs( (double)(xCart_orig[2] - xCart[2]) ) > EPS_ABS) {
      printf("Error. Cartesian disagreement: ( %.15e %.15e %.15e ) != ( %.15e %.15e %.15e )\n",
             (double)xCart_orig[0],(double)xCart_orig[1],(double)xCart_orig[2],
             (double)xCart[0],(double)xCart[1],(double)xCart[2]);
      exit(1);
    }

    if(i0_inbounds-i0 == 0 && i1_inbounds-i1 == 0 && i2_inbounds-i2 == 0) {
      bc_gz_map[IDX3(i0,i1,i2)].i0=-1;
      bc_gz_map[IDX3(i0,i1,i2)].i1=-1;
      bc_gz_map[IDX3(i0,i1,i2)].i2=-1;
      for(int which_parity=0; which_parity<10; which_parity++) {
        bc_parity_conditions[IDX3(i0,i1,i2)].parity[which_parity] = 1;
      }
    } else {
      bc_gz_map[IDX3(i0,i1,i2)].i0=i0_inbounds;
      bc_gz_map[IDX3(i0,i1,i2)].i1=i1_inbounds;
      bc_gz_map[IDX3(i0,i1,i2)].i2=i2_inbounds;
      const REAL xx0 = xx[0][i0];
      const REAL xx1 = xx[1][i1];
      const REAL xx2 = xx[2][i2];
      const REAL xx0_inbounds = xx[0][i0_inbounds];
      const REAL xx1_inbounds = xx[1][i1_inbounds];
      const REAL xx2_inbounds = xx[2][i2_inbounds];
      REAL REAL_parity_array[10];
      set_bc_parity_conditions(REAL_parity_array,  xx0,xx1,xx2, xx0_inbounds,xx1_inbounds,xx2_inbounds);
      for(int whichparity=0;whichparity<10;whichparity++) {
          //printf("Good? Parity %d evaluated to %e\n",whichparity,REAL_parity_array[whichparity]);
          // Perform sanity check on parity array output: should be +1 or -1 to within 8 significant digits:
          if( (REAL_parity_array[whichparity]  > 0 && fabs(REAL_parity_array[whichparity] - (+1)) > 1e-8) ||
              (REAL_parity_array[whichparity] <= 0 && fabs(REAL_parity_array[whichparity] - (-1)) > 1e-8) ) {
              printf("Error. Parity evaluated to %e , which is not within 8 significant digits of +1 or -1.",(double)REAL_parity_array[whichparity]);
              exit(1);
          }
          if(REAL_parity_array[whichparity] < 0.0) bc_parity_conditions[IDX3(i0,i1,i2)].parity[whichparity] = -1;
          if(REAL_parity_array[whichparity] > 0.0) bc_parity_conditions[IDX3(i0,i1,i2)].parity[whichparity] = +1;
      }
    }
  }
}
// Part P6: Declare boundary condition OB_UPDATE macro,
//          which updates a single face of the 3D grid cube
//          with
//          1. quadratic polynomial extrapolation, if the face
//             corresponds to an outer boundary, or
//          2. parity condition, if the face maps to a point
//             in the grid interior.
const int MAXFACE = -1;
const int NUL     = +0;
const int MINFACE = +1;


#define OB_UPDATE(inner,which_gf, bc_gz_map,bc_parity_conditions, i0min,i0max, i1min,i1max, i2min,i2max, FACEX0,FACEX1,FACEX2) \
  LOOP_REGION(i0min,i0max, i1min,i1max, i2min,i2max) {                  \
    const int idx3 = IDX3(i0,i1,i2);                                    \
    if(bc_gz_map[idx3].i0 == -1 && inner==0) {                          \
      gfs[IDX4(which_gf,i0,i1,i2)] =                                    \
        +3.0*gfs[IDX4(which_gf,i0+1*FACEX0,i1+1*FACEX1,i2+1*FACEX2)]    \
        -3.0*gfs[IDX4(which_gf,i0+2*FACEX0,i1+2*FACEX1,i2+2*FACEX2)]    \
        +1.0*gfs[IDX4(which_gf,i0+3*FACEX0,i1+3*FACEX1,i2+3*FACEX2)];   \
    } else if(bc_gz_map[idx3].i0 != -1 && inner==1) {                   \
     gfs[IDX4(which_gf,i0,i1,i2)] =                                    \
        ( (REAL)bc_parity_conditions[idx3].parity[gfs_parity[which_gf]] )* \
                                             gfs[IDX4(which_gf,           \
                                                    bc_gz_map[idx3].i0, \
                                                    bc_gz_map[idx3].i1, \
                                                    bc_gz_map[idx3].i2)]; \
    }                                                                   \
  }

// Part P7: Boundary condition driver routine: Apply BCs to all six
//          boundary faces of the cube, filling in the innermost
//          ghost zone first, and moving outward.
void apply_bcs(const int Nxx[3],const int Nxx_plus_2NGHOSTS[3],
               gz_map *bc_gz_map,parity_condition *bc_parity_conditions,int num_gfs,const int8_t *gfs_parity, REAL *gfs) {
#pragma omp parallel for
  for(int which_gf=0;which_gf<num_gfs;which_gf++) {
    int imin[3] = { NGHOSTS, NGHOSTS, NGHOSTS };
    int imax[3] = { Nxx_plus_2NGHOSTS[0]-NGHOSTS, Nxx_plus_2NGHOSTS[1]-NGHOSTS, Nxx_plus_2NGHOSTS[2]-NGHOSTS };
    for(int which_gz = 0; which_gz < NGHOSTS; which_gz++) {
      for(int inner=0;inner<2;inner++) {
        // After updating each face, adjust imin[] and imax[] 
        //   to reflect the newly-updated face extents.
        OB_UPDATE(inner,which_gf, bc_gz_map,bc_parity_conditions, imin[0]-1,imin[0], imin[1],imax[1], imin[2],imax[2], MINFACE,NUL,NUL); imin[0]--;
        OB_UPDATE(inner,which_gf, bc_gz_map,bc_parity_conditions, imax[0],imax[0]+1, imin[1],imax[1], imin[2],imax[2], MAXFACE,NUL,NUL); imax[0]++;

        OB_UPDATE(inner,which_gf, bc_gz_map,bc_parity_conditions, imin[0],imax[0], imin[1]-1,imin[1], imin[2],imax[2], NUL,MINFACE,NUL); imin[1]--;
        OB_UPDATE(inner,which_gf, bc_gz_map,bc_parity_conditions, imin[0],imax[0], imax[1],imax[1]+1, imin[2],imax[2], NUL,MAXFACE,NUL); imax[1]++;

        OB_UPDATE(inner,which_gf, bc_gz_map,bc_parity_conditions, imin[0],imax[0], imin[1],imax[1], imin[2]-1,imin[2], NUL,NUL,MINFACE); imin[2]--;
        OB_UPDATE(inner,which_gf, bc_gz_map,bc_parity_conditions, imin[0],imax[0], imin[1],imax[1], imax[2],imax[2]+1, NUL,NUL,MAXFACE); imax[2]++;
        if(inner==0) { for(int ii=0;ii<3;ii++) {imin[ii]++; imax[ii]--;} }
      }
    }
  }
}