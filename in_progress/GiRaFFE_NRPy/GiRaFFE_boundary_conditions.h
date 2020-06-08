
// Currently, we're using basic Cartesian boundary conditions, pending fixes by Zach.
// Part P8a: Declare boundary condition FACE_UPDATE macro,
//          which updates a single face of the 3D grid cube
//          using quadratic polynomial extrapolation.
// Basic extrapolation boundary conditions
#define  FACE_UPDATE(which_gf, i0min,i0max, i1min,i1max, i2min,i2max, FACEX0,FACEX1,FACEX2) \
  for(int i2=i2min;i2<i2max;i2++) for(int i1=i1min;i1<i1max;i1++) for(int i0=i0min;i0<i0max;i0++) { \
        gfs[IDX4(which_gf,i0,i1,i2)] =                                  \
          +2.0*gfs[IDX4(which_gf,i0+1*FACEX0,i1+1*FACEX1,i2+1*FACEX2)]  \
          -1.0*gfs[IDX4(which_gf,i0+2*FACEX0,i1+2*FACEX1,i2+2*FACEX2)]; \
      }
//          +1.0*gfs[IDX4(which_gf,i0+3*FACEX0,i1+3*FACEX1,i2+3*FACEX2)]; \
    // Basic Copy boundary conditions
#define  FACE_UPDATE_COPY(which_gf, i0min,i0max, i1min,i1max, i2min,i2max, FACEX0,FACEX1,FACEX2) \
  for(int i2=i2min;i2<i2max;i2++) for(int i1=i1min;i1<i1max;i1++) for(int i0=i0min;i0<i0max;i0++) { \
        gfs[IDX4(which_gf,i0,i1,i2)] = gfs[IDX4(which_gf,i0+1*FACEX0,i1+1*FACEX1,i2+1*FACEX2)]; \
      }

// Part P8b: Boundary condition driver routine: Apply BCs to all six
//          boundary faces of the cube, filling in the innermost
//          ghost zone first, and moving outward.
const int MAXFACE = -1;
const int NUL     = +0;
const int MINFACE = +1;// This macro acts differently in that it acts on an entire 3-vector of gfs, instead of 1.
// which_gf_0 corresponds to the zeroth component of that vector. The if statements only
// evaluate true if the velocity is directed inwards on the face in consideration.
#define  FACE_UPDATE_OUTFLOW(which_gf, i0min,i0max, i1min,i1max, i2min,i2max, FACEX0,FACEX1,FACEX2) \
  for(int i2=i2min;i2<i2max;i2++) for(int i1=i1min;i1<i1max;i1++) for(int i0=i0min;i0<i0max;i0++) { \
      aux_gfs[IDX4(which_gf,i0,i1,i2)] =                                    \
          +2.0*aux_gfs[IDX4(which_gf,i0+1*FACEX0,i1+1*FACEX1,i2+1*FACEX2)]  \
          -1.0*aux_gfs[IDX4(which_gf,i0+2*FACEX0,i1+2*FACEX1,i2+2*FACEX2)]; \
  }
/*      aux_gfs[IDX4(which_gf_0+1,i0,i1,i2)] =                                    \
          +3.0*aux_gfs[IDX4(which_gf_0+1,i0+1*FACEX0,i1+1*FACEX1,i2+1*FACEX2)]  \
          -3.0*aux_gfs[IDX4(which_gf_0+1,i0+2*FACEX0,i1+2*FACEX1,i2+2*FACEX2)]  \
          +1.0*aux_gfs[IDX4(which_gf_0+1,i0+3*FACEX0,i1+3*FACEX1,i2+3*FACEX2)]; \
      aux_gfs[IDX4(which_gf_0+2,i0,i1,i2)] =                                    \
          +3.0*aux_gfs[IDX4(which_gf_0+2,i0+1*FACEX0,i1+1*FACEX1,i2+1*FACEX2)]  \
          -3.0*aux_gfs[IDX4(which_gf_0+2,i0+2*FACEX0,i1+2*FACEX1,i2+2*FACEX2)]  \
          +1.0*aux_gfs[IDX4(which_gf_0+2,i0+3*FACEX0,i1+3*FACEX1,i2+3*FACEX2)]; \
      if(FACEX0*aux_gfs[IDX4(which_gf_0+0,i0,i1,i2)] > 0.0) {                   \
          aux_gfs[IDX4(which_gf_0+0,i0,i1,i2)] = 0.0;                           \
      }                                                                         \
      if(FACEX1*aux_gfs[IDX4(which_gf_0+1,i0,i1,i2)] > 0.0) {                   \
          aux_gfs[IDX4(which_gf_0+1,i0,i1,i2)] = 0.0;                           \
      }                                                                         \
      if(FACEX2*aux_gfs[IDX4(which_gf_0+2,i0,i1,i2)] > 0.0) {                   \
          aux_gfs[IDX4(which_gf_0+2,i0,i1,i2)] = 0.0;                           \
      }                                                                         \
*/
void apply_bcs(const int Nxx[3],const int Nxx_plus_2NGHOSTS[3],REAL *gfs,REAL *aux_gfs) {
    // First, we apply extrapolation boundary conditions to AD
#pragma omp parallel for
    for(int which_gf=0;which_gf<NUM_EVOL_GFS;which_gf++) {
    if(which_gf < STILDED0GF || which_gf > STILDED2GF) {
    int imin[3] = { NGHOSTS, NGHOSTS, NGHOSTS };
    int imax[3] = { Nxx_plus_2NGHOSTS[0]-NGHOSTS, Nxx_plus_2NGHOSTS[1]-NGHOSTS, Nxx_plus_2NGHOSTS[2]-NGHOSTS };
    for(int which_gz = 0; which_gz < NGHOSTS; which_gz++) {
      // After updating each face, adjust imin[] and imax[]
      //   to reflect the newly-updated face extents.
      FACE_UPDATE(which_gf, imin[0]-1,imin[0], imin[1],imax[1], imin[2],imax[2], MINFACE,NUL,NUL); imin[0]--;
      FACE_UPDATE(which_gf, imax[0],imax[0]+1, imin[1],imax[1], imin[2],imax[2], MAXFACE,NUL,NUL); imax[0]++;

      FACE_UPDATE(which_gf, imin[0],imax[0], imin[1]-1,imin[1], imin[2],imax[2], NUL,MINFACE,NUL); imin[1]--;
      FACE_UPDATE(which_gf, imin[0],imax[0], imax[1],imax[1]+1, imin[2],imax[2], NUL,MAXFACE,NUL); imax[1]++;

      FACE_UPDATE(which_gf, imin[0],imax[0], imin[1],imax[1], imin[2]-1,imin[2], NUL,NUL,MINFACE);
        imin[2]--;
      FACE_UPDATE(which_gf, imin[0],imax[0], imin[1],imax[1], imax[2],imax[2]+1, NUL,NUL,MAXFACE);
        imax[2]++;
    }
    }
    }    // Apply outflow/extrapolation boundary conditions to ValenciavU by passing VALENCIAVU0 as which_gf_0
    for(int which_gf=VALENCIAVU0GF;which_gf<=VALENCIAVU2GF;which_gf++) {
    int imin[3] = { NGHOSTS, NGHOSTS, NGHOSTS };
    int imax[3] = { Nxx_plus_2NGHOSTS[0]-NGHOSTS, Nxx_plus_2NGHOSTS[1]-NGHOSTS, Nxx_plus_2NGHOSTS[2]-NGHOSTS };
    for(int which_gz = 0; which_gz < NGHOSTS; which_gz++) {
      FACE_UPDATE_OUTFLOW(which_gf, imin[0]-1,imin[0], imin[1],imax[1], imin[2],imax[2], MINFACE,NUL,NUL); imin[0]--;
      FACE_UPDATE_OUTFLOW(which_gf, imax[0],imax[0]+1, imin[1],imax[1], imin[2],imax[2], MAXFACE,NUL,NUL); imax[0]++;

      FACE_UPDATE_OUTFLOW(which_gf, imin[0],imax[0], imin[1]-1,imin[1], imin[2],imax[2], NUL,MINFACE,NUL); imin[1]--;
      FACE_UPDATE_OUTFLOW(which_gf, imin[0],imax[0], imax[1],imax[1]+1, imin[2],imax[2], NUL,MAXFACE,NUL); imax[1]++;

      FACE_UPDATE_OUTFLOW(which_gf, imin[0],imax[0], imin[1],imax[1], imin[2]-1,imin[2], NUL,NUL,MINFACE);
        imin[2]--;
      FACE_UPDATE_OUTFLOW(which_gf, imin[0],imax[0], imin[1],imax[1], imax[2],imax[2]+1, NUL,NUL,MAXFACE);
        imax[2]++;
    }
    }    // Then, we apply copy boundary conditions to StildeD and psi6Phi
/*#pragma omp parallel for
    for(int which_gf=3;which_gf<NUM_EVOL_GFS;which_gf++) {
    int imin[3] = { NGHOSTS, NGHOSTS, NGHOSTS };
    int imax[3] = { Nxx_plus_2NGHOSTS[0]-NGHOSTS, Nxx_plus_2NGHOSTS[1]-NGHOSTS, Nxx_plus_2NGHOSTS[2]-NGHOSTS };
    for(int which_gz = 0; which_gz < NGHOSTS; which_gz++) {
      // After updating each face, adjust imin[] and imax[]
      //   to reflect the newly-updated face extents.
      FACE_UPDATE_COPY(which_gf, imin[0]-1,imin[0], imin[1],imax[1], imin[2],imax[2], MINFACE,NUL,NUL); imin[0]--;
      FACE_UPDATE_COPY(which_gf, imax[0],imax[0]+1, imin[1],imax[1], imin[2],imax[2], MAXFACE,NUL,NUL); imax[0]++;

      FACE_UPDATE_COPY(which_gf, imin[0],imax[0], imin[1]-1,imin[1], imin[2],imax[2], NUL,MINFACE,NUL); imin[1]--;
      FACE_UPDATE_COPY(which_gf, imin[0],imax[0], imax[1],imax[1]+1, imin[2],imax[2], NUL,MAXFACE,NUL); imax[1]++;

      FACE_UPDATE_COPY(which_gf, imin[0],imax[0], imin[1],imax[1], imin[2]-1,imin[2], NUL,NUL,MINFACE); imin[2]--;
      FACE_UPDATE_COPY(which_gf, imin[0],imax[0], imin[1],imax[1], imax[2],imax[2]+1, NUL,NUL,MAXFACE); imax[2]++;
    }
    }*/
}// A supplement to the boundary conditions for debugging. This will overwrite data with exact conditions
void FACE_UPDATE_EXACT(const int Nxx[3],const int Nxx_plus_2NGHOSTS[3],REAL *xx[3],
                       const int n, const REAL dt,REAL *out_gfs,REAL *aux_gfs,
                       const int i0min,const int i0max, const int i1min,const int i1max, const int i2min,const int i2max,
                       const int FACEX0,const int FACEX1,const int FACEX2) {
  for(int i2=i2min;i2<i2max;i2++) for(int i1=i1min;i1<i1max;i1++) for(int i0=i0min;i0<i0max;i0++) {
    REAL xx0 = xx[0][i0]-n*dt;
    REAL xx1 = xx[1][i1];
    REAL xx2 = xx[2][i2];
    if(xx0<=lbound) {
#include "GiRaFFEfood_A_v_1D_tests_left.h"
    }
    else if (xx0<rbound) {
#include "GiRaFFEfood_A_v_1D_tests_center.h"
    }
    else {
#include "GiRaFFEfood_A_v_1D_tests_right.h"
    }
    out_gfs[IDX4(PSI6PHIGF, i0,i1,i2)] = 0.0;
  }
}void apply_bcs_EXACT(const int Nxx[3],const int Nxx_plus_2NGHOSTS[3],REAL *xx[3],
                     const int n, const REAL dt,
                     REAL *out_gfs,REAL *aux_gfs) {
    int imin[3] = { NGHOSTS, NGHOSTS, NGHOSTS };
    int imax[3] = { Nxx_plus_2NGHOSTS[0]-NGHOSTS, Nxx_plus_2NGHOSTS[1]-NGHOSTS, Nxx_plus_2NGHOSTS[2]-NGHOSTS };
    for(int which_gz = 0; which_gz < NGHOSTS; which_gz++) {
      // After updating each face, adjust imin[] and imax[]
      //   to reflect the newly-updated face extents.
      // Right now, we only want to update the xmin and xmax faces with the exact data.
      FACE_UPDATE_EXACT(Nxx,Nxx_plus_2NGHOSTS,xx,n,dt,out_gfs,aux_gfs,imin[0]-1,imin[0], imin[1],imax[1], imin[2],imax[2], MINFACE,NUL,NUL);
      imin[0]--;
      FACE_UPDATE_EXACT(Nxx,Nxx_plus_2NGHOSTS,xx,n,dt,out_gfs,aux_gfs,imax[0],imax[0]+1, imin[1],imax[1], imin[2],imax[2], MAXFACE,NUL,NUL);
      imax[0]++;

      FACE_UPDATE_EXACT(Nxx,Nxx_plus_2NGHOSTS,xx,n,dt,out_gfs,aux_gfs,imin[0],imax[0], imin[1]-1,imin[1], imin[2],imax[2], NUL,MINFACE,NUL);
      imin[1]--;
      FACE_UPDATE_EXACT(Nxx,Nxx_plus_2NGHOSTS,xx,n,dt,out_gfs,aux_gfs,imin[0],imax[0], imax[1],imax[1]+1, imin[2],imax[2], NUL,MAXFACE,NUL);
      imax[1]++;

      FACE_UPDATE_EXACT(Nxx,Nxx_plus_2NGHOSTS,xx,n,dt,out_gfs,aux_gfs,imin[0],imax[0], imin[1],imax[1], imin[2]-1,imin[2], NUL,NUL,MINFACE);
      imin[2]--;
      FACE_UPDATE_EXACT(Nxx,Nxx_plus_2NGHOSTS,xx,n,dt,out_gfs,aux_gfs,imin[0],imax[0], imin[1],imax[1], imax[2],imax[2]+1, NUL,NUL,MAXFACE);
      imax[2]++;
    }
}// A supplement to the boundary conditions for debugging. This will overwrite data with exact conditions
void FACE_UPDATE_EXACT_StildeD(const int Nxx[3],const int Nxx_plus_2NGHOSTS[3],REAL *xx[3],
                               REAL *out_gfs,REAL *out_gfs_exact,
                               const int i0min,const int i0max, const int i1min,const int i1max, const int i2min,const int i2max,
                               const int FACEX0,const int FACEX1,const int FACEX2) {
    // This is currently modified to calculate more exact boundary conditions for StildeD. Rename if it works.
    /*for(int i2=i2min;i2<i2max;i2++) for(int i1=i1min;i1<i1max;i1++) for(int i0=i0min;i0<i0max;i0++) {
#include "GiRaFFEfood_HO_Stilde.h"
    }*/
      /*idx = IDX3(i0,i1,i2);
      out_gfs[IDX4pt(STILDED0GF,idx)] = out_gfs_exact[IDX4pt(STILDED0GF,idx)];
      out_gfs[IDX4pt(STILDED1GF,idx)] = out_gfs_exact[IDX4pt(STILDED1GF,idx)];
      out_gfs[IDX4pt(STILDED2GF,idx)] = out_gfs_exact[IDX4pt(STILDED2GF,idx)];*/
}void apply_bcs_EXACT_StildeD(const int Nxx[3],const int Nxx_plus_2NGHOSTS[3],REAL *xx[3],
                             REAL *out_gfs,REAL *out_gfs_exact) {
    int imin[3] = { NGHOSTS, NGHOSTS, NGHOSTS };
    int imax[3] = { Nxx_plus_2NGHOSTS[0]-NGHOSTS, Nxx_plus_2NGHOSTS[1]-NGHOSTS, Nxx_plus_2NGHOSTS[2]-NGHOSTS };
    for(int which_gz = 0; which_gz < NGHOSTS; which_gz++) {
      // After updating each face, adjust imin[] and imax[]
      //   to reflect the newly-updated face extents.
      // Right now, we only want to update the xmin and xmax faces with the exact data.
      FACE_UPDATE_EXACT_StildeD(Nxx,Nxx_plus_2NGHOSTS,xx,out_gfs,out_gfs_exact,imin[0]-1,imin[0], imin[1],imax[1], imin[2],imax[2], MINFACE,NUL,NUL);
      imin[0]--;
      FACE_UPDATE_EXACT_StildeD(Nxx,Nxx_plus_2NGHOSTS,xx,out_gfs,out_gfs_exact,imax[0],imax[0]+1, imin[1],imax[1], imin[2],imax[2], MAXFACE,NUL,NUL);
      imax[0]++;

      //FACE_UPDATE_EXACT_StildeD(Nxx,Nxx_plus_2NGHOSTS,xx,out_gfs,out_gfs_exact,imin[0],imax[0], imin[1]-1,imin[1], imin[2],imax[2], NUL,MINFACE,NUL);
      imin[1]--;
      //FACE_UPDATE_EXACT_StildeD(Nxx,Nxx_plus_2NGHOSTS,xx,out_gfs,out_gfs_exact,imin[0],imax[0], imax[1],imax[1]+1, imin[2],imax[2], NUL,MAXFACE,NUL);
      imax[1]++;

      //FACE_UPDATE_EXACT_StildeD(Nxx,Nxx_plus_2NGHOSTS,xx,out_gfs,out_gfs_exact,imin[0],imax[0], imin[1],imax[1], imin[2]-1,imin[2], NUL,NUL,MINFACE);
      imin[2]--;
      //FACE_UPDATE_EXACT_StildeD(Nxx,Nxx_plus_2NGHOSTS,xx,out_gfs,out_gfs_exact,imin[0],imax[0], imin[1],imax[1], imax[2],imax[2]+1, NUL,NUL,MAXFACE);
      imax[2]++;
    }
}
