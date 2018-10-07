#include "stdio.h"
#include "stdlib.h"
#include "math.h"

#define REAL double
#define IDX4(g,i,j,k)   ( (i) + Nxx_plus_2NGHOSTS[0] * ( (j) + Nxx_plus_2NGHOSTS[1] * ( (k) + Nxx_plus_2NGHOSTS[2] * (g) ) ) )

const REAL wavespeed = 1.0;

const REAL time = 0.0;
const REAL kk0 = 1.0;
const REAL kk1 = 1.0;
const REAL kk2 = 1.0;

#define UUGF 0
#define VVGF 1

void initial_data(const int Nxx_plus_2NGHOSTS[3],REAL *xx[3], REAL *in_gfs) {
#include "ScalarWave_InitialData.h"
}

void rhs_eval(const int Nxx[3],const int NGHOSTS,const int Nxx_plus_2NGHOSTS[3],const REAL dxx[3], const REAL *in_gfs, REAL *rhs_gfs) {
#include "ScalarWave_RHSs.h"
}

const int MAXFACE = -1;
const int NUL     = +0;
const int MINFACE = +1;
#define  FACE_UPDATE(which_gf, i0min,i0max, i1min,i1max, i2min,i2max, FACEX0,FACEX1,FACEX2) \
  for(int i2=i2min;i2<i2max;i2++) for(int i1=i1min;i1<i1max;i1++) for(int i0=i0min;i0<i0max;i0++) { \
        gfs[IDX4(which_gf,i0,i1,i2)] =                                  \
          +3.0*gfs[IDX4(which_gf,i0+1*FACEX0,i1+1*FACEX1,i2+1*FACEX2)]  \
          -3.0*gfs[IDX4(which_gf,i0+2*FACEX0,i1+2*FACEX1,i2+2*FACEX2)]  \
          +1.0*gfs[IDX4(which_gf,i0+3*FACEX0,i1+3*FACEX1,i2+3*FACEX2)]; \
      }

void apply_bcs(const int Nxx[3],const int NGHOSTS,const int Nxx_plus_2NGHOSTS[3],REAL *gfs) {
  // First xmin:
  for(int which_gf=0;which_gf<2;which_gf++) {
    int imin[3] = { NGHOSTS, NGHOSTS, NGHOSTS };
    int imax[3] = { Nxx_plus_2NGHOSTS[0]-NGHOSTS, Nxx_plus_2NGHOSTS[1]-NGHOSTS, Nxx_plus_2NGHOSTS[2]-NGHOSTS };
    for(int which_gz = 0; which_gz < NGHOSTS; which_gz++) {
      FACE_UPDATE(which_gf, imin[0]-1,imin[0], imin[1],imax[1], imin[2],imax[2], MINFACE,NUL,NUL); imin[0]--;
      FACE_UPDATE(which_gf, imax[0],imax[0]+1, imin[1],imax[1], imin[2],imax[2], MAXFACE,NUL,NUL); imax[0]++;

      FACE_UPDATE(which_gf, imin[0],imax[0], imin[1]-1,imin[1], imin[2],imax[2], NUL,MINFACE,NUL); imin[1]--;
      FACE_UPDATE(which_gf, imin[0],imax[0], imax[1],imax[1]+1, imin[2],imax[2], NUL,MAXFACE,NUL); imax[1]++;

      FACE_UPDATE(which_gf, imin[0],imax[0], imin[1],imax[1], imin[2]-1,imin[2], NUL,NUL,MINFACE); imin[2]--;
      FACE_UPDATE(which_gf, imin[0],imax[0], imin[1],imax[1], imax[2],imax[2]+1, NUL,NUL,MAXFACE); imax[2]++;
    }
  }
}

int main(int argc, const char *argv[]) {

  const int Nxx[3] = { 128, 128, 128 };
  const int NGHOSTS = 3;

  const REAL xxmin[3] = {-10.,-10.,-10. };
  const REAL xxmax[3] = { 10., 10., 10. };
  const int t_final = 20.0;
  const REAL CFL_FACTOR = 0.5;

  const int Nxx_plus_2NGHOSTS[3] = { Nxx[0]+2*NGHOSTS, Nxx[1]+2*NGHOSTS, Nxx[2]+2*NGHOSTS };
  const int Nxx_plus_2NGHOSTS_tot = Nxx_plus_2NGHOSTS[0]*Nxx_plus_2NGHOSTS[1]*Nxx_plus_2NGHOSTS[2];

  REAL *xx[3];
  for(int i=0;i<3;i++) xx[i] = (REAL *)malloc(sizeof(REAL)*Nxx_plus_2NGHOSTS[i]);
  REAL *evol_gfs = (REAL *)malloc(sizeof(REAL) * 2 /* 2 gridfunctions */ * Nxx_plus_2NGHOSTS_tot);
  REAL *rhs_gfs  = (REAL *)malloc(sizeof(REAL) * 2 /* 2 gridfunctions */ * Nxx_plus_2NGHOSTS_tot);
  REAL *aux_gfs  = (REAL *)malloc(sizeof(REAL) * 2 /* 2 gridfunctions */ * Nxx_plus_2NGHOSTS_tot);

  // xx[0][i] = xxmin[0] + (i-NGHOSTS)*dxx[0]
  REAL dxx[3];
  for(int i=0;i<3;i++) dxx[i] = (xxmax[i] - xxmin[i]) / ((REAL)Nxx[i]);
#define MIN(A, B) ( ((A) < (B)) ? (A) : (B) )
  REAL dt = CFL_FACTOR * MIN(dxx[0],MIN(dxx[1],dxx[2]));
  REAL Nt = t_final / dt + 1;

  for(int i=0;i<3;i++) {
    for(int j=0;j<Nxx_plus_2NGHOSTS[i];j++) {
      xx[i][j] = xxmin[i] + (j-NGHOSTS)*dxx[i];
    }
  }

  for(int which_gf=0;which_gf<1;which_gf++) {
    for(int i2=NGHOSTS;i2<NGHOSTS+Nxx[2];i2++) for(int i1=NGHOSTS;i1<NGHOSTS+Nxx[1];i1++) for(int i0=NGHOSTS;i0<NGHOSTS+Nxx[0];i0++) {
          evol_gfs[IDX4(which_gf,i0,i1,i2)] = 1.0;
        }
  }
  apply_bcs(Nxx,NGHOSTS,Nxx_plus_2NGHOSTS,evol_gfs);

  initial_data(Nxx_plus_2NGHOSTS,xx, evol_gfs);

  for(int n=0;n<=Nt;n++) {
    rhs_eval(Nxx,NGHOSTS,Nxx_plus_2NGHOSTS,dxx, evol_gfs, rhs_gfs);
  }
  
  free(aux_gfs);
  free(rhs_gfs);
  free(evol_gfs);
  for(int i=0;i<3;i++) free(xx[i]);
  return 0;
}
