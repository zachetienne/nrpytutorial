#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#ifndef REAL
#define REAL double
#include "NGHOSTS.h" // A NRPy+-generated file, which is set based on FD_CENTDERIVS_ORDER.
#include "../CurviBoundaryConditions/gridfunction_defines.h"
#define IDX4(g,i,j,k) \
( (i) + Nxx_plus_2NGHOSTS[0] * ( (j) + Nxx_plus_2NGHOSTS[1] * ( (k) + Nxx_plus_2NGHOSTS[2] * (g) ) ) )
#define IDX3(i,j,k) ( (i) + Nxx_plus_2NGHOSTS[0] * ( (j) + Nxx_plus_2NGHOSTS[1] * (k) ) )
// Assuming idx = IDX3(i,j,k). Much faster if idx can be reused over and over:
#define IDX4pt(g,idx)   ( (idx) + (Nxx_plus_2NGHOSTS[0]*Nxx_plus_2NGHOSTS[1]*Nxx_plus_2NGHOSTS[2]) * (g) )
#endif

// In the standalone, these are already defined by the boundary conditions module.
/*const int MAXFACE = -1;
*const int NUL     = +0;
*const int MINFACE = +1;
*/
// Declare boundary condition FACE_UPDATE function,
// which fills in the ghost zones with successively
// lower order finite differencing
void AtoB(const int ORDER, const int Nxx_plus_2NGHOSTS[3], const REAL *in_gfs, REAL *auxevol_gfs, const REAL dxx[3],
          const int i0min, const int i0max, 
          const int i1min, const int i1max, 
          const int i2min, const int i2max, 
          const int FACEX0, const int FACEX1, const int FACEX2) {
  
  const REAL invdx0 = 1.0 / dxx[0];
  const REAL invdx1 = 1.0 / dxx[1];
  const REAL invdx2 = 1.0 / dxx[2];

  if(ORDER==8) {
    //printf("Computing A to B with Order = 8...\n");
    for(int i2=i2min;i2<i2max;i2++) for(int i1=i1min;i1<i1max;i1++) for(int i0=i0min;i0<i0max;i0++) {
        #include "B_from_A_order8.h"
    }
  } else if(ORDER==6) {
    //printf("Computing A to B with Order = 6...\n");
    for(int i2=i2min;i2<i2max;i2++) for(int i1=i1min;i1<i1max;i1++) for(int i0=i0min;i0<i0max;i0++) {
        #include "B_from_A_order6.h"
    }
  } else if(ORDER==4) {
    //printf("Computing A to B with Order = 4...\n");
    for(int i2=i2min;i2<i2max;i2++) for(int i1=i1min;i1<i1max;i1++) for(int i0=i0min;i0<i0max;i0++) {
        #include "B_from_A_order4.h"
    }
  } else if(ORDER==2) {
    //printf("Computing A to B with Order = 2...\n");
    for(int i2=i2min;i2<i2max;i2++) for(int i1=i1min;i1<i1max;i1++) for(int i0=i0min;i0<i0max;i0++) {
        #include "B_from_A_order2.h"
    } 
  } else if(ORDER==0) {
    if(FACEX0==MAXFACE) {
    //printf("Computing A to B at x = max...\n");
        for(int i2=i2min;i2<i2max;i2++) for(int i1=i1min;i1<i1max;i1++) for(int i0=i0min;i0<i0max;i0++) {
            #include "B_from_A_order2_dirx0_dnwind.h"
        }
    } else if(FACEX0==MINFACE) {
    //printf("Computing A to B at x = min...\n");
        for(int i2=i2min;i2<i2max;i2++) for(int i1=i1min;i1<i1max;i1++) for(int i0=i0min;i0<i0max;i0++) {
            #include "B_from_A_order2_dirx0_upwind.h"
        }
    } else if(FACEX1==MAXFACE) {
    //printf("Computing A to B at y = max...\n");
        for(int i2=i2min;i2<i2max;i2++) for(int i1=i1min;i1<i1max;i1++) for(int i0=i0min;i0<i0max;i0++) {
            #include "B_from_A_order2_dirx1_dnwind.h"
        }
    } else if(FACEX1==MINFACE) {
    //printf("Computing A to B at y = min...\n");
        for(int i2=i2min;i2<i2max;i2++) for(int i1=i1min;i1<i1max;i1++) for(int i0=i0min;i0<i0max;i0++) {
            #include "B_from_A_order2_dirx1_upwind.h"
        }
    } else if(FACEX2==MAXFACE) {
    //printf("Computing A to B at z = max...\n");
        for(int i2=i2min;i2<i2max;i2++) for(int i1=i1min;i1<i1max;i1++) for(int i0=i0min;i0<i0max;i0++) {
            #include "B_from_A_order2_dirx2_dnwind.h"
        }
    } else if(FACEX2==MINFACE) {
    //printf("Computing A to B at z = min...\n");
        for(int i2=i2min;i2<i2max;i2++) for(int i1=i1min;i1<i1max;i1++) for(int i0=i0min;i0<i0max;i0++) {
            #include "B_from_A_order2_dirx2_upwind.h"
        }
    } else {
        printf("ERROR. FACEX parameters not set properly.\n");
        exit(1);
    }
  } else {
    printf("ERROR. ORDER = %d not supported!\n",ORDER);
    exit(1);
  }
}
void driver_A_to_B(const int Nxx[3],const int Nxx_plus_2NGHOSTS[3], const REAL dxx[3],const REAL *in_gfs,REAL *auxevol_gfs) {
  
  int ORDER = NGHOSTS*2;
  for(int ii=0;ii<Nxx_plus_2NGHOSTS[2]*Nxx_plus_2NGHOSTS[1]*Nxx_plus_2NGHOSTS[0];ii++) {
      auxevol_gfs[IDX4pt(BU0GF,ii)] = 1.0 / 0.0;
      auxevol_gfs[IDX4pt(BU1GF,ii)] = 1.0 / 0.0;
      auxevol_gfs[IDX4pt(BU2GF,ii)] = 1.0 / 0.0;
  }

  //printf("Starting A to B driver with Order = %d...\n",ORDER);
  int imin[3] = { NGHOSTS, NGHOSTS, NGHOSTS };
  int imax[3] = { Nxx[0]+NGHOSTS, Nxx[1]+NGHOSTS, Nxx[2]+NGHOSTS };
  AtoB(ORDER, Nxx_plus_2NGHOSTS, in_gfs, auxevol_gfs,dxx, imin[0],imax[0],imin[1],imax[1],imin[2],imax[2], NUL,NUL,NUL);
  while(ORDER>0) {
      // After updating each face, adjust imin[] and imax[] 
      //   to reflect the newly-updated face extents.
      ORDER -= 2;
      AtoB(ORDER, Nxx_plus_2NGHOSTS, in_gfs, auxevol_gfs,dxx, imin[0]-1,imin[0], imin[1],imax[1], imin[2],imax[2], MINFACE,NUL,NUL); 
      if(ORDER!=0) imin[0]--;
      AtoB(ORDER, Nxx_plus_2NGHOSTS, in_gfs, auxevol_gfs,dxx, imax[0],imax[0]+1, imin[1],imax[1], imin[2],imax[2], MAXFACE,NUL,NUL); 
      if(ORDER!=0) imax[0]++;

      AtoB(ORDER, Nxx_plus_2NGHOSTS, in_gfs, auxevol_gfs,dxx, imin[0],imax[0], imin[1]-1,imin[1], imin[2],imax[2], NUL,MINFACE,NUL); 
      if(ORDER!=0) imin[1]--;
      AtoB(ORDER, Nxx_plus_2NGHOSTS, in_gfs, auxevol_gfs,dxx, imin[0],imax[0], imax[1],imax[1]+1, imin[2],imax[2], NUL,MAXFACE,NUL); 
      if(ORDER!=0) imax[1]++;

      AtoB(ORDER, Nxx_plus_2NGHOSTS, in_gfs, auxevol_gfs,dxx, imin[0],imax[0], imin[1],imax[1], imin[2]-1,imin[2], NUL,NUL,MINFACE); 
      if(ORDER!=0) imin[2]--;
      AtoB(ORDER, Nxx_plus_2NGHOSTS, in_gfs, auxevol_gfs,dxx, imin[0],imax[0], imin[1],imax[1], imax[2],imax[2]+1, NUL,NUL,MAXFACE); 
      if(ORDER!=0) imax[2]++;
    }
}
