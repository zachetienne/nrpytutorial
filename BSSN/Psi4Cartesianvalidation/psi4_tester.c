#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <sys/time.h>

/* EVOLVED VARIABLES: */
#define NUM_AUX_GFS 15
#define KDD00GF 0
#define KDD01GF 1
#define KDD02GF 2
#define KDD11GF 3
#define KDD12GF 4
#define KDD22GF 5
#define GAMMADD00GF 6
#define GAMMADD01GF 7
#define GAMMADD02GF 8
#define GAMMADD11GF 9
#define GAMMADD12GF 10
#define GAMMADD22GF 11
#define XGF 12
#define YGF 13
#define ZGF 14

#define REAL double

#define NGHOSTS 2

// Cartesian coordinates parameters
const REAL xmin = -10.,xmax=10.;
const REAL ymin = -10.,ymax=10.;
const REAL zmin = -10.,zmax=10.;

#define IDX4(g,i,j,k) \
( (i) + Nxx_plus_2NGHOSTS[0] * ( (j) + Nxx_plus_2NGHOSTS[1] * ( (k) + Nxx_plus_2NGHOSTS[2] * (g) ) ) )
#define IDX3(i,j,k) ( (i) + Nxx_plus_2NGHOSTS[0] * ( (j) + Nxx_plus_2NGHOSTS[1] * (k) ) )
// Assuming idx = IDX3(i,j,k). Much faster if idx can be reused over and over:
#define IDX4pt(g,idx)   ( (idx) + (Nxx_plus_2NGHOSTS[0]*Nxx_plus_2NGHOSTS[1]*Nxx_plus_2NGHOSTS[2]) * (g) )

int main(int argc, const char *argv[]) {
  // Step 0a: Read command-line input, error out if nonconformant
  if((argc != 4) || atoi(argv[1]) < NGHOSTS || atoi(argv[2]) < NGHOSTS || atoi(argv[3]) < 2 /* FIXME; allow for axisymmetric sims */) {
    fprintf(stderr,"Error: Expected three command-line arguments: ./BrillLindquist_Playground Nx0 Nx1 Nx2,\n");
    fprintf(stderr,"where Nx[0,1,2] is the number of grid points in the 0, 1, and 2 directions.\n");
    fprintf(stderr,"Nx[] MUST BE larger than NGHOSTS (= %d)\n",NGHOSTS);
    exit(1);
  }
  // Step 0b: Set up numerical grid structure, first in space...
  const int Nxx[3] = { atoi(argv[1]), atoi(argv[2]), atoi(argv[3]) };
  if(Nxx[0]%2 != 0 || Nxx[1]%2 != 0 || Nxx[2]%2 != 0) {
    fprintf(stderr,"Error: Cannot guarantee a proper cell-centered grid if number of grid cells not set to even number.\n");
    fprintf(stderr,"       For example, in case of angular directions, proper symmetry zones will not exist.\n");
    exit(1);
  }
  const int Nxx_plus_2NGHOSTS[3] = { Nxx[0]+2*NGHOSTS, Nxx[1]+2*NGHOSTS, Nxx[2]+2*NGHOSTS };
  const int Nxx_plus_2NGHOSTS_tot = Nxx_plus_2NGHOSTS[0]*Nxx_plus_2NGHOSTS[1]*Nxx_plus_2NGHOSTS[2];

  REAL *aux_gfs_new  = (REAL *)malloc(sizeof(REAL) * NUM_AUX_GFS * Nxx_plus_2NGHOSTS_tot);
  REAL *aux_gfs_old  = (REAL *)malloc(sizeof(REAL) * NUM_AUX_GFS * Nxx_plus_2NGHOSTS_tot);

  // Step 0d: Set up space and time coordinates
  // Step 0d.i: Set \Delta x^i on uniform grids.
  REAL xxmin[3],xxmax[3];
  xxmin[0] = xmin;
  xxmin[1] = ymin;
  xxmin[2] = zmin;
  xxmax[0] = xmax;
  xxmax[1] = ymax;
  xxmax[2] = zmax;

  REAL dxx[3];
  for(int i=0;i<3;i++) dxx[i] = (xxmax[i] - xxmin[i]) / ((REAL)Nxx[i]);

  // Step 0d.ii: Set up uniform coordinate grids
  REAL *xx[3];
  for(int i=0;i<3;i++) {
    xx[i] = (REAL *)malloc(sizeof(REAL)*Nxx_plus_2NGHOSTS[i]);
    for(int j=0;j<Nxx_plus_2NGHOSTS[i];j++) {
      xx[i][j] = xxmin[i] + ((REAL)(j-NGHOSTS) + (1.0/2.0))*dxx[i]; // Cell-centered grid.
    }
  }
  struct timeval time; 
  gettimeofday(&time,NULL);
  
  // microsecond has 1 000 000
  // Assuming you did not need quite that accuracy
  // Also do not assume the system clock has that accuracy.
  srand48((time.tv_sec * 1000) + (time.tv_usec / 1000));
  /*
  time_t t;
  srand48((unsigned) time(&t));
  */
#define PERT_SIZE 1e-6
#define LOOP_REGION(i0min,i0max, i1min,i1max, i2min,i2max) \
  for(int i2=i2min;i2<i2max;i2++) for(int i1=i1min;i1<i1max;i1++) for(int i0=i0min;i0<i0max;i0++)
  
  LOOP_REGION(0,Nxx_plus_2NGHOSTS[0], 0,Nxx_plus_2NGHOSTS[1], 0,Nxx_plus_2NGHOSTS[2]) {
    for(int gf=KDD00GF;gf<=KDD22GF;gf++) aux_gfs_new[IDX4(gf,i0,i1,i2)] = aux_gfs_old[IDX4(gf,i0,i1,i2)] = PERT_SIZE*(0.5-drand48());
    for(int gf=GAMMADD00GF;gf<=GAMMADD22GF;gf++) aux_gfs_new[IDX4(gf,i0,i1,i2)] = aux_gfs_old[IDX4(gf,i0,i1,i2)] = PERT_SIZE*(0.5-drand48());
    aux_gfs_new[IDX4(GAMMADD00GF,i0,i1,i2)] = aux_gfs_old[IDX4(GAMMADD00GF,i0,i1,i2)] = 1.0 + PERT_SIZE*(0.5-drand48());
    aux_gfs_new[IDX4(GAMMADD11GF,i0,i1,i2)] = aux_gfs_old[IDX4(GAMMADD11GF,i0,i1,i2)] = 1.0 + PERT_SIZE*(0.5-drand48());
    aux_gfs_new[IDX4(GAMMADD22GF,i0,i1,i2)] = aux_gfs_old[IDX4(GAMMADD22GF,i0,i1,i2)] = 1.0 + PERT_SIZE*(0.5-drand48());
    aux_gfs_new[IDX4(XGF,i0,i1,i2)] = aux_gfs_old[IDX4(XGF,i0,i1,i2)] = xx[0][i0];
    aux_gfs_new[IDX4(YGF,i0,i1,i2)] = aux_gfs_old[IDX4(YGF,i0,i1,i2)] = xx[1][i1];
    aux_gfs_new[IDX4(ZGF,i0,i1,i2)] = aux_gfs_old[IDX4(ZGF,i0,i1,i2)] = xx[2][i2];
  }
  
  REAL psi4_real_new,psi4_real_old;
  LOOP_REGION(NGHOSTS,Nxx_plus_2NGHOSTS[0]-NGHOSTS, NGHOSTS,Nxx_plus_2NGHOSTS[1]-NGHOSTS, NGHOSTS,Nxx_plus_2NGHOSTS[2]-NGHOSTS) {
    REAL *aux_gfs = aux_gfs_new;
    const REAL xx0 = xx[0][i0], xx1 = xx[1][i1], xx2 = xx[2][i2];
    const REAL invdx0 = 1.0/dxx[0];
    const REAL invdx1 = 1.0/dxx[1];
    const REAL invdx2 = 1.0/dxx[2];

    REAL psi4_real;
#include "Psi4_new.h"
    if(i0==Nxx_plus_2NGHOSTS[0]/2 && i1==Nxx_plus_2NGHOSTS[1]/2 && i2==Nxx_plus_2NGHOSTS[2]/2) {
      //printf("%.15e\n",psi4_real);
      psi4_real_new = psi4_real;
    }
    //printf("%d %d %d %.15e\n",i0,i1,i2,psi4_real);
  }

  LOOP_REGION(NGHOSTS,Nxx_plus_2NGHOSTS[0]-NGHOSTS, NGHOSTS,Nxx_plus_2NGHOSTS[1]-NGHOSTS, NGHOSTS,Nxx_plus_2NGHOSTS[2]-NGHOSTS) {
    REAL *aux_gfs = aux_gfs_old;
    const REAL xx0 = xx[0][i0], xx1 = xx[1][i1], xx2 = xx[2][i2];
    const REAL invdx0 = 1.0/dxx[0];
    const REAL invdx1 = 1.0/dxx[1];
    const REAL invdx2 = 1.0/dxx[2];

    REAL psi4_real;
#include "Psi4_old.h"
    if(i0==Nxx_plus_2NGHOSTS[0]/2 && i1==Nxx_plus_2NGHOSTS[1]/2 && i2==Nxx_plus_2NGHOSTS[2]/2) {
      psi4_real_old = psi4_real;
      //printf("%.15e\n",psi4_real);
    }
  }

  REAL relerror = fabs(psi4_real_new - psi4_real_old) / (0.5* ( fabs(psi4_real_new) + fabs(psi4_real_old) ));
  if(relerror < 1e-13) {
    printf("GOOD: %.15e : %.15e %.15e\n",relerror, psi4_real_new,psi4_real_old);
  } else {
    printf("BAD:  %.15e : %.15e %.15e %ld\n",relerror, psi4_real_new,psi4_real_old, (time.tv_sec * 1000) + (time.tv_usec / 1000));
  }
  
  for(int i=0;i<3;i++) free(xx[i]);
  free(aux_gfs_new);
  free(aux_gfs_old);
  
  return 0;
}
