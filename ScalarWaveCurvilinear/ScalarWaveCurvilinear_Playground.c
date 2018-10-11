// Part P0: Set the number of ghost cells, from NRPy+'s FD_CENTDERIVS_ORDER
#define NGHOSTS 2

// Part P1: Import needed header files
#include "stdio.h"
#include "stdlib.h"
#include "math.h"

// Part P2: Add needed #define's to set data type, the IDX4() macro, and the gridfunctions
// Part P2a: set REAL=double, so that all gridfunctions are set to 
#define REAL double
// Part P2b: Declare the IDX4(gf,i,j,k) macro, which enables us to store 4-dimensions of
//           data in a 1D array. In this case, consecutive values of "i" 
//           (all other indices held to a fixed value) are consecutive in memory, where 
//           consecutive values of "j" (fixing all other indices) are separated by 
//           Nxx_plus_2NGHOSTS[0] elements in memory. Similarly, consecutive values of
//           "k" are separated by Nxx_plus_2NGHOSTS[0]*Nxx_plus_2NGHOSTS[1] in memory, etc.
#define IDX4(g,i,j,k) \
( (i) + Nxx_plus_2NGHOSTS[0] * ( (j) + Nxx_plus_2NGHOSTS[1] * ( (k) + Nxx_plus_2NGHOSTS[2] * (g) ) ) )
#define IDX3(i,j,k) ( (i) + Nxx_plus_2NGHOSTS[0] * ( (j) + Nxx_plus_2NGHOSTS[1] * (k) ) )
// Part P2c: Set UUGF and VVGF macros
#define NUM_GFS 2
#define UUGF 0
#define VVGF 1

// Step P3: Set free parameters for the initial data
const REAL wavespeed = 1.0;
const REAL kk0 = 1.0;
const REAL kk1 = 1.0;
const REAL kk2 = 1.0;

typedef struct ghostzone_map {
  short i0,i1,i2;
} gz_map;

#define LOOP_REGION(i0min,i0max, i1min,i1max, i2min,i2max) \
  for(int i2=i2min;i2<i2max;i2++) for(int i1=i1min;i1<i1max;i1++) for(int i0=i0min;i0<i0max;i0++)

// Part P4: Declare the function for the exact solution. time==0 corresponds to the initial data.
void exact_solution(const int Nxx_plus_2NGHOSTS[3],REAL time,REAL *xx[3], REAL *in_gfs) {
#pragma omp parallel for
  LOOP_REGION(0,Nxx_plus_2NGHOSTS[0], 1,Nxx_plus_2NGHOSTS[1], 2,Nxx_plus_2NGHOSTS[2]) {
    REAL xx_to_Cart0,xx_to_Cart1,xx_to_Cart2;
    REAL xx0 = xx[0][i0];
    REAL xx1 = xx[1][i1];
    REAL xx2 = xx[2][i2];
#include "xxCart.h"

    xx0 = xx_to_Cart0;
    xx1 = xx_to_Cart1;
    xx2 = xx_to_Cart2;
#include "ScalarWaveCartesian_ExactSolution.h"
  }
}

// Part P5: Declare the function to evaluate the scalar wave RHSs
void rhs_eval(const int Nxx[3],const int Nxx_plus_2NGHOSTS[3],const REAL dxx[3], REAL *xx[3], const REAL *in_gfs,REAL *rhs_gfs) {
#include "ScalarWaveCurvilinear_RHSs.h"
}

// Part P6: Declare boundary condition OB_UPDATE macro,
//          which updates a single face of the 3D grid cube
//          using quadratic polynomial extrapolation.
const int MAXFACE = -1;
const int NUL     = +0;
const int MINFACE = +1;

#define OB_UPDATE(inner,which_gf, bc_gz_map, i0min,i0max, i1min,i1max, i2min,i2max, FACEX0,FACEX1,FACEX2) \
  LOOP_REGION(i0min,i0max, i1min,i1max, i2min,i2max) {                  \
    const int idx3 = IDX3(i0,i1,i2);                                    \
    if(bc_gz_map[idx3].i0 == -1 && inner==0) {                          \
      gfs[IDX4(which_gf,i0,i1,i2)] =                                    \
        +3.0*gfs[IDX4(which_gf,i0+1*FACEX0,i1+1*FACEX1,i2+1*FACEX2)]    \
        -3.0*gfs[IDX4(which_gf,i0+2*FACEX0,i1+2*FACEX1,i2+2*FACEX2)]    \
        +1.0*gfs[IDX4(which_gf,i0+3*FACEX0,i1+3*FACEX1,i2+3*FACEX2)];   \
    } else if(bc_gz_map[idx3].i0 != -1 && inner==1) {                   \
      gfs[IDX4(which_gf,i0,i1,i2)] =                                    \
        gfs[IDX4(which_gf,                                              \
                 bc_gz_map[idx3].i0,                                    \
                 bc_gz_map[idx3].i1,                                    \
                 bc_gz_map[idx3].i2)];                                  \
    }                                                                   \
  }
//      printf("%e %d %d %d | %d %d %d\n",gfs[IDX4(which_gf,i0,i1,i2)],i0,i1,i2, \
  //           bc_gz_map[idx3].i0,                                        \
    //         bc_gz_map[idx3].i1,                                        \
      //       bc_gz_map[idx3].i2);                                       \

// Part P7: Boundary condition driver routine: Apply BCs to all six
//          boundary faces of the cube, filling in the innermost
//          ghost zone first, and moving outward.
void apply_bcs(const int Nxx[3],const int Nxx_plus_2NGHOSTS[3],gz_map *bc_gz_map,REAL *gfs) {
#pragma omp parallel for
    for(int which_gf=0;which_gf<NUM_GFS;which_gf++) {
    int imin[3] = { NGHOSTS, NGHOSTS, NGHOSTS };
    int imax[3] = { Nxx_plus_2NGHOSTS[0]-NGHOSTS, Nxx_plus_2NGHOSTS[1]-NGHOSTS, Nxx_plus_2NGHOSTS[2]-NGHOSTS };
    for(int which_gz = 0; which_gz < NGHOSTS; which_gz++) {
      for(int inner=0;inner<2;inner++) {
        // After updating each face, adjust imin[] and imax[] 
        //   to reflect the newly-updated face extents.
        //if(which_gz==2 && inner==0) exit(1);
        OB_UPDATE(inner,which_gf, bc_gz_map, imin[0]-1,imin[0], imin[1],imax[1], imin[2],imax[2], MINFACE,NUL,NUL); imin[0]--;
        OB_UPDATE(inner,which_gf, bc_gz_map, imax[0],imax[0]+1, imin[1],imax[1], imin[2],imax[2], MAXFACE,NUL,NUL); imax[0]++;

        OB_UPDATE(inner,which_gf, bc_gz_map, imin[0],imax[0], imin[1]-1,imin[1], imin[2],imax[2], NUL,MINFACE,NUL); imin[1]--;
        OB_UPDATE(inner,which_gf, bc_gz_map, imin[0],imax[0], imax[1],imax[1]+1, imin[2],imax[2], NUL,MAXFACE,NUL); imax[1]++;

        OB_UPDATE(inner,which_gf, bc_gz_map, imin[0],imax[0], imin[1],imax[1], imin[2]-1,imin[2], NUL,NUL,MINFACE); imin[2]--;
        OB_UPDATE(inner,which_gf, bc_gz_map, imin[0],imax[0], imin[1],imax[1], imax[2],imax[2]+1, NUL,NUL,MAXFACE); imax[2]++;
        if(inner==0) { for(int ii=0;ii<3;ii++) {imin[ii]++; imax[ii]--;} }
        //if(which_gz==1 && inner==1) {printf("hewwo %d\n",NGHOSTS);}
      }
    }
  }
}

// main() function:
// Step 0: Read command-line input, set up grid structure, allocate memory for gridfunctions, set up coordinates
// Step 1: Set up scalar wave initial data
// Step 2: Evolve scalar wave initial data forward in time using Method of Lines with RK4 algorithm,
//         applying quadratic extrapolation outer boundary conditions.
// Step 3: Output relative error between numerical and exact solution.
// Step 4: Free all allocated memory
int main(int argc, const char *argv[]) {

  // Step 0a: Read command-line input, error out if nonconformant
  if(argc != 4 || atoi(argv[1]) < NGHOSTS) {
      printf("Error: Expected one command-line argument: ./ScalarWaveCurvilinear_Playground Nx0 Nx1 Nx2,\n");
      printf("where Nx[0,1,2] is the number of grid points in the 0, 1, and 2 directions.\n");
      printf("Nx[] MUST BE larger than NGHOSTS (= %d)\n",NGHOSTS);
      exit(1);
  }
  // Step 0b: Set up numerical grid structure, first in space...
  const int Nx0 = atoi(argv[1]);
  const int Nx1 = atoi(argv[2]);
  const int Nx2 = atoi(argv[3]);
  if(Nx0%2 != 0 || Nx1%2 != 0 || Nx2%2 != 0) {
    printf("Error: Cannot guarantee a proper cell-centered grid if number of grid cells not set to even number.\n");
    printf("       For example, in case of angular directions, proper symmetry zones will not exist.\n");
    exit(1);
  }
  const int Nxx[3] = { Nx0, Nx1, Nx2 };
  const int Nxx_plus_2NGHOSTS[3] = { Nxx[0]+2*NGHOSTS, Nxx[1]+2*NGHOSTS, Nxx[2]+2*NGHOSTS };
  const int Nxx_plus_2NGHOSTS_tot = Nxx_plus_2NGHOSTS[0]*Nxx_plus_2NGHOSTS[1]*Nxx_plus_2NGHOSTS[2];
  const REAL RMAX = 10.0;
#include "xxminmax.h"

  //          ... and then set up the numerical grid structure in time:
  const int t_final = xxmax[0]*0.1; // Final time is set so that at t=t_final, 
                                    // data at the origin have not been corrupted 
                                    // by the approximate outer boundary condition
  //const REAL CFL_FACTOR = 0.015; // Set the CFL Factor
  const REAL CFL_FACTOR = 0.015; // Set the CFL Factor

  // Step 0c: Allocate memory for gridfunctions
  REAL *evol_gfs    = (REAL *)malloc(sizeof(REAL) * NUM_GFS * Nxx_plus_2NGHOSTS_tot);
  REAL *next_in_gfs = (REAL *)malloc(sizeof(REAL) * NUM_GFS * Nxx_plus_2NGHOSTS_tot);
  REAL *k1_gfs   = (REAL *)malloc(sizeof(REAL) * NUM_GFS * Nxx_plus_2NGHOSTS_tot);
  REAL *k2_gfs   = (REAL *)malloc(sizeof(REAL) * NUM_GFS * Nxx_plus_2NGHOSTS_tot);
  REAL *k3_gfs   = (REAL *)malloc(sizeof(REAL) * NUM_GFS * Nxx_plus_2NGHOSTS_tot);
  REAL *k4_gfs   = (REAL *)malloc(sizeof(REAL) * NUM_GFS * Nxx_plus_2NGHOSTS_tot);
  
  // Step 0d: Set up coordinates: Set dx, and then dt based on dx_min and CFL condition
#define MIN(A, B) ( ((A) < (B)) ? (A) : (B) )
  // xx[0][i] = xxmin[0] + (i-NGHOSTS)*dxx[0]
  REAL dxx[3];
  for(int i=0;i<3;i++) dxx[i] = (xxmax[i] - xxmin[i]) / ((REAL)Nxx[i]);
  REAL dt = CFL_FACTOR * MIN(dxx[0],MIN(dxx[1],dxx[2])); // CFL condition
  REAL Nt = t_final / dt + 1; // The number of points in time

  // Step 0e: Set up Cartesian coordinate grids
  REAL *xx[3];
  for(int i=0;i<3;i++) {
    xx[i] = (REAL *)malloc(sizeof(REAL)*Nxx_plus_2NGHOSTS[i]);
    for(int j=0;j<Nxx_plus_2NGHOSTS[i];j++) {
      xx[i][j] = xxmin[i] + ((REAL)(j-NGHOSTS) + (1.0/2.0))*dxx[i]; // Cell-centered grid.
    }
  }

  // Step 0f: Find ghostzone mappings:
  //const int num_gz_tot = Nxx_plus_2NGHOSTS_tot - Nxx[0]*Nxx[1]*Nxx[2];
  gz_map *bc_gz_map = (gz_map *)malloc(sizeof(gz_map)*Nxx_plus_2NGHOSTS_tot);
  LOOP_REGION(0,Nxx_plus_2NGHOSTS[0],0,Nxx_plus_2NGHOSTS[1],0,Nxx_plus_2NGHOSTS[2]) {
    REAL xx_to_Cart0,xx_to_Cart1,xx_to_Cart2;
    REAL xx0 = xx[0][i0];
    REAL xx1 = xx[1][i1];
    REAL xx2 = xx[2][i2];
#include "xxCart.h"
    REAL Cartx = xx_to_Cart0;
    REAL Carty = xx_to_Cart1;
    REAL Cartz = xx_to_Cart2;
    
    REAL Cart_to_xx0_inbounds,Cart_to_xx1_inbounds,Cart_to_xx2_inbounds;
#include "Cart_to_xx.h"
    int i0_inbounds = (int)( (Cart_to_xx0_inbounds - xxmin[0] - (1.0/2.0)*dxx[0] + ((REAL)NGHOSTS)*dxx[0])/dxx[0] + 0.5 ); 
    int i1_inbounds = (int)( (Cart_to_xx1_inbounds - xxmin[1] - (1.0/2.0)*dxx[1] + ((REAL)NGHOSTS)*dxx[1])/dxx[1] + 0.5 );
    int i2_inbounds = (int)( (Cart_to_xx2_inbounds - xxmin[2] - (1.0/2.0)*dxx[2] + ((REAL)NGHOSTS)*dxx[2])/dxx[2] + 0.5 );

    REAL xx_to_Cart0_orig = xx_to_Cart0;
    REAL xx_to_Cart1_orig = xx_to_Cart1;
    REAL xx_to_Cart2_orig = xx_to_Cart2;
    xx0 = xx[0][i0_inbounds];
    xx1 = xx[1][i1_inbounds];
    xx2 = xx[2][i2_inbounds];
#include "xxCart.h"
           
#define EPS_ABS 1e-8
    if(fabs( (xx_to_Cart0_orig - xx_to_Cart0) ) > EPS_ABS ||
       fabs( (xx_to_Cart1_orig - xx_to_Cart1) ) > EPS_ABS ||
       fabs( (xx_to_Cart2_orig - xx_to_Cart2) ) > EPS_ABS) {
             
      REAL r = sqrt(xx_to_Cart0*xx_to_Cart0 + xx_to_Cart1*xx_to_Cart1 + xx_to_Cart2*xx_to_Cart2);
      REAL th = acos(xx_to_Cart2/r);
      REAL ph = atan2(xx_to_Cart1,xx_to_Cart0);

      REAL rorig = sqrt(xx_to_Cart0_orig*xx_to_Cart0_orig + xx_to_Cart1_orig*xx_to_Cart1_orig + xx_to_Cart2_orig*xx_to_Cart2_orig);
      REAL thorig = acos(xx_to_Cart2_orig/rorig);
      REAL phorig = atan2(xx_to_Cart1_orig,xx_to_Cart0_orig);
      printf("Error. Cartesian disagreement: ( %.15e %.15e %.15e ) != ( %.15e %.15e %.15e )\n",xx_to_Cart0_orig,xx_to_Cart1_orig,xx_to_Cart2_orig,
             xx_to_Cart0,xx_to_Cart1,xx_to_Cart2);
      printf("Error. Cartesian disagreement2: %e %e %e ( %.15e %.15e %.15e ) != ( %.15e %.15e %.15e)\n",xx0,xx1,xx2,rorig,thorig,phorig,r,th,ph);
      exit(1);
    }

    if(i0_inbounds-i0 == 0 && i1_inbounds-i1 == 0 && i2_inbounds-i2 == 0) {
      bc_gz_map[IDX3(i0,i1,i2)].i0=-1;
      bc_gz_map[IDX3(i0,i1,i2)].i1=-1;
      bc_gz_map[IDX3(i0,i1,i2)].i2=-1;
    } else {
      bc_gz_map[IDX3(i0,i1,i2)].i0=i0_inbounds;
      bc_gz_map[IDX3(i0,i1,i2)].i1=i1_inbounds;
      bc_gz_map[IDX3(i0,i1,i2)].i2=i2_inbounds;
    }
  }
  
  // Step 1: Set up initial data to be exact solution at time=0:
  exact_solution(Nxx_plus_2NGHOSTS, 0.0, xx, evol_gfs);

  for(int n=0;n<=Nt;n++) { // Main loop to progress forward in time.
    // Step 2: Evolve scalar wave initial data forward in time using Method of Lines with RK4 algorithm,
    //         applying quadratic extrapolation outer boundary conditions.
    /***************************************************/
    /* Implement RK4 for Method of Lines timestepping: */
    /***************************************************/
    /* -= RK4: Step 1 of 4 =- */
    /* First evaluate k1 = RHSs expression             */
    rhs_eval(Nxx,Nxx_plus_2NGHOSTS,dxx, xx,evol_gfs, k1_gfs);
    /* Next k1 -> k1*dt, and then set the input for    */
    /*    the next RHS eval call to y_n+k1/2           */
    for(int i=0;i<Nxx_plus_2NGHOSTS_tot*NUM_GFS;i++) {
      k1_gfs[i] *= dt;
      next_in_gfs[i] = evol_gfs[i] + k1_gfs[i]*0.5;
    }
    /* Finally, apply boundary conditions to           */
    /* next_in_gfs, so its data are set everywhere.    */
    apply_bcs(Nxx,Nxx_plus_2NGHOSTS,bc_gz_map,next_in_gfs);

    /* -= RK4: Step 2 of 4 =- */
    rhs_eval(Nxx,Nxx_plus_2NGHOSTS,dxx, xx,next_in_gfs, k2_gfs);
    for(int i=0;i<Nxx_plus_2NGHOSTS_tot*NUM_GFS;i++) {
      k2_gfs[i] *= dt;
      next_in_gfs[i] = evol_gfs[i] + k2_gfs[i]*0.5;
    }
    apply_bcs(Nxx,Nxx_plus_2NGHOSTS,bc_gz_map,next_in_gfs);

    /* -= RK4: Step 3 of 4 =- */
    rhs_eval(Nxx,Nxx_plus_2NGHOSTS,dxx, xx,next_in_gfs, k3_gfs);
    for(int i=0;i<Nxx_plus_2NGHOSTS_tot*NUM_GFS;i++) {
      k3_gfs[i] *= dt;
      next_in_gfs[i] = evol_gfs[i] + k3_gfs[i];
    }
    apply_bcs(Nxx,Nxx_plus_2NGHOSTS,bc_gz_map,next_in_gfs);

    /* -= RK4: Step 4 of 4 =- */
    rhs_eval(Nxx,Nxx_plus_2NGHOSTS,dxx, xx,next_in_gfs, k4_gfs);
    for(int i=0;i<Nxx_plus_2NGHOSTS_tot*NUM_GFS;i++) {
      k4_gfs[i] *= dt;
      evol_gfs[i] += (1.0/6.0)*(k1_gfs[i] + 2.0*k2_gfs[i] + 2.0*k3_gfs[i] + k4_gfs[i]);
    }
    apply_bcs(Nxx,Nxx_plus_2NGHOSTS,bc_gz_map,evol_gfs);

    /* Step 3: Output relative error between numerical and exact solution, */
    const int i0mid=Nxx_plus_2NGHOSTS[0]/2;
    const int i1mid=Nxx_plus_2NGHOSTS[1]/2;
    const int i2mid=Nxx_plus_2NGHOSTS[2]/2;
    exact_solution(Nxx_plus_2NGHOSTS,(n+1)*dt, xx, k1_gfs);
    const REAL exact     = k1_gfs[IDX4(0,i0mid,i1mid,i2mid)];
    const REAL numerical = evol_gfs[IDX4(0,i0mid,i1mid,i2mid)];
    const REAL relative_error = fabs((exact-numerical)/exact);
    printf("%e %e || %e %e %e: %e %e\n",(n+1)*dt, log10(relative_error),
           xx[0][i0mid],xx[1][i1mid],xx[2][i2mid], numerical,exact);
  } // End main loop to progress forward in time.
  
  // Step 4: Free all allocated memory
  free(k4_gfs);
  free(k3_gfs);
  free(k2_gfs);
  free(k1_gfs);
  free(next_in_gfs);
  free(evol_gfs);
  for(int i=0;i<3;i++) free(xx[i]);
  return 0;
}
