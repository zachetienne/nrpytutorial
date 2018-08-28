#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

CCTK_REAL re(const CCTK_REAL x) {
  return x;
}

CCTK_REAL im(const CCTK_REAL x) {
  return x;
}

void calc_psis(const cGH* restrict const cctkGH,const int *cctk_lsh,const int *cctk_nghostzones,
	       const CCTK_REAL *x,const CCTK_REAL *y,const CCTK_REAL *z,
	       const CCTK_REAL invdx0,const CCTK_REAL invdx1,const CCTK_REAL invdx2,
	       const CCTK_REAL *gammaDD00GF,const CCTK_REAL *gammaDD01GF,const CCTK_REAL *gammaDD02GF,const CCTK_REAL *gammaDD11GF,const CCTK_REAL *gammaDD12GF,const CCTK_REAL *gammaDD22GF,
	       const CCTK_REAL     *kDD00GF,const CCTK_REAL     *kDD01GF,const CCTK_REAL     *kDD02GF,const CCTK_REAL     *kDD11GF,const CCTK_REAL     *kDD12GF,const CCTK_REAL     *kDD22GF,
	       CCTK_REAL *psi4rGF,CCTK_REAL *psi4iGF,
	       CCTK_REAL *psi3rGF,CCTK_REAL *psi3iGF,
	       CCTK_REAL *psi2rGF,CCTK_REAL *psi2iGF,
	       CCTK_REAL *psi1rGF,CCTK_REAL *psi1iGF,
	       CCTK_REAL *psi0rGF,CCTK_REAL *psi0iGF)  {

  DECLARE_CCTK_PARAMETERS;
#include "WeylScal4NRPy_psis.h"

/*#pragma omp parallel for
  //  for(int i2=cctk_nghostzones[2];i2<cctk_lsh[2]-cctk_nghostzones[2];i2++) {
  //  for(int i1=cctk_nghostzones[1];i1<cctk_lsh[1]-cctk_nghostzones[1];i1++) {
  //   for(int i0=cctk_nghostzones[0];i0<cctk_lsh[0]-cctk_nghostzones[0];i0++) {
  for(int i2=2;i2<cctk_lsh[2]-2;i2++) {
    for(int i1=2;i1<cctk_lsh[1]-2;i1++) {
      for(int i0=2;i0<cctk_lsh[0]-2;i0++) {
	const CCTK_REAL xx0 = x[CCTK_GFINDEX3D(cctkGH, i0,i1,i2)];
	const CCTK_REAL xx1 = y[CCTK_GFINDEX3D(cctkGH, i0,i1,i2)];
	const CCTK_REAL xx2 = z[CCTK_GFINDEX3D(cctkGH, i0,i1,i2)];
#include "WeylScal4NRPy_psis.h"
      }
    }
  }*/
}

void calc_invars(const cGH* restrict const cctkGH,const int *cctk_lsh,const int *cctk_nghostzones,
	         const CCTK_REAL *psi4rGF,const CCTK_REAL *psi4iGF,
	         const CCTK_REAL *psi3rGF,const CCTK_REAL *psi3iGF,
	         const CCTK_REAL *psi2rGF,const CCTK_REAL *psi2iGF,
	         const CCTK_REAL *psi1rGF,const CCTK_REAL *psi1iGF,
	         const CCTK_REAL *psi0rGF,const CCTK_REAL *psi0iGF,
		 CCTK_REAL *curvIrGF,CCTK_REAL *curvIiGF,
		 CCTK_REAL *curvJrGF,CCTK_REAL *curvJiGF,
		 CCTK_REAL *curvJ1GF,CCTK_REAL *curvJ2GF,
		 CCTK_REAL *curvJ3GF,CCTK_REAL *curvJ4GF)  {
  DECLARE_CCTK_PARAMETERS;

#include "WeylScal4NRPy_invars.h"
}

extern void weylscal4_mainfunction(CCTK_ARGUMENTS) {

  DECLARE_CCTK_PARAMETERS;
  DECLARE_CCTK_ARGUMENTS;

  if(cctk_iteration % WeylScal4NRPy_calc_every != 0) { return; }
  
  const CCTK_REAL invdx0 = 1.0 / (CCTK_DELTA_SPACE(0));
  const CCTK_REAL invdx1 = 1.0 / (CCTK_DELTA_SPACE(1));
  const CCTK_REAL invdx2 = 1.0 / (CCTK_DELTA_SPACE(2));

  /* Now, to calculate psi4: */
  calc_psis(cctkGH,cctk_lsh,cctk_nghostzones,
	    x,y,z,
	    invdx0,invdx1,invdx2,
	    gxx,gxy,gxz,gyy,gyz,gzz,
	    kxx,kxy,kxz,kyy,kyz,kzz,
	    psi4r,psi4i,
	    psi3r,psi3i,
	    psi2r,psi2i,
	    psi1r,psi1i,
	    psi0r,psi0i);

  if (CCTK_EQUALS(output_scalars, "all_psis_and_invariants")) {
    calc_invars(cctkGH,cctk_lsh,cctk_nghostzones,
      	        psi4r,psi4i,
	        psi3r,psi3i,
	        psi2r,psi2i,
	        psi1r,psi1i,
                psi0r,psi0i,
                NRPycurvIr,NRPycurvIi,
	        NRPycurvJr,NRPycurvJi,
	        NRPycurvJ1,NRPycurvJ2,
                NRPycurvJ3,NRPycurvJ4);
		}

}
