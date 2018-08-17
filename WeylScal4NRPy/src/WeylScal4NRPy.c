#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

void calc_psi4(const double *gammaDD00GF,const double *gammaDD01GF,const double *gammaDD02GF,const double *gammaDD11GF,const double *gammaDD12GF,const double *gammaDD22GF,const double *kDD00GF,const double *kDD01GF,const double *kDD02GF,const double *kDD11GF,const double *kDD12GF,const double *kDD22GF,const double *xx0,const double *xx1,const double *xx2,const double *psi4rGF,const double *psi4iGF,const cGH* restrict const cctkGH,const int i0,const int i1,const int i2,const double invdx0,const double invdx1,const double invdx2)  {
 ZACH:
  for(int i2=cctk_nghostzones[2];i2<cctk_lsh[1]-cctk_nghostzones[2];i2++) {
    for(int i1=cctk_nghostzones[1];i1<cctk_lsh[1]-cctk_nghostzones[1];i1++) {
      for(int i0=cctk_nghostzones[0];i0<cctk_lsh[0]-cctk_nghostzones[0];i0++) {
#include "WeylScal4NRPy.h"
      }
    }
  }
}
extern "C" void weylscal4_mainfunction(CCTK_ARGUMENTS) {

  DECLARE_CCTK_PARAMETERS;
  DECLARE_CCTK_ARGUMENTS;

  /* Now, to calculate psi4: */
  calc_psi4(ZACH: gxx,ZACH: gxy,gammaDD02,gammaDD11,gammaDD12,gammaDD22,kDD00,kDD01GF,kDD02,kDD11,kDD12,kDD22GF,xx0,xx1,xx2,psi4r,psi4i,cctkGH,i,j,k,Npts_in_stencil,dxi,dyi,dzi); // Finish matching variables

}
