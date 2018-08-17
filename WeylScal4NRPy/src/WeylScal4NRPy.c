#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

void calc_psi4(const double *gammaDD00GF,const double *gammaDD01GF,const double *gammaDD02GF,const double *gammaDD11GF,const double *gammaDD12GF,const double *gammaDD22GF,const double *kDD00GF,const double *kDD01GF,const double *kDD02GF,const double *kDD11GF,const double *kDD12GF,const double *kDD22GF,const double *xx0,const double *xx1,const double *xx2,const double *psi4rGF,const double *psi4iGF,const cGH* restrict const cctkGH,const int i0,const int i1,const int i2,const double invdx0,const double invdx1,const double invdx2)  {
#include "WeylScal4NRPy.h"
}

namespace WeylScal4NRPy {

extern "C" void WeylScal4_psi4_calc_Nth_SelectBCs(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
  
  if (cctk_iteration % WeylScal4NRPy_calc_every != WeylScal4NRPy_calc_offset)
    return;
  CCTK_INT ierr CCTK_ATTRIBUTE_UNUSED = 0;
  ierr = Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, GetBoundaryWidth(cctkGH), -1 /* no table */, "WeylScal4NRPy::Psi4i_group","flat");
  if (ierr < 0)
    CCTK_WARN(CCTK_WARN_ALERT, "Failed to register flat BC for WeylScal4NRPy::Psi4i_group.");
  ierr = Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, GetBoundaryWidth(cctkGH), -1 /* no table */, "WeylScal4NRPy::Psi4r_group","flat");
  if (ierr < 0)
    CCTK_WARN(CCTK_WARN_ALERT, "Failed to register flat BC for WeylScal4NRPy::Psi4r_group.");
  return;
}

static void WeylScal4NRPy_Body(const cGH* restrict const cctkGH, const int dir, const int face, const CCTK_REAL normal[3], const CCTK_REAL tangentA[3], const CCTK_REAL tangentB[3], const int imin[3], const int imax[3], const int n_subblock_gfs, CCTK_REAL* restrict const subblock_gfs[])
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  /* Initialise finite differencing variables */
  const ptrdiff_t di CCTK_ATTRIBUTE_UNUSED = 1;
  const ptrdiff_t dj CCTK_ATTRIBUTE_UNUSED = 
    CCTK_GFINDEX3D(cctkGH,0,1,0) - CCTK_GFINDEX3D(cctkGH,0,0,0);
  const ptrdiff_t dk CCTK_ATTRIBUTE_UNUSED = 
    CCTK_GFINDEX3D(cctkGH,0,0,1) - CCTK_GFINDEX3D(cctkGH,0,0,0);
  const ptrdiff_t cdi CCTK_ATTRIBUTE_UNUSED = sizeof(CCTK_REAL) * di;
  const ptrdiff_t cdj CCTK_ATTRIBUTE_UNUSED = sizeof(CCTK_REAL) * dj;
  const ptrdiff_t cdk CCTK_ATTRIBUTE_UNUSED = sizeof(CCTK_REAL) * dk;
  const ptrdiff_t cctkLbnd1 CCTK_ATTRIBUTE_UNUSED = cctk_lbnd[0];
  const ptrdiff_t cctkLbnd2 CCTK_ATTRIBUTE_UNUSED = cctk_lbnd[1];
  const ptrdiff_t cctkLbnd3 CCTK_ATTRIBUTE_UNUSED = cctk_lbnd[2];
  const CCTK_REAL_VEC t CCTK_ATTRIBUTE_UNUSED = ToReal(cctk_time);
  const CCTK_REAL_VEC cctkOriginSpace1 CCTK_ATTRIBUTE_UNUSED = 
    ToReal(CCTK_ORIGIN_SPACE(0));
  const CCTK_REAL_VEC cctkOriginSpace2 CCTK_ATTRIBUTE_UNUSED = 
    ToReal(CCTK_ORIGIN_SPACE(1));
  const CCTK_REAL_VEC cctkOriginSpace3 CCTK_ATTRIBUTE_UNUSED = 
    ToReal(CCTK_ORIGIN_SPACE(2));
  const CCTK_REAL_VEC dt CCTK_ATTRIBUTE_UNUSED = 
    ToReal(CCTK_DELTA_TIME);
  const CCTK_REAL_VEC dx CCTK_ATTRIBUTE_UNUSED = 
    ToReal(CCTK_DELTA_SPACE(0));
  const CCTK_REAL_VEC dy CCTK_ATTRIBUTE_UNUSED = 
    ToReal(CCTK_DELTA_SPACE(1));
  const CCTK_REAL_VEC dz CCTK_ATTRIBUTE_UNUSED = 
    ToReal(CCTK_DELTA_SPACE(2));
  const CCTK_REAL_VEC dxi CCTK_ATTRIBUTE_UNUSED = kdiv(ToReal(1),dx);
  const CCTK_REAL_VEC dyi CCTK_ATTRIBUTE_UNUSED = kdiv(ToReal(1),dy);
  const CCTK_REAL_VEC dzi CCTK_ATTRIBUTE_UNUSED = kdiv(ToReal(1),dz);

  /* Calculate temporaries and arrays functions */
  /* Copy local copies back to grid functions */
  /* Loop over the grid points */
  const int imin0=imin[0];
  const int imin1=imin[1];
  const int imin2=imin[2];
  const int imax0=imax[0];
  const int imax1=imax[1];
  const int imax2=imax[2];
  #pragma omp parallel
  CCTK_LOOP3STR(WeylScal4NRPy,
    i,j,k, imin0,imin1,imin2, imax0,imax1,imax2,
    cctk_ash[0],cctk_ash[1],cctk_ash[2],
    vecimin,vecimax, CCTK_REAL_VEC_SIZE)
  {
    const ptrdiff_t index CCTK_ATTRIBUTE_UNUSED = di*i + dj*j + dk*k;
    /* Assign local copies of grid functions */
    
    CCTK_REAL_VEC gammaDD00 CCTK_ATTRIBUTE_UNUSED = vec_load(gxx[index]); // gf = 1
    CCTK_REAL_VEC gammaDD01 CCTK_ATTRIBUTE_UNUSED = vec_load(gxy[index]); // gf = 2
    CCTK_REAL_VEC gammaDD02 CCTK_ATTRIBUTE_UNUSED = vec_load(gxz[index]); // gf = 3
    CCTK_REAL_VEC gammaDD11 CCTK_ATTRIBUTE_UNUSED = vec_load(gyy[index]); // gf = 4
    CCTK_REAL_VEC gammaDD12 CCTK_ATTRIBUTE_UNUSED = vec_load(gyz[index]); // gf = 5
    CCTK_REAL_VEC gammaDD22 CCTK_ATTRIBUTE_UNUSED = vec_load(gzz[index]); // gf = 6
    CCTK_REAL_VEC kDD00 CCTK_ATTRIBUTE_UNUSED = vec_load(kxx[index]); // gf = 7
    CCTK_REAL_VEC kDD01 CCTK_ATTRIBUTE_UNUSED = vec_load(kxy[index]); // gf = 8
    CCTK_REAL_VEC kDD02 CCTK_ATTRIBUTE_UNUSED = vec_load(kxz[index]); // gf = 9
    CCTK_REAL_VEC kDD11 CCTK_ATTRIBUTE_UNUSED = vec_load(kyy[index]); // gf = 10
    CCTK_REAL_VEC kDD12 CCTK_ATTRIBUTE_UNUSED = vec_load(kyz[index]); // gf = 11
    CCTK_REAL_VEC kDD22 CCTK_ATTRIBUTE_UNUSED = vec_load(kzz[index]); // gf = 12
    CCTK_REAL_VEC xx0 CCTK_ATTRIBUTE_UNUSED = vec_load(x[index]); // gf = 13
    CCTK_REAL_VEC xx1 CCTK_ATTRIBUTE_UNUSED = vec_load(y[index]); // gf = 14
    CCTK_REAL_VEC xx2 CCTK_ATTRIBUTE_UNUSED = vec_load(z[index]); // gf = 15
    CCTK_REAL_VEC psi4r CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL_VEC psi4i CCTK_ATTRIBUTE_UNUSED;
    
    /* Now, to calculate psi4: */
    calc_psi4(gammaDD00,gammaDD01,gammaDD02,gammaDD11,gammaDD12,gammaDD22,kDD00,kDD01GF,kDD02,kDD11,kDD12,kDD22GF,xx0,xx1,xx2,psi4r,psi4i,cctkGH,i,j,k,Npts_in_stencil,dxi,dyi,dzi); // Finish matching variables

    // Write psi4 to memory here.
  }
  CCTK_ENDLOOP3STR(WeylScal4NRPy);
}
extern "C" void WeylScal4NRPy(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
  
  if (verbose > 1)
  {
    CCTK_VInfo(CCTK_THORNSTRING,"Entering WeylScal4NRPy_Body");
  }
  if (cctk_iteration % WeylScal4NRPy_every != WeylScal4NRPy_calc_offset)
  {
    return;
  }
  
  const char* const groups[] = {
    "admbase::curv",
    "admbase::metric",
    "grid::coordinates",
    "WeylScal4::Psi4i_group",
    "WeylScal4::Psi4r_group"};
  AssertGroupStorage(cctkGH, "WeylScal4NRPy", 5, groups);
  
  EnsureStencilFits(cctkGH, "WeylScal4NRPy", 2, 2, 2);
  
  LoopOverInterior(cctkGH, WeylScal4NRPy_Body);
  if (verbose > 1)
  {
    CCTK_VInfo(CCTK_THORNSTRING,"Leaving WeylScal4NRPy_Body");
  }
}

} // namespace WeylScal4NRPy
