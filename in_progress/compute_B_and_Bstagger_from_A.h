// #include <cstdio>
// #include <cstdlib>
// #include <cmath>
// #include <sys/time.h>
// #include "GiRaFFE_headers.h"

#define LOOP_DEFINE_SIMPLE                      \
  _Pragma("omp parallel for")                   \
  for(int k=0;k<Nxx_plus_2NGHOSTS2;k++)                \
    for(int j=0;j<Nxx_plus_2NGHOSTS1;j++)              \
      for(int i=0;i<Nxx_plus_2NGHOSTS0;i++)

void GiRaFFE_compute_B_and_Bstagger_from_A(const paramstruct *params,
                                           const REAL *gxx, const REAL *gxy, const REAL *gxz, const REAL *gyy, const REAL *gyz,const REAL *gzz,
                                           REAL *psi3_bssn, const REAL *Ax, const REAL *Ay, const REAL *Az,
                                           REAL *Bx, REAL *By, REAL *Bz, REAL *Bx_stagger, REAL *By_stagger, REAL *Bz_stagger) {
#include "../set_Cparameters.h"

  const REAL dxi = invdx0;
  const REAL dyi = invdx1;
  const REAL dzi = invdx2;

  LOOP_DEFINE_SIMPLE {
    const int index=IDX3S(i,j,k);
    psi3_bssn[index] = sqrt(sqrt( gxx[index]*gyy[index]*gzz[index]
                               -  gxx[index]*gyz[index]*gyz[index]
                               +2*gxy[index]*gxz[index]*gyz[index]
                               -  gyy[index]*gxz[index]*gxz[index]
                               -  gzz[index]*gxy[index]*gxy[index]));
  }

  LOOP_DEFINE_SIMPLE {
    // Look Mom, no if() statements!
    const int shiftedim1 = (i-1)*(i!=0); // This way, i=0 yields shiftedim1=0 and shiftedi=1, used below for our COPY boundary condition.
    const int shiftedi   = shiftedim1+1;

    const int shiftedjm1 = (j-1)*(j!=0);
    const int shiftedj   = shiftedjm1+1;

    const int shiftedkm1 = (k-1)*(k!=0);
    const int shiftedk   = shiftedkm1+1;

    int index,indexim1,indexjm1,indexkm1;

    const int actual_index = IDX3S(i,j,k);

    const REAL Psim3 = 1.0/psi3_bssn[actual_index];

    // For the lower boundaries, the following applies a "copy"
    //    boundary condition on Bi_stagger where needed.
    //    E.g., Bx_stagger(i,jmin,k) = Bx_stagger(i,jmin+1,k)
    //    We find the copy BC works better than extrapolation.
    // For the upper boundaries, we do the following copy:
    //    E.g., Psi(imax+1,j,k)=Psi(imax,j,k)
    /**************/
    /* Bx_stagger */
    /**************/

    index    = IDX3S(i,shiftedj,shiftedk);
    indexjm1 = IDX3S(i,shiftedjm1,shiftedk);
    indexkm1 = IDX3S(i,shiftedj,shiftedkm1);
    // Set Bx_stagger = \partial_y A_z - partial_z A_y
    // "Grid" Ax(i,j,k) is actually Ax(i,j+1/2,k+1/2)
    // "Grid" Ay(i,j,k) is actually Ay(i+1/2,j,k+1/2)
    // "Grid" Az(i,j,k) is actually Ay(i+1/2,j+1/2,k)
    // Therefore, the 2nd order derivative \partial_z A_y at (i+1/2,j,k) is:
    //          ["Grid" Ay(i,j,k) - "Grid" Ay(i,j,k-1)]/dZ
    Bx_stagger[actual_index] = (Az[index]-Az[indexjm1])*dyi - (Ay[index]-Ay[indexkm1])*dzi;

    // Now multiply Bx and Bx_stagger by 1/sqrt(gamma(i+1/2,j,k)]) = 1/sqrt(1/2 [gamma + gamma_ip1]) = exp(-6 x 1/2 [phi + phi_ip1] )
    const int imax_minus_i = (Nxx_plus_2NGHOSTS0-1)-i;
    const int indexip1jk = IDX3S(i + ( (imax_minus_i > 0) - (0 > imax_minus_i) ),j,k);
    Bx_stagger[actual_index] *= Psim3/psi3_bssn[index];

    /**************/
    /* By_stagger */
    /**************/

    index    = IDX3S(shiftedi,j,shiftedk);
    indexim1 = IDX3S(shiftedim1,j,shiftedk);
    indexkm1 = IDX3S(shiftedi,j,shiftedkm1);
    // Set By_stagger = \partial_z A_x - \partial_x A_z
    By_stagger[actual_index] = (Ax[index]-Ax[indexkm1])*dzi - (Az[index]-Az[indexim1])*dxi;

    // Now multiply By and By_stagger by 1/sqrt(gamma(i,j+1/2,k)]) = 1/sqrt(1/2 [gamma + gamma_jp1]) = exp(-6 x 1/2 [phi + phi_jp1] )
    const int jmax_minus_j = (Nxx_plus_2NGHOSTS1-1)-j;
    const int indexijp1k = IDX3S(i,j + ( (jmax_minus_j > 0) - (0 > jmax_minus_j) ),k);
    By_stagger[actual_index] *= Psim3/psi3_bssn[index];

    /**************/
    /* Bz_stagger */
    /**************/

    index    = IDX3S(shiftedi,shiftedj,k);
    indexim1 = IDX3S(shiftedim1,shiftedj,k);
    indexjm1 = IDX3S(shiftedi,shiftedjm1,k);
    // Set Bz_stagger = \partial_x A_y - \partial_y A_x
    Bz_stagger[actual_index] = (Ay[index]-Ay[indexim1])*dxi - (Ax[index]-Ax[indexjm1])*dyi;

    // Now multiply Bz_stagger by 1/sqrt(gamma(i,j,k+1/2)]) = 1/sqrt(1/2 [gamma + gamma_kp1]) = exp(-6 x 1/2 [phi + phi_kp1] )
    const int kmax_minus_k = (Nxx_plus_2NGHOSTS2-1)-k;
    const int indexijkp1 = IDX3S(i,j,k + ( (kmax_minus_k > 0) - (0 > kmax_minus_k) ));
    Bz_stagger[actual_index] *= Psim3/psi3_bssn[index];

  }

  LOOP_DEFINE_SIMPLE {
    // Look Mom, no if() statements!
    const int shiftedim1 = (i-1)*(i!=0); // This way, i=0 yields shiftedim1=0 and shiftedi=1, used below for our COPY boundary condition.
    const int shiftedi   = shiftedim1+1;

    const int shiftedjm1 = (j-1)*(j!=0);
    const int shiftedj   = shiftedjm1+1;

    const int shiftedkm1 = (k-1)*(k!=0);
    const int shiftedk   = shiftedkm1+1;

    int index,indexim1,indexjm1,indexkm1;

    const int actual_index = IDX3S(i,j,k);

    // For the lower boundaries, the following applies a "copy"
    //    boundary condition on Bi and Bi_stagger where needed.
    //    E.g., Bx(imin,j,k) = Bx(imin+1,j,k)
    //    We find the copy BC works better than extrapolation.
    /******/
    /* Bx */
    /******/
    index = IDX3S(shiftedi,j,k);
    indexim1 = IDX3S(shiftedim1,j,k);
    // Set Bx = 0.5 ( Bx_stagger + Bx_stagger_im1 )
    // "Grid" Bx_stagger(i,j,k) is actually Bx_stagger(i+1/2,j,k)
    Bx[actual_index] = 0.5 * ( Bx_stagger[index] + Bx_stagger[indexim1] );

    /******/
    /* By */
    /******/
    index = IDX3S(i,shiftedj,k);
    indexjm1 = IDX3S(i,shiftedjm1,k);
    // Set By = 0.5 ( By_stagger + By_stagger_im1 )
    // "Grid" By_stagger(i,j,k) is actually By_stagger(i,j+1/2,k)
    By[actual_index] = 0.5 * ( By_stagger[index] + By_stagger[indexjm1] );

    /******/
    /* Bz */
    /******/
    index = IDX3S(i,j,shiftedk);
    indexkm1 = IDX3S(i,j,shiftedkm1);
    // Set Bz = 0.5 ( Bz_stagger + Bz_stagger_im1 )
    // "Grid" Bz_stagger(i,j,k) is actually Bz_stagger(i,j+1/2,k)
    Bz[actual_index] = 0.5 * ( Bz_stagger[index] + Bz_stagger[indexkm1] );
  }
}
