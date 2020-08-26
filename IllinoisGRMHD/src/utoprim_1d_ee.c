#ifndef __HARM_UTOPRIM_1D_EE__C__
#define __HARM_UTOPRIM_1D_EE__C__
/***********************************************************************************
    Copyright 2006 Charles F. Gammie, Jonathan C. McKinney, Scott C. Noble, 
                   Gabor Toth, and Luca Del Zanna

                        HARM  version 1.0   (released May 1, 2006)

    This file is part of HARM.  HARM is a program that solves hyperbolic 
    partial differential equations in conservative form using high-resolution
    shock-capturing techniques.  This version of HARM has been configured to 
    solve the relativistic magnetohydrodynamic equations of motion on a 
    stationary black hole spacetime in Kerr-Schild coordinates to evolve
    an accretion disk model. 

    You are morally obligated to cite the following two papers in his/her 
    scientific literature that results from use of any part of HARM:

    [1] Gammie, C. F., McKinney, J. C., \& Toth, G.\ 2003, 
        Astrophysical Journal, 589, 444.

    [2] Noble, S. C., Gammie, C. F., McKinney, J. C., \& Del Zanna, L. \ 2006, 
        Astrophysical Journal, 641, 626.

   
    Further, we strongly encourage you to obtain the latest version of 
    HARM directly from our distribution website:
    http://rainman.astro.uiuc.edu/codelib/


    HARM is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    HARM is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with HARM; if not, write to the Free Software
    Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA

***********************************************************************************/

/**********************************************************************************
 * This is a modified version of the original HARM code which is compatible with
 * IllinoisGRMHD. The modifications below are documented in pedagogical Jupyter
 * notebooks, which are available at: www.nrpyplus.net
 *
 * Author: Leo Werneck (leow155@gmail.com)
 * Date  : July/August 2020
 **********************************************************************************/

/*************************************************************************************/
/*************************************************************************************/
/*************************************************************************************

utoprim_1d_ee.c: 
---------------

  -- uses eq. (27) of Noble  et al. or the "momentum equation" and ignores
        the energy equation (29) in order to use the additional EOS, which 
        is  

             P = Sc rho^(GAMMA-1) / gamma
  
    Uses the 1D_W method: 
       -- solves for one independent variable (W) via a 1D
          Newton-Raphson method 
       -- solves for rho using Newton-Raphson using the definition of W : 
          W = Dc ( Dc + GAMMA Sc rho^(GAMMA-1) / (GAMMA-1) ) / rho 

       -- can be used (in principle) with a general equation of state. 

  -- Currently returns with an error state (>0) if a negative rest-mass
      density or internal energy density is calculated.  You may want 
      to change this aspect of the code so that it still calculates the 
      velocity and so that you can floor the densities.  If you want to 
      change this aspect of the code please comment out the "return(retval)"
      statement after "retval = 5;" statement in Utoprim_new_body();

******************************************************************************/


#include "u2p_defs.h"

#define NEWT_DIM (1)

#define LTRACE 0

// Declarations: 
static int Utoprim_new_body(CCTK_REAL U[], struct of_geom *geom, CCTK_REAL prim[], CCTK_REAL *gamma);
static void func_W(CCTK_REAL x[], CCTK_REAL dx[], CCTK_REAL resid[], 
		   CCTK_REAL jac[][NEWT_DIM], CCTK_REAL *f, CCTK_REAL *df, int n);

static int general_newton_raphson( CCTK_REAL x[], int n, 
				   void (*funcd) (CCTK_REAL [], CCTK_REAL [], CCTK_REAL [], 
						  CCTK_REAL [][NEWT_DIM], CCTK_REAL *, 
						  CCTK_REAL *, int) );
static int gnr2( CCTK_REAL x[], int n, 
			    void (*funcd) (CCTK_REAL [], CCTK_REAL [], CCTK_REAL [], 
					   CCTK_REAL [][NEWT_DIM],CCTK_REAL *,CCTK_REAL *,int) );
static void func_rho(CCTK_REAL x[], CCTK_REAL dx[], CCTK_REAL resid[], 
		     CCTK_REAL jac[][NEWT_DIM], CCTK_REAL *f, CCTK_REAL *df, int n);

/**********************************************************************/
/******************************************************************

  Utoprim_1d_ee():

  -- Driver for new prim. var. solver.  The driver just translates
     between the two sets of definitions for U and P.  The user may 
     wish to alter the translation as they see fit.  

     It assumes that on input/output:


              /  rho u^t           \
         U =  |  T^t_t   + rho u^t |  sqrt(-det(g_{\mu\nu}))
              |  T^t_\mu           |  
              \   B^i              /


             /    rho        \
	 P = |    uu         |
             | \tilde{u}^i   |
             \   B^i         /


         S = sqrt(-det(g_{a b})) u^t P / rho^(GAMMA-1)

     ala HARM. 

   Arguments:
       U[NPR]    = conserved variables (current values on input/output);
         S       = entropy density  = sqrt(-det(g_{a b})) u^t P / rho^(GAMMA-1)
       gcov[NDIM][NDIM] = covariant form of the metric ;
       gcon[NDIM][NDIM] = contravariant form of the metric ;
       gdet             = sqrt( - determinant of the metric) ;
       prim[NPR] = primitive variables (guess on input, calculated values on 
                                        output if there are no problems);
  
   -- NOTE: for those using this routine for special relativistic MHD and are
            unfamiliar with metrics, merely set 
              gcov = gcon = diag(-1,1,1,1)  and gdet = 1.  ;

******************************************************************/

int Utoprim_1d_ee(eos_struct eos, CCTK_REAL U[NPR], CCTK_REAL gcov[NDIM][NDIM], CCTK_REAL gcon[NDIM][NDIM],
		  CCTK_REAL gdet, CCTK_REAL prim[NPR], long &n_iter, CCTK_REAL S_star, CCTK_REAL *gamma )
{

  if( U[0] <= 0. ) {
    ERROR_MESSAGE("Negative U[0] found!!  We encourage you to figure out how this weird thing happened!! \n"); 
    return(-100);
  }

  /* First update the primitive B-fields */
  CCTK_REAL inv_gdet = 1.0 / gdet;
  for(int i = BCON1; i <= BCON3; i++) prim[i] = U[i] * inv_gdet ;

  /* Set the geometry variables: */
  CCTK_REAL alpha      = 1.0/sqrt(-gcon[0][0]);
  CCTK_REAL geomfactor = alpha * inv_gdet;

  /* Transform the CONSERVED variables into the new system */
  CCTK_REAL U_tmp[NPR];
  /* The input CONSERVED entropy variable is
   *
   * S_star := sqrt(-g) * S * u0
   *
   * Using u0 = gamma/alpha, with gamma the Lorentz factor,
   * and sqrt(-g), we get
   *
   * alpha / sqrt(-g) * S_star = alpha / sqrt(-g) * ( sqrt(-g) * S * u0 )
   *                           = ( gamma/alpha ) * alpha * S
   *                           = gamma * S,
   * and thus
   * .---------------------------------.
   * | gamma * S = geomfactor * S_star |
   * .---------------------------------.
   */
  CCTK_REAL gamma_times_S = geomfactor * S_star; 
  U_tmp[RHO]              = geomfactor * U[RHO];
  U_tmp[UU]               = geomfactor * (U[UU] - U[RHO]);
  U_tmp[UTCON1]           = geomfactor * U[UTCON1];
  U_tmp[UTCON2]           = geomfactor * U[UTCON2];
  U_tmp[UTCON3]           = geomfactor * U[UTCON3];
  U_tmp[BCON1 ]           = geomfactor * U[BCON1] ;
  U_tmp[BCON2 ]           = geomfactor * U[BCON2] ;
  U_tmp[BCON3 ]           = geomfactor * U[BCON3] ;

  /* Transform the PRIMITIVE variables into the new system */
  CCTK_REAL prim_tmp[NPR];
  prim_tmp[RHO   ] = prim[RHO   ];
  prim_tmp[UU    ] = prim[UU    ];
  prim_tmp[UTCON1] = prim[UTCON1];
  prim_tmp[UTCON2] = prim[UTCON2];
  prim_tmp[UTCON3] = prim[UTCON3];

  prim_tmp[BCON1] =  U_tmp[BCON1];
  prim_tmp[BCON2] =  U_tmp[BCON2];
  prim_tmp[BCON3] =  U_tmp[BCON3];

  /* Perform the C2P */
  int ret = Utoprim_new_body(eos, U_tmp, gcov, gcon, gdet, prim_tmp, n_iter, gamma_times_S, gamma);

  /* Transform new primitive variables back if there was no problem : */ 
  if( ret == 0 ) {
    prim[RHO   ] = prim_tmp[RHO   ];
    prim[UU    ] = prim_tmp[UU    ];
    prim[UTCON1] = prim_tmp[UTCON1];
    prim[UTCON2] = prim_tmp[UTCON2];
    prim[UTCON3] = prim_tmp[UTCON3];
  }

  return( ret ) ;

}


/**********************************************************************/
/**********************************************************************************

  Utoprim_new_body():

     -- Attempt an inversion from U to prim using the initial guess prim.

     -- This is the main routine that calculates auxiliary quantities for the 
        Newton-Raphson routine. 

  -- assumes that 
             /  rho gamma        \
         U = |  alpha T^t_\mu    |
             \  alpha B^i        /



               /    rho        \
	prim = |    uu         |
               | \tilde{u}^i   |
               \  alpha B^i   /


return:  (i*100 + j)  where 
         i = 0 ->  Newton-Raphson solver either was not called (yet or not used) 
                   or returned successfully;
             1 ->  Newton-Raphson solver did not converge to a solution with the 
                   given tolerances;
             2 ->  Newton-Raphson procedure encountered a numerical divergence 
                   (occurrence of "nan" or "+/-inf" ;
	     
         j = 0 -> success 
             1 -> failure: some sort of failure in Newton-Raphson; 
             2 -> failure: utsq<0 w/ initial p[] guess;
	     3 -> failure: W<0 or W>W_TOO_BIG
             4 -> failure: v^2 > 1 
             5 -> failure: rho,uu <= 0 ;

**********************************************************************************/

static int Utoprim_new_body(eos_struct eos, CCTK_REAL U[NPR], CCTK_REAL gcov[NDIM][NDIM],
                            CCTK_REAL gcon[NDIM][NDIM], CCTK_REAL gdet,  CCTK_REAL prim[NPR], long &n_iter, CCTK_REAL gamma_times_S )
{


  CCTK_REAL x_1d[1],x_rho[1],rho_g;
  CCTK_REAL QdotB,Bcon[NDIM],Qcov[NDIM],Qcon[NDIM],ncov[NDIM],ncon[NDIM],Qtcon[NDIM];
  CCTK_REAL rho0,u,p,w,gammasq,gamma,gamma_sq,W_last,W,utsq;
  CCTK_REAL g_o_WBsq, QdB_o_W;
  int i,j, retval, i_increase,ntries ;

  // Assume ok initially:
  retval = 0;

  // Calculate various scalars (Q.B, Q^2, etc)  from the conserved variables:
  Bcon[0] = 0. ;
  Bcon[1] = U[BCON1];
  Bcon[2] = U[BCON2];
  Bcon[3] = U[BCON3];
//-faster   Bsq = Bcon[1]*Bcov[1] + Bcon[2]*Bcov[2] + Bcon[3]*Bcov[3] ;
//-faster   lower(Bcon,geom,Bcov);
  Bsq =   gcov[1][1]*Bcon[1]*Bcon[1]
        + gcov[2][2]*Bcon[2]*Bcon[2]
        + gcov[3][3]*Bcon[3]*Bcon[3]
     + 2*(gcov[1][2]*Bcon[1]*Bcon[2]
	+ gcov[1][3]*Bcon[1]*Bcon[3]
	+ gcov[2][3]*Bcon[2]*Bcon[3]) ;


  Qcov[0] = U[QCOV0] ;  
  Qcov[1] = U[QCOV1] ;  
  Qcov[2] = U[QCOV2] ;  
  Qcov[3] = U[QCOV3] ;  
  raise(Qcov,geom,Qcon) ;

  QdotB   = Qcov[1]*Bcon[1] + Qcov[2]*Bcon[2] + Qcov[3]*Bcon[3] ;
  QdotBsq = QdotB*QdotB ;

//-fast   ncov[0] = -geom->alpha;
//-fast   ncov[1] = ncov[2] = ncov[3] = 0.;
//-faster   alpha = geom->alpha; 
  
//-faster  ncon[0] = -alpha * geom->gcon[0][0]; 
//-faster  ncon[1] = -alpha * geom->gcon[0][1]; 
//-faster  ncon[2] = -alpha * geom->gcon[0][2]; 
//-faster  ncon[3] = -alpha * geom->gcon[0][3]; 
  CCTK_REAL alpha = 1.0/sqrt(-gcon[0][0]);
  CCTK_REAL Qdotn = -alpha * Qcon[0];

  CCTK_REAL Qsq = 0. ;
  for(i=0;i<4;i++) Qsq += Qcov[i]*Qcon[i];

  CCTK_REAL Qtsq = Qsq + Qdotn*Qdotn ;

  CCTK_REAL D = U[RHO] ;

  /* calculate W from last timestep and use for guess */
  utsq = 0. ;
  for(i=1;i<4;i++)
    for(j=1;j<4;j++) utsq += gcov[i][j]*prim[UTCON1+i-1]*prim[UTCON1+j-1] ;
  
  /* utsq =  geom->gcov[1][1]*prim[U1]*prim[U1] */
  /*       + geom->gcov[2][2]*prim[U2]*prim[U2] */
  /*       + geom->gcov[3][3]*prim[U3]*prim[U3] */
  /*    + 2*(geom->gcov[1][2]*prim[U1]*prim[U2] */
  /* 	+ geom->gcov[1][3]*prim[U1]*prim[U3] */
  /* 	+ geom->gcov[2][3]*prim[U2]*prim[U3]) ; */

  if( (utsq < 0.) && (fabs(utsq) < 1.0e-13) ) { 
    utsq = fabs(utsq);
  }
  if(utsq < 0. || utsq > UTSQ_TOO_BIG) {
    retval = 2;
    return(retval) ;
  }

  CCTK_REAL gammasq = 1. + utsq ;    // Lorentz factor squared
  CCTK_REAL gamma   = sqrt(gammasq); // Lorentz factor, not to be confused with Gamma	
  CCTK_REAL rho0    = D / gamma ;

  /* Set basic cold polytrope-based hybrid EOS quantities */
  CCTK_REAL poly_idx   = find_polytropic_K_and_Gamma_index(eos,rho0); // Which piece of the PPEOS to use
  CCTK_REAL Gamma_poly = eos.Gamma_ppoly_tab[poly_idx];               // Set Gamma
  CCTK_REAL Gm1        = Gamma_poly - 1.0;                            // HARM auxiliary variable
  CCTK_REAL Gm1_inv    = 1.0/Gm1;                                     // HARM auxiliary variable
  CCTK_REAL rho_Gm1    = pow(rho0,Gm1);                               // HARM auxiliary variable

  /* The definition of the entropy density, S, is
   *
   * S = P / rho^(Gamma - 1)
   *
   * Thus we have
   * .-------------------------.
   * | P = rho^(Gamma - 1) * S |
   * .-------------------------.
   */
  CCTK_REAL p = gamma_times_S * rho_Gm1 / gamma;
  CCTK_REAL u = p * Gm1_inv;
  CCTK_REAL w = rho0 + u + p ;

#if(LTRACE)
  if( rho0 <= 0. ) { 
    fprintf(stderr,"beg neg rho fix1, rho,D,gamma = %26.20e %26.20e %26.20e \n", rho0, D, gamma); fflush(stderr);
    rho0 = fabs(rho0);
  }
#endif 
    
  W_last = w*gammasq ;

  // Make sure that W is large enough so that v^2 < 1 : 
  i_increase = 0;
  while( (( W_last*W_last*W_last * ( W_last + 2.*Bsq ) 
	    - QdotBsq*(2.*W_last + Bsq) ) <= W_last*W_last*(Qtsq-Bsq*Bsq))
	 && (i_increase < 10) ) {
    W_last *= 10.;
    i_increase++;
  }
  
  W_for_gnr2 = W_for_gnr2_old = W_last;
  rho_for_gnr2 = rho_for_gnr2_old = rho0;

  // Calculate W: 
  x_1d[0] = W_last;
  
  retval = general_newton_raphson( x_1d, 1, func_W );

  W = x_1d[0];

  /* Problem with solver, so return denoting error before doing anything further */
  if( (retval != 0) || (W == FAIL_VAL) ) {
    retval = retval*100+1;
#if(LTRACE)
    fprintf(stderr,"fix1: retval, W, rho = %d %26.20e %26.20e \n", retval,W, rho_for_gnr2);fflush(stderr);
#endif
    return(retval);
  }
  else{
    if(W <= 0. || W > W_TOO_BIG) {
      retval = 3;
#if(LTRACE)
      fprintf(stderr,"fix1: retval, W, rho = %d %26.20e %26.20e \n", retval,W, rho_for_gnr2);fflush(stderr);
#endif
      return(retval) ;
    }
  }

  rho_g = x_rho[0] = rho_for_gnr2; 

  ntries = 0;
  while (  (retval = gnr2( x_rho, 1, func_rho )) &&  ( ntries++ < 10 )  ) { 
    rho_g *= 10.;
    x_rho[0] = rho_g;
  }

  rho_for_gnr2 = x_rho[0];

  if( (retval != 0) ) {
    retval = 10;
    return(retval);
  }

  // Calculate v^2 : 
  rho0       = rho_for_gnr2;
  poly_idx   = find_polytropic_K_and_Gamma_index(eos,rho0);
  Gamma_poly = eos.Gamma_ppoly_tab[poly_idx];
  Gm1        = Gamma_poly - 1.0;
  rho_Gm1    = pow(rho0,Gm1);

  utsq = (D-rho0)*(D+rho0)/(rho0*rho0);

  gamma_sq = 1.+utsq;
  gamma = sqrt(gamma_sq);
  *gamma_out = gamma;
  
  if( utsq < 0. ) {
    retval = 4;
#if(LTRACE)
    fprintf(stderr,"fix1: retval, W, rho, utsq = %d %26.20e %26.20e %26.20e \n", retval,W, rho_for_gnr2,utsq);fflush(stderr);
#endif
    return(retval) ;
  }

  // Recover the primitive variables from the scalars and conserved variables:
  
  w = W / gamma_sq;
  p = gamma_times_S * rho_Gm1 / gamma; 
  u = p * Gm1_inv;

  // User may want to handle this case differently, e.g. do NOT return upon 
  // a negative rho/u, calculate v^i so that rho/u can be floored by other routine:
  if( treat_floor_as_failure && ((rho0 <= 0.) || (u <= 0.)) ) { 
    retval = 5;
#if(LTRACE)
    fprintf(stderr,"fix1: retval, W, rho,utsq,u = %d %26.20e %26.20e %26.20e %26.20e \n", retval,W, rho0,utsq,u);fflush(stderr);
#endif
    return(retval) ;
  }

  prim[RHO] = rho0 ;
  prim[UU] = u ;

  g_o_WBsq = gamma/(W+Bsq);
  QdB_o_W  = QdotB / W; 
  prim[UTCON1] = g_o_WBsq * ( Qcon[1] + geom->ncon[1] * Qdotn + QdB_o_W*Bcon[1] ) ;
  prim[UTCON2] = g_o_WBsq * ( Qcon[2] + geom->ncon[2] * Qdotn + QdB_o_W*Bcon[2] ) ;
  prim[UTCON3] = g_o_WBsq * ( Qcon[3] + geom->ncon[3] * Qdotn + QdB_o_W*Bcon[3] ) ;


  /* done! */
  return(retval) ;

}


/**********************************************************************/
/************************************************************

  general_newton_raphson(): 

    -- performs Newton-Rapshon method on an arbitrary system.

    -- inspired in part by Num. Rec.'s routine newt();

    Arguements: 

       -- x[]   = set of independent variables to solve for;
       -- n     = number of independent variables and residuals;
       -- funcd = name of function that calculates residuals, etc.;

*****************************************************************/
static int general_newton_raphson( CCTK_REAL x[], int n, 
				   void (*funcd) (CCTK_REAL [], CCTK_REAL [], CCTK_REAL [], 
						  CCTK_REAL [][NEWT_DIM], CCTK_REAL *, 
						  CCTK_REAL *, int) )
{
  CCTK_REAL f, df, dx[NEWT_DIM], x_old[NEWT_DIM], resid[NEWT_DIM], 
    jac[NEWT_DIM][NEWT_DIM];
  CCTK_REAL errx, x_orig[NEWT_DIM];
  int    n_iter, id, jd, i_extra, doing_extra;
  CCTK_REAL dW,W,W_old;

  int   keep_iterating, i_increase;


  // Initialize various parameters and variables:
  errx = 1. ; 
  df = f = 1.;
  i_extra = doing_extra = 0;
  //-fast  for( id = 0; id < n ; id++)  x_old[id] = x_orig[id] = x[id] ;
  x_old[0] = x_orig[0] = x[0] ;

  W = W_old = 0.;

  n_iter = 0;


  /* Start the Newton-Raphson iterations : */
  keep_iterating = 1;
  while( keep_iterating ) { 

    (*funcd) (x, dx, resid, jac, &f, &df, n);  /* returns with new dx, f, df */

    /* Save old values before calculating the new: */
    errx = 0.;

    //-fast    for( id = 0; id < n ; id++) { x_old[id] = x[id] ;  }
    x_old[0] = x[0] ;

    /* don't use line search : */
    //-fast    for( id = 0; id < n ; id++) { x[id] += dx[id]  ;  }
    x[0] += dx[0]  ;

//    //METHOD specific:
//    i_increase = 0;
//    while( (( x[0]*x[0]*x[0] * ( x[0] + 2.*Bsq ) - 
//	      QdotBsq*(2.*x[0] + Bsq) ) <= x[0]*x[0]*(Qtsq-Bsq*Bsq))
//	   && (i_increase < 10) ) {
//      x[0] -= (1.*i_increase) * dx[0] / 10. ;
//      i_increase++;
//    }

    /****************************************/
    /* Calculate the convergence criterion */
    /****************************************/

    /* For the new criterion, always look at error in "W" : */
    // METHOD specific:
    errx  = (x[0]==0.) ?  fabs(dx[0]) : fabs(dx[0]/x[0]);


    /****************************************/
    /* Make sure that the new x[] is physical : */
    /****************************************/
    x[0] = fabs(x[0]);


    /*****************************************************************************/
    /* If we've reached the tolerance level, then just do a few extra iterations */
    /*  before stopping                                                          */
    /*****************************************************************************/
    
    if( (fabs(errx) <= NEWT_TOL) && (doing_extra == 0) && (EXTRA_NEWT_ITER > 0) ) {
      doing_extra = 1;
    }

    if( doing_extra == 1 ) i_extra++ ;

    if( ((fabs(errx) <= NEWT_TOL)&&(doing_extra == 0)) || 
	(i_extra > EXTRA_NEWT_ITER) || (n_iter >= (MAX_NEWT_ITER-1)) ) {
      keep_iterating = 0;
    }

    n_iter++;

  }   // END of while(keep_iterating)


  /*  Check for bad untrapped divergences : */
  if( (!isfinite(f)) || (!isfinite(df)) || (!isfinite(x[0]))  ) {
#if(LTRACE)
    fprintf(stderr,"\ngnr not finite, f,df,x_o,x,W_o,W,rho_o,rho = %26.20e %26.20e %26.20e %26.20e %26.20e %26.20e %26.20e %26.20e \n",
	    f,df,x[0],x_old[0],W_for_gnr2_old,W_for_gnr2,rho_for_gnr2_old,rho_for_gnr2); fflush(stderr); 
#endif
    return(2);
  }

  if( fabs(errx) <= NEWT_TOL ){
    return(0);
  }
  else if( (fabs(errx) <= MIN_NEWT_TOL) && (fabs(errx) > NEWT_TOL) ){
    return(0);
  }
  else {
    return(1);
  } 

  return(0);

}


/***********************************************************/
/********************************************************************** 

  gnr2()

    -- used to calculate rho from W

*****************************************************************/
static int gnr2( CCTK_REAL x[], int n, 
			    void (*funcd) (CCTK_REAL [], CCTK_REAL [], CCTK_REAL [], 
					 CCTK_REAL [][NEWT_DIM],CCTK_REAL *,CCTK_REAL *,int) )
{
  CCTK_REAL f, df, dx[NEWT_DIM], x_old[NEWT_DIM], resid[NEWT_DIM], 
    jac[NEWT_DIM][NEWT_DIM];
  CCTK_REAL errx, x_orig[NEWT_DIM];
  int    n_iter, id,jd, i_extra, doing_extra;
  CCTK_REAL dW,W,W_old;

  int   keep_iterating;


  // Initialize various parameters and variables:
  errx = 1. ; 
  df = f = 1.;
  i_extra = doing_extra = 0;
  //-fast  for( id = 0; id < n ; id++)  x_old[id] = x_orig[id] = x[id] ;
  x_old[0] = x_orig[0] = x[0] ;

  n_iter = 0;


  /* Start the Newton-Raphson iterations : */
  keep_iterating = 1;
  while( keep_iterating ) { 

    (*funcd) (x, dx, resid, jac, &f, &df, n);  /* returns with new dx, f, df */

    /* Save old values before calculating the new: */
    //-fast  errx = 0.;
    //-fast    for( id = 0; id < n ; id++) { x_old[id] = x[id] ;  }
    x_old[0] = x[0] ;  

    /* Make the newton step: */
    //-fast    for( id = 0; id < n ; id++) { x[id] += dx[id] ;  }
    x[0] += dx[0] ;

    /* Calculate the convergence criterion */
    //-fast    for( id = 0; id < n ; id++) { errx  += (x[id]==0.) ?  fabs(dx[id]) : fabs(dx[id]/x[id]); }
    //-fast    errx /= 1.*n;
    errx  = (x[0]==0.) ?  fabs(dx[0]) : fabs(dx[0]/x[0]); 

    /* Make sure that the new x[] is physical : */
    // METHOD specific:
    x[0] = fabs(x[0]);


    /* If we've reached the tolerance level, then just do a few extra iterations */
    /*   before stopping                                                         */
    if( (fabs(errx) <= NEWT_TOL2) && (doing_extra == 0) && (EXTRA_NEWT_ITER > 0) ) {
      doing_extra = 0;
    }

    if( doing_extra == 1 ) i_extra++ ;

    // See if we've done the extra iterations, or have done too many iterations:
    if( ((fabs(errx) <= NEWT_TOL2)&&(doing_extra == 0)) || 
	(i_extra > EXTRA_NEWT_ITER) || (n_iter >= (MAX_NEWT_ITER-1)) ) {
      keep_iterating = 0;
    }

    n_iter++;

  }  


  /*  Check for bad untrapped divergences : */
  if( (!isfinite(f)) || (!isfinite(df)) || (!isfinite(x[0]))  ) {
#if(LTRACE)
    fprintf(stderr,"\ngnr2 not finite, f,df,x_o,x,W_o,W,rho_o,rho = %26.20e %26.20e %26.20e %26.20e %26.20e %26.20e %26.20e %26.20e \n",
	    f,df,x[0],x_old[0],W_for_gnr2_old,W_for_gnr2,rho_for_gnr2_old,rho_for_gnr2); fflush(stderr); 
#endif
    return(2);
  }

  // Return in different ways depending on whether a solution was found:
  if( fabs(errx) <= NEWT_TOL ){
    return(0);
  }
  else if( (fabs(errx) <= MIN_NEWT_TOL) && (fabs(errx) > NEWT_TOL) ){
    return(0);
  }
  else {
    return(1);
  } 

  return(0);

}



/**********************************************************************/
/*********************************************************************************
   func_W()

        -- calculates the residuals, and Newton step for general_newton_raphson();
        -- for this method, x=W here;

     Arguments:
          x   = current value of independent var's (on input & output);
         dx   = Newton-Raphson step (on output);
        resid = residuals based on x (on output);
         jac  = Jacobian matrix based on x (on output);
         f    =  resid.resid/2  (on output)
        df    = -2*f;  (on output)
         n    = dimension of x[];
 *********************************************************************************/
static void func_W(CCTK_REAL x[], CCTK_REAL dx[], CCTK_REAL resid[], 
			 CCTK_REAL jac[][NEWT_DIM], CCTK_REAL *f, CCTK_REAL *df, int n)
{
  int retval, ntries;
  CCTK_REAL W, x_rho[1], rho, rho_g ;
  CCTK_REAL t2,t4,t7,t15,t24,t200,t300,t400,t1000;
  CCTK_REAL s100,s200;

  /* HEYHEY */
  t2 = D*D;
  t4 = QdotBsq*t2;
  t7 = Bsq*Bsq;
  t24 = 1/t2;
  two_Bsq = Bsq + Bsq;
  t300 = QdotBsq*Bsq*t2;
  t400 = Qtsq*t2;

  s200 = D*GAMMA*Sc;
  s100 = g_o_gm1*Sc;

  W  = x[0];
  W_for_gnr2_old = W_for_gnr2;
  W_for_gnr2 = W;

  // get rho from NR:
  rho_g = x_rho[0] = rho_for_gnr2;
  
  ntries = 0;
  while (  (retval = gnr2( x_rho, 1, func_rho)) &&  ( ntries++ < 10 )  ) { 
    rho_g *= 10.;
    x_rho[0] = rho_g;
  }

#if(LTRACE)
  if( x_rho[0] <= 0. ) { 
    fprintf(stderr,"gnr2 neg rho = %d ,rho_n,rho,rho_o,W,W_o = %26.20e %26.20e %26.20e %26.20e %26.20e \n", retval, x_rho[0], rho_for_gnr2, rho_for_gnr2_old, x[0], W_for_gnr2_old);
    fflush(stderr);
  }

  if( retval ) { 
    fprintf(stderr,"gnr2 retval = %d ,rho_n,rho,rho_o,W,W_o = %26.20e %26.20e %26.20e %26.20e %26.20e \n", retval, x_rho[0], rho_for_gnr2, rho_for_gnr2_old, x[0], W_for_gnr2_old);
    fflush(stderr);
  }
#endif 

  rho_for_gnr2_old = rho_for_gnr2; 
  rho = rho_for_gnr2 = x_rho[0];


  rho_Gm1 = pow(rho,gm1);
  drho_dW = -rho*rho/( -rho_Gm1*s200 + W*rho);

  //  t6 = rho*rho;
  //  t100 = (D-rho)*(D+rho);  // t2 - t6;
  t15 = -(D-rho)*(D+rho);  // t6-t2
  t200 = W + two_Bsq;
  t1000 = rho*drho_dW;
  resid[0] = (t300+(t4+t4+(t400+t15*(t7+(t200)*W))*W)*W)*t24;
  jac[0][0] = 2*(t4+(t400+t15*t7+(3.0*t15*Bsq+t7*t1000+(t15+t15+t1000*(t200))*W)*W)*W)*t24;

  dx[0] = -resid[0]/jac[0][0];

  *df = - resid[0]*resid[0];
  *f = -0.5*(*df);


//  fprintf(stdout,"QdotBsq = %28.18e ; \n",QdotBsq );
//  fprintf(stdout,"Sc      = %28.18e ; \n",Sc     );
//  fprintf(stdout,"Bsq     = %28.18e ; \n",Bsq     );
//  fprintf(stdout,"Qtsq    = %28.18e ; \n",Qtsq    );
//  fprintf(stdout,"Dc      = %28.18e ; \n",D       );
//  fprintf(stdout,"drhodW  = %28.18e ; \n",drho_dW  );
//  fprintf(stdout,"W       = %28.18e ; \n",W       );
//  fprintf(stdout,"rho     = %28.18e ; \n",rho     );
//  fprintf(stdout,"resid_W = %28.18e ; \n",resid[0] );
//  fprintf(stdout,"jac_W   = %28.18e ; \n",jac[0][0]);
//  fprintf(stdout,"deriv1 %g %g %g %g \n",W,resid[0],jac[0][0],dx[0]);

  return;

}


/**********************************************************************/
/*********************************************************************************
   func_rho():

        -- residual/jacobian routine to calculate rho from W via the polytrope:

        W  =  ( 1 + GAMMA * K_atm * rho^(GAMMA-1)/(GAMMA-1) ) D^2 / rho

     Arguments:
          x   = current value of independent var's (on input & output);
         dx   = Newton-Raphson step (on output);
        resid = residuals based on x (on output);
         jac  = Jacobian matrix based on x (on output);
         f    =  resid.resid/2  (on output)
        df    = -2*f;  (on output)
         n    = dimension of x[];
 *********************************************************************************/
// for the isentropic version:   eq.  (27)
static void func_rho(CCTK_REAL x[], CCTK_REAL dx[], CCTK_REAL resid[], 
			 CCTK_REAL jac[][NEWT_DIM], CCTK_REAL *f, CCTK_REAL *df, int n)
{

  CCTK_REAL A, B, C, rho, W, B0;
  CCTK_REAL t40,t14;
  
  rho = x[0];
  W = W_for_gnr2;

  rho_Gm1 = t40 = pow(rho,gm1);

  resid[0] = (rho*W+(-t40*s100-D)*D);
  t14 = t40/rho;  // rho^(g-2)
  jac[0][0] = -t14*s200 + W;
  //  drho_dW = -rho/jac[0][0];

  dx[0] = -resid[0]/jac[0][0];
  *df = - resid[0]*resid[0];
  *f = -0.5*(*df);

  //  fprintf(stdout,"deriv3 %g %g %g %g %g \n",rho,W,resid[0],jac[0][0],dx[0]);
//  fprintf(stdout,"Dc        := %28.18e ; \n",D);
//  fprintf(stdout,"GAMMA     := %28.18e ; \n",GAMMA);
//  fprintf(stdout,"Sc        := %28.18e ; \n",Sc);
//  fprintf(stdout,"rho       := %28.18e ; \n",rho);
//  fprintf(stdout,"W         := %28.18e ; \n",W);
//  fprintf(stdout,"resid_rho := %28.18e ; \n",resid[0] );
//  fprintf(stdout,"jac_rho   := %28.18e ; \n",jac[0][0] );

  return;

}

/****************************************************************************** 
             END   OF   UTOPRIM_1D_EE.C
 ******************************************************************************/
#endif
