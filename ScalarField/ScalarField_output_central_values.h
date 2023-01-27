// This simple function computes the central values of the scalar field,
// the lapse function, and energy-density. If the lapse falls below
// a certain threshold, the function returns 1, indicating that the
// lapse has collapsed, otherwise it returns 0.
//
// Author: Leonardo R. Werneck
//         wernecklr **at** gmail **dot** com

int output_central_values( const REAL t, const paramstruct *restrict params, REAL *restrict in_gfs ) {
#include "set_Cparameters.h"
      
  /* Set indices */
  const int i0 = NGHOSTS;
  const int i1 = NGHOSTS;
  const int i2 = NGHOSTS;
  /* Set needed values of the scalar field */
  const REAL sf_i0p1     = in_gfs[IDX4S(SFGF, i0,i1,i2)];
  const REAL sf_i0p2     = in_gfs[IDX4S(SFGF, i0+1,i1,i2)];
  const REAL sf_i0p3     = in_gfs[IDX4S(SFGF, i0+2,i1,i2)];
  /* Set needed values of alpha */
  const REAL alpha_i0p1  = in_gfs[IDX4S(ALPHAGF, i0,i1,i2)];
  const REAL alpha_i0p2  = in_gfs[IDX4S(ALPHAGF, i0+1,i1,i2)];
  const REAL alpha_i0p3  = in_gfs[IDX4S(ALPHAGF, i0+2,i1,i2)];
  /* Compute the central values of the scalar field, alpha, and rho */
  const REAL sf_c    = 3.0*sf_i0p1    - 3.0*sf_i0p2    + sf_i0p3;
  const REAL alpha_c = 3.0*alpha_i0p1 - 3.0*alpha_i0p2 + alpha_i0p3;

  /* Set the output file */
  FILE *outfile;
  if( t > 0.0 ) {
    outfile = fopen("out_central_values.dat","a");
  }
  else {
    outfile = fopen("out_central_values.dat","w");
  }

  /* Output the central values of the scalar field and alpha */
  fprintf(outfile,"%.15e %.15e %.15e\n",t,alpha_c,sf_c);

  /* Close the file */
  fclose(outfile);

  if( alpha_c < alpha_threshold )
    return 1;
  else
    return 0;

}
