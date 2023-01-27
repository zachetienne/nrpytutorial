// This C header file reads a Scalar Field solution from a data file
// and performs 1D interpolation of the solution to a desired radius.

// Author: Leonardo Werneck
//         werneck **at** if **dot** usp **dot** br
//
// Based on TOV_interp.h by Zachariah B. Etienne

#include "stdio.h"
#include "stdlib.h"
#include "math.h"
#include "string.h"

#define REAL double

//#define STANDALONE_UNIT_TEST

int count_num_lines_in_file(FILE *in1DSFID) {
  int numlines_in_file = 0;
  char * line = NULL;

  size_t len = 0;
  ssize_t read;
  while ((read = getline(&line, &len, in1DSFID)) != -1) {
    numlines_in_file++;
  }
  rewind(in1DSFID);

  free(line);
  return numlines_in_file;
}

int read_checkpoint_datafile(FILE *in1DSFID, REAL *in_gfs) {

  char * line = NULL;

  size_t len = 0;
  ssize_t read;

  int which_line = 0;

  while ((read = getline(&line, &len, in1DSFID)) != -1) {
    in_gfs[which_line] = atof(line);
    /* fprintf(stderr,"%d %.15e %.15e\n",which_line,atof(line),in_gfs[which_line]); */
    which_line++;
  }
  
  
  return 0;
}

int read_datafile__set_arrays(FILE *in1DSFID,
                              REAL *r_arr,
                              REAL *sf_arr,
                              REAL *gammaDD00_or_psi4_arr,
                              REAL *alpha_arr) {
  char * line = NULL;

  size_t len = 0;
  ssize_t read;

  int which_line = 0;
  while ((read = getline(&line, &len, in1DSFID)) != -1) {
    // Define the line delimiters (i.e., the stuff that goes between the data on a given
    //     line of data.  Here, we define both spaces " " and tabs "\t" as data delimiters.
    const char delimiters[] = " \t";

    //Now we define "token", a pointer to the first column of data
    char *token;

    //Each successive time we call strtok(NULL,blah), we read in a new column of data from
    //     the originally defined character array, as pointed to by token.

    token                             = strtok(line, delimiters); if(token==NULL) { printf("BADDDD\n"); return 1; }
    r_arr[which_line]                 = strtod(token, NULL); token = strtok( NULL, delimiters );
    sf_arr[which_line]                = strtod(token, NULL); token = strtok( NULL, delimiters );
    gammaDD00_or_psi4_arr[which_line] = strtod(token, NULL); token = strtok( NULL, delimiters );
    alpha_arr[which_line]             = strtod(token, NULL); token = strtok( NULL, delimiters );

    which_line++;
  }
  free(line);
  return 0;
}

int bisection_idx_finder(const REAL r_star, const int numlines_in_file, const REAL *r_arr) {
  int x1 = 0;
  int x2 = numlines_in_file-1;
  REAL y1 = r_star-r_arr[x1];
  REAL y2 = r_star-r_arr[x2];
  if(y1*y2 >= 0) {
    fprintf(stderr,"%d %d | %lf %lf | %lf %lf | %lf \n", x1, x2, y1, y2, r_arr[x1], r_arr[x2], r_star);
    fprintf(stderr,"INTERPOLATION BRACKETING ERROR %e | %e %e \n",r_star,y1,y2);
    exit(1);
  }
  for(int i=0;i<numlines_in_file;i++) {
    int x_midpoint = (x1+x2)/2;
    REAL y_midpoint = r_star-r_arr[x_midpoint];
    if(y_midpoint*y1 < 0) {
      x2 = x_midpoint;
      y2 = y_midpoint;
    } else {
      x1 = x_midpoint;
      y1 = y_midpoint;
    }
    if( abs(x2-x1) == 1 ) {
      // If r_arr[x1] is closer to r_star than r_arr[x2] then return x1:
      if(fabs(r_star-r_arr[x1]) < fabs(r_star-r_arr[x2])){
          return x1;
      }
      // Otherwiser return x2:
      return x2;
    }
  }
  fprintf(stderr,"INTERPOLATION BRACKETING ERROR: DID NOT CONVERGE.\n");
  exit(1);
}

void scalarfield_interpolate_1D( REAL const  r_star_in,                       // Point of interest for interpolation
                                 int  const interp_stencil_size,              // Number of points used to interpolate
                                 int  const numlines_in_file,                 // Number of lines in initial data file
                                 REAL const *restrict r_arr,                  // Initial data array for r
                                 REAL const *restrict sf_arr,                 // Initial data array for varphi
                                 REAL const *restrict gammaDD00_or_psi4_arr,  // Initial data array for gamma_{rr}
                                 REAL const *restrict alpha_arr,              // Initial data array for alpha
                                 REAL       *restrict sf_star,                // "Output": varphi(t=0,r_star)
                                 REAL       *restrict gammaDD00_or_psi4_star, // "Output": gamma_{rr}(t=0,r_star) or psi^{4}(t=0,r_star)
                                 REAL       *restrict alpha_star) {           // "Output": alpha(t=0,r_star)
    
    
  // For this case, we know that for all functions, f(r) = f(-r)
  REAL r_star = r_star_in;
  if ( r_star < 0 ) r_star = -r_star;
    
  // First we find the central interpolation stencil index:
  int idx = bisection_idx_finder(r_star,numlines_in_file,r_arr);

#ifdef MAX
#undef MAX
#endif
#define MAX(A, B) ( ((A) > (B)) ? (A) : (B) )

  int idxmin = MAX(0,idx-interp_stencil_size/2-1);

  // Now perform the Lagrange polynomial interpolation:

  // First set the interpolation coefficients:
  REAL r_sample[interp_stencil_size];
  for(int i=idxmin;i<idxmin+interp_stencil_size;i++) {
    r_sample[i-idxmin] = r_arr[i];
  }
  REAL l_i_of_r[interp_stencil_size];
  for(int i=0;i<interp_stencil_size;i++) {
    REAL numer = 1.0;
    REAL denom = 1.0;
    for(int j=0;j<i;j++) {
      numer *= r_star - r_sample[j];
      denom *= r_sample[i] - r_sample[j];
    }
    for(int j=i+1;j<interp_stencil_size;j++) {
      numer *= r_star - r_sample[j];
      denom *= r_sample[i] - r_sample[j];
    }
    l_i_of_r[i] = numer/denom;
  }

  // Then perform the interpolation:
  *sf_star                = 0.0;
  *gammaDD00_or_psi4_star = 0.0;
  *alpha_star             = 0.0;

  for(int i=idxmin;i<idxmin+interp_stencil_size;i++) {
    *sf_star                += l_i_of_r[i-idxmin] * sf_arr[i];
    *gammaDD00_or_psi4_star += l_i_of_r[i-idxmin] * gammaDD00_or_psi4_arr[i];
    *alpha_star             += l_i_of_r[i-idxmin] * alpha_arr[i];
  }
 
}
