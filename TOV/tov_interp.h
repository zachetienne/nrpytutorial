// This C header file reads a TOV solution from data file and performs
//    1D interpolation of the solution to a desired radius.

// Author: Zachariah B. Etienne
//         zachetie **at** gmail **dot* com

#include "stdio.h"
#include "stdlib.h"
#include "math.h"
#include "string.h"

#define REAL double

//#define STANDALONE_UNIT_TEST

int count_num_lines_in_file(FILE *in1Dpolytrope) {
  int numlines_in_file = 0;
  char * line = NULL;

  size_t len = 0;
  ssize_t read;
  while ((read = getline(&line, &len, in1Dpolytrope)) != -1) {
    numlines_in_file++;
  }
  rewind(in1Dpolytrope);

  free(line);
  return numlines_in_file;
}

int read_datafile__set_arrays(FILE *in1Dpolytrope, REAL *r_Schw_arr,REAL *rho_arr,REAL *P_arr,REAL *M_arr,REAL *expnu_arr,REAL *exp4phi_arr,REAL *rbar_arr) {
  char * line = NULL;

  size_t len = 0;
  ssize_t read;

  int which_line = 0;
  while ((read = getline(&line, &len, in1Dpolytrope)) != -1) {
    // Define the line delimiters (i.e., the stuff that goes between the data on a given
    //     line of data.  Here, we define both spaces " " and tabs "\t" as data delimiters.
    const char delimiters[] = " \t";

    //Now we define "token", a pointer to the first column of data
    char *token;

    //Each successive time we call strtok(NULL,blah), we read in a new column of data from
    //     the originally defined character array, as pointed to by token.

    token=strtok(line, delimiters); if(token==NULL) { printf("BADDDD\n"); return 1; }
    r_Schw_arr[which_line]  = strtod(token, NULL); token = strtok( NULL, delimiters );
    rho_arr[which_line]     = strtod(token, NULL); token = strtok( NULL, delimiters );
    P_arr[which_line]       = strtod(token, NULL); token = strtok( NULL, delimiters );
    M_arr[which_line]       = strtod(token, NULL); token = strtok( NULL, delimiters );
    expnu_arr[which_line]   = strtod(token, NULL); token = strtok( NULL, delimiters );
    exp4phi_arr[which_line] = strtod(token, NULL); token = strtok( NULL, delimiters );
    rbar_arr[which_line]    = strtod(token, NULL);

    which_line++;
  }
  free(line);
  return 0;
}


void TOV_interpolate_1D(REAL rrbar,const REAL Rbar,const int Rbar_idx,const int interp_stencil_size,
                        const int numlines_in_file,const REAL *r_Schw_arr,const REAL *rho_arr,const REAL *P_arr,const REAL *M_arr,const REAL *expnu_arr,const REAL *exp4phi_arr,const REAL *rbar_arr,
                        REAL *rho,REAL *P,REAL *M,REAL *expnu,REAL *exp4phi) {
  // Find interpolation index using Bisection root-finding algorithm:
  int bisection_idx_finder(const REAL rr, const int numlines_in_file, const REAL *rbar_arr) {
    int x1 = 0;
    int x2 = numlines_in_file-1;
    REAL y1 = rrbar-rbar_arr[x1];
    REAL y2 = rrbar-rbar_arr[x2];
    if(y1*y2 >= 0) {
      fprintf(stderr,"INTERPOLATION BRACKETING ERROR %e | %e %e\n",rr,y1,y2);
      exit(1);
    }
    for(int i=0;i<numlines_in_file;i++) {
      int x_midpoint = (x1+x2)/2;
      REAL y_midpoint = rrbar-rbar_arr[x_midpoint];
      if(y_midpoint*y1 < 0) {
        x2 = x_midpoint;
        y2 = y_midpoint;
      } else {
        x1 = x_midpoint;
        y1 = y_midpoint;
      }
      if( abs(x2-x1) == 1 ) {
        // If rbar_arr[x1] is closer to rrbar than rbar_arr[x2] then return x1:
        if(fabs(rrbar-rbar_arr[x1]) < fabs(rrbar-rbar_arr[x2])) return x1;
        // Otherwiser return x2:
        return x2;
      }
    }
    fprintf(stderr,"INTERPOLATION BRACKETING ERROR: DID NOT CONVERGE.\n");
    exit(1);
  }

  // For this case, we know that for all functions, f(r) = f(-r)
  if(rrbar < 0) rrbar = -rrbar;

  // First find the central interpolation stencil index:
  int idx = bisection_idx_finder(rrbar,numlines_in_file,rbar_arr);

  
#ifdef MAX
#undef MAX
#endif
#define MAX(A, B) ( ((A) > (B)) ? (A) : (B) )

  int idxmin = MAX(0,idx-interp_stencil_size/2-1);

#ifdef MIN
#undef MIN
#endif
#define MIN(A, B) ( ((A) < (B)) ? (A) : (B) )

  // -= Do not allow the interpolation stencil to cross the star's surface =-
  // max index is when idxmin + (interp_stencil_size-1) = Rbar_idx
  //  -> idxmin at most can be Rbar_idx - interp_stencil_size + 1
  if(rrbar < Rbar) {
    idxmin = MIN(idxmin,Rbar_idx - interp_stencil_size + 1);
  } else {
    idxmin = MAX(idxmin,Rbar_idx+1);
    idxmin = MIN(idxmin,numlines_in_file - interp_stencil_size + 1);
  }
  // Now perform the Lagrange polynomial interpolation:

  // First set the interpolation coefficients:
  REAL rbar_sample[interp_stencil_size];
  for(int i=idxmin;i<idxmin+interp_stencil_size;i++) {
    rbar_sample[i-idxmin] = rbar_arr[i];
  }
  REAL l_i_of_r[interp_stencil_size];
  for(int i=0;i<interp_stencil_size;i++) {
    REAL numer = 1.0;
    REAL denom = 1.0;
    for(int j=0;j<i;j++) {
      numer *= rrbar - rbar_sample[j];
      denom *= rbar_sample[i] - rbar_sample[j];
    }
    for(int j=i+1;j<interp_stencil_size;j++) {
      numer *= rrbar - rbar_sample[j];
      denom *= rbar_sample[i] - rbar_sample[j];
    }
    l_i_of_r[i] = numer/denom;
  }

  // Then perform the interpolation:
  *rho = 0.0;
  *P = 0.0;
  *M = 0.0;
  *expnu = 0.0;
  *exp4phi = 0.0;

  REAL r_Schw = 0.0;
  for(int i=idxmin;i<idxmin+interp_stencil_size;i++) {
    r_Schw   += l_i_of_r[i-idxmin] * r_Schw_arr[i];
    *rho     += l_i_of_r[i-idxmin] * rho_arr[i];
    *P       += l_i_of_r[i-idxmin] * P_arr[i];
    *M       += l_i_of_r[i-idxmin] * M_arr[i];
    *expnu   += l_i_of_r[i-idxmin] * expnu_arr[i];
    *exp4phi += l_i_of_r[i-idxmin] * exp4phi_arr[i];
  }

  if(rrbar > Rbar) {
    *rho   = 0;
    *P     = 0;
    *M     = M_arr[Rbar_idx+1];
    *expnu = 1. - 2.*(*M) / r_Schw;
    *exp4phi = pow(r_Schw / rrbar,2.0);
  }
}

#ifdef STANDALONE_UNIT_TEST
int main() {

  // Open the data file:
  char filename[100];
  sprintf(filename,"../outputTOVpolytrope.txt");
  FILE *in1Dpolytrope = fopen(filename, "r");
  if (in1Dpolytrope == NULL) {
    fprintf(stderr,"ERROR: could not open file %s\n",filename);
    exit(1);
  }

  // Count the number of lines in the data file:
  int numlines_in_file = count_num_lines_in_file(in1Dpolytrope);

  // Allocate space for all data arrays:
  REAL *r_Schw_arr  = (REAL *)malloc(sizeof(REAL)*numlines_in_file);
  REAL *rho_arr     = (REAL *)malloc(sizeof(REAL)*numlines_in_file);
  REAL *P_arr       = (REAL *)malloc(sizeof(REAL)*numlines_in_file);
  REAL *M_arr       = (REAL *)malloc(sizeof(REAL)*numlines_in_file);
  REAL *expnu_arr   = (REAL *)malloc(sizeof(REAL)*numlines_in_file);
  REAL *exp4phi_arr = (REAL *)malloc(sizeof(REAL)*numlines_in_file);
  REAL *rbar_arr    = (REAL *)malloc(sizeof(REAL)*numlines_in_file);

  // Read from the data file, filling in arrays
  if(read_datafile__set_arrays(in1Dpolytrope, r_Schw_arr,rho_arr,P_arr,M_arr,expnu_arr,exp4phi_arr,rbar_arr) == 1) {
    fprintf(stderr,"ERROR WHEN READING FILE %s!\n",filename);
    exit(1);
  }
  fclose(in1Dpolytrope);

  REAL Rbar = -100;
  int Rbar_idx = -100;
  for(int i=1;i<numlines_in_file;i++) {
    if(rho_arr[i-1]>0 && rho_arr[i]==0) { Rbar = rbar_arr[i-1]; Rbar_idx = i-1; }
  }
  if(Rbar<0) {
    fprintf(stderr,"Error: could not find r=R from data file.\n");
    exit(1);
  }
      
  // Next, interpolate!
  // Create trial radius array:
  int num_r_pts = 100000;
  //REAL *r_out_arr = (REAL *)malloc(sizeof(REAL)*num_r_pts);
  struct drand48_data randBuffer;
  srand48_r(1313, &randBuffer);
#pragma omp parallel for
  for(int i=0;i<num_r_pts;i++) {
    REAL rrbar;
    drand48_r(&randBuffer,&rrbar);
    //rrbar *= 10.; //rbar_arr[numlines_in_file-1];
    rrbar = rrbar*0.1 + 0.8; //rbar_arr[numlines_in_file-1];
    REAL rho,P,M,expnu,exp4phi;
    TOV_interpolate_1D(rrbar,Rbar,Rbar_idx,4,  numlines_in_file,r_Schw_arr,rho_arr,P_arr,M_arr,expnu_arr,exp4phi_arr,rbar_arr,  &rho,&P,&M,&expnu,&exp4phi);
    printf("%e %e %e %e %e %e\n",rrbar,rho,P,M,expnu,exp4phi);
  }

  // Free the malloc()'s!
  free(r_Schw_arr);
  free(rho_arr);
  free(P_arr);
  free(M_arr);
  free(expnu_arr);
  free(exp4phi_arr);
  free(rbar_arr);
  
  return 0;
}
#endif
