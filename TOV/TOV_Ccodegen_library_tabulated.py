## TOV C codegen library, for reading/interpolating TOV solution data
## Documented in Tutorial-ADM_Initial_Data-TOV.ipynb
## Author: Zachariah B. Etienne

# Import needed Python/NRPy+ modules
from outputC import outputC,add_to_Cfunction_dict  # NRPy+: Core C code output module
import sympy as sp                  # SymPy: The Python computer algebra package upon which NRPy+ depends
import NRPy_param_funcs as par      # NRPy+: Parameter interface
import indexedexp as ixp            # NRPy+: Symbolic indexed expression (e.g., tensors, vectors, etc.) support
import reference_metric as rfm      # NRPy+: Reference metric support


def ID_persist_str():
    return r"""
  REAL Rbar;      // rbar corresponding to the outermost radius at which density is nonzero
  int Rbar_idx;   // Index (line of data file) corresponding to Rbar
  int interp_stencil_size;  // Lagrange polynomial interpolation stencil size
  int numlines_in_file;     // Number of radii stored in the file
  // The following arrays store stellar information at all numlines_in_file radii:
  REAL *restrict r_Schw_arr; // Stellar radial coordinate in units of Schwarzschild radius
  REAL *restrict rho_arr;    // Mass-energy density
  REAL *restrict rho_baryon_arr; // Baryonic mass density
  REAL *restrict P_arr;  // Pressure
  REAL *restrict Y_e_arr;  // electron fraction
  //REAL *restrict entropy_arr;  // entropy
  REAL *restrict M_arr;  // Integrated rest mass
  REAL *restrict expnu_arr;    // Metric quantity
  REAL *restrict exp4phi_arr;  // Metric quantity
  REAL *restrict rbar_arr;  // Rbar coordinate
"""

# TOV_read_data_file_set_ID_persist(): Read TOV data file and store data to the ID_persist struct
def add_to_Cfunction_dict_TOV_read_data_file_set_ID_persist(interp_stencil_size=12):
    includes = ["NRPy_basic_defines.h"]
    desc = "Returns the number of lines in a TOV data file."
    c_type = "void"
    name = "TOV_read_data_file_set_ID_persist"
    params = "const char *input_filename, ID_persist_struct *ID_persist"
    body = r"""
  char filename[100];
  snprintf(filename, 100, input_filename);
  FILE *TOV_solution_datafile = fopen(filename, "r");
  if(TOV_solution_datafile == NULL) {
    fprintf(stderr,"ERROR: could not open TOV solution data file %s\n",filename);
    exit(1);
  }

  // First set the interpolation stencil size. Note that the interpolation is set up to
  //   avoid interpolations across the star's surface. Hence we can use a super high
  //   order interpolant.
  ID_persist->interp_stencil_size = """+str(interp_stencil_size)+r""";

  int numlines_in_file = 0;
  {
    char * line = NULL;

    size_t len = 0;
    ssize_t read;
    while ((read = getline(&line, &len, TOV_solution_datafile)) != -1) {
      numlines_in_file++;
    }
    rewind(TOV_solution_datafile);

    free(line);
  }
  ID_persist->numlines_in_file = numlines_in_file;

  // Now that numlines_in_file is set, we can now allocate memory for all arrays.
  {
    ID_persist->r_Schw_arr     = (REAL *restrict)malloc(sizeof(REAL)*numlines_in_file);
    ID_persist->rho_arr        = (REAL *restrict)malloc(sizeof(REAL)*numlines_in_file);
    ID_persist->rho_baryon_arr = (REAL *restrict)malloc(sizeof(REAL)*numlines_in_file);
    ID_persist->P_arr          = (REAL *restrict)malloc(sizeof(REAL)*numlines_in_file);
    ID_persist->Y_e_arr        = (REAL *restrict)malloc(sizeof(REAL)*numlines_in_file);
    //ID_persist->entropy_arr    = (REAL *restrict)malloc(sizeof(REAL)*numlines_in_file);
    ID_persist->M_arr          = (REAL *restrict)malloc(sizeof(REAL)*numlines_in_file);
    ID_persist->expnu_arr      = (REAL *restrict)malloc(sizeof(REAL)*numlines_in_file);
    ID_persist->exp4phi_arr    = (REAL *restrict)malloc(sizeof(REAL)*numlines_in_file);
    ID_persist->rbar_arr       = (REAL *restrict)malloc(sizeof(REAL)*numlines_in_file);
  }

  {
    char * line = NULL;

    size_t len = 0;
    ssize_t read;

    int which_line = 0;
    while ((read = getline(&line, &len, TOV_solution_datafile)) != -1) {
      // Define the line delimiters (i.e., the stuff that goes between the data on a given
      //     line of data.  Here, we define both spaces " " and tabs "\t" as data delimiters.
      const char delimiters[] = " \t";

      // Now we define "token", a pointer to the first column of data
      char *token;

      // Each successive time we call strtok(NULL,blah), we read in a new column of data from
      //     the originally defined character array, as pointed to by token.

      token=strtok(line, delimiters); if(token==NULL) { fprintf(stderr, "Error reading %s\n", filename); exit(1); }
      ID_persist->r_Schw_arr[which_line]     = strtod(token, NULL); token = strtok( NULL, delimiters );
      ID_persist->rho_arr[which_line]        = strtod(token, NULL); token = strtok( NULL, delimiters );
      ID_persist->rho_baryon_arr[which_line] = strtod(token, NULL); token = strtok( NULL, delimiters );
      ID_persist->P_arr[which_line]          = strtod(token, NULL); token = strtok( NULL, delimiters );
      ID_persist->Y_e_arr[which_line]          = strtod(token, NULL); token = strtok( NULL, delimiters );
      //ID_persist->entropy_arr[which_line]          = strtod(token, NULL); token = strtok( NULL, delimiters );
      ID_persist->M_arr[which_line]          = strtod(token, NULL); token = strtok( NULL, delimiters );
      ID_persist->expnu_arr[which_line]      = strtod(token, NULL); token = strtok( NULL, delimiters );
      ID_persist->exp4phi_arr[which_line]    = strtod(token, NULL); token = strtok( NULL, delimiters );
      ID_persist->rbar_arr[which_line]       = strtod(token, NULL);

      which_line++;
    }
    free(line);

    fclose(TOV_solution_datafile);
  }

  {
    // Finally set Rbar and Rbar_idx
    ID_persist->Rbar     = -100.0;
    ID_persist->Rbar_idx = -100;
    for(int i=1;i<numlines_in_file;i++) {
      if(ID_persist->rho_arr[i-1] > 0  &&  ID_persist->rho_arr[i] == 0) {
        ID_persist->Rbar = ID_persist->rbar_arr[i-1];
        ID_persist->Rbar_idx = i-1;
      }
    }
    if(ID_persist->Rbar < 0) {
      fprintf(stderr,"Error: could not find rbar=Rbar (i.e., the surface of the star) from data file.\n");
      exit(1);
    }
  }
"""
    add_to_Cfunction_dict(
        includes=includes,
        desc=desc,
        c_type=c_type, name=name, params=params,
        body=body,
        enableCparameters=False)

# TOV_interpolate_1D(): Interpolate TOV data to any desired distance from the center of the star
def add_to_Cfunction_dict_TOV_interpolate_1D():
    includes=["NRPy_basic_defines.h"]
    prefunc = r"""
// Find interpolation index using Bisection root-finding algorithm:
static inline int bisection_idx_finder(const REAL rrbar, const int numlines_in_file, const REAL *restrict rbar_arr) {
  int x1 = 0;
  int x2 = numlines_in_file-1;
  REAL y1 = rrbar-rbar_arr[x1];
  REAL y2 = rrbar-rbar_arr[x2];
  if(y1*y2 >= 0) {
    fprintf(stderr,"INTERPOLATION BRACKETING ERROR %e | %e %e\n",rrbar,y1,y2);
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
"""
    desc = """Read a TOV solution from data file and perform
1D interpolation of the solution to a desired radius.

Author: Zachariah B. Etienne
        zachetie **at** gmail **dot* com
"""
    c_type="void"
    name="TOV_interpolate_1D"
    params="""REAL rrbar,const ID_persist_struct *ID_persist,
                        REAL *restrict rho,REAL *restrict rho_baryon,REAL *restrict P, REAL *restrict M,REAL *restrict expnu,REAL *restrict exp4phi"""
    body = r"""
  // First unpack ID_persist struct (this should be pretty quick relative to the rest of the routine):
  const REAL Rbar               = ID_persist->Rbar;
  const int Rbar_idx            = ID_persist->Rbar_idx;
  const int interp_stencil_size = ID_persist->interp_stencil_size;
  const int numlines_in_file    = ID_persist->numlines_in_file;
  const REAL *restrict rbar_arr       = ID_persist->rbar_arr;
  const REAL *restrict r_Schw_arr     = ID_persist->r_Schw_arr;
  const REAL *restrict rho_arr        = ID_persist->rho_arr;
  const REAL *restrict rho_baryon_arr = ID_persist->rho_baryon_arr;
  const REAL *restrict P_arr          = ID_persist->P_arr;
  const REAL *restrict Y_e_arr        = ID_persist->Y_e_arr;
  //const REAL *restrict entropy_arr    = ID_persist->entropy_arr;
  const REAL *restrict M_arr          = ID_persist->M_arr;
  const REAL *restrict expnu_arr      = ID_persist->expnu_arr;
  const REAL *restrict exp4phi_arr    = ID_persist->exp4phi_arr;

  // For this case, we know that for all functions, f(r) = f(-r)
  if(rrbar < 0) rrbar = -rrbar;

  // First find the central interpolation stencil index:
  int idx = bisection_idx_finder(rrbar,numlines_in_file,rbar_arr);


  int idxmin = MAX(0,idx-interp_stencil_size/2-1);

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
  *rho_baryon = 0.0;
  *P = 0.0;
  //*Y_e = 0.0;
  //*entropy = 0.0;
  *M = 0.0;
  *expnu = 0.0;
  *exp4phi = 0.0;

  REAL r_Schw = 0.0;
  for(int i=idxmin;i<idxmin+interp_stencil_size;i++) {
    r_Schw      += l_i_of_r[i-idxmin] * r_Schw_arr[i];
    *rho        += l_i_of_r[i-idxmin] * rho_arr[i];
    *rho_baryon += l_i_of_r[i-idxmin] * rho_baryon_arr[i];
    *P          += l_i_of_r[i-idxmin] * P_arr[i];
    //*Y_e        += l_i_of_r[i-idxmin] * Y_e_arr[i];
    //*entropy    += l_i_of_r[i-idxmin] * entropy_arr[i];
    *M          += l_i_of_r[i-idxmin] * M_arr[i];
    *expnu      += l_i_of_r[i-idxmin] * expnu_arr[i];
    *exp4phi    += l_i_of_r[i-idxmin] * exp4phi_arr[i];
  }

  if(rrbar > Rbar) {
    *rho        = 0;
    *rho_baryon = 0;
    *P          = 0;
    //*Y_e        = 0;
    //*entropy    = 0;
    *M          = M_arr[Rbar_idx+1];
    *expnu      = 1. - 2.*(*M) / r_Schw;
    *exp4phi    = pow(r_Schw / rrbar,2.0);
  }
"""
    add_to_Cfunction_dict(
        includes=includes,
        prefunc=prefunc,
        desc=desc,
        c_type=c_type, name=name, params=params,
        body=body,
        enableCparameters=False)


def ADM_quantities_ito_TOV_soln(rbar, theta, expnu, exp4phi):
    # in TOV ID, betaU=BU=KDD=0
    alpha = sp.sqrt(expnu)

    betaU = ixp.zerorank1()
    BU = ixp.zerorank1()

    gammaDD = ixp.zerorank2(DIM=3)
    gammaDD[0][0] = exp4phi
    gammaDD[1][1] = exp4phi * rbar ** 2
    gammaDD[2][2] = exp4phi * rbar ** 2 * sp.sin(theta) ** 2

    KDD = ixp.zerorank2()

    return alpha, betaU, BU, gammaDD, KDD


def T4UU_ito_TOV_soln(rbar, theta, rho, P, expnu, exp4phi):

    T4UU = ixp.zerorank2(DIM=4)
    T4UU[0][0] = rho/expnu
    T4UU[1][1] = P/exp4phi
    T4UU[2][2] = P/(exp4phi*rbar**2)
    T4UU[3][3] = P/(exp4phi*rbar**2*sp.sin(theta)**2)

    return T4UU


# TOV_ID_function(): Driver routine for setting up ADM metric quantities and T4UU from interpolated TOV data
def add_to_Cfunction_dict_TOV_ID_function():
    includes = ["NRPy_basic_defines.h", "NRPy_function_prototypes.h"]
    desc = """This function sets TOV initial data at a single point xCart[3] = (x,y,z),
outputting to initial_data_struct *restrict initial_data"""
    c_type = "void"
    name = "TOV_ID_function"
    params = "const paramstruct *params, const REAL xCart[3], const ID_persist_struct *restrict ID_persist, initial_data_struct *restrict initial_data"
    body  = """  // First set r(=rbar), theta, phi in terms of xCart[3]:
  const REAL Cartx=xCart[0], Carty=xCart[1], Cartz=xCart[2];"""
    desired_rfm_coord = par.parval_from_str("reference_metric::CoordSystem")
    # "Spherical" is the native basis for TOV initial data.
    par.set_parval_from_str("reference_metric::CoordSystem", "Spherical")
    rfm.reference_metric()
    body += r"""
  REAL rbar, theta, phi;
  {
""" + outputC(rfm.Cart_to_xx[:3], ["rbar", "theta", "phi"], filename="returnstring",
              params="outCverbose=False,includebraces=False,preindent=2")
    body += r"""
  }

  // Next set gamma_{ij} in spherical basis
  REAL rho,rho_baryon,P,M,expnu,exp4phi;
  TOV_interpolate_1D(rbar,ID_persist, &rho,&rho_baryon,&P,&M,&expnu,&exp4phi);
"""
    rbar, theta, rho, P, expnu, exp4phi = sp.symbols('rbar theta rho P expnu exp4phi', real=True)
    alpha, betaU, BU, gammaDD, KDD = ADM_quantities_ito_TOV_soln(rbar, theta, expnu, exp4phi)
    T4UU = T4UU_ito_TOV_soln(rbar, theta, rho, P, expnu, exp4phi)

    list_of_output_exprs = [alpha]
    list_of_output_varnames = ["initial_data->alpha"]
    for i in range(3):
        list_of_output_exprs += [betaU[i]]
        list_of_output_varnames += ["initial_data->betaSphorCartU" + str(i)]
        list_of_output_exprs += [BU[i]]
        list_of_output_varnames += ["initial_data->BSphorCartU" + str(i)]
        for j in range(i, 3):
            list_of_output_exprs += [gammaDD[i][j]]
            list_of_output_varnames += ["initial_data->gammaSphorCartDD" + str(i) + str(j)]
            list_of_output_exprs += [KDD[i][j]]
            list_of_output_varnames += ["initial_data->KSphorCartDD" + str(i) + str(j)]
    for mu in range(4):
        for nu in range(mu, 4):
            list_of_output_exprs += [T4UU[mu][nu]]
            list_of_output_varnames += ["initial_data->T4SphorCartUU" + str(mu) + str(nu)]

    # Sort the outputs before calling outputC()
    # https://stackoverflow.com/questions/9764298/is-it-possible-to-sort-two-listswhich-reference-each-other-in-the-exact-same-w
    list_of_output_varnames, list_of_output_exprs = (list(t) for t in zip(*sorted(zip(list_of_output_varnames, list_of_output_exprs))))

    body += outputC(list_of_output_exprs, list_of_output_varnames,
                    filename="returnstring", params="outCverbose=False,includebraces=False,preindent=1")

    # Restore CoordSystem:
    par.set_parval_from_str("reference_metric::CoordSystem", desired_rfm_coord)
    rfm.reference_metric()

    add_to_Cfunction_dict(
        includes=includes,
        desc=desc,
        c_type=c_type, name=name, params=params,
        body=body,
        enableCparameters=False)
