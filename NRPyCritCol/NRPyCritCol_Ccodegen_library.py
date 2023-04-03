## Title    : NRPyCritCol C codegen library
## Author   : Leo Werneck
## Tutorial :
## Desc     : C functions for reading/interpolating initial data

from sympy import symbols, exp, tanh, sympify
from outputC import outputC, add_to_Cfunction_dict, lhrh
from outputC import outC_NRPy_basic_defines_h_dict
from indexedexp import zerorank1, zerorank2
from NRPy_param_funcs import Cparameters
from NRPy_param_funcs import parval_from_str
from NRPy_param_funcs import set_parval_from_str
from reference_metric import reference_metric
from reference_metric import Cart_to_xx
import BSSN.ADM_Initial_Data_Reader__BSSN_Converter as IDread
import ScalarField.ScalarField_quantities as SFq

def ID_persist_str():
    return r"""
  /* The following store the scalar field and spacetime initial data */
  int Nr;             /* Number of points in ID radial grid */
  REAL *restrict r;   /* ID radial grid                     */
  REAL *restrict phi; /* Initial scalar field profile       */
  REAL *restrict psi; /* Initial conformal factor           */

  /* Interpolation parameters */
  int interp_stencil_size; /* Size of intepolation stencil  */
"""

def add_to_Cfunction_dict_tridiagonal_solver_thomas_algorithm():
    desc     = """(c) 2023 Leo Werneck

Tridiagonal system solver using Thomas algorithm
Reference: https://en.wikipedia.org/wiki/Tridiagonal_matrix_algorithm"""
    includes = ["NRPy_basic_defines.h"]
    name     = "tridiagonal_solver_thomas_algorithm"
    params   = r"""
      const int n,            /* Number of points */
      const REAL *restrict a, /* Lower diagonal   */
      REAL *restrict b,       /* Main  diagonal   */
      const REAL *restrict c, /* Upper diagonal   */
      REAL *restrict d,       /* Source   vector  */
      REAL *restrict x        /* Solution vector  */"""
    body     = r"""  // Step 1: First pass prepares b and d
  if( b[0] == 0.0 ) {
    fprintf(stderr, "(%s) Error: first diagonal element is zero\n", __func__);
    exit(1);
  }
  for(int i=1;i<n;i++) {
    const REAL w = a[i-1]/b[i-1];
    b[i] -= w * c[i-1];
    d[i] -= w * d[i-1];
  }

  // Step 2: Second pass performs the back-substitution
  const int N = n-1;
  x[N] = d[N] / b[N];
  for(int i=N;i>0;i--)
    x[i-1] = (d[i-1]-c[i-1]*x[i])/b[i-1];
"""
    add_to_Cfunction_dict(includes=includes, desc=desc, name=name,
                          params=params, body=body, enableCparameters=False)

# Step 3: Initial scalar field profile
def add_to_Cfunction_dict_NRPyCritCol_initial_profiles():
    # Step 3.a: Implement symbolic expressions for phi(0,r)
    Nr = Cparameters("int", "NRPyCritCol", "Nr", 320)
    r = symbols("r", real=True)
    eta, r0, sigma = Cparameters("REAL", "NRPyCritCol", ["eta", "r0", "sigma"], [1, 0, 1])
    global sf_profiles
    sf_profiles = {}
    sf_profiles["gaussian_pulse"   ] = eta * exp( -(r-r0)**2 / sigma**2 )
    sf_profiles["gaussian_pulse_v2"] = eta * r**3 * exp( -(r-r0)**2 / sigma**2 )
    sf_profiles["tanh_pulse"       ] = eta * ( 1 - tanh((r-r0)**2/sigma**2) )

    # Step 3.b: Get maximum length of keys (Optional)
    max_length = 0
    for key in sf_profiles.keys():
        if len(key) > max_length:
            max_length = len(key)

    # Step 3.c: Python function that generates C functions for the initial 
    #           scalar field profile and its radial derivative
    def generate_phi_or_Phi_function(var, expr, profile):
        desc     = """(c) 2023 Leo Werneck

NRPyCritCol: Initial """+var+" profile - "+profile
        outCparams = "includebraces=False,outCverbose=False,preindent=1"
        includes   = ["NRPy_basic_defines.h"]
        c_type     = "REAL"
        name       = "NRPyCritCol_initial_"+var+"_"+profile
        params     = """
      const REAL r,    /* Radial coordinate    */
      const REAL eta,  /* Pulse amplitude      */
      const REAL r0,   /* Pulse center         */
      const REAL sigma /* Pulse width          */"""
        body = outputC(expr, "return", "returnstring",
                       params=outCparams).replace("_ =","")
        add_to_Cfunction_dict(includes=includes, prefunc=prefunc, desc=desc,
                              name=name, params=params, body=body, c_type=c_type,
                              enableCparameters=False)

    # Step 3.d: Finally, generate all C code for computing the initial
    #           scalar field profile and its radial derivative, including
    #           "driver" functions that select which profile to use
    includes = ["NRPy_basic_defines.h", "NRPy_function_prototypes.h"]
    prefunc  = ""
    c_type   = "REAL"
    name     = "NRPyCritCol_initial_var"
    params   = """
      const int id_type, /* Type of initial data */
      const REAL r,      /* Radial coordinate    */
      const REAL eta,    /* Pulse amplitude      */
      const REAL r0,     /* Pulse center         */
      const REAL sigma   /* Pulse width          */"""
    body     = """
  switch (id_type) {
"""
    counter = 0
    defines = ""
    for profile, phi in sf_profiles.items():
        Phi = phi.diff(r)
        generate_phi_or_Phi_function("phi", phi, profile)
        generate_phi_or_Phi_function("phi_radial_deriv", Phi, profile)
        defines += "#define "+profile.upper()
        for i in range(max_length - len(profile) + 1):
            defines += " "
        defines += str(counter)+"\n"
        body    += "  case "+profile.upper()+""":
    return NRPyCritCol_initial_var_"""+profile+"(r, eta, r0, sigma);\n    break;\n"
        counter += 1
    body += r"""  default:
    fprintf(stderr, "Unknown initial data key (%d)\n", id_type);
    exit(1);
  }
"""
    outC_NRPy_basic_defines_h_dict["ScalarFieldID"] = defines
    for var in ["phi", "phi_radial_deriv"]:
        newname = name.replace("var", var)
        newbody = body.replace("var", var)
        add_to_Cfunction_dict(includes=includes, name=newname, params=params,
                              c_type=c_type, body=newbody, enableCparameters=False)

def add_to_Cfunction_dict_NRPyCritCol_initialize_ID_persist():
    
    includes = ["NRPy_basic_defines.h"]
    desc     = "(c) 2023 Leo Werneck\n\nNRPyCritCol - Solve the Hamiltonian constraint and initialize ID persist"
    c_type   = "void"
    name     = "NRPyCritCol_solve_Hamiltonian_constraint_set_ID_persist"
    params   = """
      const paramstruct *restrict params,
      ID_persist_struct *restrict ID_persist"""
    body     = """
  // Step 1: Set auxiliary quantities
  ID_persist->Nr  = Nr;                                      /* Number of points */
  ID_persist->r   = (double *)malloc(sizeof(double)*Nr);     /* Radial grid      */
  ID_persist->phi = (double *)malloc(sizeof(double)*Nr);     /* Scalar field     */
  ID_persist->psi = (double *)malloc(sizeof(double)*Nr);     /* Conformal factor */
  double *a       = (double *)malloc(sizeof(double)*(Nr-1)); /* Lower diagonal   */
  double *b       = (double *)malloc(sizeof(double)*Nr);     /* Main  diagonal   */
  double *c       = (double *)malloc(sizeof(double)*(Nr-1)); /* Upper diagonal   */
  double *d       = (double *)malloc(sizeof(double)*Nr);     /* Sources          */

  // Step 2: Set ID type
  const int id_type = GAUSSIAN_PULSE;

  // Step 3: Compute step sizes
  const double dr  = rmax/Nr;
  const double dr2 = dr*dr;

  // Step 4: Build tridiagonal system
  for(int i=0;i<Nr-1;i++) {
    // Step 4.a: Set helper variables
    const double r    = (i-0.5)*dr;
    const double Phi  = NRPyCritCol_initial_phi_radial_deriv(id_type, r, eta, r0, sigma);
    const double tmp  = dr/(r+dr);
  
    // Step 4.b: Set arrays
    ID_persist->r[i]   = r;
    ID_persist->phi[i] = NRPyCritCol_initial_phi(id_type, r, eta, r0, sigma);
    a[i]               = 1 - tmp + (i==Nr-2)*(1 + tmp);
    b[i]               = M_PI * dr2 * Phi * Phi - 2;
    c[i]               = 1 + dr/r;
    d[i]               = 0;
  }

  // Step 5: Adjust remaining entries
  const int N        = Nr-1;
  const double r_N   = (N-0.5)*dr;
  const double tmp0  = dr/r_N;
  const double tmp1  = -2 * tmp0 * (1 + tmp0);
  const double Phi   = NRPyCritCol_initial_phi_radial_deriv(id_type, r_N, eta, r0, sigma);
  ID_persist->r[N]   = r_N;
  ID_persist->phi[N] = NRPyCritCol_initial_phi(id_type, r_N, eta, r0, sigma);
  b[0]              += 1 - dr/rr[0];
  b[N]               = M_PI * dr2 * Phi * Phi - 2 + tmp1;
  d[N]               = tmp1;

  // Step 6: Solve the system
  tridiagonal_solver_thomas_algorithm(Nr, a, b, c, d, ID_persist->psi);

  // Step 7: Free memory of auxiliary variables
  free(a);
  free(b);
  free(c);
  free(d);
"""

def add_to_Cfunction_dict_NRPyCritCol_interpolate_1D():
    includes = ["NRPy_basic_defines.h"]
    desc     = "(c) 2023 Leo Werneck\n\nNRPyCritCol - Interpolate ID to numerical grid"
    prefunc  = r"""
// Bisection routine
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
"""
    c_type   = "void"
    name     = "NRPyCritCol_interpolate_1D"
    params   = """
      const REAL r_star,
      const ID_persist_struct *restrict ID_persist,
      REAL *restrict phi_star,
      REAL *restrict psi_star"""

    body     = r"""
  // Step 1: Unpack the ID_persist struct
  const int interp_stencil_size = ID_persist->interp_stencil_size;
  const int Nr                  = ID_persist->Nr;
  const REAL *restrict r_arr    = ID_persist->r;
  const REAL *restrict phi_arr  = ID_persist->phi;
  const REAL *restrict psi_arr  = ID_persist->psi;

  // Step 2: Set interpolation radius
  const REAL r_interp = r_star*( 1 - 2*(r_star<0) );

  // Step 3: Find the central interpolation stencil index
  int idx    = bisection_idx_finder(r_interp, Nr, r_arr);
  int idxmin = MAX(0, idx - interp_stencil_size/2 - 1);

  // Step 4: Find interpolation coefficients
  REAL r_sample[interp_stencil_size];
  for(int i=idxmin;i<idxmin+interp_stencil_size;i++)
    r_sample[i-idxmin] = r_arr[i];

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

  // Step 5: Now perform the interpolation
  // Then perform the interpolation:
  *phi_star = 0.0;
  *psi_star = 0.0;
  for(int i=idxmin;i<idxmin+interp_stencil_size;i++) {
    *phi_star += l_i_of_r[i-idxmin] * phi_arr[i];
    *psi_star += l_i_of_r[i-idxmin] * psi_arr[i];
  }
"""
    add_to_Cfunction_dict(
        includes=includes, prefunc=prefunc, desc=desc,
        c_type=c_type, name=name, params=params, body=body,
        enableCparameters=False)

def ADM_quantities_ito_NRPyCritCol_soln(phi, psi, lapse_id_type="pre-collapsed"):

    # Step 1: Lapse initial data
    if lapse_id_type == "pre-collapsed":
        alpha = psi**(-2)
    elif lapse_id_type == "unity":
        alpha = sympify(1)
    else:
        raise ValueError("Unsupported lapse ID type \""+lapse_id_type+"\"")

    # Step 2: Shift initial data
    betaU = zerorank1()
    BU    = zerorank1()

    # Step 3: Physical metric
    gammaDD = zerorank2()
    psi4    = psi**4
    for i in range(3):
        gammaDD[i][i] = psi4

    # Step 4: Extrinsic curvature
    KDD = zerorank2()

    return alpha, betaU, BU, gammaDD, KDD

def add_to_Cfunction_dict_NRPyCritCol_ID_function(lapse_id_type="pre-collapsed"):
    includes = ["NRPy_basic_defines.h", "NRPy_function_prototypes.h"]
    desc     = "(c) 2023 Leo Werneck\n\nNRPyCritCol ID driver"
    c_type   = "void"
    name     = "NRPyCritCol_ID_function"
    params   = """
      const paramstruct *params,
      const REAL xCart[3],
      const ID_persist_struct *restrict ID_persist,
      initial_data_struct *restrict initial_data"""
    body     = """  // Step 1: Set r, theta, phi in terms of xCart[3]
  const REAL Cartx=xCart[0], Carty=xCart[1], Cartz=xCart[2];
  REAL r_star, th, ph;\n"""
    desired_rfm_coord = parval_from_str("reference_metric::CoordSystem")
    # "Spherical" is the native basis for NRPyCritCol initial data.
    set_parval_from_str("reference_metric::CoordSystem", "Spherical")
    reference_metric()
    body += outputC(Cart_to_xx[:3],
                    ["r_star", "const REAL th", "const REAL ph"],
                    filename="returnstring",
                    params="outCverbose=False,preindent=1")
    body += r"""
  // Step 2: Interpolate the initial data to r_star
  REAL sf_star, psi_star;
  NRPyCritCol_interpolate_1D(r_star, ID_persist, &sf_star, &psi_star);
"""
    # Add scalar field quantities to the initial data structs
    
    # Now compute metric and scalar field quantities
    SFq.ScalarField_quantities()
    sf_star, psi_star = symbols("sf_star psi_star", real=True)
    alpha, betaU, BU, gammaDD, KDD = ADM_quantities_ito_NRPyCritCol_soln(sf_star, psi_star, lapse_id_type)

    list_of_output_exprs = [alpha, sf_star, sympify(0)]
    list_of_output_varnames = ["initial_data->alpha", "initial_data->sf", "initial_data->sfM"]
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

    # Sort the outputs before calling outputC()
    # https://stackoverflow.com/questions/9764298/is-it-possible-to-sort-two-listswhich-reference-each-other-in-the-exact-same-w
    list_of_output_varnames, list_of_output_exprs = (list(t) for t in zip(*sorted(zip(list_of_output_varnames, list_of_output_exprs))))

    body += outputC(list_of_output_exprs,
                    list_of_output_varnames,
                    filename="returnstring",
                    params="outCverbose=False,includebraces=False,preindent=1")

    # Restore CoordSystem
    set_parval_from_str("reference_metric::CoordSystem", desired_rfm_coord)
    reference_metric()
    
    add_to_Cfunction_dict(
        includes=includes, desc=desc, c_type=c_type, name=name,
        params=params, body=body, enableCparameters=False)

def NRPyCritCol_register_C_functions_and_NRPy_basic_defines(lapse_id_type="pre-collapsed"):
    IDread.add_to_Cfunction_dict_initial_data_reader__convert_ADM_Sph_or_Cart_to_BSSN(input_Coord="Spherical")
    IDread.register_NRPy_basic_defines(ID_persist_struct_contents_str=ID_persist_str())
    add_to_Cfunction_dict_tridiagonal_solver_thomas_algorithm()
    add_to_Cfunction_dict_NRPyCritCol_initial_profiles()
    add_to_Cfunction_dict_NRPyCritCol_initialize_ID_persist()
    add_to_Cfunction_dict_NRPyCritCol_interpolate_1D()
    add_to_Cfunction_dict_NRPyCritCol_ID_function(lapse_id_type)
