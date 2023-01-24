# This module provides functions for setting up Scalar Field initial data
#     as documented in Tutorial-ADM_Initial_Data-ScalarField.ipynb

# Authors: Leonardo R. Werneck
#          wernecklr **at** gmail **dot* com
#          Zachariah B. Etienne

# First we import needed core NRPy+ modules
import os,sys                           # Standard Python modules for multiplatform OS-level functions
import sympy as sp                      # SymPy: The Python computer algebra package upon which NRPy+ depends
import numpy as np                      # NumPy: A large collection of mathematical functions for Python
from scipy.sparse import spdiags        # SciPy: Sparse, tri-diagonal matrix setup function
from scipy.sparse import csc_matrix     # SciPy: Sparse matrix optimization function
from scipy.sparse.linalg import spsolve # SciPy: Solver of linear systems involving sparse matrices
import outputC as outC                  # NRPy+: Core C code output module
import reference_metric as rfm          # NRPy+: Reference metric support

def ScalarField_InitialData(outputname, ID_Family, Nr, rmax,
                            pulse_amplitude, pulse_center, pulse_width):

    # Set the step size and the r coordinate
    dr = rmax/Nr
    r  = np.empty(Nr)
    for i in range(Nr):
        r[i] = (i-0.5)*dr

    # Set the step size squared
    dr2   = dr**2

    # Let's begin by setting the parameters involved in the initial data
    phi0,rr,rr0,sigma = sp.symbols("phi0 rr rr0 sigma",real=True)

    # Now set the initial profile of the scalar field
    if ID_Family == "Gaussian_pulse":
        phiID = phi0 * sp.exp( -rr**2/sigma**2 )
    elif ID_Family == "Gaussian_pulsev2":
        phiID = phi0 * rr**3 * sp.exp( -(rr-rr0)**2/sigma**2 )
    elif ID_Family == "Tanh_pulse":
        phiID = phi0 * ( 1 - sp.tanh( (rr-rr0)**2/sigma**2 ) )
    else:
        print("Unkown initial data family: ",ID_Family)
        print("Available options are: Gaussian_pulse, Gaussian_pulsev2, and Tanh_pulse")
        sys.exit(1)

    # Now compute Phi := \partial_{r}phi
    PhiID = sp.diff(phiID,rr)

    # Now set numpy functions for phi and Phi
    phi = sp.lambdify((phi0,rr,rr0,sigma),phiID)
    Phi = sp.lambdify((phi0,rr,rr0,sigma),PhiID)

    # ## Part A.1c: populating the varphi(0,r) array
    phi0  = pulse_amplitude
    r0    = pulse_center
    sigma = pulse_width
    ID_sf = phi(phi0,r,r0,sigma)

    # Set the main diagonal
    main_diag = np.pi * dr2 * Phi(phi0,r,r0,sigma)**2 - 2

    # Update the first element of the main diagonal
    main_diag[0] += 1 - dr/r[0]

    # Update the last element of the main diagonal
    main_diag[Nr-1] += - (2 * dr / r[Nr-1])*(1 + dr / r[Nr-1])

    # Set the upper diagonal, ignoring the last point in the r array
    upper_diag = np.zeros(Nr)
    upper_diag[1:] = 1 + dr/r[:-1]

    # Set the lower diagonal, start counting the r array at the second element
    lower_diag = np.zeros(Nr)
    lower_diag[:-1] = 1 - dr/r[1:]

    # Change the last term in the lower diagonal to its correct value
    lower_diag[Nr-2] = 2

    # Set the sparse matrix A by adding up the three diagonals
    # https://docs.scipy.org/doc/scipy/reference/generated/scipy.sparse.spdiags.html
    A = spdiags([main_diag,upper_diag,lower_diag],[0,1,-1],Nr,Nr)

    # Then compress the sparse matrix A column wise, so that SciPy can invert it later
    # https://docs.scipy.org/doc/scipy/reference/generated/scipy.sparse.csc_matrix.html
    A = csc_matrix(A)

    # Set up the right-hand side of the linear system: s
    s = np.zeros(Nr)

    # Update the last entry of the vector s
    s[Nr-1] = - (2 * dr / r[Nr-1])*(1 + dr / r[Nr-1])

    # Compress the vector s column-wise
    # https://docs.scipy.org/doc/scipy/reference/generated/scipy.sparse.csc_matrix.html
    s = csc_matrix(s)

    # Solve the sparse linear system using scipy
    # https://docs.scipy.org/doc/scipy-0.14.0/reference/generated/scipy.sparse.linalg.spsolve.html
    psi = spsolve(A, s.T)

    if sys.version_info[0] == 3:
        np.savetxt(outputname, list(zip( r, ID_sf, psi )),fmt="%.15e")
    elif sys.version_info[0] == 2:
        np.savetxt(outputname, zip( r, ID_sf, psi ),fmt="%.15e")

    print("Generated the ADM initial data for the gravitational collapse \n" \
          "of a massless scalar field in spherical coordinates.\n")
    print("Type of initial condition: Scalar field: \"Gaussian\" Shell\n"\
          "                         ADM quantities: Time-symmetric\n")
    print("Parameters: amplitude         = "+str(phi0)+",\n" \
          "            center            = "+str(r0)+",\n"   \
          "            width             = "+str(sigma)+",\n"   \
          "            domain size       = "+str(rmax)+",\n"   \
          "            number of points  = "+str(Nr)+",\n"
          "            Initial data file = "+str(outputname)+".\n")

def ID_scalarfield_ADM_quantities(Ccodesdir=".",new_way=False):
    includes = ["NRPy_basic_defines.h", "NRPy_function_prototypes.h"]
    desc = """(c) 2021 Leo Werneck
This function takes as input either (x,y,z) or (r,th,ph) and outputs
all ADM quantities in the Cartesian or Spherical basis, respectively.
"""
    c_type = "void"
    name   = "ID_scalarfield_ADM_quantities"
    params = """const REAL xyz_or_rthph[3],const ID_inputs other_inputs,
                REAL *restrict gammaDD00,REAL *restrict gammaDD01,REAL *restrict gammaDD02,
                REAL *restrict gammaDD11,REAL *restrict gammaDD12,REAL *restrict gammaDD22,
                REAL *restrict KDD00,REAL *restrict KDD01,REAL *restrict KDD02,
                REAL *restrict KDD11,REAL *restrict KDD12,REAL *restrict KDD22,
                REAL *restrict alpha,
                REAL *restrict betaU0,REAL *restrict betaU1,REAL *restrict betaU2,
                REAL *restrict BU0,REAL *restrict BU1,REAL *restrict BU2"""
    body   = """
  const REAL r  = xyz_or_rthph[0];
  const REAL th = xyz_or_rthph[1];
  const REAL ph = xyz_or_rthph[2];

  REAL sf_star,psi4_star,alpha_star;

  scalarfield_interpolate_1D(r,
                             other_inputs.interp_stencil_size,
                             other_inputs.numlines_in_file,
                             other_inputs.r_arr,
                             other_inputs.sf_arr,
                             other_inputs.psi4_arr,
                             other_inputs.alpha_arr,
                             &sf_star,&psi4_star,&alpha_star);

  // Update alpha
  *alpha = alpha_star;
  // gamma_{rr} = psi^4
  *gammaDD00 = psi4_star;
  // gamma_{thth} = psi^4 r^2
  *gammaDD11 = psi4_star*r*r;
  // gamma_{phph} = psi^4 r^2 sin^2(th)
  *gammaDD22 = psi4_star*r*r*sin(th)*sin(th);

  // All other quantities ARE ZERO:
  *gammaDD01 = 0.0; *gammaDD02 = 0.0;
  /**/              *gammaDD12 = 0.0;

  *KDD00 = 0.0; *KDD01 = 0.0; *KDD02 = 0.0;
  /**/          *KDD11 = 0.0; *KDD12 = 0.0;
  /**/                        *KDD22 = 0.0;

  *betaU0 = 0.0; *betaU1 = 0.0; *betaU2 = 0.0;

  *BU0 = 0.0; *BU1 = 0.0; *BU2 = 0.0;
"""
    if new_way == True:
        outC.add_to_Cfunction_dict(includes=includes,desc=desc,c_type=c_type,name=name,
                                   params=params,body=body,enableCparameters=False)
    else:
        outfile = os.path.join(Ccodesdir,"ID_scalarfield_ADM_quantities.h")
        outC.outCfunction(outfile=outfile,
                          includes=None,desc=desc,c_type=c_type,name=name,
                          params=params,body=body,enableCparameters=False)

def ID_scalarfield_spherical(Ccodesdir=".",new_way=False):
    includes = ["NRPy_basic_defines.h", "NRPy_function_prototypes.h"]
    desc = """(c) 2021 Leo Werneck
This function takes as input either (x,y,z) or (r,th,ph) and outputs all
scalar field quantities in the Cartesian or Spherical basis, respectively.
"""
    c_type = "void"
    name   = "ID_scalarfield_spherical"
    params = "const REAL xyz_or_rthph[3],const ID_inputs other_inputs,REAL *restrict sf,REAL *restrict sfM"
    body   = """
  const REAL r  = xyz_or_rthph[0];
  const REAL th = xyz_or_rthph[1];
  const REAL ph = xyz_or_rthph[2];

  REAL sf_star,psi4_star,alpha_star;

  scalarfield_interpolate_1D(r,
                             other_inputs.interp_stencil_size,
                             other_inputs.numlines_in_file,
                             other_inputs.r_arr,
                             other_inputs.sf_arr,
                             other_inputs.psi4_arr,
                             other_inputs.alpha_arr,
                             &sf_star,&psi4_star,&alpha_star);

  // Update varphi
  *sf  = sf_star;
  // Update Pi
  *sfM = 0;
"""
    if new_way == True:
        outC.add_to_Cfunction_dict(includes=includes,desc=desc,c_type=c_type,name=name,
                                   params=params,body=body,enableCparameters=False)
    else:
        outfile = os.path.join(Ccodesdir,"ID_scalarfield_spherical.h")
        outC.outCfunction(outfile=outfile,
                          includes=None,desc=desc,c_type=c_type,name=name,
                          params=params,body=body,enableCparameters=False)

def ID_scalarfield_xx0xx1xx2_to_BSSN_xx0xx1xx2(Ccodesdir=".",pointer_to_ID_inputs=False,new_way=False):

    rfm.reference_metric()

    rthph = outC.outputC(rfm.xxSph[0:3],["rthph[0]", "rthph[1]", "rthph[2]"],
                         "returnstring", "includebraces=False,outCverbose=False,preindent=1")

    includes = ["NRPy_basic_defines.h", "NRPy_function_prototypes.h"]
    desc = """(c) 2021 Leo Werneck
This function takes as input either (x,y,z) or (r,th,ph) and outputs all
scalar field quantities in the Cartesian or Spherical basis, respectively.
"""
    c_type = "void"
    name   = "ID_scalarfield_xx0xx1xx2_to_BSSN_xx0xx1xx2"
    params = "const paramstruct *restrict params,const REAL xx0xx1xx2[3],\n"
    if pointer_to_ID_inputs == True:
        params += "ID_inputs *other_inputs,\n"
    else:
        params += "ID_inputs other_inputs,\n"
    params += "REAL *restrict sf, REAL *restrict sfM"
    body   = """
  const REAL xx0 = xx0xx1xx2[0];
  const REAL xx1 = xx0xx1xx2[1];
  const REAL xx2 = xx0xx1xx2[2];
  REAL rthph[3];
"""+rthph+"""
  ID_scalarfield_spherical(rthph,other_inputs,sf,sfM);
"""

    if new_way == True:
        outC.add_to_Cfunction_dict(includes=includes,desc=desc,c_type=c_type,name=name,
                                   params=params,body=body)
    else:
        outfile = os.path.join(Ccodesdir,"ID_scalarfield_xx0xx1xx2_to_BSSN_xx0xx1xx2.h")
        outC.outCfunction(outfile=outfile,
                          includes=None,desc=desc,c_type=c_type,name=name,
                          params=params,body=body)

def ID_scalarfield(Ccodesdir=".",new_way=False):
    includes = ["NRPy_basic_defines.h", "NRPy_function_prototypes.h"]
    desc = """(c) 2021 Leo Werneck
This is the scalar field initial data driver functiono.
"""
    c_type = "void"
    name   = "ID_scalarfield"
    params = """const paramstruct *restrict params,REAL *restrict xx[3],
                ID_inputs other_inputs,REAL *restrict in_gfs"""
    body   = """
  const int idx = IDX3S(i0,i1,i2);
  const REAL xx0xx1xx2[3] = {xx0,xx1,xx2};
  ID_scalarfield_xx0xx1xx2_to_BSSN_xx0xx1xx2(params,xx0xx1xx2,other_inputs,
                                             &in_gfs[IDX4ptS(SFGF,idx)],
                                             &in_gfs[IDX4ptS(SFMGF,idx)]);
"""
    loopopts = "AllPoints,Read_xxs"
    if new_way == True:
        outC.add_to_Cfunction_dict(includes=includes,desc=desc,c_type=c_type,name=name,
                                   params=params,body=body,loopopts=loopopts)
    else:
        outfile = os.path.join(Ccodesdir,"ID_scalarfield.h")
        outC.outCfunction(outfile=outfile,
                          includes=None,desc=desc,c_type=c_type,name=name,
                          params=params,body=body,loopopts=loopopts)

def NRPy_param_funcs_register_C_functions_and_NRPy_basic_defines(Ccodesdir=".",pointer_to_ID_inputs=False,new_way=False):
    ID_scalarfield_ADM_quantities(Ccodesdir=Ccodesdir,new_way=new_way)
    ID_scalarfield_spherical(Ccodesdir=Ccodesdir,new_way=new_way)
    ID_scalarfield_xx0xx1xx2_to_BSSN_xx0xx1xx2(Ccodesdir=Ccodesdir,pointer_to_ID_inputs=pointer_to_ID_inputs,new_way=new_way)
    ID_scalarfield(Ccodesdir=Ccodesdir,new_way=new_way)
