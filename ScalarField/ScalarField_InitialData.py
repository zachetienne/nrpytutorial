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
from outputC import outputC             # NRPy+: Core C code output module
import loop as lp                       # NRPy+: C loops module
import reference_metric as rfm          # NRPy+: Reference metric support

def ScalarField_InitialData(outputname,Ccodesdir,ID_Family,
                            pulse_amplitude,pulse_center,pulse_width,NR,rmax,
                            lapse_condition="Pre-collapsed",CoordSystem="Spherical",
                            sinhA=None,sinhW=None):
    
    if CoordSystem == "Spherical":
        r  = np.linspace(0,rmax,NR+1) # Set the r array
        dr = np.zeros(NR)
        for i in range(NR):
            dr[i] = r[1]-r[0]
        r  = np.delete(r-dr[0]/2,0)      # Shift the vector by -dr/2 and remove the negative entry
    elif CoordSystem == "SinhSpherical":
        if sinhA is None or sinhW is None:
            print("Error: SinhSpherical coordinates require initialization of both sinhA and sinhW")
            sys.exit(1)
        else:
            x  = np.linspace(0,1.0,NR+1)
            dx = 1.0/(NR+1)
            x  = np.delete(x-dx/2,0)      # Shift the vector by -dx/2 and remove the negative entry
            r  = sinhA * np.sinh( x/sinhW ) / np.sinh( 1.0/sinhW )
            dr = sinhA * np.cosh( x/sinhW ) / np.sinh( 1.0/sinhW ) * dx
    else:
        print("Error: Unknown coordinate system")
        sys.exit(1)

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
    main_diag[0] += 1 - dr[0]/r[0]

    # Update the last element of the main diagonal
    main_diag[NR-1] += - (2 * dr[NR-1] / r[NR-1])*(1 + dr[NR-1] / r[NR-1])

    # Set the upper diagonal, ignoring the last point in the r array
    upper_diag = np.zeros(NR)
    upper_diag[1:] = 1 + dr[:-1]/r[:-1]

    # Set the lower diagonal, start counting the r array at the second element
    lower_diag = np.zeros(NR)
    lower_diag[:-1] = 1 - dr[1:]/r[1:]

    # Change the last term in the lower diagonal to its correct value
    lower_diag[NR-2] = 2

    # Set the sparse matrix A by adding up the three diagonals
    # https://docs.scipy.org/doc/scipy/reference/generated/scipy.sparse.spdiags.html
    A = spdiags([main_diag,upper_diag,lower_diag],[0,1,-1],NR,NR)

    # Then compress the sparse matrix A column wise, so that SciPy can invert it later
    # https://docs.scipy.org/doc/scipy/reference/generated/scipy.sparse.csc_matrix.html
    A = csc_matrix(A)

    # Set up the right-hand side of the linear system: s
    s = np.zeros(NR)

    # Update the last entry of the vector s
    s[NR-1] = - (2 * dr[NR-1] / r[NR-1])*(1 + dr[NR-1] / r[NR-1])

    # Compress the vector s column-wise
    # https://docs.scipy.org/doc/scipy/reference/generated/scipy.sparse.csc_matrix.html
    s = csc_matrix(s)

    # Solve the sparse linear system using scipy
    # https://docs.scipy.org/doc/scipy-0.14.0/reference/generated/scipy.sparse.linalg.spsolve.html
    psi = spsolve(A, s.T)

    if lapse_condition == "Pre-collapsed":
        ID_alpha = psi**(-2)

        if sys.version_info[0] == 3:
            np.savetxt(outputname, list(zip( r, ID_sf, psi**4, ID_alpha )),fmt="%.15e")

        elif sys.version_info[0] == 2:
            np.savetxt(outputname, zip( r, ID_sf, psi**4, ID_alpha ),fmt="%.15e")

    elif lapse_condition == "Unity":
        ID_alpha = np.ones(NR)

        if sys.version_info[0] == 3:
            np.savetxt(outputname, list(zip( r, ID_sf, psi**4, ID_alpha )),fmt="%.15e")

        elif sys.version_info[0] == 2:
            np.savetxt(outputname, zip( r, ID_sf, psi**4, ID_alpha ),fmt="%.15e")

    else:
        print("Error: unknown lapse condition. Available options are: \"Pre-collapsed\" and \"Unity\"")
        return

    print("Generated the ADM initial data for the gravitational collapse \n" \
          "of a massless scalar field in %s coordinates.\n"%CoordSystem)
    print("Type of initial condition: Scalar field: \"Gaussian\" Shell\n"\
          "                         ADM quantities: Time-symmetric\n"\
          "                        Lapse condition: "+lapse_condition)
    print("Parameters: amplitude         = "+str(phi0)+",\n" \
          "            center            = "+str(r0)+",\n"   \
          "            width             = "+str(sigma)+",\n"   \
          "            domain size       = "+str(rmax)+",\n"   \
          "            number of points  = "+str(NR)+",\n"
          "            Initial data file = "+str(outputname)+".\n")

    with open(os.path.join(Ccodesdir,"ID_scalar_field_ADM_quantities.h"), "w") as file:
        file.write("""
// This function takes as input either (x,y,z) or (r,th,ph) and outputs
//   all ADM quantities in the Cartesian or Spherical basis, respectively.
void ID_scalar_field_ADM_quantities(
                     const REAL xyz_or_rthph[3],

                     const ID_inputs other_inputs,

                     REAL *gammaDD00,REAL *gammaDD01,REAL *gammaDD02,REAL *gammaDD11,REAL *gammaDD12,REAL *gammaDD22,
                     REAL *KDD00,REAL *KDD01,REAL *KDD02,REAL *KDD11,REAL *KDD12,REAL *KDD22,
                     REAL *alpha,
                     REAL *betaU0,REAL *betaU1,REAL *betaU2,
                     REAL *BU0,REAL *BU1,REAL *BU2) {

      const REAL r  = xyz_or_rthph[0];
      const REAL th = xyz_or_rthph[1];
      const REAL ph = xyz_or_rthph[2];

      REAL sf_star,psi4_star,alpha_star;

      scalar_field_interpolate_1D(r,
                                  other_inputs.interp_stencil_size,
                                  other_inputs.numlines_in_file,
                                  other_inputs.r_arr,
                                  other_inputs.sf_arr,
                                  other_inputs.psi4_arr,
                                  other_inputs.alpha_arr,
                                  &sf_star,&psi4_star,&alpha_star);

      // Update alpha
      *alpha = alpha_star;
      // \gamma_{rr} = psi^4
      *gammaDD00 = psi4_star;
      // \gamma_{thth} = psi^4 r^2
      *gammaDD11 = psi4_star*r*r;
      // \gamma_{phph} = psi^4 r^2 sin^2(th)
      *gammaDD22 = psi4_star*r*r*sin(th)*sin(th);

      // All other quantities ARE ZERO:
      *gammaDD01 = 0.0; *gammaDD02 = 0.0;
      /**/              *gammaDD12 = 0.0;

      *KDD00 = 0.0; *KDD01 = 0.0; *KDD02 = 0.0;
      /**/          *KDD11 = 0.0; *KDD12 = 0.0;
      /**/                        *KDD22 = 0.0;

      *betaU0 = 0.0; *betaU1 = 0.0; *betaU2 = 0.0;

      *BU0 = 0.0; *BU1 = 0.0; *BU2 = 0.0;
}\n""")

    print("Wrote to file "+os.path.join(Ccodesdir,"ID_scalar_field_ADM_quantities.h"))

    with open(os.path.join(Ccodesdir,"ID_scalar_field_spherical.h"), "w") as file:
        file.write("""
    void ID_scalarfield_spherical(
                     const REAL xyz_or_rthph[3],
                     const ID_inputs other_inputs,
                     REAL *sf, REAL *sfM) {

      const REAL r  = xyz_or_rthph[0];
      const REAL th = xyz_or_rthph[1];
      const REAL ph = xyz_or_rthph[2];

      REAL sf_star,psi4_star,alpha_star;

      scalar_field_interpolate_1D(r,
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

}\n""")

    print("Wrote to file "+os.path.join(Ccodesdir,"ID_scalar_field_spherical.h"))

    # Make sure that rfm.reference_metric() has been called.
    #    We'll need the variables it defines throughout this module.
    rfm.reference_metric()

    CoordType_in = "Spherical"
    pointer_to_ID_inputs=False

    sf, sfM = sp.symbols("sfSphorCart sfMSphorCart")

    r_th_ph_or_Cart_xyz_oID_xx = []
    if CoordType_in == "Spherical":
        r_th_ph_or_Cart_xyz_oID_xx = rfm.xxSph
    else:
        print("Error: Can only convert scalar field Spherical initial data to BSSN Curvilinear coords.")
        exit(1)

    with open(os.path.join(Ccodesdir,"ID_scalarfield_xx0xx1xx2_to_BSSN_xx0xx1xx2.h"), "w") as file:
        file.write("void ID_scalarfield_xx0xx1xx2_to_BSSN_xx0xx1xx2(const paramstruct *restrict params,const REAL xx0xx1xx2[3],")
        if pointer_to_ID_inputs == True:
            file.write("ID_inputs *other_inputs,")
        else:
            file.write("ID_inputs other_inputs,")
        file.write("""
                    REAL *restrict sf, REAL *restrict sfM ) {
#include \"set_Cparameters.h\"

      REAL sfSphorCart,sfMSphorCart;
      const REAL xx0 = xx0xx1xx2[0];
      const REAL xx1 = xx0xx1xx2[1];
      const REAL xx2 = xx0xx1xx2[2];
      REAL xyz_or_rthph[3];\n""")
    outCparams = "preindent=1,outCfileaccess=a,outCverbose=False,includebraces=False"
    outputC(r_th_ph_or_Cart_xyz_oID_xx[0:3], ["xyz_or_rthph[0]", "xyz_or_rthph[1]", "xyz_or_rthph[2]"],
            os.path.join(Ccodesdir,"ID_scalarfield_xx0xx1xx2_to_BSSN_xx0xx1xx2.h"), outCparams + ",CSE_enable=False")

    with open(os.path.join(Ccodesdir,"ID_scalarfield_xx0xx1xx2_to_BSSN_xx0xx1xx2.h"), "a") as file:
        file.write("""ID_scalarfield_spherical(xyz_or_rthph, other_inputs,
                      &sfSphorCart, &sfMSphorCart);
        // Next compute all rescaled BSSN curvilinear quantities:\n""")
    outCparams = "preindent=1,outCfileaccess=a,outCverbose=False,includebraces=False"
    outputC([sf, sfM], ["*sf", "*sfM"],
            os.path.join(Ccodesdir,"ID_scalarfield_xx0xx1xx2_to_BSSN_xx0xx1xx2.h"), params=outCparams)

    with open(os.path.join(Ccodesdir,"ID_scalarfield_xx0xx1xx2_to_BSSN_xx0xx1xx2.h"), "a") as file:
        file.write("}\n")

    # Driver
    with open(os.path.join(Ccodesdir,"ID_scalarfield.h"), "w") as file:
        file.write("""void ID_scalarfield(const paramstruct *restrict params,REAL *restrict xx[3],
                                          ID_inputs other_inputs,REAL *restrict in_gfs) {
#include \"set_Cparameters.h\"\n""")
        file.write(lp.loop(["i2", "i1", "i0"], ["0", "0", "0"],
                           ["Nxx_plus_2NGHOSTS2", "Nxx_plus_2NGHOSTS1", "Nxx_plus_2NGHOSTS0"],
                           ["1", "1", "1"], ["#pragma omp parallel for",
                                             "    const REAL xx2 = xx[2][i2];",
                                             "        const REAL xx1 = xx[1][i1];"], "",
                           """const REAL xx0 = xx[0][i0];
const int idx = IDX3S(i0,i1,i2);
const REAL xx0xx1xx2[3] = {xx0,xx1,xx2};
ID_scalarfield_xx0xx1xx2_to_BSSN_xx0xx1xx2(params,xx0xx1xx2,other_inputs,
                    &in_gfs[IDX4ptS(SFGF,idx)],&in_gfs[IDX4ptS(SFMGF,idx)]);\n}\n"""))
