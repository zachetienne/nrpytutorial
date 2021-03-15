# Step 0: Add NRPy's directory to the path
# https://stackoverflow.com/questions/16780014/import-file-from-parent-directory
import os, sys                   # Standard Python modules for multiplatform OS-level functions
nrpy_dir_path = os.path.join("..")
if nrpy_dir_path not in sys.path:
    sys.path.append(nrpy_dir_path)

from outputC import outCfunction, outputC # NRPy+: Core C code output module
import sympy as sp               # SymPy: The Python computer algebra package upon which NRPy+ depends
import NRPy_param_funcs as par   # NRPy+: Parameter interface
import indexedexp as ixp         # NRPy+: Symbolic indexed expression (e.g., tensors, vectors, etc.) support

thismodule = "GiRaFFE_NRPy-Induction_Equation"

import GRHD.equations as GRHD
# import GRFFE.equations as GRFFE

# We'll write this as a function so that we can calculate the expressions on-demand for any choice of i
def find_cp_cm(lapse,shifti,gammaUUii):
    # Inputs:  u0,vi,lapse,shift,gammadet,gupii
    # Outputs: cplus,cminus

    # a = 1/(alpha^2)
    a = sp.sympify(1)/(lapse*lapse)
    # b = 2 beta^i / alpha^2
    b = sp.sympify(2) * shifti /(lapse*lapse)
    # c = -g^{ii} + (beta^i)^2 / alpha^2
    c = - gammaUUii + shifti*shifti/(lapse*lapse)

    # Now, we are free to solve the quadratic equation as usual. We take care to avoid passing a
    # negative value to the sqrt function.
    detm = b*b - sp.sympify(4)*a*c

    import Min_Max_and_Piecewise_Expressions as noif
    detm = sp.sqrt(noif.max_noif(sp.sympify(0),detm))
    global cplus,cminus
    cplus  = sp.Rational(1,2)*(-b/a + detm/a)
    cminus = sp.Rational(1,2)*(-b/a - detm/a)

# We'll write this as a function, and call it within HLLE_solver, below.
def find_cmax_cmin(field_comp,gamma_faceDD,beta_faceU,alpha_face):
    # Inputs:  flux direction field_comp, Inverse metric gamma_faceUU, shift beta_faceU,
    #          lapse alpha_face, metric determinant gammadet_face
    # Outputs: maximum and minimum characteristic speeds cmax and cmin
    # First, we need to find the characteristic speeds on each face
    gamma_faceUU,unusedgammaDET = ixp.generic_matrix_inverter3x3(gamma_faceDD)
    find_cp_cm(alpha_face,beta_faceU[field_comp],gamma_faceUU[field_comp][field_comp])
    cp = cplus
    cm = cminus

    # The following algorithms have been verified with random floats:

    global cmax,cmin
    # Now, we need to set cmax to the larger of cpr,cpl, and 0

    import Min_Max_and_Piecewise_Expressions as noif
    cmax = noif.max_noif(cp,sp.sympify(0))

    # And then, set cmin to the smaller of cmr,cml, and 0
    cmin = -noif.min_noif(cm,sp.sympify(0))

def calculate_flux_and_state_for_Induction(field_comp,flux_dirn, gammaDD,betaU,alpha,ValenciavU,BU):
    # Define Levi-Civita symbol
    def define_LeviCivitaSymbol_rank3(DIM=-1):
        if DIM == -1:
            DIM = par.parval_from_str("DIM")

        LeviCivitaSymbol = ixp.zerorank3()

        for i in range(DIM):
            for j in range(DIM):
                for k in range(DIM):
                    # From https://codegolf.stackexchange.com/questions/160359/levi-civita-symbol :
                    LeviCivitaSymbol[i][j][k] = (i - j) * (j - k) * (k - i) * sp.Rational(1,2)
        return LeviCivitaSymbol
    GRHD.compute_sqrtgammaDET(gammaDD)
    # Here, we import the Levi-Civita tensor and compute the tensor with lower indices
    LeviCivitaDDD = define_LeviCivitaSymbol_rank3()
    for i in range(3):
        for j in range(3):
            for k in range(3):
                LeviCivitaDDD[i][j][k] *= GRHD.sqrtgammaDET

    global U,F
    # Flux F = \epsilon_{ijk} v^j B^k
    F = sp.sympify(0)
    for j in range(3):
        for k in range(3):
            F += LeviCivitaDDD[field_comp][j][k] * (alpha*ValenciavU[j]-betaU[j]) * BU[k]
    # U = B^i
    U = BU[flux_dirn]

def HLLE_solver(cmax, cmin, Fr, Fl, Ur, Ul):
    # This solves the Riemann problem for the flux of E_i in one direction

    # F^HLL = (c_\min f_R + c_\max f_L - c_\min c_\max ( st_j_r - st_j_l )) / (c_\min + c_\max)
#     return (- cmin*cmax*(Ur-Ul) )/(cmax + cmin)
    return (cmin*Fr + cmax*Fl - cmin*cmax*(Ur-Ul) )/(cmax + cmin)

def calculate_E_i_flux(flux_dirn,alpha_face=None,gamma_faceDD=None,beta_faceU=None,\
                       Valenciav_rU=None,B_rU=None,Valenciav_lU=None,B_lU=None):
    global E_fluxD
    E_fluxD = ixp.zerorank1()
    for field_comp in range(3):
        find_cmax_cmin(field_comp,gamma_faceDD,beta_faceU,alpha_face)
        calculate_flux_and_state_for_Induction(field_comp,flux_dirn, gamma_faceDD,beta_faceU,alpha_face,\
                                               Valenciav_rU,B_rU)
        Fr = F
        Ur = U
        calculate_flux_and_state_for_Induction(field_comp,flux_dirn, gamma_faceDD,beta_faceU,alpha_face,\
                                               Valenciav_lU,B_lU)
        Fl = F
        Ul = U
        E_fluxD[field_comp] += HLLE_solver(cmax, cmin, Fr, Fl, Ur, Ul)

def generate_Afield_flux_function_files(out_dir,subdir,alpha_face,gamma_faceDD,beta_faceU,\
                                        Valenciav_rU,B_rU,Valenciav_lU,B_lU,inputs_provided=True):
    if not inputs_provided:
        # declare all variables
        alpha_face = sp.symbols(alpha_face)
        beta_faceU = ixp.declarerank1("beta_faceU")
        gamma_faceDD = ixp.declarerank2("gamma_faceDD","sym01")
        Valenciav_rU = ixp.declarerank1("Valenciav_rU")
        B_rU = ixp.declarerank1("B_rU")
        Valenciav_lU = ixp.declarerank1("Valenciav_lU")
        B_lU = ixp.declarerank1("B_lU")

    Memory_Read = """const double alpha_face = auxevol_gfs[IDX4S(ALPHA_FACEGF, i0,i1,i2)];
    const double gamma_faceDD00 = auxevol_gfs[IDX4S(GAMMA_FACEDD00GF, i0,i1,i2)];
    const double gamma_faceDD01 = auxevol_gfs[IDX4S(GAMMA_FACEDD01GF, i0,i1,i2)];
    const double gamma_faceDD02 = auxevol_gfs[IDX4S(GAMMA_FACEDD02GF, i0,i1,i2)];
    const double gamma_faceDD11 = auxevol_gfs[IDX4S(GAMMA_FACEDD11GF, i0,i1,i2)];
    const double gamma_faceDD12 = auxevol_gfs[IDX4S(GAMMA_FACEDD12GF, i0,i1,i2)];
    const double gamma_faceDD22 = auxevol_gfs[IDX4S(GAMMA_FACEDD22GF, i0,i1,i2)];
    const double beta_faceU0 = auxevol_gfs[IDX4S(BETA_FACEU0GF, i0,i1,i2)];
    const double beta_faceU1 = auxevol_gfs[IDX4S(BETA_FACEU1GF, i0,i1,i2)];
    const double beta_faceU2 = auxevol_gfs[IDX4S(BETA_FACEU2GF, i0,i1,i2)];
    const double Valenciav_rU0 = auxevol_gfs[IDX4S(VALENCIAV_RU0GF, i0,i1,i2)];
    const double Valenciav_rU1 = auxevol_gfs[IDX4S(VALENCIAV_RU1GF, i0,i1,i2)];
    const double Valenciav_rU2 = auxevol_gfs[IDX4S(VALENCIAV_RU2GF, i0,i1,i2)];
    const double B_rU0 = auxevol_gfs[IDX4S(B_RU0GF, i0,i1,i2)];
    const double B_rU1 = auxevol_gfs[IDX4S(B_RU1GF, i0,i1,i2)];
    const double B_rU2 = auxevol_gfs[IDX4S(B_RU2GF, i0,i1,i2)];
    const double Valenciav_lU0 = auxevol_gfs[IDX4S(VALENCIAV_LU0GF, i0,i1,i2)];
    const double Valenciav_lU1 = auxevol_gfs[IDX4S(VALENCIAV_LU1GF, i0,i1,i2)];
    const double Valenciav_lU2 = auxevol_gfs[IDX4S(VALENCIAV_LU2GF, i0,i1,i2)];
    const double B_lU0 = auxevol_gfs[IDX4S(B_LU0GF, i0,i1,i2)];
    const double B_lU1 = auxevol_gfs[IDX4S(B_LU1GF, i0,i1,i2)];
    const double B_lU2 = auxevol_gfs[IDX4S(B_LU2GF, i0,i1,i2)];
    REAL A_rhsD0 = 0; REAL A_rhsD1 = 0; REAL A_rhsD2 = 0;
    """
    Memory_Write = """rhs_gfs[IDX4S(AD0GF,i0,i1,i2)] += A_rhsD0;
    rhs_gfs[IDX4S(AD1GF,i0,i1,i2)] += A_rhsD1;
    rhs_gfs[IDX4S(AD2GF,i0,i1,i2)] += A_rhsD2;
    """

    indices = ["i0","i1","i2"]
    indicesp1 = ["i0+1","i1+1","i2+1"]

    for flux_dirn in range(3):
        calculate_E_i_flux(flux_dirn,alpha_face,gamma_faceDD,beta_faceU,\
                              Valenciav_rU,B_rU,Valenciav_lU,B_lU)

        E_field_to_print = [\
                            sp.Rational(1,4)*E_fluxD[(flux_dirn+1)%3],
                            sp.Rational(1,4)*E_fluxD[(flux_dirn+2)%3],
                           ]
        E_field_names = [\
                         "A_rhsD"+str((flux_dirn+1)%3),
                         "A_rhsD"+str((flux_dirn+2)%3),
                        ]

        desc = "Calculate the electric flux on the left face in direction " + str(flux_dirn) + "."
        name = "calculate_E_field_D" + str(flux_dirn) + "_right"
        outCfunction(
            outfile  = os.path.join(out_dir,subdir,name+".h"), desc=desc, name=name,
            params   ="const paramstruct *params,const REAL *auxevol_gfs,REAL *rhs_gfs",
            body     =  Memory_Read \
                       +outputC(E_field_to_print,E_field_names,"returnstring",params="outCverbose=False").replace("IDX4","IDX4S")\
                       +Memory_Write,
            loopopts ="InteriorPoints",
            rel_path_to_Cparams=os.path.join("../"))

        desc = "Calculate the electric flux on the left face in direction " + str(flux_dirn) + "."
        name = "calculate_E_field_D" + str(flux_dirn) + "_left"
        outCfunction(
            outfile  = os.path.join(out_dir,subdir,name+".h"), desc=desc, name=name,
            params   ="const paramstruct *params,const REAL *auxevol_gfs,REAL *rhs_gfs",
            body     =  Memory_Read.replace(indices[flux_dirn],indicesp1[flux_dirn]) \
                       +outputC(E_field_to_print,E_field_names,"returnstring",params="outCverbose=False").replace("IDX4","IDX4S")\
                       +Memory_Write,
            loopopts ="InteriorPoints",
            rel_path_to_Cparams=os.path.join("../"))

