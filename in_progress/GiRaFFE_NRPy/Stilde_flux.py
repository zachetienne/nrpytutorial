from outputC import outCfunction, outputC # NRPy+: Core C code output module
import sympy as sp               # SymPy: The Python computer algebra package upon which NRPy+ depends
import NRPy_param_funcs as par   # NRPy+: Parameter interface
import grid as gri               # NRPy+: Functions having to do with numerical grids
import indexedexp as ixp         # NRPy+: Symbolic indexed expression (e.g., tensors, vectors, etc.) support
import os, sys           # Standard Python modules for multiplatform OS-level functions

par.set_parval_from_str("grid::DIM", 3)
DIM = par.parval_from_str("grid::DIM")

thismodule = __name__

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
def find_cmax_cmin(flux_dirn,gamma_faceDD,beta_faceU,alpha_face):
    # Inputs:  flux direction flux_dirn, Inverse metric gamma_faceUU, shift beta_faceU,
    #          lapse alpha_face, metric determinant gammadet_face
    # Outputs: maximum and minimum characteristic speeds cmax and cmin
    # First, we need to find the characteristic speeds on each face
    gamma_faceUU,unusedgammaDET = ixp.generic_matrix_inverter3x3(gamma_faceDD)
    find_cp_cm(alpha_face,beta_faceU[flux_dirn],gamma_faceUU[flux_dirn][flux_dirn])
    cp = cplus
    cm = cminus

    # The following algorithms have been verified with random floats:

    global cmax,cmin
    # Now, we need to set cmax to the larger of cpr,cpl, and 0

    import Min_Max_and_Piecewise_Expressions as noif
    cmax = noif.max_noif(cp,sp.sympify(0))

    # And then, set cmin to the smaller of cmr,cml, and 0
    cmin = -noif.min_noif(cm,sp.sympify(0))
# We'll rewrite this assuming that we've passed the entire reconstructed
# gridfunctions. You could also do this with only one point, but then you'd
# need to declare everything as a Cparam in NRPy+
nrpy_dir_path = os.path.join("..")
if nrpy_dir_path not in sys.path:
    sys.path.append(nrpy_dir_path)

# We'll rewrite this assuming that we've passed the entire reconstructed
# gridfunctions. You could also do this with only one point, but then you'd
# need to declare everything as a Cparam in NRPy+

import GRHD.equations as GRHD
import GRFFE.equations as GRFFE

def calculate_GRFFE_Tmunu_and_contractions(flux_dirn, mom_comp, gammaDD,betaU,alpha,ValenciavU,BU,sqrt4pi):
    GRHD.compute_sqrtgammaDET(gammaDD)

    GRHD.u4U_in_terms_of_ValenciavU__rescale_ValenciavU_by_applying_speed_limit(alpha, betaU, gammaDD, ValenciavU)
    GRFFE.compute_smallb4U(gammaDD, betaU, alpha, GRHD.u4U_ito_ValenciavU, BU, sqrt4pi)
    GRFFE.compute_smallbsquared(gammaDD, betaU, alpha, GRFFE.smallb4U)

    GRFFE.compute_TEM4UU(gammaDD, betaU, alpha, GRFFE.smallb4U, GRFFE.smallbsquared, GRHD.u4U_ito_ValenciavU)
    GRFFE.compute_TEM4UD(gammaDD, betaU, alpha, GRFFE.TEM4UU)

    # Compute conservative variables in terms of primitive variables
    GRHD.compute_S_tildeD(alpha, GRHD.sqrtgammaDET, GRFFE.TEM4UD)

    global U,F
    # Flux F = alpha*sqrt{gamma}*T^i_j
    F = alpha*GRHD.sqrtgammaDET*GRFFE.TEM4UD[flux_dirn+1][mom_comp+1]
    # U = alpha*sqrt{gamma}*T^0_j = Stilde_j
    U = GRHD.S_tildeD[mom_comp]


def HLLE_solver(cmax, cmin, Fr, Fl, Ur, Ul):
    # This solves the Riemann problem for the mom_comp component of the momentum
    # flux StildeD in the flux_dirn direction.

    # st_j_flux = (c_\min f_R + c_\max f_L - c_\min c_\max ( st_j_r - st_j_l )) / (c_\min + c_\max)
    return (cmin*Fr + cmax*Fl - cmin*cmax*(Ur-Ul) )/(cmax + cmin)

def calculate_Stilde_flux(flux_dirn,alpha_face,gamma_faceDD,beta_faceU,\
                          Valenciav_rU,B_rU,Valenciav_lU,B_lU,sqrt4pi):

    find_cmax_cmin(flux_dirn,gamma_faceDD,beta_faceU,alpha_face)

    global Stilde_fluxD
    Stilde_fluxD = ixp.zerorank3()
    for mom_comp in range(3):
        calculate_GRFFE_Tmunu_and_contractions(flux_dirn, mom_comp, gamma_faceDD,beta_faceU,alpha_face,\
                                               Valenciav_rU,B_rU,sqrt4pi)
        Fr = F
        Ur = U
        calculate_GRFFE_Tmunu_and_contractions(flux_dirn, mom_comp, gamma_faceDD,beta_faceU,alpha_face,\
                                               Valenciav_lU,B_lU,sqrt4pi)
        Fl = F
        Ul = U
        Stilde_fluxD[mom_comp] = HLLE_solver(cmax, cmin, Fr, Fl, Ur, Ul)

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
REAL Stilde_fluxD0 = 0; REAL Stilde_fluxD1 = 0; REAL Stilde_fluxD2 = 0;
"""
Memory_Write = """rhs_gfs[IDX4S(STILDED0GF, i0, i1, i2)] += invdx0*Stilde_fluxD0;
rhs_gfs[IDX4S(STILDED1GF, i0, i1, i2)] += invdx0*Stilde_fluxD1;
rhs_gfs[IDX4S(STILDED2GF, i0, i1, i2)] += invdx0*Stilde_fluxD2;
"""

indices = ["i0","i1","i2"]
indicesp1 = ["i0+1","i1+1","i2+1"]
assignment = "+="
assignmentp1 = "-="
invdx = ["invdx0","invdx1","invdx2"]

def generate_C_code_for_Stilde_flux(out_dir,inputs_provided = False, alpha_face=None, gamma_faceDD=None, beta_faceU=None,
                                    Valenciav_rU=None, B_rU=None, Valenciav_lU=None, B_lU=None, sqrt4pi=None,
                                    outCparams = "outCverbose=False,CSE_sorting=none"):
    if not inputs_provided:
        # We will pass values of the gridfunction on the cell faces into the function. This requires us
        # to declare them as C parameters in NRPy+. We will denote this with the _face infix/suffix.
        alpha_face = gri.register_gridfunctions("AUXEVOL","alpha_face")
        gamma_faceDD = ixp.register_gridfunctions_for_single_rank2("AUXEVOL","gamma_faceDD","sym01")
        beta_faceU = ixp.register_gridfunctions_for_single_rank1("AUXEVOL","beta_faceU")

        # We'll need some more gridfunctions, now, to represent the reconstructions of BU and ValenciavU
        # on the right and left faces
        Valenciav_rU = ixp.register_gridfunctions_for_single_rank1("AUXEVOL","Valenciav_rU",DIM=3)
        B_rU = ixp.register_gridfunctions_for_single_rank1("AUXEVOL","B_rU",DIM=3)
        Valenciav_lU = ixp.register_gridfunctions_for_single_rank1("AUXEVOL","Valenciav_lU",DIM=3)
        B_lU = ixp.register_gridfunctions_for_single_rank1("AUXEVOL","B_lU",DIM=3)
        sqrt4pi = par.Cparameters("REAL",thismodule,"sqrt4pi","sqrt(4.0*M_PI)")

    for flux_dirn in range(3):
        calculate_Stilde_flux(flux_dirn,alpha_face,gamma_faceDD,beta_faceU,\
                              Valenciav_rU,B_rU,Valenciav_lU,B_lU,sqrt4pi)
        Stilde_flux_to_print = [\
                                Stilde_fluxD[0],\
                                Stilde_fluxD[1],\
                                Stilde_fluxD[2],\
                               ]
        Stilde_flux_names = [\
                             "Stilde_fluxD0",\
                             "Stilde_fluxD1",\
                             "Stilde_fluxD2",\
                            ]

        desc = "Compute the flux of all 3 components of tilde{S}_i on the right face in the " + str(flux_dirn) + "."
        name = "calculate_Stilde_flux_D" + str(flux_dirn) + "_right"
        outCfunction(
            outfile  = os.path.join(out_dir,name+".h"), desc=desc, name=name,
            params   ="const paramstruct *params,const REAL *auxevol_gfs,REAL *rhs_gfs",
            body     =  Memory_Read \
                       +outputC(Stilde_flux_to_print,Stilde_flux_names,"returnstring",params=outCparams).replace("IDX4","IDX4S")\
                       +Memory_Write.replace(invdx[0],invdx[flux_dirn]),
            loopopts ="InteriorPoints",
            rel_path_for_Cparams=os.path.join("../"))

        desc = "Compute the flux of all 3 components of tilde{S}_i on the left face in the " + str(flux_dirn) + "."
        name = "calculate_Stilde_flux_D" + str(flux_dirn) + "_left"
        outCfunction(
            outfile  = os.path.join(out_dir,name+".h"), desc=desc, name=name,
            params   ="const paramstruct *params,const REAL *auxevol_gfs,REAL *rhs_gfs",
            body     =  Memory_Read.replace(indices[flux_dirn],indicesp1[flux_dirn]) \
                       +outputC(Stilde_flux_to_print,Stilde_flux_names,"returnstring",params=outCparams).replace("IDX4","IDX4S")\
                       +Memory_Write.replace(invdx[0],invdx[flux_dirn]).replace(assignment,assignmentp1),
            loopopts ="InteriorPoints",
            rel_path_for_Cparams=os.path.join("../"))
