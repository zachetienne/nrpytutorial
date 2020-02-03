# Step 0: Add NRPy's directory to the path
# https://stackoverflow.com/questions/16780014/import-file-from-parent-directory
import os,sys
nrpy_dir_path = os.path.join("..")
if nrpy_dir_path not in sys.path:
    sys.path.append(nrpy_dir_path)

from outputC import *            # NRPy+: Core C code output module
import finite_difference as fin  # NRPy+: Finite difference C code generation module
import NRPy_param_funcs as par   # NRPy+: Parameter interface
import grid as gri               # NRPy+: Functions having to do with numerical grids
import loop as lp                # NRPy+: Generate C code loops
import indexedexp as ixp         # NRPy+: Symbolic indexed expression (e.g., tensors, vectors, etc.) support
import reference_metric as rfm   # NRPy+: Reference metric support
import cmdline_helper as cmd     # NRPy+: Multi-platform Python command-line interface
import shutil, os, sys           # Standard Python modules for multiplatform OS-level functions

thismodule = "GiRaFFE_NRPy-Induction_Equation"

import GRHD.equations as GRHD
import GRFFE.equations as GRFFE

# Import the Levi-Civita symbol and build the corresponding tensor.
# We already have a handy function to define the Levi-Civita symbol in WeylScalars
import WeylScal4NRPy.WeylScalars_Cartesian as weyl

# We'll write this as a function so that we can calculate the expressions on-demand for any choice of i
def find_cp_cm(lapse,shifti,gupii):
    # Inputs:  u0,vi,lapse,shift,gammadet,gupii
    # Outputs: cplus,cminus 
    
    # a = 1/(alpha^2)
    a = 1/(lapse*lapse)
    # b = 2 beta^i / alpha^2
    b = 2 * shifti /(lapse*lapse)
    # c = -g^{ii} + (beta^i)^2 / alpha^2
    c = - gupii + shifti*shifti/(lapse*lapse)
    
    # Now, we are free to solve the quadratic equation as usual. We take care to avoid passing a
    # negative value to the sqrt function.
    detm = b*b - 4*a*c
    detm = sp.sqrt(sp.Rational(1,2)*(detm + nrpyAbs(detm)))
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
    cpr = cplus
    cmr = cminus
    find_cp_cm(alpha_face,beta_faceU[flux_dirn],gamma_faceUU[flux_dirn][flux_dirn])
    cpl = cplus
    cml = cminus
    
    # The following algorithms have been verified with random floats:
    
    global cmax,cmin
    # Now, we need to set cmax to the larger of cpr,cpl, and 0
    cmax = sp.Rational(1,2)*(cpr+cpl+nrpyAbs(cpr-cpl))
    cmax = sp.Rational(1,2)*(cmax+nrpyAbs(cmax))
    
    # And then, set cmin to the smaller of cmr,cml, and 0
    cmin =  sp.Rational(1,2)*(cmr+cml-nrpyAbs(cmr-cml))
    cmin = -sp.Rational(1,2)*(cmin-nrpyAbs(cmin))

def calculate_flux_and_state_for_Induction(flux_dirn, gammaDD,betaU,alpha,ValenciavU,BU):
    GRHD.compute_sqrtgammaDET(gammaDD)
    # Here, we import the Levi-Civita tensor and compute the tensor with lower indices
    LeviCivitaDDD = weyl.define_LeviCivitaSymbol_rank3()
    for i in range(3):
        for j in range(3):
            for k in range(3):
                LCijk = LeviCivitaDDD[i][j][k]
                LeviCivitaDDD[i][j][k] = LCijk * GRHD.sqrtgammaDET
    #             LeviCivitaUUU[i][j][k] = LCijk / sp.sqrt(gammadet)

    global U,F
    # Flux F = \epsilon_{ijk} v^j B^k
    F = sp.sympify(0)
    for j in range(3):
        for k in range(3):
            F += LeviCivitaDDD[flux_dirn][j][k] * (alpha*ValenciavU[j]-betaU[j]) * BU[k]
    # U = B^i
    U = BU[flux_dirn]
    
def HLLE_solver(cmax, cmin, Fr, Fl, Ur, Ul): 
    # This solves the Riemann problem for the flux of E_i in one direction
    
    # F^HLL = (c_\min f_R + c_\max f_L - c_\min c_\max ( st_j_r - st_j_l )) / (c_\min + c_\max)
    return (cmin*Fr + cmax*Fl - cmin*cmax*(Ur-Ul) )/(cmax + cmin)

def calculate_E_i_flux(inputs_provided=True,alpha_face=None,gamma_faceDD=None,beta_faceU=None,\
                       Valenciav_rU=None,B_rU=None,Valenciav_lU=None,B_lU=None):
    if not inputs_provided:
        # declare all variables
        alpha_face = sp.symbols(alpha_face)
        beta_faceU = ixp.declarerank1("beta_faceU")
        gamma_faceDD = ixp.declarerank2("gamma_faceDD","sym01")
        Valenciav_rU = ixp.declarerank1("Valenciav_rU")
        B_rU = ixp.declarerank1("B_rU")
        Valenciav_lU = ixp.declarerank1("Valenciav_lU")
        B_lU = ixp.declarerank1("B_lU")
    global E_fluxD
    E_fluxD = ixp.zerorank1()
    for flux_dirn in range(3):
        find_cmax_cmin(flux_dirn,gamma_faceDD,beta_faceU,alpha_face)
        calculate_flux_and_state_for_Induction(flux_dirn, gamma_faceDD,beta_faceU,alpha_face,\
                                               Valenciav_rU,B_rU)
        Fr = F
        Ur = U
        calculate_flux_and_state_for_Induction(flux_dirn, gamma_faceDD,beta_faceU,alpha_face,\
                                               Valenciav_lU,B_lU)
        Fl = F
        Ul = U
        E_fluxD[flux_dirn] += HLLE_solver(cmax, cmin, Fr, Fl, Ur, Ul)
