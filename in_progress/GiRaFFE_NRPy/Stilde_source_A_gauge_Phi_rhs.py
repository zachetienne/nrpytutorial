# Step 1: Add NRPy's directory to the path
# https://stackoverflow.com/questions/16780014/import-file-from-parent-directory
import shutil, os, sys           # Standard Python modules for multiplatform OS-level functions
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

thismodule = "GiRaFFE_NRPy-Stilde_source_A_gauge_Phi_rhs"

# Here, we'll import the GRHD and GRFFE modules. 
import GRHD.equations as GRHD
import GRFFE.equations as GRFFE

def compute_StildeD_source_term(gammaDD,betaU,alpha,gammaDD_dD,betaU_dD,alpha_dD,sqrt4pi,ValenciavU,BU):
    # First, we'll use the GRHD module to generate the metric derivatives, four-velocity, and small b
    GRHD.u4U_in_terms_of_ValenciavU__rescale_ValenciavU_by_applying_speed_limit(alpha,betaU,gammaDD, ValenciavU)

    # Next sqrt(gamma)
    GRHD.compute_sqrtgammaDET(gammaDD)

    # Then compute g4DD_zerotimederiv_dD
    GRHD.compute_g4DD_zerotimederiv_dD(gammaDD,betaU,alpha, gammaDD_dD,betaU_dD,alpha_dD)

    # small b
    GRFFE.compute_smallb4U(gammaDD,betaU,alpha, GRHD.u4U_ito_ValenciavU, BU, sqrt4pi)
    GRFFE.compute_smallbsquared(gammaDD, betaU, alpha, GRFFE.smallb4U)
    # Electromagnetic stress-energy tensor
    GRFFE.compute_TEM4UU(gammaDD,betaU,alpha, GRFFE.smallb4U, GRFFE.smallbsquared,GRHD.u4U_ito_ValenciavU)

    # Add the source term to the RHS
    global Stilde_rhsD
    Stilde_rhsD = ixp.zerorank1(DIM=3)
    for i in range(3):
        for mu in range(4):
            for nu in range(4):
                # \frac{1}{2} \alpha \sqrt{\gamma} T^{\mu \nu}_{\rm EM} \partial_i g_{\mu \nu}
                Stilde_rhsD[i] += sp.Rational(1,2) * alpha * GRHD.sqrtgammaDET * \
                                  GRFFE.TEM4UU[mu][nu] * GRHD.g4DD_zerotimederiv_dD[mu][nu][i+1]

def compute_AD_gauge_term(gammaDD,betaU,alpha,psi6Phi,AD):
    GRHD.compute_sqrtgammaDET(gammaDD)
    Phi = psi6Phi/GRHD.sqrtgammaDET
    global AevolParen
    # \alpha \Phi - \beta^j A_j
    AevolParen = alpha * Phi
    for j in range(3):
        AevolParen += -betaU[j] * AD[j]

    # Take the gradient of the parenthetical and subtract it from the RHS
    AevolParen_dD = ixp.declarerank1("AevolParen_dD",DIM=3)
    global A_rhsD
    A_rhsD = ixp.declarerank1("A_rhsD",DIM=3)
    for i in range(3):
        A_rhsD[i] += -AevolParen_dD[i]
        
def compute_psi6Phi_rhs(gammaDD,betaU,alpha,AD,psi6Phi,xi_damping):
    GRHD.compute_sqrtgammaDET(gammaDD)
    gammaUU,unusedgammaDET = ixp.symm_matrix_inverter3x3(gammaDD)
    AU = ixp.zerorank1()
    # Raise the index on A in the usual way:
    for i in range(3):
        for j in range(3):
            AU[i] = gammaUU[i][j] * AD[j]
    
    global PhievolParenU
    PhievolParenU = ixp.zerorank1(DIM=3)
    
    for i in range(3):
        # \alpha\sqrt{\gamma}A^j - \beta^j [\sqrt{\gamma} \Phi]
        PhievolParenU[j] += alpha*GRHD.sqrtgammaDET*AU[j] - betaU[j]*psi6Phi
    
    # Tell NRPy+ to take the derivative numerically:
    PhievolParenU_dD = ixp.declarerank2("PhievolParenU_dD","nosym",DIM=3)
    # -\partial_j (\alpha\sqrt{\gamma}A^j - \beta^j [\sqrt{\gamma} \Phi]) - \xi \alpha [\sqrt{\gamma} \Phi]
    PhievolParen_innerProduct = sp.sympify(0)
    # Divergence of the parenthetical term
    for j in range(3):
        PhievolParen_innerProduct += PhievolParenU_dD[j][j]
    # Combine the divergence and the damping term
    global psi6Phi_rhs
    psi6Phi_rhs = -PhievolParen_innerProduct - xi_damping * alpha * psi6Phi

