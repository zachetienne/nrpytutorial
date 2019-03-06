# This module sets up Shifted Kerr-Schild initial data in terms of
# the variables used in BSSN_RHSs.py

# Authors: George Vopal, gvopal **at** gmail **dot** com
#          Zachariah B. Etienne, zachetie **at** gmail **dot** com


# ### NRPy+ Source Code for this module: [BSSN/ShiftedKerrSchild.py](../edit/BSSN/BrillLindquist.py)
# 
# WARNING: This module has not yet undergone code testing.

# **Inputs for initial data**:
# 
# * The black hole mass, M
# * The black hole spin parameter, a
 

# Step P0: Load needed modules
import sympy as sp
import NRPy_param_funcs as par
from outputC import *
import indexedexp as ixp
import reference_metric as rfm
import BSSN.ADM_Exact_Spherical_or_Cartesian_to_BSSNCurvilinear as AtoB
import BSSN.BSSN_ID_function_string as bIDf

# ComputeADMGlobalsOnly == True will only set up the ADM global quantities. 
#                       == False will perform the full ADM SphorCart->BSSN Curvi conversion
def ShiftedKerrSchild(ComputeADMGlobalsOnly = False):
    global Sph_r_th_ph,r,th,ph, rho2, gammaSphDD, KSphDD, alphaSph, betaSphU, BSphU
    
    # All gridfunctions will be written in terms of spherical coordinates (r, th, ph):
    r,th,ph = sp.symbols('r th ph', real=True)

    thismodule = "ShiftedKerrSchild"

    DIM = 3
    par.set_parval_from_str("grid::DIM",DIM)

    # Input parameters:
    M, a = par.Cparameters("REAL", thismodule, ["M", "a"])

    # Auxiliary variables:
    rho2 = sp.symbols('rho2', real=True)
    
    #rho^2 = r^2 + a^2*cos^2(theta)
    rho2 = r*r + a*a*sp.cos(th)**2

    # alpha = 1/sqrt{1 + Mr/rho^2}
    alphaSph = 1/(sp.sqrt(1 + 2*M*r/rho2))

    # Initialize the shift vector, \beta^i, to zero.
    betaSphU = ixp.zerorank1()
    #beta^r = alpha^2*2Mr/rho^2
    betaSphU[0] = alphaSph*alphaSph*2*M*r/rho2

    # Time derivative of shift vector beta^i, B^i, is zero.
    BSphU = ixp.zerorank1()

    # Initialize \gamma_{ij} to zero.
    gammaSphDD = ixp.zerorank2()
    
    #gamma_{rr} = 1 +2Mr/rho^2
    gammaSphDD[0][0] = 1 + 2*M*r/rho2

    # gamma_{r phi} = -a*gammaDD{r r}*sin^2(theta)
    gammaSphDD[0][2] = -a*gammaSphDD[0][0]*sp.sin(th)**2

    #gamma_{thetatheta} = rho^2
    gammaSphDD[1][1] = rho2

    # gamma_{phi phi} = (r^2 + a^2 + 2Mr/rho^2*a^2*sin^2(theta))*sin^2(theta)
    gammaSphDD[2][2] = (r*r + a*a + 2*M*r*a*a*sp.sin(th)**2/rho2)*sp.sin(th)**2
    
    # *** Define Useful Quantities A, B, D ***
    # A = (a^2*cos^2(2theta) + a^2 + 2r^2)
    A = (a*a*sp.cos(2*th) + a*a + 2*r*r)

    # B = A + 4Mr
    B = A + 4*M*r

    # D = \sqrt(2Mr/(a^2cos^2(theta) + r^2) + 1)
    D = sp.sqrt(2*M*r/(a*a*sp.cos(th)**2 + r*r) + 1)
                
    
    # *** The extrinsic curvature in spherical polar coordinates ***
    
    # Establish the 3x3 zero-matrix
    KSphDD = ixp.zerorank2()

    # *** Fill in the nonzero components ***
    # *** This will create an upper-triangular matrix ***
    # K_{r r} = D(A+2Mr)/(A^2*B)[4M(a^2*cos(2theta) + a^2 - 2r^2)]
    KSphDD[0][0] = D*(A+2*M*r)/(A*A*B)*(4*M*(a*a*sp.cos(2*th)+a*a-2*r*r))


    # K_{r theta} = D/(AB)[8a^2*Mr*sin(theta)cos(theta)]
    KSphDD[0][1] = D/(A*B)*(8*a*a*M*r*sp.sin(th)*sp.cos(th))

    # K_{r phi} = D/A^2[-2aMsin^2(theta)(a^2cos(2theta)+a^2-2r^2)]
    KSphDD[0][2] = D/(A*A)*(-2*a*M*sp.sin(th)**2*(a*a*sp.cos(2*th)+a*a-2*r*r))

    # K_{theta theta} = D/B[4Mr^2]
    KSphDD[1][1] = D/B*(4*M*r*r)

    # K_{theta phi} = D/(AB)*(-8*a^3*Mr*sin^3(theta)cos(theta))
    KSphDD[1][2] = D/(A*B)*(-8*a**3*M*r*sp.sin(th)**3*sp.cos(th))

    # K_{phi phi} = D/(A^2*B)[2Mr*sin^2(theta)(a^4(M+3r)
    #   +4a^2r^2(2r-M)+4a^2r*cos(2theta)(a^2+r(M+2r))+8r^5)]
    KSphDD[2][2] = D/(A*A*B)*(2*M*r*sp.sin(th)**2*(a**4*(r-M)*sp.cos(4*th)\
                            + a**4*(M+3*r)+4*a*a*r*r*(2*r-M)\
                            + 4*a*a*r*sp.cos(2*th)*(a*a + r*(M + 2*r)) + 8*r**5))           
    
    
    if ComputeADMGlobalsOnly == True:
        return
    
    # Validated against original SENR: 
    #print(sp.mathematica_code(gammaSphDD[1][1]))

    Sph_r_th_ph = [r,th,ph]
    cf,hDD,lambdaU,aDD,trK,alpha,vetU,betU = \
        AtoB.Convert_Spherical_or_Cartesian_ADM_to_BSSN_curvilinear("Spherical", Sph_r_th_ph, 
                                                                    gammaSphDD,KSphDD,alphaSph,betaSphU,BSphU)

    global returnfunction
    returnfunction = bIDf.BSSN_ID_function_string(cf,hDD,lambdaU,aDD,trK,alpha,vetU,betU)
