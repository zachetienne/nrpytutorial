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


# Step 1: Initialize core Python/NRPy+ modules
import sympy as sp             # SymPy: The Python computer algebra package upon which NRPy+ depends
import NRPy_param_funcs as par # NRPy+: Parameter interface
import indexedexp as ixp       # NRPy+: Symbolic indexed expression (e.g., tensors, vectors, etc.) support
import BSSN.ADM_Exact_Spherical_or_Cartesian_to_BSSNCurvilinear as AtoB

thismodule = __name__

# Input parameters:
M, a, r0 = par.Cparameters("REAL", thismodule,
                           ["M", "a", "r0"],
                           [1.0, 0.9, 1.0])

# ComputeADMGlobalsOnly == True will only set up the ADM global quantities.
#                       == False will perform the full ADM SphorCart->BSSN Curvi conversion
def ShiftedKerrSchild(ComputeADMGlobalsOnly = False):
    global Sph_r_th_ph,r,th,ph, rho2, gammaSphDD, KSphDD, alphaSph, betaSphU, BSphU

    # All gridfunctions will be written in terms of spherical coordinates (r, th, ph):
    r,th,ph = sp.symbols('r th ph', real=True)

    DIM = 3
    par.set_parval_from_str("grid::DIM",DIM)

    # Auxiliary variables:
    rho2 = sp.symbols('rho2', real=True)

    # Step 1: Define rho^2, alpha, beta^(r_{KS}), beta^(theta), beta^(phi), gamma_{r_{KS}theta}, gamma_{theta\phi}

    # r_{KS} = r + r0
    rKS = r+r0

    # rho^2 = rKS^2 + a^2*cos^2(theta)
    rho2 = rKS*rKS + a*a*sp.cos(th)**2

    # alpha = 1/sqrt{1 + M*rKS/rho^2}
    alphaSph = 1/(sp.sqrt(1 + 2*M*rKS/rho2))

    # Initialize the shift vector, \beta^i, to zero.
    betaSphU = ixp.zerorank1()
    # beta^r = alpha^2*2Mr/rho^2
    betaSphU[0] = alphaSph*alphaSph*2*M*rKS/rho2

    # Time derivative of shift vector beta^i, B^i, is zero.
    BSphU = ixp.zerorank1()

    # Step 2: Define and construct nonzero components gamma_{r_{KS}r_{KS}}$, gamma_{r_{KS}phi},
    #         gamma_{thetatheta}, gamma_{phiphi}

    # Initialize \gamma_{ij} to zero.
    gammaSphDD = ixp.zerorank2()

    # gammaDD{rKS rKS} = 1 +2M*rKS/rho^2
    gammaSphDD[0][0] = 1 + 2*M*rKS/rho2

    # gammaDD{rKS phi} = -a*gammaDD{r r}*sin^2(theta)
    gammaSphDD[0][2] = gammaSphDD[2][0] = -a*gammaSphDD[0][0]*sp.sin(th)**2

    # gammaDD{theta theta} = rho^2
    gammaSphDD[1][1] = rho2

    # gammaDD{phi phi} = (rKS^2 + a^2 + 2Mr/rho^2*a^2*sin^2(theta))*sin^2(theta)
    gammaSphDD[2][2] = (rKS*rKS + a*a + 2*M*rKS*a*a*sp.sin(th)**2/rho2)*sp.sin(th)**2

    # Step 3: Define useful quantities A, B, C
    # A = (a^2*cos^2(2theta) + a^2 + 2r^2)
    A = (a*a*sp.cos(2*th) + a*a + 2*rKS*rKS)

    # B = A + 4M*rKS
    B = A + 4*M*rKS

    # D = \sqrt(2M*rKS/(a^2cos^2(theta) + rKS^2) + 1)
    D = sp.sqrt(2*M*rKS/(a*a*sp.cos(th)**2 + rKS*rKS) + 1)


    # Step 4: Define the extrinsic curvature in spherical polar coordinates

    # Establish the 3x3 zero-matrix
    KSphDD = ixp.zerorank2()

    # *** Fill in the nonzero components ***
    # *** This will create an upper-triangular matrix ***
    # K_{r r} = D(A+2Mr)/(A^2*B)[4M(a^2*cos(2theta) + a^2 - 2r^2)]
    KSphDD[0][0] = D*(A+2*M*rKS)/(A*A*B)*(4*M*(a*a*sp.cos(2*th)+a*a-2*rKS*rKS))

    # K_{r theta} = D/(AB)[8a^2*Mr*sin(theta)cos(theta)]
    KSphDD[0][1] = KSphDD[1][0] = D/(A*B)*(8*a*a*M*rKS*sp.sin(th)*sp.cos(th))

    # K_{r phi} = D/A^2[-2aMsin^2(theta)(a^2cos(2theta)+a^2-2r^2)]
    KSphDD[0][2] = KSphDD[2][0] =  D/(A*A)*(-2*a*M*sp.sin(th)**2*(a*a*sp.cos(2*th)+a*a-2*rKS*rKS))

    # K_{theta theta} = D/B[4Mr^2]
    KSphDD[1][1] = D/B*(4*M*rKS*rKS)

    # K_{theta phi} = D/(AB)*(-8*a^3*Mr*sin^3(theta)cos(theta))
    KSphDD[1][2] = KSphDD[2][1] = D/(A*B)*(-8*a**3*M*rKS*sp.sin(th)**3*sp.cos(th))

    # K_{phi phi} = D/(A^2*B)[2Mr*sin^2(theta)(a^4(M+3r)
    #   +4a^2r^2(2r-M)+4a^2r*cos(2theta)(a^2+r(M+2r))+8r^5)]
    KSphDD[2][2] = D/(A*A*B)*(2*M*rKS*sp.sin(th)**2*(a**4*(rKS-M)*sp.cos(4*th)\
                            + a**4*(M+3*rKS)+4*a*a*rKS*rKS*(2*rKS-M)\
                            + 4*a*a*rKS*sp.cos(2*th)*(a*a + rKS*(M + 2*rKS)) + 8*rKS**5))


    if ComputeADMGlobalsOnly == True:
        return

    # Validated against original SENR:
    #print(sp.mathematica_code(gammaSphDD[1][1]))

    Sph_r_th_ph = [r,th,ph]
    cf,hDD,lambdaU,aDD,trK,alpha,vetU,betU = \
        AtoB.Convert_Spherical_or_Cartesian_ADM_to_BSSN_curvilinear("Spherical", Sph_r_th_ph,
                                                                    gammaSphDD,KSphDD,alphaSph,betaSphU,BSphU)

    import BSSN.BSSN_ID_function_string as bIDf
    # Generates initial_data() C function & stores to outC_function_dict["initial_data"]
    bIDf.BSSN_ID_function_string(cf, hDD, lambdaU, aDD, trK, alpha, vetU, betU)
