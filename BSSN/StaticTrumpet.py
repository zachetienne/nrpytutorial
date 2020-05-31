# This module sets up Static Trumpet initial data in terms of
# the variables used in BSSN_RHSs.py

# Authors: Terrence Pierre Jacques, terrencepierrej **at** gmail **dot** com
#          Zachariah B. Etienne, zachetie **at** gmail **dot** com
#          Ian Ruchlin

# ## This module sets up initial data for a static trumpet geometry in spherical coordinates. We can convert from spherical to any coordinate system defined in [reference_metric.py](../edit/reference_metric.py) (e.g., SinhSpherical, Cylindrical, Cartesian, etc.) using the [Exact ADM Spherical-to-BSSNCurvilinear converter module](Tutorial-ADM_Initial_Data-Converting_Exact_ADM_Spherical_or_Cartesian_to_BSSNCurvilinear.ipynb)
#
# ### NRPy+ Source Code for this module: [BSSN/BrillLindquist.py](../edit/BSSN/BrillLindquist.py)
#
# <font color='green'>**All quantities have been validated against the [original SENR code](https://bitbucket.org/zach_etienne/nrpy).**</font>

# ### Here we set up Static Trumpet initial data ([Dennison and Baumgarte, 2014](https://arxiv.org/abs/1403.5484)):
#
# Description of Static Trumpet geometry.
#
# **Inputs for initial data**:
#
# * The black hole mass $M$.
#
# **Additional variables needed for spacetime evolution**:
#
# * Desired coordinate system
# * Desired initial lapse $\alpha$ and shift $\beta^i$. We will choose our gauge conditions as $\alpha=1$ and $\beta^i=B^i=0$. $\alpha = \psi^{-2}$ will yield much better behavior, but the conformal factor $\psi$ depends on the desired *destination* coordinate system (which may not be spherical coordinates).

# Step 1: Initialize core Python/NRPy+ modules
import sympy as sp             # SymPy: The Python computer algebra package upon which NRPy+ depends
import NRPy_param_funcs as par # NRPy+: Parameter interface
import indexedexp as ixp       # NRPy+: Symbolic indexed expression (e.g., tensors, vectors, etc.) support
import BSSN.ADM_Exact_Spherical_or_Cartesian_to_BSSNCurvilinear as AtoB

thismodule = __name__

# Input parameters:
M = par.Cparameters("REAL", thismodule, ["M"], [1.0])

# ComputeADMGlobalsOnly == True will only set up the ADM global quantities.
#                       == False will perform the full ADM SphorCart->BSSN Curvi conversion
def StaticTrumpet(ComputeADMGlobalsOnly = False):
    global Sph_r_th_ph,r,th,ph, gammaSphDD, KSphDD, alphaSph, betaSphU, BSphU

    # All gridfunctions will be written in terms of spherical coordinates (r, th, ph):
    r,th,ph = sp.symbols('r th ph', real=True)

    # Step 0: Set spatial dimension (must be 3 for BSSN)
    DIM = 3
    par.set_parval_from_str("grid::DIM",DIM)

    # Step 1: Set psi, the StaticTrumpet conformal factor
    # Dennison and Baumgarte (2014) Eq. 13
    # https://arxiv.org/pdf/1403.5484.pdf

    # psi = sqrt{1 + M/r }
    psi0 = sp.sqrt(1 + M/r)

    # *** The physical spatial metric in spherical basis ***
    # Set the upper-triangle of the matrix...
    # Eq. 15
    # gamma_{ij} = psi^4 * eta_{ij}
    # eta_00 = 1, eta_11 = r^2, eta_22 = r^2 * sin^2 (theta)
    gammaSphDD = ixp.zerorank2()
    gammaSphDD[0][0] = psi0**4
    gammaSphDD[1][1] = psi0**4 * r**2
    gammaSphDD[2][2] = psi0**4 * r**2*sp.sin(th)**2
    # ... then apply symmetries to get the other components

    # *** The physical trace-free extrinsic curvature in spherical basis ***
    # Set the upper-triangle of the matrix...

    # Eq.19 and 20
    KSphDD = ixp.zerorank2()

    # K_{rr} = M / r^2
    KSphDD[0][0] = -M / r**2

    # K_{theta theta} = K_{phi phi} / sin^2 theta = M
    KSphDD[1][1] = M

    KSphDD[2][2] = M * sp.sin(th)**2
    # ... then apply symmetries to get the other components

    # Lapse function and shift vector
    # Eq. 15
    # alpha = r / (r+M)
    alphaSph = r / (r + M)

    betaSphU = ixp.zerorank1()
    # beta^r = Mr / (r + M)^2
    betaSphU[0] = M*r / (r + M)**2

    BSphU    = ixp.zerorank1()

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
