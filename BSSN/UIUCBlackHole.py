# This module sets up UIUC Black Hole initial data in terms of
# the variables used in BSSN_RHSs.py

# Author: Zachariah B. Etienne
#         zachetie **at** gmail **dot* com

# ## This module sets up initial data for a merging black hole system in spherical coordinates. We can convert from spherical to any coordinate system defined in [reference_metric.py](../edit/reference_metric.py) (e.g., SinhSpherical, Cylindrical, Cartesian, etc.) using the [ADMSpherical-to-BSSNCurvilinear converter module](Tutorial-ADM_Initial_Data-Converting_ADMSpherical_to_BSSNCurvilinear)
# 
# ### NRPy+ Source Code for this module: [BSSN/BrillLindquist.py](../edit/BSSN/BrillLindquist.py)
# 
# <font color='green'>**All quantities have been validated against the [original SENR code](https://bitbucket.org/zach_etienne/nrpy).**</font>

# ### Here we set up UIUC Black Hole initial data ([Liu, Etienne, & Shapiro, PRD 80 121503, 2009](https://arxiv.org/abs/1001.4077)):
# 
# UIUC black holes have the advantage of finite coordinate radius in the maximal spin limit. It is therefore excellent for studying very highly spinning black holes. This module sets the UIUC black hole at the origin. 
# 
# **Inputs for initial data**:
# 
# * The black hole mass $M$.
# * The dimensionless spin parameter $\chi = a/M$
# 
# **Additional variables needed for spacetime evolution**:
# 
# * Desired coordinate system
# * Desired initial lapse $\alpha$ and shift $\beta^i$. We will choose our gauge conditions as $\alpha=1$ and $\beta^i=B^i=0$. $\alpha = \psi^{-2}$ will yield much better behavior, but the conformal factor $\psi$ depends on the desired *destination* coordinate system (which may not be spherical coordinates).

# Step P0: Load needed modules
import sympy as sp
import NRPy_param_funcs as par
from outputC import *
import indexedexp as ixp
import reference_metric as rfm
#import BSSN.SphericalADMID_to_BSSNCurvilinearID as ctob
#import BSSN.BSSN_ID_function_string as bIDf

def UIUCBlackHole():
    global r,th,ph, gammaSphDD, KSphDD, alphaSph, betaSphU, BSphU
    
    # All gridfunctions will be written in terms of spherical coordinates (r, th, ph):
    r,th,ph = sp.symbols('r th ph', real=True)

    thismodule = "UIUCBlackHole"

    # Step 0: Set spatial dimension (must be 3 for BSSN)
    DIM = 3
    par.set_parval_from_str("grid::DIM",DIM)
    
    # Step 1: Set psi, the conformal factor:

    # The UIUC initial data represent a Kerr black hole with mass M
    #  and dimensionless spin chi in UIUC quasi-isotropic coordinates,
    #   see https://arxiv.org/abs/1001.4077
    # Input parameters:
    M,chi = par.Cparameters("REAL", thismodule, ["M","chi"])

    # Auxiliary variables:
    a,rp,rm,rBL,SIG,DEL,AA = sp.symbols('a rp rm rBL SIG DEL AA', real=True)
    # Spin per unit mass
    a = M*chi

    # Boyer - Lindquist outer horizon
    rp = M + sp.sqrt(M**2 - a**2)
    # Boyer - Lindquist inner horizon
    rm = M - sp.sqrt(M**2 - a**2)

    # Boyer - Lindquist radius in terms of UIUC radius
    rBL = r*(1 + rp / (4*r))**2

    # UIUC definitions (Just below Eq 2)
    SIG = rBL**2 + a**2*sp.cos(th)**2
    DEL = rBL**2 - 2*M*rBL + a**2
    AA = (rBL**2 + a**2)**2 - DEL*a**2*sp.sin(th)**2

    # *** The ADM 3-metric in spherical basis ***
    gammaSphDD = ixp.zerorank2()
    # Declare the nonzero components of the 3-metric (Eq 2):
    gammaSphDD[0][0] = ((SIG*(r + rp/4)**2)/(r**3*(rBL - rm)))
    gammaSphDD[1][1] = SIG
    gammaSphDD[2][2] = AA/SIG*sp.sin(th)**2

    # *** The physical trace-free extrinsic curvature in spherical basis ***
    # Declare the nonzero components of the extrinsic curvature (Eqs 14-15):
    KSphDD     = ixp.zerorank2() # K_{ij} = 0 for these initial data
    KSphDD[0][2] = KSphDD[2][0] = (M*a*sp.sin(th)**2)/(SIG*sp.sqrt(AA*SIG))*\
                    (3*rBL**4 + 2*a**2*rBL**2 - a**4- a**2*(rBL**2 - a**2)*sp.sin(th)**2)*(1 + rp/(4*r))*1/sp.sqrt(r*(rBL - rm))
    KSphDD[1][2] = KSphDD[2][1] = -((2*a**3*M*rBL*sp.cos(th)*sp.sin(th)**3)/(SIG*sp.sqrt(AA*SIG)))*(r - rp/4)*sp.sqrt((rBL - rm)/r)

    alphaSph = sp.sympify(1)   # We generally choose alpha = 1/psi**2 (psi = BSSN conformal factor) for these initial data
    betaSphU = ixp.zerorank1() # We generally choose \beta^i = 0 for these initial data
    BSphU    = ixp.zerorank1() # We generally choose B^i = 0 for these initial data

    # Validated against original SENR: KSphDD[0][2], KSphDD[1][2], gammaSphDD[2][2], gammaSphDD[0][0], gammaSphDD[1][1]
    #print(sp.mathematica_code(gammaSphDD[1][1]))

    cf,hDD,lambdaU,aDD,trK,alpha,vetU,betU = \
        ctob.Convert_Spherical_ADM_to_BSSN_curvilinear(Sphxyz, gammaSphDD,KSphDD,alphaSph,betaSphU,BSphU)

    global returnfunction
    returnfunction = bIDf.BSSN_ID_function_string(cf,hDD,lambdaU,aDD,trK,alpha,vetU,betU)
