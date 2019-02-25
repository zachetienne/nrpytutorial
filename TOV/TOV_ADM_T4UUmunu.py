# This module sets up TOV initial data in terms of
# the variables used in BSSN_RHSs.py

# Author: Zachariah B. Etienne
#         zachetie **at** gmail **dot* com

# For full documentation, see Tutorial-ADM_and_T4UUmunu_Initial_Data-TOV.ipynb

# Step 4.1: Load NRPy+ modules and set basic parameters

import sympy as sp
import NRPy_param_funcs as par
import indexedexp as ixp
import grid as gri
from outputC import *
import BSSN.ADM_Numerical_Spherical_or_Cartesian_to_BSSNCurvilinear as AtoB
import BSSN.T4UUmunu_Numerical_Spherical_or_Cartesian_to_BSSNCurvilinear as TmunutoB

import reference_metric as rfm

def TOV_ADM_T4UUmunu(ComputeADMT4UUmunuGlobalsOnly = False):
    global Sph_r_th_ph,r,th,ph, gammaSphDD, KSphDD, alphaSph, betaSphU, BSphU, T4UU

    # All gridfunctions will be written in terms of spherical coordinates (r, th, ph):
    r,th,ph = sp.symbols('r th ph', real=True)

    thismodule = "TOV"

    DIM = 3
    par.set_parval_from_str("grid::DIM",DIM)

    # Input parameters read in from the TOV data file:
    rbar,expnu,exp4phi,P,rho = par.Cparameters("REAL", thismodule, ["rbar","expnu","exp4phi","P","rho"])

    # Step 4.2: Construct ADM quantities:

    # *** The physical spatial metric in spherical basis ***
    # In isotropic coordinates,
    #  gamma_{ij} = e^{4 phi} eta_{ij},
    # where eta is the flat-space 3-metric in spherical coordinates
    gammaSphDD = ixp.zerorank2()
    gammaSphDD[0][0] = exp4phi
    gammaSphDD[1][1] = exp4phi * rbar**2
    gammaSphDD[2][2] = exp4phi * rbar**2*sp.sin(th)**2

    # *** The extrinsic curvature in spherical basis ***
    # K_{ij} = 0 for the TOV solution
    KSphDD = ixp.zerorank2()

    # *** The lapse and shift in spherical basis ***
    # alpha = exp^{nu/2} for the TOV solution
    # \beta^i = 0 for the TOV solution
    alphaSph = sp.sqrt(expnu)
    betaSphU = ixp.zerorank1()
    BSphU = ixp.zerorank1()

    # Step 4.3: Construct T^{mu nu}:
    T4UU = ixp.zerorank2(DIM=4)

    # T^tt = e^(-nu) * rho
    T4UU[0][0] = rho / expnu
    # T^{ii} = P / gamma_{ii}
    for i in range(3):
        T4UU[i+1][i+1] = P / gammaSphDD[i][i]

    if ComputeADMT4UUmunuGlobalsOnly == True:
        return

    Sph_r_th_ph = [r, th, ph]
    cf, hDD, lambdaU, aDD, trK, alpha, vetU, betU = \
        AtoB.Convert_Spherical_or_Cartesian_ADM_to_BSSN_curvilinear("Spherical", Sph_r_th_ph,
                                                                    gammaSphDD, KSphDD, alphaSph, betaSphU, BSphU)
    sDD, sD, S, rho = \
        AtoB.Convert_Spherical_or_Cartesian_T4UUmunu_to_BSSN_curvilinear("Spherical", Sph_r_th_ph, T4UU)
    
    global returnfunction
    returnfunction = bIDf.BSSN_ID_function_string(cf, hDD, lambdaU, aDD, trK, alpha, vetU, betU)