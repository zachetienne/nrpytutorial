# This module performs the conversion between ADM
# spacetime variables in Spherical or Cartesian coordinates
# given as *exact* SymPy expressions (i.e., direct functions
# of r,th,ph or x,y,z), to rescaled BSSN-in-curvilinear
# coordinate quantities, as defined in BSSN_RHSs.py

# Author: Zachariah B. Etienne
#         zachetie **at** gmail **dot* com

import sympy as sp
import NRPy_param_funcs as par
from outputC import *
import indexedexp as ixp
import reference_metric as rfm
import BSSN.BSSN_quantities as Bq # The EvolvedConformalFactor_cf parameter is used below

def Convert_Spherical_or_Cartesian_ADM_to_BSSN_curvilinear(CoordType_in, Sph_r_th_ph_or_Cart_xyz,
        gammaDD_inSphorCart, KDD_inSphorCart, alpha_inSphorCart, betaU_inSphorCart, BU_inSphorCart):
    # This routine converts the ADM variables
    #    $$\left\{\gamma_{ij}, K_{ij}, \alpha, \beta^i\right\}$$
    #    in Spherical or Cartesian basis+coordinates, first to the BSSN variables
    #    in the chosen reference_metric::CoordSystem coordinate system+basis:
    # $$\left\{\bar{\gamma}_{i j},\bar{A}_{i j},\phi, K, \bar{\Lambda}^{i}, \alpha, \beta^i, B^i\right\},$$
    #    and then to the rescaled variables:
    # $$\left\{h_{i j},a_{i j},\phi, K, \lambda^{i}, \alpha, \mathcal{V}^i, \mathcal{B}^i\right\}.$$

    # The ADM & BSSN formalisms only work in 3D; they are 3+1 decompositions of Einstein's equations.
    #    To implement axisymmetry or spherical symmetry, simply set all spatial derivatives in
    #    the relevant angular directions to zero; DO NOT SET DIM TO ANYTHING BUT 3.
    # Step 0: Set spatial dimension (must be 3 for BSSN)
    DIM = 3

    # Step 0: Copy gammaSphDD_in to gammaSphDD, KSphDD_in to KSphDD, etc.
    #    This ensures that the input arrays are not modified below;
    #    modifying them would result in unexpected effects outside
    #    this function.
    alphaSphorCart = alpha_inSphorCart
    betaSphorCartU   = ixp.zerorank1()
    BSphorCartU      = ixp.zerorank1()
    gammaSphorCartDD = ixp.zerorank2()
    KSphorCartDD     = ixp.zerorank2()
    for i in range(DIM):
        betaSphorCartU[i] = betaU_inSphorCart[i]
        BSphorCartU[i]    = BU_inSphorCart[i]
        for j in range(DIM):
            gammaSphorCartDD[i][j] = gammaDD_inSphorCart[i][j]
            KSphorCartDD[i][j]     = KDD_inSphorCart[i][j]

    # Make sure that rfm.reference_metric() has been called.
    #    We'll need the variables it defines throughout this module.
    if rfm.have_already_called_reference_metric_function == False:
        print("Error. Called Convert_Spherical_ADM_to_BSSN_curvilinear() without")
        print("       first setting up reference metric, by calling rfm.reference_metric().")
        exit(1)
    
    # Step 1: All input quantities are in terms of r,th,ph or x,y,z. We want them in terms 
    #         of xx0,xx1,xx2, so here we call sympify_integers__replace_rthph() to replace
    #         r,th,ph or x,y,z, respectively, with the appropriate functions of xx0,xx1,xx2
    #         as defined for this particular reference metric in reference_metric.py's 
    #         xxSph[] or xxCart[], respectively:
    # Note that substitution only works when the variable is not an integer. Hence the 
    #         if isinstance(...,...) stuff:
    def sympify_integers__replace_rthph_or_Cartxyz(obj, rthph_or_xyz, rthph_or_xyz_of_xx):
        if isinstance(obj, int):
            return sp.sympify(obj)
        else:
            return obj.subs(rthph_or_xyz[0], rthph_or_xyz_of_xx[0]).\
                subs(rthph_or_xyz[1], rthph_or_xyz_of_xx[1]).\
                subs(rthph_or_xyz[2], rthph_or_xyz_of_xx[2])

    r_th_ph_or_Cart_xyz_of_xx = []
    if CoordType_in == "Spherical":
        r_th_ph_or_Cart_xyz_of_xx = rfm.xxSph
    elif CoordType_in == "Cartesian":
        r_th_ph_or_Cart_xyz_of_xx = rfm.xxCart
    else:
        print("Error: Can only convert ADM Cartesian or Spherical initial data to BSSN Curvilinear coords.")
        exit(1)

    alphaSphorCart = sympify_integers__replace_rthph_or_Cartxyz(
        alphaSphorCart, Sph_r_th_ph_or_Cart_xyz, r_th_ph_or_Cart_xyz_of_xx)
    for i in range(DIM):
        betaSphorCartU[i] = sympify_integers__replace_rthph_or_Cartxyz(
            betaSphorCartU[i], Sph_r_th_ph_or_Cart_xyz, r_th_ph_or_Cart_xyz_of_xx)
        BSphorCartU[i]    = sympify_integers__replace_rthph_or_Cartxyz(
            BSphorCartU[i],    Sph_r_th_ph_or_Cart_xyz, r_th_ph_or_Cart_xyz_of_xx)
        for j in range(DIM):
            gammaSphorCartDD[i][j] = sympify_integers__replace_rthph_or_Cartxyz(
                gammaSphorCartDD[i][j], Sph_r_th_ph_or_Cart_xyz, r_th_ph_or_Cart_xyz_of_xx)
            KSphorCartDD[i][j]     = sympify_integers__replace_rthph_or_Cartxyz(
                KSphorCartDD[i][j],     Sph_r_th_ph_or_Cart_xyz, r_th_ph_or_Cart_xyz_of_xx)

    # Step 2: All ADM initial data quantities are now functions of xx0,xx1,xx2, but
    #         they are still in the Spherical or Cartesian basis. We can now directly apply
    #         Jacobian transformations to get them in the correct xx0,xx1,xx2 basis:

    # alpha is a scalar, so no Jacobian transformation is necessary.
    alpha = alphaSphorCart

    Jac_dUSphorCart_dDrfmUD = ixp.zerorank2()
    for i in range(DIM):
        for j in range(DIM):
            Jac_dUSphorCart_dDrfmUD[i][j] = sp.diff(r_th_ph_or_Cart_xyz_of_xx[i],rfm.xx[j])

    Jac_dUrfm_dDSphorCartUD, dummyDET = ixp.generic_matrix_inverter3x3(Jac_dUSphorCart_dDrfmUD)

    betaU   = ixp.zerorank1()
    BU      = ixp.zerorank1()
    gammaDD = ixp.zerorank2()
    KDD     = ixp.zerorank2()
    for i in range(DIM):
        for j in range(DIM):
            betaU[i] += Jac_dUrfm_dDSphorCartUD[i][j] * betaSphorCartU[j]
            BU[i]    += Jac_dUrfm_dDSphorCartUD[i][j] * BSphorCartU[j]
            for k in range(DIM):
                for l in range(DIM):
                    gammaDD[i][j] += Jac_dUSphorCart_dDrfmUD[k][i]*Jac_dUSphorCart_dDrfmUD[l][j] * gammaSphorCartDD[k][l]
                    KDD[i][j]     += Jac_dUSphorCart_dDrfmUD[k][i]*Jac_dUSphorCart_dDrfmUD[l][j] *     KSphorCartDD[k][l]


    # Step 3: All ADM quantities were input into this function in the Spherical or Cartesian
    #         basis, as functions of r,th,ph or x,y,z, respectively. In Steps 1 and 2 above,
    #         we converted them to the xx0,xx1,xx2 basis, and as functions of xx0,xx1,xx2.
    #         Here we convert ADM quantities to their BSSN Curvilinear counterparts:
    
    # Step 3.1: Convert ADM $\gamma_{ij}$ to BSSN $\bar{\gamma}_{ij}$:
    #   We have (Eqs. 2 and 3 of [Ruchlin *et al.*](https://arxiv.org/pdf/1712.07658.pdf)):
    # \bar{\gamma}_{i j} = \left(\frac{\bar{\gamma}}{\gamma}\right)^{1/3} \gamma_{ij}.
    gammaUU, gammaDET = ixp.symm_matrix_inverter3x3(gammaDD)
    gammabarDD = ixp.zerorank2()
    for i in range(DIM):
        for j in range(DIM):
            gammabarDD[i][j] = (rfm.detgammahat/gammaDET)**(sp.Rational(1,3))*gammaDD[i][j]

    # Step 3.2: Convert the extrinsic curvature $K_{ij}$ to the trace-free extrinsic 
    #           curvature $\bar{A}_{ij}$, plus the trace of the extrinsic curvature $K$, 
    #           where (Eq. 3 of [Baumgarte *et al.*](https://arxiv.org/pdf/1211.6632.pdf)):

    # K = \gamma^{ij} K_{ij}, and
    # \bar{A}_{ij} &= \left(\frac{\bar{\gamma}}{\gamma}\right)^{1/3} \left(K_{ij} - \frac{1}{3} \gamma_{ij} K \right)
    trK = sp.sympify(0)
    for i in range(DIM):
        for j in range(DIM):
            trK += gammaUU[i][j]*KDD[i][j]

    AbarDD = ixp.zerorank2()
    for i in range(DIM):
        for j in range(DIM):
            AbarDD[i][j] = (rfm.detgammahat/gammaDET)**(sp.Rational(1,3))*(KDD[i][j] - sp.Rational(1,3)*gammaDD[i][j]*trK)


    # Step 3.3: Define $\bar{\Lambda}^i$ (Eqs. 4 and 5 of [Baumgarte *et al.*](https://arxiv.org/pdf/1211.6632.pdf)):

    # \bar{\Lambda}^i = \bar{\gamma}^{jk}\left(\bar{\Gamma}^i_{jk} - \hat{\Gamma}^i_{jk}\right).
    gammabarUU, gammabarDET = ixp.symm_matrix_inverter3x3(gammabarDD)

    # First compute \bar{\Gamma}^i_{jk}:
    GammabarUDD = ixp.zerorank3()
    for i in range(DIM):
        for j in range(DIM):
            for k in range(DIM):
                for l in range(DIM):
                    GammabarUDD[i][j][k] += sp.Rational(1,2)*gammabarUU[i][l]*( sp.diff(gammabarDD[l][j],rfm.xx[k]) +
                                                                                sp.diff(gammabarDD[l][k],rfm.xx[j]) -
                                                                                sp.diff(gammabarDD[j][k],rfm.xx[l]) )
    # Next evaluate \bar{\Lambda}^i, based on GammabarUDD above and GammahatUDD 
    #       (from the reference metric):
    LambdabarU = ixp.zerorank1()
    for i in range(DIM):
        for j in range(DIM):
            for k in range(DIM):
                LambdabarU[i] += gammabarUU[j][k] * (GammabarUDD[i][j][k] - rfm.GammahatUDD[i][j][k])


    # Step 3.4: Set the conformal factor variable $\texttt{cf}$, which is set 
    #           by the "BSSN_quantities::EvolvedConformalFactor_cf" parameter. For example if 
    #           "EvolvedConformalFactor_cf" is set to "phi", we can use Eq. 3 of 
    #           [Ruchlin *et al.*](https://arxiv.org/pdf/1712.07658.pdf), 
    #           which in arbitrary coordinates is written:

    # \phi = \frac{1}{12} \log\left(\frac{\gamma}{\bar{\gamma}}\right).

    # Alternatively if "BSSN_quantities::EvolvedConformalFactor_cf" is set to "chi", then

    # \chi = e^{-4 \phi} = \exp\left(-4 \frac{1}{12} \left(\frac{\gamma}{\bar{\gamma}}\right)\right)
    #      = \exp\left(-\frac{1}{3} \log\left(\frac{\gamma}{\bar{\gamma}}\right)\right) = \left(\frac{\gamma}{\bar{\gamma}}\right)^{-1/3}.
    #
    # Finally if "BSSN_quantities::EvolvedConformalFactor_cf" is set to "W", then

    # W = e^{-2 \phi} = \exp\left(-2 \frac{1}{12} \log\left(\frac{\gamma}{\bar{\gamma}}\right)\right) =
    # \exp\left(-\frac{1}{6} \log\left(\frac{\gamma}{\bar{\gamma}}\right)\right) =
    # \left(\frac{\gamma}{\bar{\gamma}}\right)^{-1/6}.

    cf = sp.sympify(0)

    if par.parval_from_str("EvolvedConformalFactor_cf") == "phi":
        cf = sp.Rational(1,12)*sp.log(gammaDET/gammabarDET)
    elif par.parval_from_str("EvolvedConformalFactor_cf") == "chi":
        cf = (gammaDET/gammabarDET)**(-sp.Rational(1,3))
    elif par.parval_from_str("EvolvedConformalFactor_cf") == "W":
        cf = (gammaDET/gammabarDET)**(-sp.Rational(1,6))
    else:
        print("Error EvolvedConformalFactor_cf type = \""+par.parval_from_str("EvolvedConformalFactor_cf")+"\" unknown.")
        exit(1)


    # Step 4: Rescale tensorial quantities according to the prescription described in 
    #         the [BSSN in curvilinear coordinates tutorial module](Tutorial-BSSNCurvilinear.ipynb) 
    #         (also [Ruchlin *et al.*](https://arxiv.org/pdf/1712.07658.pdf)):
    #
    # h_{ij} &= (\bar{\gamma}_{ij} - \hat{\gamma}_{ij})/\text{ReDD[i][j]}\\
    # a_{ij} &= \bar{A}_{ij}/\text{ReDD[i][j]}\\
    # \lambda^i &= \bar{\Lambda}^i/\text{ReU[i]}\\
    # \mathcal{V}^i &= \beta^i/\text{ReU[i]}\\
    # \mathcal{B}^i &= B^i/\text{ReU[i]}\\
    hDD     = ixp.zerorank2()
    aDD     = ixp.zerorank2()
    lambdaU = ixp.zerorank1()
    vetU    = ixp.zerorank1()
    betU    = ixp.zerorank1()
    for i in range(DIM):
        lambdaU[i] = LambdabarU[i] / rfm.ReU[i]
        vetU[i]    =      betaU[i] / rfm.ReU[i]
        betU[i]    =         BU[i] / rfm.ReU[i]
        for j in range(DIM):
            hDD[i][j] = (gammabarDD[i][j] - rfm.ghatDD[i][j]) / rfm.ReDD[i][j]
            aDD[i][j] =                          AbarDD[i][j] / rfm.ReDD[i][j]
    #print(sp.mathematica_code(hDD[0][0]))

    # Step 5: Return the BSSN Curvilinear variables in the desired xx0,xx1,xx2
    #         basis, and as functions of the consistent xx0,xx1,xx2 coordinates.
    return cf, hDD, lambdaU, aDD, trK, alpha, vetU, betU