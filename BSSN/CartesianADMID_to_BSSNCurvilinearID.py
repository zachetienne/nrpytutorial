# This module converts Cartesian ADM initial data to BSSN
# curvilinear initial data in terms of the variables used
# in BSSN_RHSs.py

# Author: Zachariah B. Etienne
#         zachetie **at** gmail **dot* com

import sympy as sp
import NRPy_param_funcs as par
from outputC import *
import indexedexp as ixp
import reference_metric as rfm
import BSSN.BSSN_RHSs as bssn # needed for parameters

def Convert_Cartesian_ADM_to_BSSN_curvilinear(Cartxyz, gammaCartDD, KCartDD, alphaCart, betaCartU, BCartU):
    # This routine convert the ADM variables
    # $$\left\{\gamma_{ij}, K_{ij}, \alpha, \beta^i\right\}$$
    # in Cartesian coordinates to the BSSN variables
    # $$\left\{\bar{\gamma}_{i j},\bar{A}_{i j},\phi, K, \bar{\Lambda}^{i}, \alpha, \beta^i, B^i\right\},$$
    # to the rescaled variables
    # $$\left\{h_{i j},a_{i j},\phi, K, \lambda^{i}, \alpha, \mathcal{V}^i, \mathcal{B}^i\right\}.$$

    # The ADM & BSSN formalisms only work in 3D; they are 3+1 decompositions of Einstein's equations.
    #    To implement axisymmetry or spherical symmetry, simply set all spatial derivatives in
    #    the relevant angular directions to zero; DO NOT SET DIM TO ANYTHING BUT 3.
    # Step 0: Set spatial dimension (must be 3 for BSSN)
    DIM = 3

    # First we have (Eqs. 2 and 3 of [Ruchlin *et al.*](https://arxiv.org/pdf/1712.07658.pdf)):

    # \bar{\gamma}_{i j} = \left(\frac{\bar{\gamma}}{\gamma}\right)^{1/3} \gamma_{ij}.

    # Brill-Linquist initial data are conformally flat, meaning that $\bar{\gamma}=\hat{\gamma}=1$
    # in Cartesian coordinates. So for these initial data in Cartesian coordinates we have

    # \bar{\gamma}_{i j} = \gamma^{-1/3} \gamma_{ij}.

    # Step 1: Convert ADM $\gamma_{ij}$ to BSSN $\bar{\gamma}_{ij}$
    gammaCartUU, gammaCartDET = ixp.symm_matrix_inverter3x3(gammaCartDD)
    gammaCartbarDD = ixp.zerorank2()
    for i in range(DIM):
        for j in range(DIM):
            gammaCartbarDD[i][j] = gammaCartDET ** (-sp.Rational(1, 3)) * gammaCartDD[i][j]

    # Second, we convert $K_{ij}$ to $\bar{A}_{ij}$ and $K$, where (Eq. 3 of [Baumgarte *et al.*](https://arxiv.org/pdf/1211.6632.pdf)):
    #
    # \begin{align}
    # K &= \gamma^{ij} K_{ij} \\
    # \bar{A}_{ij} &= \left(\frac{\bar{\gamma}}{\gamma}\right)^{1/3} \left(K_{ij} - \frac{1}{3} \gamma_{ij} K \right)\\
    # &= \gamma^{-1/3} \left(K_{ij} - \frac{1}{3} \gamma_{ij} K \right)
    # \end{align}

    # Step 2: Convert ADM $K_{ij}$ to $\bar{A}_{ij}$ and $K$,
    #          where (Eq. 3 of Baumgarte et al.: https://arxiv.org/pdf/1211.6632.pdf)
    trKCart = sp.sympify(0)
    for i in range(DIM):
        for j in range(DIM):
            trKCart += gammaCartUU[i][j] * KCartDD[i][j]

    ACartbarDD = ixp.zerorank2()
    for i in range(DIM):
        for j in range(DIM):
            ACartbarDD[i][j] = gammaCartDET ** (-sp.Rational(1, 3)) * (
                    KCartDD[i][j] - sp.Rational(1, 3) * gammaCartDD[i][j] * trKCart)

    # Third, we define $\bar{\Lambda}^i$ (Eqs. 4 and 5 of [Baumgarte *et al.*](https://arxiv.org/pdf/1211.6632.pdf)):
    #
    # $$
    # \bar{\Lambda}^i = \bar{\gamma}^{jk}\left(\bar{\Gamma}^i_{jk} - \hat{\Gamma}^i_{jk}\right).
    # $$
    #
    # In Cartesian coordinates, $\hat{\Gamma}^i_{jk}=0$, so we get in Cartesian coordinates
    # \begin{align}
    # \bar{\Lambda}^i &= \bar{\gamma}^{jk} \bar{\Gamma}^i_{jk} \\
    # &= \bar{\gamma}^{jk} \frac{1}{2} \bar{\gamma}^{il} \left(
    # \bar{\gamma}_{lj,k} + \bar{\gamma}_{lk,j} - \bar{\gamma}_{jk,l} \right)
    # \end{align}

    # Step 3: Define $\bar{\Lambda}^i$ from Eqs. 4 and 5 of Baumgarte et al.: https://arxiv.org/pdf/1211.6632.pdf
    LambdaCartbarU = ixp.zerorank1()
    # Need to compute \bar{\gamma}^{ij} from \bar{\gamma}_{ij}:
    gammaCartbarUU, gammaCartbarDET = ixp.symm_matrix_inverter3x3(gammaCartbarDD)

    for i in range(DIM):
        for j in range(DIM):
            for k in range(DIM):
                for l in range(DIM):
                    LambdaCartbarU[i] += sp.Rational(1, 2) * gammaCartbarUU[j][k] * gammaCartbarUU[i][l] * (
                            sp.diff(gammaCartbarUU[l][j], Cartxyz[k]) +
                            sp.diff(gammaCartbarUU[l][k], Cartxyz[j]) -
                            sp.diff(gammaCartbarUU[j][k], Cartxyz[l]))

    # All BSSN tensors and vectors are in Cartesian coordinates $x^i_{\rm Cart} = (x,y,z)$, but we need them in the curvilinear coordinate system $x^i_{\rm rfm}= ({\rm xx0},{\rm xx1},{\rm xx2})$ set by the "reference_metric::CoordSystem" variable. Empirically speaking, it is far easier to write $(x({\rm xx0},{\rm xx1},{\rm xx2}),y({\rm xx0},{\rm xx1},{\rm xx2}),z({\rm xx0},{\rm xx1},{\rm xx2}))$ than the inverse, so we will compute the Jacobian matrix
    #
    # $$
    # {\rm Jac\_dUCart\_dDrfmUD[i][j]} = \frac{\partial x^i_{\rm Cart}}{\partial x^j_{\rm rfm}},
    # $$
    #
    # via exact differentiation (courtesy SymPy), and the inverse Jacobian
    # $$
    # {\rm Jac\_dUrfm\_dDCartUD[i][j]} = \frac{\partial x^i_{\rm rfm}}{\partial x^j_{\rm Cart}},
    # $$
    #
    # using NRPy+'s ${\rm generic\_matrix\_inverter3x3()}$ function. In terms of these, the transformation of BSSN tensors from Cartesian to "reference_metric::CoordSystem" coordinates may be written:

    # \bar{\gamma}^{\rm rfm}_{ij} &=
    # \frac{\partial x^\ell_{\rm Cart}}{\partial x^i_{\rm rfm}}
    # \frac{\partial x^m_{\rm Cart}}{\partial x^j_{\rm rfm}} \bar{\gamma}^{\rm Cart}_{\ell m}\\
    # \bar{A}^{\rm rfm}_{ij} &=
    # \frac{\partial x^\ell_{\rm Cart}}{\partial x^i_{\rm rfm}}
    # \frac{\partial x^m_{\rm Cart}}{\partial x^j_{\rm rfm}} \bar{A}^{\rm Cart}_{\ell m}\\
    # \bar{\Lambda}^i_{\rm rfm} &= \frac{\partial x^i_{\rm rfm}}{\partial x^\ell_{\rm Cart}} \bar{\Lambda}^\ell_{\rm Cart} \\
    # \beta^i_{\rm rfm} &= \frac{\partial x^i_{\rm rfm}}{\partial x^\ell_{\rm Cart}} \beta^\ell_{\rm Cart}

    # Step 4: Convert BSSN tensors to basis specified by CoordSystem variable above.

    # trK and alpha are scalars, so no Jacobian transformations are necessary.
    trK = trKCart
    alpha = alphaCart

    Jac_dUCart_dDrfmUD = ixp.zerorank2()
    for i in range(DIM):
        for j in range(DIM):
            Jac_dUCart_dDrfmUD[i][j] = sp.diff(rfm.xxCart[i], rfm.xx[j])

    Jac_dUrfm_dDCartUD, dummyDET = ixp.generic_matrix_inverter3x3(Jac_dUCart_dDrfmUD)

    gammaDD = ixp.zerorank2()
    gammabarDD = ixp.zerorank2()
    AbarDD = ixp.zerorank2()
    LambdabarU = ixp.zerorank1()
    betaU = ixp.zerorank1()
    BU = ixp.zerorank1()
    for i in range(DIM):
        for j in range(DIM):
            LambdabarU[i] += Jac_dUrfm_dDCartUD[i][j] * LambdaCartbarU[j]
            betaU[i]      += Jac_dUrfm_dDCartUD[i][j] * betaCartU[j]
            BU[i]         += Jac_dUrfm_dDCartUD[i][j] * BCartU[j]
            for k in range(DIM):
                for l in range(DIM):
                    gammaDD[i][j]    += Jac_dUCart_dDrfmUD[k][i] * Jac_dUCart_dDrfmUD[l][j] * gammaCartDD[k][l]
                    gammabarDD[i][j] += Jac_dUCart_dDrfmUD[k][i] * Jac_dUCart_dDrfmUD[l][j] * gammaCartbarDD[k][l]
                    AbarDD[i][j]     += Jac_dUCart_dDrfmUD[k][i] * Jac_dUCart_dDrfmUD[l][j] * ACartbarDD[k][l]

    # Next we set the conformal factor variable $\texttt{cf}$, which is set by the "BSSN_RHSs::ConformalFactor" parameter. For example if "ConformalFactor" is set to "phi", we can use Eq. 3 of [Ruchlin *et al.*](https://arxiv.org/pdf/1712.07658.pdf), which in arbitrary coordinates is written:

    # \phi = \frac{1}{12} \log\left(\frac{\gamma}{\bar{\gamma}}\right).

    # Alternatively if "BSSN_RHSs::ConformalFactor" is set to "chi", then

    # \chi = e^{-4 \phi} = \exp\left(-4 \frac{1}{12} \left(\frac{\gamma}{\bar{\gamma}}\right)\right)
    # = \exp\left(-\frac{1}{3} \log\left(\frac{\gamma}{\bar{\gamma}}\right)\right) = \left(\frac{\gamma}{\bar{\gamma}}\right)^{-1/3}.

    # Finally if "BSSN_RHSs::ConformalFactor" is set to "W", then

    # W = e^{-2 \phi} = \exp\left(-2 \frac{1}{12} \log\left(\frac{\gamma}{\bar{\gamma}}\right)\right) =
    # \exp\left(-\frac{1}{6} \log\left(\frac{\gamma}{\bar{\gamma}}\right)\right) =
    # \left(\frac{\gamma}{\bar{\gamma}}\right)^{-1/6}.

    # Step 5: Set the conformal factor variable according to the parameter BSSN_RHSs::ConformalFactor

    cf = sp.sympify(0)

    gammaUU, gammaDET = ixp.symm_matrix_inverter3x3(gammaDD)
    gammabarUU, gammabarDET = ixp.symm_matrix_inverter3x3(gammabarDD)

    if par.parval_from_str("ConformalFactor") == "phi":
        cf = sp.Rational(1, 12) * sp.log(gammaDET / gammabarDET)
    elif par.parval_from_str("ConformalFactor") == "chi":
        cf = (gammaDET / gammabarDET) ** (-sp.Rational(1, 3))
    elif par.parval_from_str("ConformalFactor") == "W":
        cf = (gammaDET / gammabarDET) ** (-sp.Rational(1, 6))
    else:
        print("Error ConformalFactor type = \"" + par.parval_from_str("ConformalFactor") + "\" unknown.")
        exit(1)

    # Finally, we rescale tensorial quantities according to the prescription described in the [BSSN in curvilinear coordinates tutorial module](Tutorial-BSSNCurvilinear.ipynb) (also [Ruchlin *et al.*](https://arxiv.org/pdf/1712.07658.pdf)):
    # h_{ij} &= (\bar{\gamma}_{ij} - \hat{\gamma}_{ij})/\text{ReDD[i][j]}\\
    # a_{ij} &= \bar{A}_{ij}/\text{ReDD[i][j]}\\
    # \lambda^i &= \bar{\Lambda}^i/\text{ReU[i]}\\
    # \mathcal{V}^i &= \beta^i/\text{ReU[i]}\\
    # \mathcal{B}^i &= B^i/\text{ReU[i]}\\

    # Step 6: Rescale all tensorial quantities.
    
    if rfm.have_already_called_reference_metric_function == False:
        print("Error. Called Convert_Cartesian_ADM_to_BSSN_curvilinear() without")
        print("       first setting up reference metric, by calling rfm.reference_metric().")
        exit(1)
    
    hDD = ixp.zerorank2()
    aDD = ixp.zerorank2()
    lambdaU = ixp.zerorank1()
    vetU = ixp.zerorank1()
    betU = ixp.zerorank1()
    for i in range(DIM):
        lambdaU[i] = LambdabarU[i] / rfm.ReU[i]
        vetU[i] = betaU[i] / rfm.ReU[i]
        betU[i] = BU[i] / rfm.ReU[i]
        for j in range(DIM):
            hDD[i][j] = (gammabarDD[i][j] - rfm.ghatDD[i][j]) / rfm.ReDD[i][j]
            aDD[i][j] = AbarDD[i][j] / rfm.ReDD[i][j]
    # print(sp.mathematica_code(hDD[0][0]))
    return cf, hDD, lambdaU, aDD, trK, alpha, vetU, betU
