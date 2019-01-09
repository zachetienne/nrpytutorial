#!/usr/bin/env python
# coding: utf-8

# # Converting Spherical ADM Initial Data to BSSN Curvilinear Initial Data
# 
# Here we use NRPy+ to first convert the ADM variables
# 
# $$\left\{\gamma_{ij}, K_{ij}, \alpha, \beta^i\right\}$$
# 
# to the BSSN variables
# 
# $$\left\{\bar{\gamma}_{i j},\bar{A}_{i j},\phi, K, \bar{\Lambda}^{i}, \alpha, \beta^i, B^i\right\},$$ 
# 
# and the BSSN variables to the rescaled BSSNCurvilinear variables (as defined in [the BSSN Curvilinear tutorial](Tutorial-BSSNCurvilinear.ipynb)):
# 
# $$\left\{h_{i j},a_{i j},\phi, K, \lambda^{i}, \alpha, \mathcal{V}^i, \mathcal{B}^i\right\}.$$ 

# First we have (Eqs. 2 and 3 of [Ruchlin *et al.*](https://arxiv.org/pdf/1712.07658.pdf)):
# $$
# \bar{\gamma}_{i j} = \left(\frac{\bar{\gamma}}{\gamma}\right)^{1/3} \gamma_{ij}.
# $$
# 
# We choose $\bar{\gamma}=\hat{\gamma}$, which in spherical coordinates implies $\bar{\gamma}=\hat{\gamma}=r^4 \sin^2\theta$. Instead of manually setting the determinant, we instead call rfm.reference_metric(), which automatically sets up the metric based on the chosen reference_metric::CoordSystem parameter.

# Step P1: Import needed modules
import sympy as sp
import NRPy_param_funcs as par
from outputC import *
import indexedexp as ixp
import reference_metric as rfm
import BSSN.UIUCBlackHole as uibh
import BSSN.BSSN_RHSs as bssnrhs # The ConformalFactor parameter is used below


# Temporarily set the source reference metric to Spherical, so 
#     we can use some of the CoordSystem==Spherical reference_metric 
#     functionality.
# First backup the desired destination coordinate system for
#     BSSNCurvilinear output:
CoordSystem_dest = par.parval_from_str("reference_metric::CoordSystem")

par.set_parval_from_str("reference_metric::CoordSystem","Spherical")
rfm.reference_metric()

# Import UIUC Black Hole initial data
uibh.UIUCBlackHole()
Sph_r_th_ph = [uibh.r,uibh.th,uibh.ph]
gammaSphDD  = uibh.gammaSphDD
KSphDD   = uibh.KSphDD
alphaSph = uibh.alphaSph
betaSphU = uibh.betaSphU
BSphU    = uibh.BSphU


# The ADM & BSSN formalisms only work in 3D; they are 3+1 decompositions of Einstein's equations.
#    To implement axisymmetry or spherical symmetry, simply set all spatial derivatives in
#    the relevant angular directions to zero; DO NOT SET DIM TO ANYTHING BUT 3.
# Step 0: Set spatial dimension (must be 3 for BSSN)
DIM = 3

# Step 1: Convert ADM $\gamma_{ij}$ to BSSN $\bar{\gamma}_{ij}$
gammaSphUU, gammaSphDET = ixp.symm_matrix_inverter3x3(gammaSphDD)
gammaSphbarDD = ixp.zerorank2()
for i in range(DIM):
    for j in range(DIM):
        gammaSphbarDD[i][j] = (rfm.detgammahat/gammaSphDET)**(sp.Rational(1,3))*gammaSphDD[i][j]


# Second, we convert the extrinsic curvature $K_{ij}$ to the trace-free extrinsic curvature $\bar{A}_{ij}$, plus the trace of the extrinsic curvature $K$, where (Eq. 3 of [Baumgarte *et al.*](https://arxiv.org/pdf/1211.6632.pdf)):
# 
# \begin{align}
# K &= \gamma^{ij} K_{ij} \\
# \bar{A}_{ij} &= \left(\frac{\bar{\gamma}}{\gamma}\right)^{1/3} \left(K_{ij} - \frac{1}{3} \gamma_{ij} K \right)
# \end{align}

# Step 2: Convert ADM $K_{ij}$ to $\bar{A}_{ij}$ and $K$, 
#          where (Eq. 3 of Baumgarte et al.: https://arxiv.org/pdf/1211.6632.pdf)
trKSph = sp.sympify(0)
for i in range(DIM):
    for j in range(DIM):
        trKSph += gammaSphUU[i][j]*KSphDD[i][j]

ASphbarDD = ixp.zerorank2()
for i in range(DIM):
    for j in range(DIM):
        ASphbarDD[i][j] = (rfm.detgammahat/gammaSphDET)**(sp.Rational(1,3))*          (KSphDD[i][j] - sp.Rational(1,3)*gammaSphDD[i][j]*trKSph)


# Third, we define $\bar{\Lambda}^i$ (Eqs. 4 and 5 of [Baumgarte *et al.*](https://arxiv.org/pdf/1211.6632.pdf)):
# 
# $$
# \bar{\Lambda}^i = \bar{\gamma}^{jk}\left(\bar{\Gamma}^i_{jk} - \hat{\Gamma}^i_{jk}\right).
# $$

# All BSSN tensors and vectors are in Spherical coordinates $x^i_{\rm Sph} = (r,\theta,\phi)$, but we need them in the curvilinear coordinate system $x^i_{\rm rfm}= ({\rm xx0},{\rm xx1},{\rm xx2})$ set by the "reference_metric::CoordSystem" variable. Empirically speaking, it is far easier to write $(x({\rm xx0},{\rm xx1},{\rm xx2}),y({\rm xx0},{\rm xx1},{\rm xx2}),z({\rm xx0},{\rm xx1},{\rm xx2}))$ than the inverse, so we will compute the Jacobian matrix
# 
# $$
# {\rm Jac\_dUSph\_dDrfmUD[i][j]} = \frac{\partial x^i_{\rm Sph}}{\partial x^j_{\rm rfm}},
# $$
# 
# via exact differentiation (courtesy SymPy), and the inverse Jacobian
# $$
# {\rm Jac\_dUrfm\_dDSphUD[i][j]} = \frac{\partial x^i_{\rm rfm}}{\partial x^j_{\rm Sph}},
# $$
# 
# using NRPy+'s ${\rm generic\_matrix\_inverter3x3()}$ function. In terms of these, the transformation of BSSN tensors from Spherical to "reference_metric::CoordSystem" coordinates may be written:
# 
# \begin{align}
# \bar{\gamma}^{\rm rfm}_{ij} &= 
# \frac{\partial x^\ell_{\rm Sph}}{\partial x^i_{\rm rfm}}
# \frac{\partial x^m_{\rm Sph}}{\partial x^j_{\rm rfm}} \bar{\gamma}^{\rm Sph}_{\ell m}\\
# \bar{A}^{\rm rfm}_{ij} &= 
# \frac{\partial x^\ell_{\rm Sph}}{\partial x^i_{\rm rfm}}
# \frac{\partial x^m_{\rm Sph}}{\partial x^j_{\rm rfm}} \bar{A}^{\rm Sph}_{\ell m}\\
# \bar{\Lambda}^i_{\rm rfm} &= \frac{\partial x^i_{\rm rfm}}{\partial x^\ell_{\rm Sph}} \bar{\Lambda}^\ell_{\rm Sph} \\
# \beta^i_{\rm rfm} &= \frac{\partial x^i_{\rm rfm}}{\partial x^\ell_{\rm Sph}} \beta^\ell_{\rm Sph}
# \end{align}

# Step 3: Define $\bar{\Lambda}^i$ from Eqs. 4 and 5 of Baumgarte et al.: https://arxiv.org/pdf/1211.6632.pdf
LambdaSphbarU = ixp.zerorank1()
# Need to compute \bar{\gamma}^{ij} from \bar{\gamma}_{ij}:
gammaSphbarUU, gammaSphbarDET = ixp.symm_matrix_inverter3x3(gammaSphbarDD)

GammabarUDD = ixp.zerorank3()
for i in range(DIM):
    for j in range(DIM):
        for k in range(DIM):
            for l in range(DIM):
                GammabarUDD[i][j][k] +=                     sp.Rational(1,2)*gammaSphbarUU[i][l]*( sp.diff(gammaSphbarUU[l][j],Sph_r_th_ph[k]) + 
                                                           sp.diff(gammaSphbarUU[l][k],Sph_r_th_ph[j]) -
                                                           sp.diff(gammaSphbarUU[j][k],Sph_r_th_ph[l]) )

for i in range(DIM):
    for j in range(DIM):
        for k in range(DIM):
                LambdaSphbarU[i] += gammaSphbarUU[j][k]* (GammabarUDD[i][j][k] - rfm.GammahatUDD[i][j][k])


# Step 4: Convert BSSN tensors to destination basis specified by CoordSystem_dest variable above.

# Step 4.1: First restore reference_metric::CoordSystem to 
#           CoordSystem_dest, and call rfm.reference_metric(),
#           so that all hatted quantities are updated to
#           CoordSystem_dest.
par.set_parval_from_str("reference_metric::CoordSystem",CoordSystem_dest)
rfm.reference_metric()

# trK and alpha are scalars, so no Jacobian transformations are necessary.
trK   = trKSph
alpha = alphaSph

Jac_dUSph_dDrfmUD = ixp.zerorank2()
for i in range(DIM):
    for j in range(DIM):
        Jac_dUSph_dDrfmUD[i][j] = sp.diff(rfm.xxSph[i],rfm.xx[j])

Jac_dUrfm_dDSphUD, dummyDET = ixp.generic_matrix_inverter3x3(Jac_dUSph_dDrfmUD)

gammaDD    = ixp.zerorank2()
gammabarDD = ixp.zerorank2()
AbarDD     = ixp.zerorank2()
LambdabarU = ixp.zerorank1()
betaU      = ixp.zerorank1()
BU         = ixp.zerorank1()
for i in range(DIM):
    for j in range(DIM):
        LambdabarU[i] += Jac_dUrfm_dDSphUD[i][j]*LambdaSphbarU[j]
        betaU[i]      += Jac_dUrfm_dDSphUD[i][j]*betaSphU[j]
        BU[i]         += Jac_dUrfm_dDSphUD[i][j]*BSphU[j]
        for k in range(DIM):
            for l in range(DIM):
                gammaDD[i][j]    += Jac_dUSph_dDrfmUD[k][i]*Jac_dUSph_dDrfmUD[l][j]*   gammaSphDD[k][l]
                gammabarDD[i][j] += Jac_dUSph_dDrfmUD[k][i]*Jac_dUSph_dDrfmUD[l][j]*gammaSphbarDD[k][l]
                AbarDD[i][j]     += Jac_dUSph_dDrfmUD[k][i]*Jac_dUSph_dDrfmUD[l][j]*    ASphbarDD[k][l]


# Next we set the conformal factor variable $\texttt{cf}$, which is set by the "BSSN_RHSs::ConformalFactor" parameter. For example if "ConformalFactor" is set to "phi", we can use Eq. 3 of [Ruchlin *et al.*](https://arxiv.org/pdf/1712.07658.pdf), which in arbitrary coordinates is written:
# 
# $$
# \phi = \frac{1}{12} \log\left(\frac{\gamma}{\bar{\gamma}}\right).
# $$
# 
# Alternatively if "BSSN_RHSs::ConformalFactor" is set to "chi", then
# $$
# \chi = e^{-4 \phi} = \exp\left(-4 \frac{1}{12} \left(\frac{\gamma}{\bar{\gamma}}\right)\right) 
# = \exp\left(-\frac{1}{3} \log\left(\frac{\gamma}{\bar{\gamma}}\right)\right) = \left(\frac{\gamma}{\bar{\gamma}}\right)^{-1/3}.
# $$
# 
# Finally if "BSSN_RHSs::ConformalFactor" is set to "W", then
# $$
# W = e^{-2 \phi} = \exp\left(-2 \frac{1}{12} \log\left(\frac{\gamma}{\bar{\gamma}}\right)\right) = 
# \exp\left(-\frac{1}{6} \log\left(\frac{\gamma}{\bar{\gamma}}\right)\right) = 
# \left(\frac{\gamma}{\bar{\gamma}}\right)^{-1/6}.
# $$

# Step 5: Set the conformal factor variable according to the parameter BSSN_RHSs::ConformalFactor
cf = sp.sympify(0)

gammaUU, gammaDET = ixp.symm_matrix_inverter3x3(gammaDD)
gammabarUU, gammabarDET = ixp.symm_matrix_inverter3x3(gammabarDD)

if par.parval_from_str("ConformalFactor") == "phi":
    cf = sp.Rational(1,12)*sp.log(gammaDET/gammabarDET)
elif par.parval_from_str("ConformalFactor") == "chi":
    cf = (gammaDET/gammabarDET)**(-sp.Rational(1,3))
elif par.parval_from_str("ConformalFactor") == "W":
    cf = (gammaDET/gammabarDET)**(-sp.Rational(1,6))
else:
    print("Error ConformalFactor type = \""+par.parval_from_str("ConformalFactor")+"\" unknown.")
    exit(1)


# Finally, we rescale tensorial quantities according to the prescription described in the [BSSN in curvilinear coordinates tutorial module](Tutorial-BSSNCurvilinear.ipynb) (also [Ruchlin *et al.*](https://arxiv.org/pdf/1712.07658.pdf)):
# 
# \begin{align}
# h_{ij} &= (\bar{\gamma}_{ij} - \hat{\gamma}_{ij})/\text{ReDD[i][j]}\\
# a_{ij} &= \bar{A}_{ij}/\text{ReDD[i][j]}\\
# \lambda^i &= \bar{\Lambda}^i/\text{ReU[i]}\\
# \mathcal{V}^i &= \beta^i/\text{ReU[i]}\\
# \mathcal{B}^i &= B^i/\text{ReU[i]}\\
# \end{align}

# Step 6: Rescale all tensorial quantities.

if rfm.have_already_called_reference_metric_function == False:
    print("Error. Called Convert_Spherical_ADM_to_BSSN_curvilinear() without")
    print("       first setting up reference metric, by calling rfm.reference_metric().")
    exit(1)

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
