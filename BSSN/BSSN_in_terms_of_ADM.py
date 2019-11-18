# As documented in the NRPy+ tutorial module
#   Tutorial-BSSN_in_terms_of_ADM.ipynb,
#   this module will construct expressions for BSSN
#   quantities in terms of ADM quantities.

# Author: Zachariah B. Etienne
#         zachetie **at** gmail **dot* com

# Step 1: Import needed core NRPy+ modules
from outputC import *             # NRPy+: Core C code output module
import NRPy_param_funcs as par    # NRPy+: Parameter interface
import indexedexp as ixp          # NRPy+: Symbolic indexed expression (e.g., tensors, vectors, etc.) support
import reference_metric as rfm    # NRPy+: Reference metric support
import sys                        # Standard Python modules for multiplatform OS-level functions
import BSSN.BSSN_quantities as Bq # NRPy+: This module depends on the parameter EvolvedConformalFactor_cf,
                                  #        which is defined in BSSN.BSSN_quantities

# Step 1.a: Set DIM=3, as we're using a 3+1 decomposition of Einstein's equations
DIM=3

# Step 2: All ADM quantities were input into this function in the Spherical or Cartesian
#         basis, as functions of r,th,ph or x,y,z, respectively. In Steps 1 and 2 above,
#         we converted them to the xx0,xx1,xx2 basis, and as functions of xx0,xx1,xx2.
#         Here we convert ADM quantities to their BSSN Curvilinear counterparts:

# Step 2.a: Convert ADM $\gamma_{ij}$ to BSSN $\bar{gamma}_{ij}$:
#           We have (Eqs. 2 and 3 of [Ruchlin *et al.*](https://arxiv.org/pdf/1712.07658.pdf)):
def gammabarDD_hDD(gammaDD):
    global gammabarDD,hDD

    if gammaDD == None:
        gammaDD = ixp.declarerank2("gammaDD","sym01")

    if rfm.have_already_called_reference_metric_function == False:
        print("BSSN.BSSN_in_terms_of_ADM.hDD_given_ADM(): Must call reference_metric() first!")
        sys.exit(1)
    # \bar{gamma}_{ij} = (\frac{\bar{gamma}}{gamma})^{1/3}*gamma_{ij}.
    gammaUU, gammaDET = ixp.symm_matrix_inverter3x3(gammaDD)
    gammabarDD = ixp.zerorank2()
    hDD        = ixp.zerorank2()
    for i in range(DIM):
        for j in range(DIM):
            gammabarDD[i][j] = (rfm.detgammahat/gammaDET)**(sp.Rational(1,3))*gammaDD[i][j]
            hDD[i][j] = (gammabarDD[i][j] - rfm.ghatDD[i][j]) / rfm.ReDD[i][j]

# Step 2.b: Convert the extrinsic curvature K_{ij} to the trace-free extrinsic
#           curvature \bar{A}_{ij}, plus the trace of the extrinsic curvature K,
#           where (Eq. 3 of [Baumgarte *et al.*](https://arxiv.org/pdf/1211.6632.pdf)):
def trK_AbarDD_aDD(gammaDD,KDD):
    global trK,AbarDD,aDD

    if gammaDD == None:
        gammaDD = ixp.declarerank2("gammaDD","sym01")
    if KDD == None:
        KDD = ixp.declarerank2("KDD","sym01")

    if rfm.have_already_called_reference_metric_function == False:
        print("BSSN.BSSN_in_terms_of_ADM.trK_AbarDD(): Must call reference_metric() first!")
        sys.exit(1)
    # \bar{gamma}_{ij} = (\frac{\bar{gamma}}{gamma})^{1/3}*gamma_{ij}.
    gammaUU, gammaDET = ixp.symm_matrix_inverter3x3(gammaDD)
    # K = gamma^{ij} K_{ij}, and
    # \bar{A}_{ij} &= (\frac{\bar{gamma}}{gamma})^{1/3}*(K_{ij} - \frac{1}{3}*gamma_{ij}*K)
    trK = sp.sympify(0)
    for i in range(DIM):
        for j in range(DIM):
            trK += gammaUU[i][j]*KDD[i][j]

    AbarDD = ixp.zerorank2()
    aDD    = ixp.zerorank2()
    for i in range(DIM):
        for j in range(DIM):
            AbarDD[i][j] = (rfm.detgammahat/gammaDET)**(sp.Rational(1,3))*(KDD[i][j] - sp.Rational(1,3)*gammaDD[i][j]*trK)
            aDD[i][j]    = AbarDD[i][j] / rfm.ReDD[i][j]


# Step 2.c: Define \bar{Lambda}^i (Eqs. 4 and 5 of [Baumgarte *et al.*](https://arxiv.org/pdf/1211.6632.pdf)):
def LambdabarU_lambdaU__exact_gammaDD(gammaDD):
    global LambdabarU, lambdaU

    if gammaDD == None:
        gammaDD = ixp.declarerank2("gammaDD","sym01")

    # \bar{Lambda}^i = \bar{gamma}^{jk}(\bar{Gamma}^i_{jk} - \hat{Gamma}^i_{jk}).
    gammabarDD_hDD(gammaDD)
    gammabarUU, gammabarDET = ixp.symm_matrix_inverter3x3(gammabarDD)

    # First compute Christoffel symbols \bar{Gamma}^i_{jk}, with respect to barred metric:
    GammabarUDD = ixp.zerorank3()
    for i in range(DIM):
        for j in range(DIM):
            for k in range(DIM):
                for l in range(DIM):
                    GammabarUDD[i][j][k] += sp.Rational(1, 2) * gammabarUU[i][l] * (
                                sp.diff(gammabarDD[l][j], rfm.xx[k]) +
                                sp.diff(gammabarDD[l][k], rfm.xx[j]) -
                                sp.diff(gammabarDD[j][k], rfm.xx[l]))
    # Next evaluate \bar{Lambda}^i, based on GammabarUDD above and GammahatUDD
    #       (from the reference metric):
    LambdabarU = ixp.zerorank1()
    for i in range(DIM):
        for j in range(DIM):
            for k in range(DIM):
                LambdabarU[i] += gammabarUU[j][k] * (GammabarUDD[i][j][k] - rfm.GammahatUDD[i][j][k])
    lambdaU    = ixp.zerorank1()
    for i in range(DIM):
        lambdaU[i] = LambdabarU[i] / rfm.ReU[i]


# Step 2.d: Set the conformal factor variable cf, which is set
#           by the "BSSN_quantities::EvolvedConformalFactor_cf" parameter. For example if
#           "EvolvedConformalFactor_cf" is set to "phi", we can use Eq. 3 of
#           [Ruchlin *et al.*](https://arxiv.org/pdf/1712.07658.pdf),
#           which in arbitrary coordinates is written:
def cf_from_gammaDD(gammaDD):
    global cf

    if gammaDD == None:
        gammaDD = ixp.declarerank2("gammaDD","sym01")

    # \bar{Lambda}^i = \bar{gamma}^{jk}(\bar{Gamma}^i_{jk} - \hat{Gamma}^i_{jk}).
    gammabarDD_hDD(gammaDD)
    gammabarUU, gammabarDET = ixp.symm_matrix_inverter3x3(gammabarDD)
    gammaUU, gammaDET       = ixp.symm_matrix_inverter3x3(gammaDD)

    cf = sp.sympify(0)

    if par.parval_from_str("EvolvedConformalFactor_cf") == "phi":
        # phi = \frac{1}{12} log(\frac{gamma}{\bar{gamma}}).
        cf = sp.Rational(1, 12) * sp.log(gammaDET / gammabarDET)
    elif par.parval_from_str("EvolvedConformalFactor_cf") == "chi":
        # chi = exp(-4*phi) = exp(-4*\frac{1}{12}*(\frac{gamma}{\bar{gamma}}))
        #      = exp(-\frac{1}{3}*log(\frac{gamma}{\bar{gamma}})) = (\frac{gamma}{\bar{gamma}})^{-1/3}.
        #
        cf = (gammaDET / gammabarDET) ** (-sp.Rational(1, 3))
    elif par.parval_from_str("EvolvedConformalFactor_cf") == "W":
        # W = exp(-2*phi) = exp(-2*\frac{1}{12}*log(\frac{gamma}{\bar{gamma}}))
        #   = exp(-\frac{1}{6}*log(\frac{gamma}{\bar{gamma}})) = (\frac{gamma}{bar{gamma}})^{-1/6}.
        cf = (gammaDET / gammabarDET) ** (-sp.Rational(1, 6))
    else:
        print("Error EvolvedConformalFactor_cf type = \"" + par.parval_from_str(
            "EvolvedConformalFactor_cf") + "\" unknown.")
        sys.exit(1)

# Step 2.e: Rescale beta^i and B^i according to the prescription described in
#         the [BSSN in curvilinear coordinates tutorial module](Tutorial-BSSNCurvilinear.ipynb)
#         (also [Ruchlin *et al.*](https://arxiv.org/pdf/1712.07658.pdf)):
#
# \mathcal{V}^i &= beta^i/(ReU[i])
# \mathcal{B}^i &= B^i/(ReU[i])
def betU_vetU(betaU,BU):
    global vetU,betU

    if betaU == None:
        betaU = ixp.declarerank1("betaU")
    if BU == None:
        BU = ixp.declarerank1("BU")

    if rfm.have_already_called_reference_metric_function == False:
        print("BSSN.BSSN_in_terms_of_ADM.bet_vet(): Must call reference_metric() first!")
        sys.exit(1)
    vetU    = ixp.zerorank1()
    betU    = ixp.zerorank1()
    for i in range(DIM):
        vetU[i]    =      betaU[i] / rfm.ReU[i]
        betU[i]    =         BU[i] / rfm.ReU[i]
