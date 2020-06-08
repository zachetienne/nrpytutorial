# As documented in the NRPy+ tutorial module
#   Tutorial-BSSN_in_terms_of_ADM.ipynb,
#   this module will construct expressions for BSSN
#   quantities in terms of ADM quantities.

# Author: Zachariah B. Etienne
#         zachetie **at** gmail **dot* com

# Step 1: Import needed core NRPy+ modules
import sympy as sp                # SymPy: The Python computer algebra package upon which NRPy+ depends
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

    if gammaDD is None:
        gammaDD = ixp.declarerank2("gammaDD","sym01")

    if rfm.have_already_called_reference_metric_function == False:
        print("BSSN.BSSN_in_terms_of_ADM.hDD_given_ADM(): Must call reference_metric() first!")
        sys.exit(1)
    # \bar{gamma}_{ij} = (\frac{\bar{gamma}}{gamma})^{1/3}*gamma_{ij}.
    _gammaUU, gammaDET = ixp.symm_matrix_inverter3x3(gammaDD) # _gammaUU unused.
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

    if gammaDD is None: # Use "is None" instead of "==None", as the former is more correct.
        gammaDD = ixp.declarerank2("gammaDD","sym01")
    if KDD is None: # Use "is None" instead of "==None", as the former is more correct.
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

    if gammaDD is None: # Use "is None" instead of "==None", as the former is more correct.
        gammaDD = ixp.declarerank2("gammaDD","sym01")

    # \bar{Lambda}^i = \bar{gamma}^{jk}(\bar{Gamma}^i_{jk} - \hat{Gamma}^i_{jk}).
    gammabarDD_hDD(gammaDD)
    gammabarUU, _gammabarDET = ixp.symm_matrix_inverter3x3(gammabarDD) # _gammabarDET unused.

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
    for i in range(DIM):
        # We evaluate LambdabarU[i] here to ensure proper cancellations. If these cancellations
        #   are not applied, certain expressions (e.g., lambdaU[0] in StaticTrumpet) will
        #   cause SymPy's (v1.5+) CSE algorithm to hang
        LambdabarU[i] = LambdabarU[i].doit()
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

    if gammaDD is None: # Use "is None" instead of "==None", as the former is more correct.
        gammaDD = ixp.declarerank2("gammaDD","sym01")

    # \bar{Lambda}^i = \bar{gamma}^{jk}(\bar{Gamma}^i_{jk} - \hat{Gamma}^i_{jk}).
    gammabarDD_hDD(gammaDD)
    _gammabarUU, gammabarDET = ixp.symm_matrix_inverter3x3(gammabarDD) # _gammabarUU unused.
    _gammaUU, gammaDET       = ixp.symm_matrix_inverter3x3(gammaDD)    # _gammaUU unused.

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

    if betaU is None: # Use "is None" instead of "==None", as the former is more correct.
        betaU = ixp.declarerank1("betaU")
    if BU is None: # Use "is None" instead of "==None", as the former is more correct.
        BU = ixp.declarerank1("BU")

    if rfm.have_already_called_reference_metric_function == False:
        print("BSSN.BSSN_in_terms_of_ADM.bet_vet(): Must call reference_metric() first!")
        sys.exit(1)
    vetU    = ixp.zerorank1()
    betU    = ixp.zerorank1()
    for i in range(DIM):
        vetU[i]    =      betaU[i] / rfm.ReU[i]
        betU[i]    =         BU[i] / rfm.ReU[i]
