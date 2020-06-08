# As documented in the NRPy+ tutorial module
#   Tutorial-ADM_in_terms_of_BSSN.ipynb,
#   this module will construct expressions for ADM
#   quantities in terms of BSSN quantities.

# Author: Zachariah B. Etienne
#         zachetie **at** gmail **dot* com

# Step 1.a: import all needed modules from NRPy+:
import NRPy_param_funcs as par    # NRPy+: parameter interface
import sympy as sp                # SymPy: The Python computer algebra package upon which NRPy+ depends
import indexedexp as ixp          # NRPy+: Symbolic indexed expression (e.g., tensors, vectors, etc.) support
import reference_metric as rfm    # NRPy+: Reference metric support
import sys                        # Standard Python module for multiplatform OS-level functions

def ADM_in_terms_of_BSSN():
    global gammaDD, gammaDDdD, gammaDDdDD, gammaUU, detgamma, GammaUDD, KDD, KDDdD
    # Step 1.c: Given the chosen coordinate system, set up
    #           corresponding reference metric and needed
    #           reference metric quantities
    # The following function call sets up the reference metric
    #    and related quantities, including rescaling matrices ReDD,
    #    ReU, and hatted quantities.
    rfm.reference_metric()

    # Step 1.d: Set spatial dimension (must be 3 for BSSN, as BSSN is
    #           a 3+1-dimensional decomposition of the general
    #           relativistic field equations)
    DIM = 3

    # Step 1.e: Import all basic (unrescaled) BSSN scalars & tensors
    import BSSN.BSSN_quantities as Bq
    Bq.BSSN_basic_tensors()
    gammabarDD = Bq.gammabarDD
    cf         = Bq.cf
    AbarDD     = Bq.AbarDD
    trK        = Bq.trK

    Bq.gammabar__inverse_and_derivs()
    gammabarDD_dD  = Bq.gammabarDD_dD
    gammabarDD_dDD = Bq.gammabarDD_dDD

    Bq.AbarUU_AbarUD_trAbar_AbarDD_dD()
    AbarDD_dD = Bq.AbarDD_dD

    # Step 2: The ADM three-metric gammaDD and its
    #         derivatives in terms of BSSN quantities.
    gammaDD = ixp.zerorank2()

    exp4phi = sp.sympify(0)
    if par.parval_from_str("EvolvedConformalFactor_cf") == "phi":
        exp4phi = sp.exp(4 * cf)
    elif par.parval_from_str("EvolvedConformalFactor_cf") == "chi":
        exp4phi = (1 / cf)
    elif par.parval_from_str("EvolvedConformalFactor_cf") == "W":
        exp4phi = (1 / cf ** 2)
    else:
        print("Error EvolvedConformalFactor_cf type = \"" + par.parval_from_str("EvolvedConformalFactor_cf") + "\" unknown.")
        sys.exit(1)

    for i in range(DIM):
        for j in range(DIM):
            gammaDD[i][j] = exp4phi * gammabarDD[i][j]

    # Step 2.a: Derivatives of $e^{4\phi}$
    phidD = ixp.zerorank1()
    phidDD = ixp.zerorank2()
    cf_dD  = ixp.declarerank1("cf_dD")
    cf_dDD = ixp.declarerank2("cf_dDD","sym01")
    if par.parval_from_str("EvolvedConformalFactor_cf") == "phi":
        for i in range(DIM):
            phidD[i]  = cf_dD[i]
            for j in range(DIM):
                phidDD[i][j] = cf_dDD[i][j]
    elif par.parval_from_str("EvolvedConformalFactor_cf") == "chi":
        for i in range(DIM):
            phidD[i]  = -sp.Rational(1,4)*exp4phi*cf_dD[i]
            for j in range(DIM):
                phidDD[i][j] = sp.Rational(1,4)*( exp4phi**2*cf_dD[i]*cf_dD[j] - exp4phi*cf_dDD[i][j] )
    elif par.parval_from_str("EvolvedConformalFactor_cf") == "W":
        exp2phi = (1 / cf)
        for i in range(DIM):
            phidD[i]  = -sp.Rational(1,2)*exp2phi*cf_dD[i]
            for j in range(DIM):
                phidDD[i][j] = sp.Rational(1,2)*( exp4phi*cf_dD[i]*cf_dD[j] - exp2phi*cf_dDD[i][j] )
    else:
        print("Error EvolvedConformalFactor_cf type = \""+par.parval_from_str("EvolvedConformalFactor_cf")+"\" unknown.")
        sys.exit(1)

    exp4phidD  = ixp.zerorank1()
    exp4phidDD = ixp.zerorank2()
    for i in range(DIM):
        exp4phidD[i] = 4*exp4phi*phidD[i]
        for j in range(DIM):
            exp4phidDD[i][j] = 16*exp4phi*phidD[i]*phidD[j] + 4*exp4phi*phidDD[i][j]

    # Step 2.b: Derivatives of gammaDD, the ADM three-metric
    gammaDDdD = ixp.zerorank3()
    gammaDDdDD = ixp.zerorank4()

    for i in range(DIM):
        for j in range(DIM):
            for k in range(DIM):
                gammaDDdD[i][j][k] = exp4phidD[k] * gammabarDD[i][j] + exp4phi * gammabarDD_dD[i][j][k]
                for l in range(DIM):
                    gammaDDdDD[i][j][k][l] = exp4phidDD[k][l] * gammabarDD[i][j] + \
                                             exp4phidD[k] * gammabarDD_dD[i][j][l] + \
                                             exp4phidD[l] * gammabarDD_dD[i][j][k] + \
                                             exp4phi * gammabarDD_dDD[i][j][k][l]

    # Step 2.c: 3-Christoffel symbols associated with ADM 3-metric gammaDD
    # Step 2.c.i: First compute the inverse 3-metric gammaUU:
    gammaUU, detgamma = ixp.symm_matrix_inverter3x3(gammaDD)

    GammaUDD = ixp.zerorank3()

    for i in range(DIM):
        for j in range(DIM):
            for k in range(DIM):
                for l in range(DIM):
                    GammaUDD[i][j][k] += sp.Rational(1,2)*gammaUU[i][l]* \
                                    (gammaDDdD[l][j][k] + gammaDDdD[l][k][j] - gammaDDdD[j][k][l])

    # Step 3: Define ADM extrinsic curvature KDD and
    #         its first spatial derivatives KDDdD
    #         in terms of BSSN quantities
    KDD = ixp.zerorank2()

    for i in range(DIM):
        for j in range(DIM):
            KDD[i][j] = exp4phi * AbarDD[i][j] + sp.Rational(1, 3) * gammaDD[i][j] * trK

    KDDdD = ixp.zerorank3()
    trK_dD = ixp.declarerank1("trK_dD")
    for i in range(DIM):
        for j in range(DIM):
            for k in range(DIM):
                KDDdD[i][j][k] = exp4phidD[k] * AbarDD[i][j] + exp4phi * AbarDD_dD[i][j][k] + \
                                 sp.Rational(1, 3) * (gammaDDdD[i][j][k] * trK + gammaDD[i][j] * trK_dD[k])
