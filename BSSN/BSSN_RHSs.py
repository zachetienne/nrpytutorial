# As documented in the NRPy+ tutorial module
#   Tutorial-BSSN_time_evolution-BSSN_RHSs.ipynb,
#   this module will construct the right-hand sides (RHSs)
#   expressions of the BSSN time evolution equations.
#
# Time-evolution equations for the BSSN gauge conditions are
#   specified in the BSSN_gauge_RHSs module and documented in
#   the Tutorial-BSSN_time_evolution-BSSN_gauge_RHSs.ipynb
#   NRPy+ tutorial module.

# Author: Zachariah B. Etienne
#         zachetie **at** gmail **dot* com

# Step 1.a: import all needed modules from NRPy+:
import sympy as sp                # SymPy: The Python computer algebra package upon which NRPy+ depends
import NRPy_param_funcs as par    # NRPy+: Parameter interface
import indexedexp as ixp          # NRPy+: Symbolic indexed expression (e.g., tensors, vectors, etc.) support
import reference_metric as rfm    # NRPy+: Reference metric support
import sys                        # Standard Python module for multiplatform OS-level functions
import BSSN.BSSN_quantities as Bq # NRPy+: Basic BSSN quantities

have_already_called_BSSN_RHSs_function = False

# Step 1.b: Set the coordinate system for the numerical grid:
#  DO NOT SET IN STANDALONE PYTHON MODULE
# par.set_parval_from_str("reference_metric::CoordSystem","Spherical")

def BSSN_RHSs():
    # Step 1.c: Given the chosen coordinate system, set up
    #           corresponding reference metric and needed
    #           reference metric quantities
    # The following function call sets up the reference metric
    #    and related quantities, including rescaling matrices ReDD,
    #    ReU, and hatted quantities.
    rfm.reference_metric()

    global have_already_called_BSSN_RHSs_function # setting to global enables other modules to see updated value.
    have_already_called_BSSN_RHSs_function = True

    # Step 1.d: Set spatial dimension (must be 3 for BSSN, as BSSN is
    #           a 3+1-dimensional decomposition of the general
    #           relativistic field equations)
    DIM = 3

    # Step 1.e: Import all basic (unrescaled) BSSN scalars & tensors
    Bq.BSSN_basic_tensors()
    gammabarDD = Bq.gammabarDD
    AbarDD     = Bq.AbarDD
    LambdabarU = Bq.LambdabarU
    trK   = Bq.trK
    alpha = Bq.alpha
    betaU = Bq.betaU

    # Step 1.f: Import all needed rescaled BSSN tensors:
    cf  = Bq.cf
    lambdaU = Bq.lambdaU

    # Step 2.a.i: Import derivative expressions for betaU defined in the BSSN.BSSN_quantities module:
    Bq.betaU_derivs()
    betaU_dD = Bq.betaU_dD
    betaU_dDD = Bq.betaU_dDD
    # Step 2.a.ii: Import derivative expression for gammabarDD
    Bq.gammabar__inverse_and_derivs()
    gammabarDD_dupD = Bq.gammabarDD_dupD

    # Step 2.a.iii: First term of \partial_t \bar{\gamma}_{i j} right-hand side:
    # \beta^k \bar{\gamma}_{ij,k} + \beta^k_{,i} \bar{\gamma}_{kj} + \beta^k_{,j} \bar{\gamma}_{ik}
    gammabar_rhsDD = ixp.zerorank2()
    for i in range(DIM):
        for j in range(DIM):
            for k in range(DIM):
                gammabar_rhsDD[i][j] += betaU[k] * gammabarDD_dupD[i][j][k] + betaU_dD[k][i] * gammabarDD[k][j] \
                                        + betaU_dD[k][j] * gammabarDD[i][k]

    # Step 2.b.i: First import \bar{A}_{ij} = AbarDD[i][j], and its contraction trAbar = \bar{A}^k_k
    #           from BSSN.BSSN_quantities
    Bq.AbarUU_AbarUD_trAbar_AbarDD_dD()
    trAbar = Bq.trAbar

    # Step 2.b.ii: Import detgammabar quantities from BSSN.BSSN_quantities:
    Bq.detgammabar_and_derivs()
    detgammabar = Bq.detgammabar
    detgammabar_dD = Bq.detgammabar_dD

    # Step 2.b.ii: Compute the contraction \bar{D}_k \beta^k = \beta^k_{,k} + \frac{\beta^k \bar{\gamma}_{,k}}{2 \bar{\gamma}}
    Dbarbetacontraction = sp.sympify(0)
    for k in range(DIM):
        Dbarbetacontraction += betaU_dD[k][k] + betaU[k] * detgammabar_dD[k] / (2 * detgammabar)

    # Step 2.b.iii: Second term of \partial_t \bar{\gamma}_{i j} right-hand side:
    # \frac{2}{3} \bar{\gamma}_{i j} \left (\alpha \bar{A}_{k}^{k} - \bar{D}_{k} \beta^{k}\right )
    for i in range(DIM):
        for j in range(DIM):
            gammabar_rhsDD[i][j] += sp.Rational(2, 3) * gammabarDD[i][j] * (alpha * trAbar - Dbarbetacontraction)

    # Step 2.c: Third term of \partial_t \bar{\gamma}_{i j} right-hand side:
    # -2 \alpha \bar{A}_{ij}
    for i in range(DIM):
        for j in range(DIM):
            gammabar_rhsDD[i][j] += -2 * alpha * AbarDD[i][j]

    # Step 3.a: First term of \partial_t \bar{A}_{i j}:
    # \beta^k \partial_k \bar{A}_{ij} + \partial_i \beta^k \bar{A}_{kj} + \partial_j \beta^k \bar{A}_{ik}

    # First define AbarDD_dupD:
    AbarDD_dupD = Bq.AbarDD_dupD # From Bq.AbarUU_AbarUD_trAbar_AbarDD_dD()

    Abar_rhsDD = ixp.zerorank2()
    for i in range(DIM):
        for j in range(DIM):
            for k in range(DIM):
                Abar_rhsDD[i][j] += betaU[k] * AbarDD_dupD[i][j][k] + betaU_dD[k][i] * AbarDD[k][j] \
                                    + betaU_dD[k][j] * AbarDD[i][k]

    # Step 3.b: Second term of \partial_t \bar{A}_{i j}:
    # - (2/3) \bar{A}_{i j} \bar{D}_{k} \beta^{k} - 2 \alpha \bar{A}_{i k} {\bar{A}^{k}}_{j} + \alpha \bar{A}_{i j} K
    gammabarUU = Bq.gammabarUU  # From Bq.gammabar__inverse_and_derivs()
    AbarUD = Bq.AbarUD  # From Bq.AbarUU_AbarUD_trAbar()
    for i in range(DIM):
        for j in range(DIM):
            Abar_rhsDD[i][j] += -sp.Rational(2, 3) * AbarDD[i][j] * Dbarbetacontraction + alpha * AbarDD[i][j] * trK
            for k in range(DIM):
                Abar_rhsDD[i][j] += -2 * alpha * AbarDD[i][k] * AbarUD[k][j]

    # Step 3.c.i: Define partial derivatives of \phi in terms of evolved quantity "cf":
    Bq.phi_and_derivs()
    phi_dD = Bq.phi_dD
    phi_dupD = Bq.phi_dupD
    exp_m4phi = Bq.exp_m4phi
    phi_dBarD = Bq.phi_dBarD  # phi_dBarD = Dbar_i phi = phi_dD (since phi is a scalar)
    phi_dBarDD = Bq.phi_dBarDD  # phi_dBarDD = Dbar_i Dbar_j phi (covariant derivative)

    # Step 3.c.ii: Define RbarDD
    Bq.RicciBar__gammabarDD_dHatD__DGammaUDD__DGammaU()
    RbarDD = Bq.RbarDD

    # Step 3.c.iii: Define first and second derivatives of \alpha, as well as
    #         \bar{D}_i \bar{D}_j \alpha, which is defined just like phi
    alpha_dD = ixp.declarerank1("alpha_dD")
    alpha_dDD = ixp.declarerank2("alpha_dDD", "sym01")
    alpha_dBarD = alpha_dD
    alpha_dBarDD = ixp.zerorank2()
    GammabarUDD = Bq.GammabarUDD  # Defined in Bq.gammabar__inverse_and_derivs()
    for i in range(DIM):
        for j in range(DIM):
            alpha_dBarDD[i][j] = alpha_dDD[i][j]
            for k in range(DIM):
                alpha_dBarDD[i][j] += - GammabarUDD[k][i][j] * alpha_dD[k]

    # Step 3.c.iv: Define the terms in curly braces:
    curlybrackettermsDD = ixp.zerorank2()
    for i in range(DIM):
        for j in range(DIM):
            curlybrackettermsDD[i][j] = -2 * alpha * phi_dBarDD[i][j] + 4 * alpha * phi_dBarD[i] * phi_dBarD[j] \
                                        + 2 * alpha_dBarD[i] * phi_dBarD[j] \
                                        + 2 * alpha_dBarD[j] * phi_dBarD[i] \
                                        - alpha_dBarDD[i][j] + alpha * RbarDD[i][j]

    # Step 3.c.v: Compute the trace:
    curlybracketterms_trace = sp.sympify(0)
    for i in range(DIM):
        for j in range(DIM):
            curlybracketterms_trace += gammabarUU[i][j] * curlybrackettermsDD[i][j]

    # Step 3.c.vi: Third and final term of Abar_rhsDD[i][j]:
    for i in range(DIM):
        for j in range(DIM):
            Abar_rhsDD[i][j] += exp_m4phi * (curlybrackettermsDD[i][j] -
                                             sp.Rational(1, 3) * gammabarDD[i][j] * curlybracketterms_trace)

    # Step 4: Right-hand side of conformal factor variable "cf". Supported
    #          options include: cf=phi, cf=W=e^(-2*phi) (default), and cf=chi=e^(-4*phi)
    # \partial_t phi = \left[\beta^k \partial_k \phi \right] <- TERM 1
    #                  + \frac{1}{6} \left (\bar{D}_{k} \beta^{k} - \alpha K \right ) <- TERM 2
    global cf_rhs
    cf_rhs = sp.Rational(1, 6) * (Dbarbetacontraction - alpha * trK)  # Term 2
    for k in range(DIM):
        cf_rhs += betaU[k] * phi_dupD[k]  # Term 1

    # Next multiply to convert phi_rhs to cf_rhs.
    if par.parval_from_str("BSSN.BSSN_quantities::EvolvedConformalFactor_cf") == "phi":
        pass  # do nothing; cf_rhs = phi_rhs
    elif par.parval_from_str("BSSN.BSSN_quantities::EvolvedConformalFactor_cf") == "W":
        cf_rhs *= -2 * cf  # cf_rhs = -2*cf*phi_rhs
    elif par.parval_from_str("BSSN.BSSN_quantities::EvolvedConformalFactor_cf") == "chi":
        cf_rhs *= -4 * cf  # cf_rhs = -4*cf*phi_rhs
    else:
        print("Error: EvolvedConformalFactor_cf == " +
              par.parval_from_str("BSSN.BSSN_quantities::EvolvedConformalFactor_cf") + " unsupported!")
        sys.exit(1)

    # Step 5: right-hand side of trK (trace of extrinsic curvature):
    # \partial_t K = \beta^k \partial_k K <- TERM 1
    #           + \frac{1}{3} \alpha K^{2} <- TERM 2
    #           + \alpha \bar{A}_{i j} \bar{A}^{i j} <- TERM 3
    #           - - e^{-4 \phi} (\bar{D}_{i} \bar{D}^{i} \alpha + 2 \bar{D}^{i} \alpha \bar{D}_{i} \phi ) <- TERM 4
    global trK_rhs
    # TERM 2:
    trK_rhs = sp.Rational(1, 3) * alpha * trK * trK
    trK_dupD = ixp.declarerank1("trK_dupD")
    for i in range(DIM):
        # TERM 1:
        trK_rhs += betaU[i] * trK_dupD[i]
    for i in range(DIM):
        for j in range(DIM):
            # TERM 4:
            trK_rhs += -exp_m4phi * gammabarUU[i][j] * (alpha_dBarDD[i][j] + 2 * alpha_dBarD[j] * phi_dBarD[i])
    AbarUU = Bq.AbarUU  # From Bq.AbarUU_AbarUD_trAbar()
    for i in range(DIM):
        for j in range(DIM):
            # TERM 3:
            trK_rhs += alpha * AbarDD[i][j] * AbarUU[i][j]

    # Step 6: right-hand side of \partial_t \bar{\Lambda}^i:
    # \partial_t \bar{\Lambda}^i = \beta^k \partial_k \bar{\Lambda}^i - \partial_k \beta^i \bar{\Lambda}^k <- TERM 1
    #                            + \bar{\gamma}^{j k} \hat{D}_{j} \hat{D}_{k} \beta^{i} <- TERM 2
    #                            + \frac{2}{3} \Delta^{i} \bar{D}_{j} \beta^{j} <- TERM 3
    #                            + \frac{1}{3} \bar{D}^{i} \bar{D}_{j} \beta^{j} <- TERM 4
    #                            - 2 \bar{A}^{i j} (\partial_{j} \alpha - 6 \partial_{j} \phi) <- TERM 5
    #                            + 2 \alpha \bar{A}^{j k} \Delta_{j k}^{i} <- TERM 6
    #                            - \frac{4}{3} \alpha \bar{\gamma}^{i j} \partial_{j} K <- TERM 7

    # Step 6.a: Term 1 of \partial_t \bar{\Lambda}^i: \beta^k \partial_k \bar{\Lambda}^i - \partial_k \beta^i \bar{\Lambda}^k
    # First we declare \bar{\Lambda}^i and \bar{\Lambda}^i_{,j} in terms of \lambda^i and \lambda^i_{,j}
    global LambdabarU_dupD # Used on the RHS of the Gamma-driving shift conditions
    LambdabarU_dupD = ixp.zerorank2()
    lambdaU_dupD = ixp.declarerank2("lambdaU_dupD", "nosym")
    for i in range(DIM):
        for j in range(DIM):
            LambdabarU_dupD[i][j] = lambdaU_dupD[i][j] * rfm.ReU[i] + lambdaU[i] * rfm.ReUdD[i][j]

    global Lambdabar_rhsU # Used on the RHS of the Gamma-driving shift conditions
    Lambdabar_rhsU = ixp.zerorank1()
    for i in range(DIM):
        for k in range(DIM):
            Lambdabar_rhsU[i] += betaU[k] * LambdabarU_dupD[i][k] - betaU_dD[i][k] * LambdabarU[k]  # Term 1

    # Step 6.b: Term 2 of \partial_t \bar{\Lambda}^i = \bar{\gamma}^{jk} (Term 2a + Term 2b + Term 2c)
    # Term 2a: \bar{\gamma}^{jk} \beta^i_{,kj}
    Term2aUDD = ixp.zerorank3()
    for i in range(DIM):
        for j in range(DIM):
            for k in range(DIM):
                Term2aUDD[i][j][k] += betaU_dDD[i][k][j]
    # Term 2b: \hat{\Gamma}^i_{mk,j} \beta^m + \hat{\Gamma}^i_{mk} \beta^m_{,j}
    #          + \hat{\Gamma}^i_{dj}\beta^d_{,k} - \hat{\Gamma}^d_{kj} \beta^i_{,d}
    Term2bUDD = ixp.zerorank3()
    for i in range(DIM):
        for j in range(DIM):
            for k in range(DIM):
                for m in range(DIM):
                    Term2bUDD[i][j][k] += rfm.GammahatUDDdD[i][m][k][j] * betaU[m] \
                                          + rfm.GammahatUDD[i][m][k] * betaU_dD[m][j] \
                                          + rfm.GammahatUDD[i][m][j] * betaU_dD[m][k] \
                                          - rfm.GammahatUDD[m][k][j] * betaU_dD[i][m]
    # Term 2c: \hat{\Gamma}^i_{dj}\hat{\Gamma}^d_{mk} \beta^m - \hat{\Gamma}^d_{kj} \hat{\Gamma}^i_{md} \beta^m
    Term2cUDD = ixp.zerorank3()
    for i in range(DIM):
        for j in range(DIM):
            for k in range(DIM):
                for m in range(DIM):
                    for d in range(DIM):
                        Term2cUDD[i][j][k] += (rfm.GammahatUDD[i][d][j] * rfm.GammahatUDD[d][m][k] \
                                               - rfm.GammahatUDD[d][k][j] * rfm.GammahatUDD[i][m][d]) * betaU[m]

    Lambdabar_rhsUpieceU = ixp.zerorank1()

    # Put it all together to get Term 2:
    for i in range(DIM):
        for j in range(DIM):
            for k in range(DIM):
                Lambdabar_rhsU[i] += gammabarUU[j][k] * (Term2aUDD[i][j][k] + Term2bUDD[i][j][k] + Term2cUDD[i][j][k])
                Lambdabar_rhsUpieceU[i] += gammabarUU[j][k] * (
                            Term2aUDD[i][j][k] + Term2bUDD[i][j][k] + Term2cUDD[i][j][k])

    # Step 6.c: Term 3 of \partial_t \bar{\Lambda}^i:
    #    \frac{2}{3} \Delta^{i} \bar{D}_{j} \beta^{j}
    DGammaU = Bq.DGammaU  # From Bq.RicciBar__gammabarDD_dHatD__DGammaUDD__DGammaU()
    for i in range(DIM):
        Lambdabar_rhsU[i] += sp.Rational(2, 3) * DGammaU[i] * Dbarbetacontraction  # Term 3

    # Step 6.d: Term 4 of \partial_t \bar{\Lambda}^i:
    #           \frac{1}{3} \bar{D}^{i} \bar{D}_{j} \beta^{j}
    detgammabar_dDD = Bq.detgammabar_dDD  # From Bq.detgammabar_and_derivs()
    Dbarbetacontraction_dBarD = ixp.zerorank1()
    for k in range(DIM):
        for m in range(DIM):
            Dbarbetacontraction_dBarD[m] += betaU_dDD[k][k][m] + \
                                            (betaU_dD[k][m] * detgammabar_dD[k] +
                                             betaU[k] * detgammabar_dDD[k][m]) / (2 * detgammabar) \
                                            - betaU[k] * detgammabar_dD[k] * detgammabar_dD[m] / (
                                                        2 * detgammabar * detgammabar)
    for i in range(DIM):
        for m in range(DIM):
            Lambdabar_rhsU[i] += sp.Rational(1, 3) * gammabarUU[i][m] * Dbarbetacontraction_dBarD[m]

    # Step 6.e: Term 5 of \partial_t \bar{\Lambda}^i:
    #           - 2 \bar{A}^{i j} (\partial_{j} \alpha - 6 \alpha \partial_{j} \phi)
    for i in range(DIM):
        for j in range(DIM):
            Lambdabar_rhsU[i] += -2 * AbarUU[i][j] * (alpha_dD[j] - 6 * alpha * phi_dD[j])

    # Step 6.f: Term 6 of \partial_t \bar{\Lambda}^i:
    #           2 \alpha \bar{A}^{j k} \Delta^{i}_{j k}
    DGammaUDD = Bq.DGammaUDD  # From RicciBar__gammabarDD_dHatD__DGammaUDD__DGammaU()
    for i in range(DIM):
        for j in range(DIM):
            for k in range(DIM):
                Lambdabar_rhsU[i] += 2 * alpha * AbarUU[j][k] * DGammaUDD[i][j][k]

    # Step 6.g: Term 7 of \partial_t \bar{\Lambda}^i:
    #           -\frac{4}{3} \alpha \bar{\gamma}^{i j} \partial_{j} K
    trK_dD = ixp.declarerank1("trK_dD")
    for i in range(DIM):
        for j in range(DIM):
            Lambdabar_rhsU[i] += -sp.Rational(4, 3) * alpha * gammabarUU[i][j] * trK_dD[j]

    # Step 7: Rescale the RHS quantities so that the evolved
    #         variables are smooth across coord singularities
    global h_rhsDD,a_rhsDD,lambda_rhsU
    h_rhsDD = ixp.zerorank2()
    a_rhsDD = ixp.zerorank2()
    lambda_rhsU = ixp.zerorank1()
    for i in range(DIM):
        lambda_rhsU[i] = Lambdabar_rhsU[i] / rfm.ReU[i]
        for j in range(DIM):
            h_rhsDD[i][j] = gammabar_rhsDD[i][j] / rfm.ReDD[i][j]
            a_rhsDD[i][j] = Abar_rhsDD[i][j] / rfm.ReDD[i][j]
    # print(str(Abar_rhsDD[2][2]).replace("**","^").replace("_","").replace("xx","x").replace("sin(x2)","Sin[x2]").replace("sin(2*x2)","Sin[2*x2]").replace("cos(x2)","Cos[x2]").replace("detgbaroverdetghat","detg"))
    # print(str(Dbarbetacontraction).replace("**","^").replace("_","").replace("xx","x").replace("sin(x2)","Sin[x2]").replace("detgbaroverdetghat","detg"))
    # print(betaU_dD)
    # print(str(trK_rhs).replace("xx2","xx3").replace("xx1","xx2").replace("xx0","xx1").replace("**","^").replace("_","").replace("sin(xx2)","Sinx2").replace("xx","x").replace("sin(2*x2)","Sin2x2").replace("cos(x2)","Cosx2").replace("detgbaroverdetghat","detg"))
    # print(str(bet_rhsU[0]).replace("xx2","xx3").replace("xx1","xx2").replace("xx0","xx1").replace("**","^").replace("_","").replace("sin(xx2)","Sinx2").replace("xx","x").replace("sin(2*x2)","Sin2x2").replace("cos(x2)","Cosx2").replace("detgbaroverdetghat","detg"))
