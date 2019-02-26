# This module will declare and set the right-hand sides (RHSs) of the
#    time evolution equations in the BSSN formulation of Einstein's
#    equations.
# This includes the Ricci tensor, metric tensor, and Christoffel Symbols

# Author: Zachariah B. Etienne
#         zachetie **at** gmail **dot* com

# Step P1: import all needed modules from NRPy+:
import NRPy_param_funcs as par
import indexedexp as ixp
import grid as gri
import finite_difference as fin
import reference_metric as rfm
from outputC import *
import BSSN.BSSN_rescaled_vars as Brv
import BSSN.BSSN_unrescaled_and_barred_vars as Bubv

# Calculates the useful tensors and tensor-like quantities in the BSSN metric.
# Step P2: Initialize BSSN_RHS parameters
thismodule = __name__
par.initialize_param(par.glb_param("char", thismodule, "LapseEvolutionOption", "UsualOnePluslog"))
par.initialize_param(par.glb_param("char", thismodule, "ShiftEvolutionOption", "GammaDriving2ndOrder_NoCovariant"))

def BSSN_RHSs():
    # Step 1: Set up reference metric
    rfm.reference_metric()

    # Step 2: Set spatial dimension (must be 3 for BSSN, a 3+1 decomposition of Einstein's equations of general relativity)
    DIM = 3
    par.set_parval_from_str("grid::DIM",DIM)

    # Step 3: All barred quantities are defined in terms of BSSN rescaled gridfunctions,
    #         which we declare here in case they haven't yet been declared elsewhere.
    Brv.declare_BSSN_rescaled_gridfunctions_if_not_declared_already()
    aDD     = Brv.aDD
    lambdaU = Brv.lambdaU
    betU    = Brv.betU
    alpha   = Brv.alpha
    trK     = Brv.trK
    cf      = Brv.cf

    Bubv.BSSN_barred_variables()
    gammabarDD      = Bubv.gammabarDD
    gammabarUU      = Bubv.gammabarUU
    detgammabar     = Bubv.detgammabar
    detgammabar_dD  = Bubv.detgammabar_dD
    detgammabar_dDD = Bubv.detgammabar_dDD
    gammabarDD_dupD = Bubv.gammabarDD_dupD
    GammabarUDD     = Bubv.GammabarUDD
    AbarDD          = Bubv.AbarDD
    AbarUU          = Bubv.AbarUU

    Bubv.RicciBar__gammabarDD_dHatD__DGammaUDD__DGammaU()
    RbarDD    = Bubv.RbarDD
    DGammaUDD = Bubv.DGammaUDD
    DGammaU   = Bubv.DGammaU

    Bubv.betaUbar_and_derivs()
    betaU      = Bubv.betaU
    betaU_dD   = Bubv.betaU_dD
    betaU_dupD = Bubv.betaU_dupD
    betaU_dDD  = Bubv.betaU_dDD

    # phi_dD, phi_dupD, phi_dDD, exp_m4phi, phi_dBarD, phi_dBarDD =
    Bubv.phi_and_derivs()
    exp_m4phi = Bubv.exp_m4phi
    phi_dD = Bubv.phi_dD
    phi_dupD = Bubv.phi_dupD
    phi_dBarD = Bubv.phi_dBarD
    phi_dBarDD = Bubv.phi_dBarDD

    # Step 8b: First term of \partial_t \bar{\gamma}_{i j} right-hand side:
    # \beta^k \bar{\gamma}_{ij,k} + \beta^k_{,i} \bar{\gamma}_{kj} + \beta^k_{,j} \bar{\gamma}_{ik}
    gammabar_rhsDD = ixp.zerorank2()
    for i in range(DIM):
        for j in range(DIM):
            for k in range(DIM):
                gammabar_rhsDD[i][j] += betaU[k] * gammabarDD_dupD[i][j][k] + betaU_dD[k][i] * gammabarDD[k][j] \
                                        + betaU_dD[k][j] * gammabarDD[i][k]

    # Step 8c: Define \bar{A}_{ij} = a_{ij} \text{ReDD[i][j]} = AbarDD[i][j], and its contraction trAbar = \bar{A}^k_k
    trAbar = sp.sympify(0)
    for i in range(DIM):
        for j in range(DIM):
            trAbar += gammabarUU[i][j] * AbarDD[i][j]

    # Step 8e: Compute the contraction \bar{D}_k \beta^k = \beta^k_{,k} + \frac{\beta^k \bar{\gamma}_{,k}}{2 \bar{\gamma}}
    Dbarbetacontraction = sp.sympify(0)
    for k in range(DIM):
        Dbarbetacontraction += betaU_dD[k][k] + betaU[k] * detgammabar_dD[k] / (2 * detgammabar)

    # Step 8f: Second term of \partial_t \bar{\gamma}_{i j} right-hand side:
    # \frac{2}{3} \bar{\gamma}_{i j} \left (\alpha \bar{A}_{k}^{k} - \bar{D}_{k} \beta^{k}\right )
    for i in range(DIM):
        for j in range(DIM):
            gammabar_rhsDD[i][j] += sp.Rational(2, 3) * gammabarDD[i][j] * (alpha * trAbar - Dbarbetacontraction)

    # Step 8g: Third term of \partial_t \bar{\gamma}_{i j} right-hand side:
    # -2 \alpha \bar{A}_{ij}
    for i in range(DIM):
        for j in range(DIM):
            gammabar_rhsDD[i][j] += -2 * alpha * AbarDD[i][j]

    # Step 9a: First term of \partial_t \bar{A}_{i j}:
    # \beta^k \partial_k \bar{A}_{ij} + \partial_i \beta^k \bar{A}_{kj} + \partial_j \beta^k \bar{A}_{ik}

    # First define aDD_dupD:
    AbarDD_dupD = ixp.zerorank3()
    aDD_dupD = ixp.declarerank3("aDD_dupD","sym01")
    for i in range(DIM):
        for j in range(DIM):
            for k in range(DIM):
                AbarDD_dupD[i][j][k] += aDD_dupD[i][j][k]*rfm.ReDD[i][j] + aDD[i][j]*rfm.ReDDdD[i][j][k]

    Abar_rhsDD = ixp.zerorank2()
    for i in range(DIM):
        for j in range(DIM):
            for k in range(DIM):
                Abar_rhsDD[i][j] += betaU[k]*AbarDD_dupD[i][j][k] + betaU_dD[k][i]*AbarDD[k][j] \
                                                                  + betaU_dD[k][j]*AbarDD[i][k]

    # Step 9b: Second term of \partial_t \bar{A}_{i j}:
    # - (2/3) \bar{A}_{i j} \bar{D}_{k} \beta^{k} - 2 \alpha \bar{A}_{i k} {\bar{A}^{k}}_{j} + \alpha \bar{A}_{i j} K
    AbarUD = ixp.zerorank2()
    for i in range(DIM):
        for j in range(DIM):
            for k in range(DIM):
                AbarUD[i][j] += gammabarUU[i][k]*AbarDD[k][j]

    for i in range(DIM):
        for j in range(DIM):
            Abar_rhsDD[i][j] += -sp.Rational(2,3)*AbarDD[i][j]*Dbarbetacontraction + alpha*AbarDD[i][j]*trK
            for k in range(DIM):
                Abar_rhsDD[i][j] += -2*alpha * AbarDD[i][k]*AbarUD[k][j]

    # Step 9e: Define first and second derivatives of \alpha, as well as
    #         \bar{D}_i \bar{D}_j \alpha, which is defined just like phi
    alpha_dD = ixp.declarerank1("alpha_dD")
    alpha_dDD = ixp.declarerank2("alpha_dDD", "sym01")
    alpha_dBarD = alpha_dD
    alpha_dBarDD = ixp.zerorank2()
    for i in range(DIM):
        for j in range(DIM):
            alpha_dBarDD[i][j] = alpha_dDD[i][j]
            for k in range(DIM):
                alpha_dBarDD[i][j] += - GammabarUDD[k][i][j] * alpha_dD[k]

    # Step 9f: Define the terms in curly braces:
    curlybrackettermsDD = ixp.zerorank2()
    for i in range(DIM):
        for j in range(DIM):
            curlybrackettermsDD[i][j] = -2 * alpha * phi_dBarDD[i][j] + 4 * alpha * phi_dBarD[i] * phi_dBarD[j] \
                                        + 2 * alpha_dBarD[i] * phi_dBarD[j] \
                                        + 2 * alpha_dBarD[j] * phi_dBarD[i] \
                                        - alpha_dBarDD[i][j] + alpha * RbarDD[i][j]

    # Step 9g: Compute the trace:
    curlybracketterms_trace = sp.sympify(0)
    for i in range(DIM):
        for j in range(DIM):
            curlybracketterms_trace += gammabarUU[i][j] * curlybrackettermsDD[i][j]

    # Step 9h: Third and final term of Abar_rhsDD[i][j]:
    for i in range(DIM):
        for j in range(DIM):
            Abar_rhsDD[i][j] += exp_m4phi * (curlybrackettermsDD[i][j] - \
                                             sp.Rational(1, 3) * gammabarDD[i][j] * curlybracketterms_trace)

    # Step 10: right-hand side of conformal factor variable "cf". Supported
    #          options include: cf=phi, cf=W=e^(-2*phi) (default), and cf=chi=e^(-4*phi)
    # \partial_t phi = \left[\beta^k \partial_k \phi \right] <- TERM 1
    #                  + \frac{1}{6} \left (\bar{D}_{k} \beta^{k} - \alpha K \right ) <- TERM 2
    global cf_rhs
    cf_rhs = sp.Rational(1,6) * (Dbarbetacontraction - alpha*trK) # Term 2
    for k in range(DIM):
        cf_rhs += betaU[k]*phi_dupD[k] # Term 1

    # Next multiply to convert phi_rhs to cf_rhs.
    if par.parval_from_str("BSSN.BSSN_unrescaled_and_barred_vars::ConformalFactor") == "phi":
        pass # do nothing; cf_rhs = phi_rhs
    elif par.parval_from_str("BSSN.BSSN_unrescaled_and_barred_vars::ConformalFactor") == "W":
        cf_rhs *= -2*cf # cf_rhs = -2*cf*phi_rhs
    elif par.parval_from_str("BSSN.BSSN_unrescaled_and_barred_vars::ConformalFactor") == "chi":
        cf_rhs *= -4*cf # cf_rhs = -4*cf*phi_rhs
    else:
        print("Error: ConformalFactor == "+par.parval_from_str("BSSN_unrescaled_and_barred_vars::ConformalFactor")+" unsupported!")
        exit(1)

    # Step 11: right-hand side of trK (trace of extrinsic curvature):
    # \partial_t K = \beta^k \partial_k K <- TERM 1
    #           + \frac{1}{3} \alpha K^{2} <- TERM 2
    #           + \alpha \bar{A}_{i j} \bar{A}^{i j} <- TERM 3
    #           - - e^{-4 \phi} (\bar{D}_{i} \bar{D}^{i} \alpha + 2 \bar{D}^{i} \alpha \bar{D}_{i} \phi ) <- TERM 4
    global trK_rhs
    trK_rhs = sp.Rational(1,3)*alpha*trK*trK # Term 2
    trK_dupD = ixp.declarerank1("trK_dupD")
    for k in range(DIM):
        trK_rhs += betaU[k]*trK_dupD[k] # Term 1
    for i in range(DIM):
        for j in range(DIM):
            trK_rhs += -exp_m4phi*gammabarUU[i][j]*(alpha_dBarDD[i][j] + 2*alpha_dBarD[j]*phi_dBarD[i]) # Term 4
    # global AbarUU # Needed for BSSN constraints, etc.
    # AbarUU = ixp.zerorank2() # Needed also for \partial_t \bar{\Lambda}^i
    # for i in range(DIM):
    #     for j in range(DIM):
    #         for k in range(DIM):
    #             for l in range(DIM):
    #                 AbarUU[i][j] += gammabarUU[i][k]*gammabarUU[j][l]*AbarDD[k][l]
    for i in range(DIM):
        for j in range(DIM):
            trK_rhs += alpha*AbarDD[i][j]*AbarUU[i][j] # Term 3

    # Step 12: right-hand side of \partial_t \bar{\Lambda}^i:
    # \partial_t \bar{\Lambda}^i = \beta^k \partial_k \bar{\Lambda}^i - \partial_k \beta^i \bar{\Lambda}^k <- TERM 1
    #                            + \bar{\gamma}^{j k} \hat{D}_{j} \hat{D}_{k} \beta^{i} <- TERM 2
    #                            + \frac{2}{3} \Delta^{i} \bar{D}_{j} \beta^{j} <- TERM 3
    #                            + \frac{1}{3} \bar{D}^{i} \bar{D}_{j} \beta^{j} <- TERM 4
    #                            - 2 \bar{A}^{i j} (\partial_{j} \alpha - 6 \partial_{j} \phi) <- TERM 5
    #                            + 2 \alpha \bar{A}^{j k} \Delta_{j k}^{i} <- TERM 6
    #                            - \frac{4}{3} \alpha \bar{\gamma}^{i j} \partial_{j} K <- TERM 7

    # Step 12a: Term 1 of \partial_t \bar{\Lambda}^i: \beta^k \partial_k \bar{\Lambda}^i - \partial_k \beta^i \bar{\Lambda}^k
    # First we declare \bar{\Lambda}^i and \bar{\Lambda}^i_{,j} in terms of \lambda^i and \lambda^i_{,j}
    LambarU = ixp.zerorank1()
    LambarU_dupD = ixp.zerorank2()
    lambdaU_dupD = ixp.declarerank2("lambdaU_dupD","nosym")
    for i in range(DIM):
        LambarU[i] = lambdaU[i]*rfm.ReU[i]
        for j in range(DIM):
            LambarU_dupD[i][j] = lambdaU_dupD[i][j]*rfm.ReU[i] + lambdaU[i]*rfm.ReUdD[i][j]

    Lambar_rhsU = ixp.zerorank1()
    for i in range(DIM):
        for k in range(DIM):
            Lambar_rhsU[i] += betaU[k]*LambarU_dupD[i][k] - betaU_dD[i][k]*LambarU[k] # Term 1

    # Step 12b: Term 2 of \partial_t \bar{\Lambda}^i = \bar{\gamma}^{jk} (Term 2a + Term 2b + Term 2c)
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
                    Term2bUDD[i][j][k] += rfm.GammahatUDDdD[i][m][k][j]*betaU[m]    \
                                          + rfm.GammahatUDD[i][m][k]*betaU_dD[m][j] \
                                          + rfm.GammahatUDD[i][m][j]*betaU_dD[m][k] \
                                          - rfm.GammahatUDD[m][k][j]*betaU_dD[i][m]
    # Term 2c: \hat{\Gamma}^i_{dj}\hat{\Gamma}^d_{mk} \beta^m - \hat{\Gamma}^d_{kj} \hat{\Gamma}^i_{md} \beta^m
    Term2cUDD = ixp.zerorank3()
    for i in range(DIM):
        for j in range(DIM):
            for k in range(DIM):
                for m in range(DIM):
                    for d in range(DIM):
                        Term2cUDD[i][j][k] += ( rfm.GammahatUDD[i][d][j]*rfm.GammahatUDD[d][m][k] \
                                               -rfm.GammahatUDD[d][k][j]*rfm.GammahatUDD[i][m][d])*betaU[m]

    Lambar_rhsUpieceU = ixp.zerorank1()

    # Put it all together to get Term 2:
    for i in range(DIM):
        for j in range(DIM):
            for k in range(DIM):
                Lambar_rhsU[i] += gammabarUU[j][k] * (Term2aUDD[i][j][k] + Term2bUDD[i][j][k] + Term2cUDD[i][j][k])
                Lambar_rhsUpieceU[i] += gammabarUU[j][k] * (Term2aUDD[i][j][k] + Term2bUDD[i][j][k] + Term2cUDD[i][j][k])

    # Step 12c: Term 3 of \partial_t \bar{\Lambda}^i:
    #    \frac{2}{3} \Delta^{i} \bar{D}_{j} \beta^{j}
    for i in range(DIM):
        Lambar_rhsU[i] += sp.Rational(2,3)*DGammaU[i]*Dbarbetacontraction # Term 3

    # Step 12d: Term 4 of \partial_t \bar{\Lambda}^i:
    #           \frac{1}{3} \bar{D}^{i} \bar{D}_{j} \beta^{j}
    Dbarbetacontraction_dBarD = ixp.zerorank1()
    for k in range(DIM):
        for m in range(DIM):
            Dbarbetacontraction_dBarD[m] += betaU_dDD[k][k][m] + \
                                            (betaU_dD[k][m]*detgammabar_dD[k] +
                                             betaU[k]*detgammabar_dDD[k][m])/(2*detgammabar) \
                                            -betaU[k]*detgammabar_dD[k]*detgammabar_dD[m]/(2*detgammabar*detgammabar)
    for i in range(DIM):
        for m in range(DIM):
            Lambar_rhsU[i] += sp.Rational(1,3)*gammabarUU[i][m]*Dbarbetacontraction_dBarD[m]

    # Step 12e: Term 5 of \partial_t \bar{\Lambda}^i:
    #           - 2 \bar{A}^{i j} (\partial_{j} \alpha - 6 \alpha \partial_{j} \phi)
    for i in range(DIM):
        for j in range(DIM):
            Lambar_rhsU[i] += -2*AbarUU[i][j]*(alpha_dD[j] - 6*alpha*phi_dD[j])

    # Step 12f: Term 6 of \partial_t \bar{\Lambda}^i:
    #           2 \alpha \bar{A}^{j k} \Delta^{i}_{j k}
    for i in range(DIM):
        for j in range(DIM):
            for k in range(DIM):
                Lambar_rhsU[i] += 2*alpha*AbarUU[j][k]*DGammaUDD[i][j][k]

    # Step 12g: Term 7 of \partial_t \bar{\Lambda}^i:
    #           -\frac{4}{3} \alpha \bar{\gamma}^{i j} \partial_{j} K
    trK_dD = ixp.declarerank1("trK_dD")
    for i in range(DIM):
        for j in range(DIM):
            Lambar_rhsU[i] += -sp.Rational(4,3)*alpha*gammabarUU[i][j]*trK_dD[j]

    # Step 13: \partial_t \alpha = \beta^i \alpha_{,i} - 2*\alpha*K
    global alpha_rhs
    if par.parval_from_str("BSSN.BSSN_RHSs::LapseEvolutionOption") == "UsualOnePluslog":
        alpha_rhs = -2*alpha*trK
        alpha_dupD = ixp.declarerank1("alpha_dupD")
        for i in range(DIM):
            alpha_rhs += betaU[i]*alpha_dupD[i]
    if par.parval_from_str("BSSN.BSSN_RHSs::LapseEvolutionOption") == "MaximalSlicing":
        # As defined on Pg 2 of https://arxiv.org/pdf/gr-qc/9902024.pdf , this is given by
        #   \partial_t \alpha = \partial_t e^{6 \phi} = 6 e^{6 \phi} \partial_t \phi
        # If cf = W = e^{-2 phi}, then
        #  6 e^{6 \phi} \partial_t \phi = 6 W^(-3) \partial_t \phi
        # But \partial_t phi = -\partial_t cf / (2 cf)  (as described above), so
        #   alpha_rhs = 6 e^{6 \phi} \partial_t \phi
        #             = 6 W^(-3) (-\partial_t W / (2 W))
        #             = -3 cf^(-4) cf_rhs
        if par.parval_from_str("BSSN.BSSN_unrescaled_and_barred_vars::ConformalFactor") == "W":
            alpha_rhs = -3*cf_rhs/(cf*cf*cf*cf)
        else:
            print("Error LapseEvolutionOption==MaximalSlicing unsupported for ConformalFactor!=W")
            exit(1)
    if par.parval_from_str("BSSN.BSSN_RHSs::LapseEvolutionOption") == "Frozen":
        alpha_rhs = sp.sympify(0)
            
    # Step 14: Set \partial_t \beta^i
    beta_rhsU = ixp.zerorank1()
    if par.parval_from_str("BSSN.BSSN_RHSs::ShiftEvolutionOption") == "GammaDriving2ndOrder_NoCovariant":
        # Step 14 Option 1: \partial_t \beta^i = \beta^j \beta^i_{,j} + B^i
        # First define BU, in terms of rescaled variable \bet^i
        BU = ixp.zerorank1()
        for i in range(DIM):
            BU[i] = betU[i]*rfm.ReU[i]

        # Then compute right-hand side:
        for i in range(DIM):
            beta_rhsU[i] += BU[i]
            for j in range(DIM):
                beta_rhsU[i] += betaU[j]*betaU_dupD[i][j]
    if par.parval_from_str("BSSN.BSSN_RHSs::ShiftEvolutionOption") == "GammaDriving2ndOrder_Covariant":
        # Step 14 Option 2: \partial_t \beta^i = \left[\beta^j \bar{D}_j \beta^i\right] + B^{i}
        # First define BU, in terms of rescaled variable \bet^i
        BU = ixp.zerorank1()
        for i in range(DIM):
            BU[i] = betU[i]*rfm.ReU[i]

        # Then compute right-hand side:
        # Term 1: \beta^j \beta^i_{,j}
        for i in range(DIM):
            for j in range(DIM):
                beta_rhsU[i] += betaU[j]*betaU_dupD[i][j]

        # Term 2: \beta^j \bar{\Gamma}^i_{mj} \beta^m
        for i in range(DIM):
            for j in range(DIM):
                for m in range(DIM):
                    beta_rhsU[i] += betaU[j]*GammabarUDD[i][m][j]*betaU[m]
        # Term 3: B^i
        for i in range(DIM):
            beta_rhsU[i] += BU[i]
    if par.parval_from_str("BSSN.BSSN_RHSs::ShiftEvolutionOption") == "Frozen":
        # Step 14 Option 3: \partial_t \beta^i = 0
        for i in range(DIM):
            beta_rhsU[i] = sp.sympify(0)
        
    # Step 15: Evaluate \partial_t B^i

    # Step 15a: Declare free parameter eta and B_rhsU
    eta = par.Cparameters("REAL",thismodule,["eta"])
    B_rhsU = ixp.zerorank1()

    if par.parval_from_str("BSSN.BSSN_RHSs::ShiftEvolutionOption") == "GammaDriving2ndOrder_NoCovariant":
        # Step 15: Non-covariant option:
        #  \partial_t B^i = \beta^j \partial_j B^i
        #                 + \frac{3}{4} \partial_{0} \bar{\Lambda}^{i} - \eta B^{i}
        # Step 15b: Define BU_dupD, in terms of derivative of rescaled variable \bet^i
        BU_dupD = ixp.zerorank2()
        betU_dupD = ixp.declarerank2("betU_dupD","nosym")
        for i in range(DIM):
            for j in range(DIM):
                BU_dupD[i][j] = betU_dupD[i][j]*rfm.ReU[i] + betU[i]*rfm.ReUdD[i][j]

        # Step 15c: Compute \partial_0 \bar{\Lambda}^i = (\partial_t - \beta^i \partial_i) \bar{\Lambda}^j 
        Lambar_partial0 = ixp.zerorank1()
        for i in range(DIM):
            Lambar_partial0[i] = Lambar_rhsU[i]
        for i in range(DIM):
            for j in range(DIM):
                Lambar_partial0[j] += -betaU[i]*LambarU_dupD[j][i]

        # Step 15d: Evaluate RHS of B^i:
        for i in range(DIM):
            B_rhsU[i] += sp.Rational(3,4)*Lambar_partial0[i] - eta*BU[i]
            for j in range(DIM):
                B_rhsU[i] += betaU[j]*BU_dupD[i][j]

    if par.parval_from_str("BSSN.BSSN_RHSs::ShiftEvolutionOption") == "GammaDriving2ndOrder_Covariant":
        # Step 15: Covariant option:
        #  \partial_t B^i = \beta^j \bar{D}_j B^i
        #               + \frac{3}{4} ( \partial_t \bar{\Lambda}^{i} - \beta^j \bar{D}_j \bar{\Lambda}^{i} ) 
        #               - \eta B^{i}
        #                 = \beta^j B^i_{,j} + \beta^j \bar{\Gamma}^i_{mj} B^m
        #               + \frac{3}{4}[ \partial_t \bar{\Lambda}^{i} 
        #                            - \beta^j (\bar{\Lambda}^i_{,j} + \bar{\Gamma}^i_{mj} \bar{\Lambda}^m)] 
        #               - \eta B^{i}
        # Term 1, part a: First compute B^i_{,j} using upwinded derivative
        BU_dupD = ixp.zerorank2()
        betU_dupD = ixp.declarerank2("betU_dupD","nosym")
        for i in range(DIM):
            for j in range(DIM):
                BU_dupD[i][j] = betU_dupD[i][j]*rfm.ReU[i] + betU[i]*rfm.ReUdD[i][j]
        # Term 1: \beta^j B^i_{,j}
        for i in range(DIM):
            for j in range(DIM):
                B_rhsU[i] += betaU[j]*BU_dupD[i][j]
        # Term 2: \beta^j \bar{\Gamma}^i_{mj} B^m
        for i in range(DIM):
            for j in range(DIM):
                for m in range(DIM):
                    B_rhsU[i] += betaU[j]*GammabarUDD[i][m][j]*BU[m]
        # Term 3: \frac{3}{4}\partial_t \bar{\Lambda}^{i}
        for i in range(DIM):
            B_rhsU[i] += sp.Rational(3,4)*Lambar_rhsU[i]
        # Term 4: -\frac{3}{4}\beta^j \bar{\Lambda}^i_{,j}
        for i in range(DIM):
            for j in range(DIM):
                B_rhsU[i] += -sp.Rational(3,4)*betaU[j]*LambarU_dupD[i][j]
        # Term 5: -\frac{3}{4}\beta^j \bar{\Gamma}^i_{mj} \bar{\Lambda}^m
        for i in range(DIM):
            for j in range(DIM):
                for m in range(DIM):
                    B_rhsU[i] += -sp.Rational(3,4)*betaU[j]*GammabarUDD[i][m][j]*LambarU[m]
        # Term 6: - \eta B^i
        for i in range(DIM):
            B_rhsU[i] += -eta*BU[i]

    if par.parval_from_str("BSSN.BSSN_RHSs::ShiftEvolutionOption") == "Frozen":
        # Step 14 Option 3: \partial_t B^i = 0
        for i in range(DIM):
            B_rhsU[i] = sp.sympify(0)

    # Step 16: Rescale the RHS quantities so that the evolved
    #          variables are smooth across coord singularities
    global h_rhsDD,a_rhsDD,lambda_rhsU,vet_rhsU,bet_rhsU
    h_rhsDD     = ixp.zerorank2()
    a_rhsDD     = ixp.zerorank2()
    lambda_rhsU = ixp.zerorank1()
    vet_rhsU    = ixp.zerorank1()
    bet_rhsU    = ixp.zerorank1()
    for i in range(DIM):
        lambda_rhsU[i] = Lambar_rhsU[i] / rfm.ReU[i]
        vet_rhsU[i]    =   beta_rhsU[i] / rfm.ReU[i]
        bet_rhsU[i]    =      B_rhsU[i] / rfm.ReU[i]
        for j in range(DIM):
            h_rhsDD[i][j] = gammabar_rhsDD[i][j] / rfm.ReDD[i][j]
            a_rhsDD[i][j] =     Abar_rhsDD[i][j] / rfm.ReDD[i][j]
