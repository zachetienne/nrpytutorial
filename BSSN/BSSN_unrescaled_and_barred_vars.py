# This module defines barred BSSN variables in terms of rescaled quantities

# Author: Zachariah B. Etienne
#         zachetie **at** gmail **dot* com

# Step P1: import all needed modules from NRPy+:
import NRPy_param_funcs as par
import sympy as sp
import indexedexp as ixp
import grid as gri
import reference_metric as rfm
import BSSN.BSSN_rescaled_vars as Brv

thismodule = __name__
par.initialize_param(par.glb_param("char", thismodule, "ConformalFactor", "W"))
par.initialize_param(par.glb_param("bool", thismodule, "detgbarOverdetghat_equals_one", "True"))

def BSSN_barred_variables():
    global gammabarDD, gammabarUU, gammabarDD_dD, gammabarDD_dupD
    global detgbarOverdetghat,detgammabar, detgammabar_dD, detgammabar_dDD
    global GammabarUDD
    global AbarDD, AbarUU

    # Step 1: All barred quantities are defined in terms of BSSN rescaled gridfunctions,
    #         which we declare here in case they haven't yet been declared elsewhere.
    Brv.declare_BSSN_rescaled_gridfunctions_if_not_declared_already()
    hDD = Brv.hDD
    aDD = Brv.aDD

    # Step 2: Set spatial dimension (must be 3 for BSSN, a 3+1 decomposition of Einstein's equations of general relativity)
    DIM = 3

    # Step 3: gammabarDD and gammabarUU:
    gammabarDD = ixp.zerorank2()
    for i in range(DIM):
        for j in range(DIM):
            gammabarDD[i][j] = hDD[i][j]*rfm.ReDD[i][j] + rfm.ghatDD[i][j]
    gammabarUU, dummydet = ixp.symm_matrix_inverter3x3(gammabarDD)

    if par.parval_from_str("BSSN.BSSN_unrescaled_and_barred_vars::detgbarOverdetghat_equals_one") == "False":
        detgbarOverdetghat = gri.register_gridfunctions("AUX", ["detgbarOverdetghat"])
        detgbarOverdetghatInitial = gri.register_gridfunctions("AUX", ["detgbarOverdetghatInitial"])
        print("Error: detgbarOverdetghat_equals_one=\"False\" is not fully implemented yet.")
        exit(1)
    else:
        detgbarOverdetghat = sp.sympify(1)

    # Step 8d: Define detgammabar, detgammabar_dD, and detgammabar_dDD (needed for \partial_t \bar{\Lambda}^i below)
    detgammabar = detgbarOverdetghat * rfm.detgammahat
    detgammabar_dD = ixp.zerorank1()
    if par.parval_from_str("BSSN.BSSN_unrescaled_and_barred_vars::detgbarOverdetghat_equals_one") == "False":
        detgbarOverdetghat_dD = ixp.declarerank1("detgbarOverdetghat_dD")
    else:
        detgbarOverdetghat_dD = ixp.zerorank1()
    for i in range(DIM):
        detgammabar_dD[i] = detgbarOverdetghat_dD[i] * rfm.detgammahat + detgbarOverdetghat * rfm.detgammahatdD[i]

    detgammabar_dDD = ixp.zerorank2()
    if par.parval_from_str("BSSN.BSSN_unrescaled_and_barred_vars::detgbarOverdetghat_equals_one") == "False":
        detgbarOverdetghat_dDD = ixp.declarerank2("detgbarOverdetghat_dDD", "sym01")
    else:
        detgbarOverdetghat_dDD = ixp.zerorank2()

    for i in range(DIM):
        for j in range(DIM):
            detgammabar_dDD[i][j] = detgbarOverdetghat_dDD[i][j] * rfm.detgammahat + \
                                    detgbarOverdetghat_dD[i] * rfm.detgammahatdD[j] + \
                                    detgbarOverdetghat_dD[j] * rfm.detgammahatdD[i] + \
                                    detgbarOverdetghat * rfm.detgammahatdDD[i][j]

    # Step 4: gammabarDD_dD = gammabarDDdD[i][j][k]
    #          = h_{ij,k} \text{ReDD[i][j]} + h_{ij} \text{ReDDdD[i][j][k]} + \hat{\gamma}_{ij,k}.
    gammabarDD_dD = ixp.zerorank3()
    hDD_dD = ixp.declarerank3("hDD_dD","sym01") # Needed for \bar{\gamma}_{ij} RHS
    hDD_dupD = ixp.declarerank3("hDD_dupD","sym01") # Needed for \bar{\gamma}_{ij} RHS
    gammabarDD_dupD = ixp.zerorank3()  # Needed for \bar{\gamma}_{ij} RHS
    for i in range(DIM):
        for j in range(DIM):
            for k in range(DIM):
                gammabarDD_dD[i][j][k] = hDD_dD[i][j][k]*rfm.ReDD[i][j] + hDD[i][j]*rfm.ReDDdD[i][j][k] \
                                       + rfm.ghatDDdD[i][j][k]
                gammabarDD_dupD[i][j][k] = hDD_dupD[i][j][k]*rfm.ReDD[i][j] + hDD[i][j]*rfm.ReDDdD[i][j][k] \
                                         + rfm.ghatDDdD[i][j][k]

    # Step 7b: Define barred Christoffel symbol \bar{\Gamma}^{i}_{jk} = GammabarUDD[i][j][k]
    GammabarUDD = ixp.zerorank3()
    for i in range(DIM):
        for k in range(DIM):
            for l in range(DIM):
                for m in range(DIM):
                    GammabarUDD[i][k][l] += (sp.Rational(1,2))*gammabarUU[i][m]* \
                                            (gammabarDD_dD[m][k][l] + gammabarDD_dD[m][l][k] - gammabarDD_dD[k][l][m])

    AbarDD = ixp.zerorank2()
    for i in range(DIM):
        for j in range(DIM):
            AbarDD[i][j] = aDD[i][j] * rfm.ReDD[i][j]

    AbarUU = ixp.zerorank2() # Needed also for \partial_t \bar{\Lambda}^i
    for i in range(DIM):
        for j in range(DIM):
            for k in range(DIM):
                for l in range(DIM):
                    AbarUU[i][j] += gammabarUU[i][k]*gammabarUU[j][l]*AbarDD[k][l]

def RicciBar__gammabarDD_dHatD__DGammaUDD__DGammaU():
    # Step 0: Declare as globals all quantities computed in this function.
    global RbarDD,gammabarDD_dHatD,DGammaUDD,DGammaU

    # Step 1: All barred quantities are defined in terms of BSSN rescaled gridfunctions,
    #         which we declare here in case they haven't yet been declared elsewhere.
    Brv.declare_BSSN_rescaled_gridfunctions_if_not_declared_already()
    lambdaU = Brv.lambdaU
    hDD     = Brv.hDD

    # Step 2: Set spatial dimension (must be 3 for BSSN, a 3+1 decomposition of Einstein's equations of general relativity)
    DIM = 3

    # Step 5a: Define \varepsilon_{ij} = epsDD[i][j]
    epsDD = ixp.zerorank3()
    for i in range(DIM):
        for j in range(DIM):
            epsDD[i][j] = hDD[i][j]*rfm.ReDD[i][j]

    # Step 5b: Define epsDD_dD[i][j][k]
    hDD_dD = ixp.declarerank3("hDD_dD","sym01")
    epsDD_dD = ixp.zerorank3()
    for i in range(DIM):
        for j in range(DIM):
            for k in range(DIM):
                epsDD_dD[i][j][k] = hDD_dD[i][j][k]*rfm.ReDD[i][j] + hDD[i][j]*rfm.ReDDdD[i][j][k]

    # Step 5c: Define epsDD_dDD[i][j][k][l]
    hDD_dDD = ixp.declarerank4("hDD_dDD","sym01_sym23")
    epsDD_dDD = ixp.zerorank4()
    for i in range(DIM):
        for j in range(DIM):
            for k in range(DIM):
                for l in range(DIM):
                    epsDD_dDD[i][j][k][l] = hDD_dDD[i][j][k][l]*rfm.ReDD[i][j] + \
                                            hDD_dD[i][j][k]*rfm.ReDDdD[i][j][l] + \
                                            hDD_dD[i][j][l]*rfm.ReDDdD[i][j][k] + \
                                            hDD[i][j]*rfm.ReDDdDD[i][j][k][l]

    # Step 5d: Define DhatgammabarDDdD[i][j][l] = \bar{\gamma}_{ij;\hat{l}}
    # \bar{\gamma}_{ij;\hat{l}} = \varepsilon_{i j,l}
    #                           - \hat{\Gamma}^m_{i l} \varepsilon_{m j}
    #                           - \hat{\Gamma}^m_{j l} \varepsilon_{i m}
    gammabarDD_dHatD = ixp.zerorank3()
    for i in range(DIM):
        for j in range(DIM):
            for l in range(DIM):
                gammabarDD_dHatD[i][j][l] = epsDD_dD[i][j][l]
                for m in range(DIM):
                    gammabarDD_dHatD[i][j][l] += - rfm.GammahatUDD[m][i][l]*epsDD[m][j] \
                                                 - rfm.GammahatUDD[m][j][l]*epsDD[i][m]

    # Step 5e: Define \bar{\gamma}_{ij;\hat{l},k} = DhatgammabarDD_dHatD_dD[i][j][l][k]:
    #        \bar{\gamma}_{ij;\hat{l},k} = \varepsilon_{ij,lk}
    #                                      - \hat{\Gamma}^m_{i l,k} \varepsilon_{m j}
    #                                      - \hat{\Gamma}^m_{i l} \varepsilon_{m j,k}
    #                                      - \hat{\Gamma}^m_{j l,k} \varepsilon_{i m}
    #                                      - \hat{\Gamma}^m_{j l} \varepsilon_{i m,k}
    gammabarDD_dHatD_dD = ixp.zerorank4()
    for i in range(DIM):
        for j in range(DIM):
            for l in range(DIM):
                for k in range(DIM):
                    gammabarDD_dHatD_dD[i][j][l][k] = epsDD_dDD[i][j][l][k]
                    for m in range(DIM):
                        gammabarDD_dHatD_dD[i][j][l][k] += -rfm.GammahatUDDdD[m][i][l][k]*epsDD[m][j]  \
                                                           -rfm.GammahatUDD[m][i][l]*epsDD_dD[m][j][k] \
                                                           -rfm.GammahatUDDdD[m][j][l][k]*epsDD[i][m]  \
                                                           -rfm.GammahatUDD[m][j][l]*epsDD_dD[i][m][k]

    # Step 5f: Define \bar{\gamma}_{ij;\hat{l}\hat{k}} = DhatgammabarDD_dHatDD[i][j][l][k]
    #          \bar{\gamma}_{ij;\hat{l}\hat{k}} = \partial_k \hat{D}_{l} \varepsilon_{i j}
    #                                           - \hat{\Gamma}^m_{lk} \left(\hat{D}_{m} \varepsilon_{i j}\right)
    #                                           - \hat{\Gamma}^m_{ik} \left(\hat{D}_{l} \varepsilon_{m j}\right)
    #                                           - \hat{\Gamma}^m_{jk} \left(\hat{D}_{l} \varepsilon_{i m}\right)
    gammabarDD_dHatDD = ixp.zerorank4()
    for i in range(DIM):
        for j in range(DIM):
            for l in range(DIM):
                for k in range(DIM):
                    gammabarDD_dHatDD[i][j][l][k] = gammabarDD_dHatD_dD[i][j][l][k]
                    for m in range(DIM):
                        gammabarDD_dHatDD[i][j][l][k] += - rfm.GammahatUDD[m][l][k]*gammabarDD_dHatD[i][j][m] \
                                                         - rfm.GammahatUDD[m][i][k]*gammabarDD_dHatD[m][j][l] \
                                                         - rfm.GammahatUDD[m][j][k]*gammabarDD_dHatD[i][m][l]

    # Step 5h: Add the first term to RbarDD:
    #         - \frac{1}{2} \bar{\gamma}^{k l} \hat{D}_{k} \hat{D}_{l} \bar{\gamma}_{i j}
    RbarDD = ixp.zerorank2()
    RbarDDpiece = ixp.zerorank2()
    for i in range(DIM):
        for j in range(DIM):
            for k in range(DIM):
                for l in range(DIM):
                    RbarDD[i][j] += -sp.Rational(1,2) * gammabarUU[k][l]*gammabarDD_dHatDD[i][j][l][k]
                    RbarDDpiece[i][j] += -sp.Rational(1,2) * gammabarUU[k][l]*gammabarDD_dHatDD[i][j][l][k]

    # Step 6a: Second term of RhatDD: compute \hat{D}_{j} \bar{\Lambda}^{k} = LambarU_dHatD[k][j]
    lambdaU_dD = ixp.declarerank2("lambdaU_dD","nosym")
    LambarU_dHatD = ixp.zerorank2()
    for j in range(DIM):
        for k in range(DIM):
            LambarU_dHatD[k][j] = lambdaU_dD[k][j]*rfm.ReU[k] + lambdaU[k]*rfm.ReUdD[k][j]
            for m in range(DIM):
                LambarU_dHatD[k][j] += rfm.GammahatUDD[k][m][j]*lambdaU[m]*rfm.ReU[m]

    # Step 6b: Add the second term to the Ricci tensor
    for i in range(DIM):
        for j in range(DIM):
            for k in range(DIM):
                RbarDD[i][j] += sp.Rational(1,2) * (gammabarDD[k][i]*LambarU_dHatD[k][j] + \
                                                    gammabarDD[k][j]*LambarU_dHatD[k][i])


    # Step 7c: Define \Delta^i_{jk} = \bar{\Gamma}^i_{jk} - \hat{\Gamma}^i_{jk} = DGammaUDD[i][j][k]
    DGammaUDD = ixp.zerorank3()
    for i in range(DIM):
        for j in range(DIM):
            for k in range(DIM):
                DGammaUDD[i][j][k] = GammabarUDD[i][j][k] - rfm.GammahatUDD[i][j][k]

    # Step 7d: Define \Delta^i = \bar{\gamma}^{jk} \Delta^i_{jk}
    DGammaU = ixp.zerorank1()
    for i in range(DIM):
        for j in range(DIM):
            for k in range(DIM):
                DGammaU[i] += gammabarUU[j][k] * DGammaUDD[i][j][k]

    # Step 7e: Define \Delta_{ijk} = \bar{\gamma}_{im} \Delta^m_{jk}
    DGammaDDD = ixp.zerorank3()
    for i in range(DIM):
        for j in range(DIM):
            for k in range(DIM):
                for m in range(DIM):
                    DGammaDDD[i][j][k] += gammabarDD[i][m] * DGammaUDD[m][j][k]

    # Step 7e: Add third term to Ricci tensor:
    #      \Delta^{k} \Delta_{(i j) k} = 1/2 \Delta^{k} (\Delta_{i j k} + \Delta_{j i k})
    for i in range(DIM):
        for j in range(DIM):
            for k in range(DIM):
                RbarDD[i][j] += sp.Rational(1,2) * DGammaU[k] * (DGammaDDD[i][j][k] + DGammaDDD[j][i][k])

    # Step 7f: Add remaining terms to Ricci tensor:
    # \bar{\gamma}^{k l} (\Delta^{m}_{k i} \Delta_{j m l}
    #                   + \Delta^{m}_{k j} \Delta_{i m l}
    #                   + \Delta^{m}_{i k} \Delta_{m j l})
    for i in range(DIM):
        for j in range(DIM):
            for k in range(DIM):
                for l in range(DIM):
                    for m in range(DIM):
                        RbarDD[i][j] += gammabarUU[k][l] * (DGammaUDD[m][k][i]*DGammaDDD[j][m][l] +
                                                            DGammaUDD[m][k][j]*DGammaDDD[i][m][l] +
                                                            DGammaUDD[m][i][k]*DGammaDDD[m][j][l])

    #return RbarDD, gammabarDD_dHatD, DGammaUDD, DGammaU

def betaUbar_and_derivs():
    # Step 0: Declare as globals all quantities computed in this function.
    global betaU,betaU_dD,betaU_dupD,betaU_dDD

    # Step 1: All barred quantities are defined in terms of BSSN rescaled gridfunctions,
    #         which we declare here in case they haven't yet been declared elsewhere.
    Brv.declare_BSSN_rescaled_gridfunctions_if_not_declared_already()
    vetU  = Brv.vetU

    # Step 2: Set spatial dimension (must be 3 for BSSN, a 3+1 decomposition of Einstein's equations of general relativity)
    DIM = 3

    # Step 8a: Define \beta^i and \beta^i_{,k} in terms of rescaled quantity vetU[i] and vetU_dD[i][j]:
    betaU = ixp.zerorank1()
    for i in range(DIM):
        betaU[i] = vetU[i] * rfm.ReU[i]

    vetU_dD = ixp.declarerank2("vetU_dD", "nosym")
    vetU_dupD = ixp.declarerank2("vetU_dupD", "nosym")  # Needed for \beta^i RHS
    vetU_dDD = ixp.declarerank3("vetU_dDD", "sym12")  # Needed for \bar{\Lambda}^i RHS
    betaU_dD = ixp.zerorank2()
    betaU_dupD = ixp.zerorank2()  # Needed for \beta^i RHS
    betaU_dDD = ixp.zerorank3()  # Needed for \bar{\Lambda}^i RHS
    for i in range(DIM):
        for j in range(DIM):
            betaU_dD[i][j] = vetU_dD[i][j] * rfm.ReU[i] + vetU[i] * rfm.ReUdD[i][j]
            betaU_dupD[i][j] = vetU_dupD[i][j] * rfm.ReU[i] + vetU[i] * rfm.ReUdD[i][j]  # Needed for \beta^i RHS
            for k in range(DIM):
                # Needed for \bar{\Lambda}^i RHS:
                betaU_dDD[i][j][k] = vetU_dDD[i][j][k] * rfm.ReU[i] + vetU_dD[i][j] * rfm.ReUdD[i][k] + \
                                     vetU_dD[i][k] * rfm.ReUdD[i][j] + vetU[i] * rfm.ReUdDD[i][j][k]

def phi_and_derivs():
    # Step 0: Declare as globals all quantities computed in this function.
    global phi_dD, phi_dupD, phi_dDD, exp_m4phi, phi_dBarD, phi_dBarDD

    # Step 1: All barred quantities are defined in terms of BSSN rescaled gridfunctions,
    #         which we declare here in case they haven't yet been declared elsewhere.
    Brv.declare_BSSN_rescaled_gridfunctions_if_not_declared_already()
    cf = Brv.cf

    # Step 2: Set spatial dimension (must be 3 for BSSN, a 3+1 decomposition of Einstein's equations of general relativity)
    DIM = 3

    # Step 9c: Define partial derivatives of \phi in terms of evolved quantity "cf":
    cf_dD = ixp.declarerank1("cf_dD")
    cf_dupD = ixp.declarerank1("cf_dupD")  # Needed for \partial_t \phi next.
    cf_dDD = ixp.declarerank2("cf_dDD", "sym01")
    phi_dD = ixp.zerorank1()
    phi_dupD = ixp.zerorank1()
    phi_dDD = ixp.zerorank2()

    exp_m4phi = sp.sympify(0)
    if par.parval_from_str("BSSN.BSSN_unrescaled_and_barred_vars::ConformalFactor") == "phi":
        for i in range(DIM):
            phi_dD[i] = cf_dD[i]
            phi_dupD[i] = cf_dupD[i]
            for j in range(DIM):
                phi_dDD[i][j] = cf_dDD[i][j]
        exp_m4phi = sp.exp(-4 * cf)
    elif par.parval_from_str("BSSN.BSSN_unrescaled_and_barred_vars::ConformalFactor") == "W":
        # \partial_i W = \partial_i (e^{-2 phi}) = -2 e^{-2 phi} \partial_i phi
        # -> \partial_i phi = -\partial_i cf / (2 cf)
        for i in range(DIM):
            phi_dD[i] = - cf_dD[i] / (2 * cf)
            phi_dupD[i] = - cf_dupD[i] / (2 * cf)
            for j in range(DIM):
                # \partial_j \partial_i phi = - \partial_j [\partial_i cf / (2 cf)]
                #                           = - cf_{,ij} / (2 cf) + \partial_i cf \partial_j cf / (2 cf^2)
                phi_dDD[i][j] = (- cf_dDD[i][j] + cf_dD[i] * cf_dD[j] / cf) / (2 * cf)
        exp_m4phi = cf * cf
    elif par.parval_from_str("BSSN.BSSN_unrescaled_and_barred_vars::ConformalFactor") == "chi":
        # \partial_i chi = \partial_i (e^{-4 phi}) = -4 e^{-4 phi} \partial_i phi
        # -> \partial_i phi = -\partial_i cf / (4 cf)
        for i in range(DIM):
            phi_dD[i] = - cf_dD[i] / (4 * cf)
            phi_dupD[i] = - cf_dupD[i] / (4 * cf)
            for j in range(DIM):
                # \partial_j \partial_i phi = - \partial_j [\partial_i cf / (4 cf)]
                #                           = - cf_{,ij} / (4 cf) + \partial_i cf \partial_j cf / (4 cf^2)
                phi_dDD[i][j] = (- cf_dDD[i][j] + cf_dD[i] * cf_dD[j] / cf) / (4 * cf)
        exp_m4phi = cf
    else:
        print("Error: ConformalFactor == " + par.parval_from_str("BSSN_unrescaled_and_barred_vars::ConformalFactor") + " unsupported!")
        exit(1)

    # Step 9d: Define phi_dBarD = phi_dD (since phi is a scalar) and phi_dBarDD (covariant derivative)
    #          \bar{D}_i \bar{D}_j \phi = \phi_{;\bar{i}\bar{j}} = \bar{D}_i \phi_{,j}
    #                                   = \phi_{,ij} - \bar{\Gamma}^k_{ij} \phi_{,k}
    phi_dBarD = phi_dD
    phi_dBarDD = ixp.zerorank2()
    for i in range(DIM):
        for j in range(DIM):
            phi_dBarDD[i][j] = phi_dDD[i][j]
            for k in range(DIM):
                phi_dBarDD[i][j] += - GammabarUDD[k][i][j] * phi_dD[k]
