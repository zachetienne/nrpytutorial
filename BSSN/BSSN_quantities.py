# This module provides functions that declare and define useful BSSN quantities

# Author: Zachariah B. Etienne
#         zachetie **at** gmail **dot* com

# Step 1: Import all needed modules from NRPy+:
import NRPy_param_funcs as par
import sympy as sp
import indexedexp as ixp
import grid as gri
import reference_metric as rfm

# Step 1.a: Set the coordinate system for the numerical grid
#  DO NOT SET IN STANDALONE PYTHON MODULE
# par.set_parval_from_str("reference_metric::CoordSystem","Spherical")

# Step 1.b: Given the chosen coordinate system, set up 
#           corresponding reference metric and needed
#           reference metric quantities
# The following function call sets up the reference metric
#    and related quantities, including rescaling matrices ReDD,
#    ReU, and hatted quantities.
#  DO NOT CALL IN STANDALONE PYTHON MODULE
# rfm.reference_metric()

# Step 1.c: Set spatial dimension (must be 3 for BSSN, as BSSN is 
#           a 3+1-dimensional decomposition of the general 
#           relativistic field equations)
#  DO NOT CALL IN STANDALONE PYTHON MODULE
# DIM = 3
# par.set_parval_from_str("grid::DIM",DIM)

# Step 1.d: Declare/initialize parameters for this module
thismodule = __name__
par.initialize_param(par.glb_param("char", thismodule, "EvolvedConformalFactor_cf", "W"))
par.initialize_param(par.glb_param("bool", thismodule, "detgbarOverdetghat_equals_one", "True"))

def declare_BSSN_gridfunctions_if_not_declared_already():
    # Step 2: Register all needed BSSN gridfunctions.

    # Declare as globals all variables that may be
    # used outside this function
    global hDD,aDD,lambdaU,vetU,betU,trK,cf,alpha

    #   Check to see if this function has already been called.
    #   If so, do not register the gridfunctions again!
    for i in range(len(gri.glb_gridfcs_list)):
        if "hDD00" in gri.glb_gridfcs_list[i].name:
            hDD = ixp.declarerank2("hDD", "sym01")
            aDD = ixp.declarerank2("aDD", "sym01")
            lambdaU = ixp.declarerank1("lambdaU")
            vetU = ixp.declarerank1("vetU")
            betU = ixp.declarerank1("betU")
            trK, cf, alpha = sp.symbols('trK cf alpha', real=True)
            return hDD, aDD, lambdaU, vetU, betU, trK, cf, alpha

    # Step 2.a: Register indexed quantities, using ixp.register_... functions
    hDD = ixp.register_gridfunctions_for_single_rank2("EVOL", "hDD", "sym01")
    aDD = ixp.register_gridfunctions_for_single_rank2("EVOL", "aDD", "sym01")
    lambdaU = ixp.register_gridfunctions_for_single_rank1("EVOL", "lambdaU")
    vetU = ixp.register_gridfunctions_for_single_rank1("EVOL", "vetU")
    betU = ixp.register_gridfunctions_for_single_rank1("EVOL", "betU")

    # Step 2.b: Register scalar quantities, using gri.register_gridfunctions()
    trK, cf, alpha = gri.register_gridfunctions("EVOL", ["trK", "cf", "alpha"])

    return hDD, aDD, lambdaU, vetU, betU, trK, cf, alpha

# Step 3: Define all basic conformal BSSN tensors
#        gammabarDD,AbarDD,LambdabarU,betaU,BU
#        in terms of BSSN gridfunctions.
def BSSN_basic_tensors():

    # Step 3.a: Declare as globals all variables that may be used
    #           outside this function, declare BSSN gridfunctions
    #           if not defined already, and set DIM=3.
    global gammabarDD,AbarDD,LambdabarU,betaU,BU
    hDD, aDD, lambdaU, vetU, betU, trK, cf, alpha = declare_BSSN_gridfunctions_if_not_declared_already()
    DIM = 3

    # Step 3.a.i: gammabarDD and AbarDD:
    gammabarDD = ixp.zerorank2()
    AbarDD = ixp.zerorank2()
    for i in range(DIM):
        for j in range(DIM):
            # gammabar_{ij}  = h_{ij}*ReDD[i][j] + gammahat_{ij}
            gammabarDD[i][j] = hDD[i][j] * rfm.ReDD[i][j] + rfm.ghatDD[i][j]
            # Abar_{ij}      = a_{ij}*ReDD[i][j]
            AbarDD[i][j] = aDD[i][j] * rfm.ReDD[i][j]

    # Step 3.a.ii: LambdabarU, betaU, and BU:
    LambdabarU = ixp.zerorank1()
    betaU = ixp.zerorank1()
    BU = ixp.zerorank1()
    for i in range(DIM):
        LambdabarU[i] = lambdaU[i] * rfm.ReU[i]
        betaU[i] = vetU[i] * rfm.ReU[i]
        BU[i] = betU[i] * rfm.ReU[i]

# Step 4: gammabarUU and spatial derivatives of gammabarDD,
#         including GammabarUDD
def gammabar__inverse_and_derivs():
    # Step 4.a: Declare as globals all expressions that may be used
    #           outside this function, declare BSSN gridfunctions
    #           if not defined already, and set DIM=3.
    global gammabarUU, gammabarDD_dD, gammabarDD_dupD, gammabarDD_dDD, GammabarUDD
    hDD, aDD, lambdaU, vetU, betU, trK, cf, alpha = declare_BSSN_gridfunctions_if_not_declared_already()
    DIM = 3
    # This function needs gammabarDD, defined in BSSN_basic_tensors()
    BSSN_basic_tensors()

    # Step 4.a.i: gammabarUU:
    gammabarUU, dummydet = ixp.symm_matrix_inverter3x3(gammabarDD)

    # Step 4.b.i: gammabarDDdD[i][j][k]
    #               = \hat{\gamma}_{ij,k} + h_{ij,k} \text{ReDD[i][j]} + h_{ij} \text{ReDDdD[i][j][k]}.
    gammabarDD_dD = ixp.zerorank3()
    gammabarDD_dupD = ixp.zerorank3()
    hDD_dD = ixp.declarerank3("hDD_dD", "sym01")
    hDD_dupD = ixp.declarerank3("hDD_dupD", "sym01")
    for i in range(DIM):
        for j in range(DIM):
            for k in range(DIM):
                gammabarDD_dD[i][j][k] = rfm.ghatDDdD[i][j][k] + \
                                         hDD_dD[i][j][k] * rfm.ReDD[i][j] + hDD[i][j] * rfm.ReDDdD[i][j][k]

                # Compute associated upwinded derivative, needed for the \bar{\gamma}_{ij} RHS
                gammabarDD_dupD[i][j][k] = rfm.ghatDDdD[i][j][k] + \
                                           hDD_dupD[i][j][k] * rfm.ReDD[i][j] + hDD[i][j] * rfm.ReDDdD[i][j][k]

    # Step 4.b.ii: Compute gammabarDD_dDD in terms of the rescaled BSSN quantity hDD
    #      and its derivatives, as well as the reference metric and rescaling
    #      matrix, and its derivatives (expression given below):
    hDD_dDD = ixp.declarerank4("hDD_dDD", "sym01_sym23")
    gammabarDD_dDD = ixp.zerorank4()
    for i in range(DIM):
        for j in range(DIM):
            for k in range(DIM):
                for l in range(DIM):
                    # gammabar_{ij,kl} = gammahat_{ij,kl}
                    #                  + h_{ij,kl} ReDD[i][j]
                    #                  + h_{ij,k} ReDDdD[i][j][l] + h_{ij,l} ReDDdD[i][j][k]
                    #                  + h_{ij} ReDDdDD[i][j][k][l]
                    gammabarDD_dDD[i][j][k][l] = rfm.ghatDDdDD[i][j][k][l]
                    gammabarDD_dDD[i][j][k][l] += hDD_dDD[i][j][k][l] * rfm.ReDD[i][j]
                    gammabarDD_dDD[i][j][k][l] += hDD_dD[i][j][k] * rfm.ReDDdD[i][j][l] + \
                                                  hDD_dD[i][j][l] * rfm.ReDDdD[i][j][k]
                    gammabarDD_dDD[i][j][k][l] += hDD[i][j] * rfm.ReDDdDD[i][j][k][l]

    # Step 4.b.iii: Define barred Christoffel symbol \bar{\Gamma}^{i}_{kl} = GammabarUDD[i][k][l] (see expression below)
    GammabarUDD = ixp.zerorank3()
    for i in range(DIM):
        for k in range(DIM):
            for l in range(DIM):
                for m in range(DIM):
                    # Gammabar^i_{kl} = 1/2 * gammabar^{im} ( gammabar_{mk,l} + gammabar_{ml,k} - gammabar_{kl,m}):
                    GammabarUDD[i][k][l] += sp.Rational(1, 2) * gammabarUU[i][m] * \
                                            (gammabarDD_dD[m][k][l] + gammabarDD_dD[m][l][k] - gammabarDD_dD[k][l][m])

# Step 5: det(gammabarDD) and its derivatives
def detgammabar_and_derivs():
    # Step 5.a: Declare as globals all expressions that may be used
    #           outside this function, declare BSSN gridfunctions
    #           if not defined already, and set DIM=3.
    global detgammabar,detgammabar_dD,detgammabar_dDD
    hDD, aDD, lambdaU, vetU, betU, trK, cf, alpha = declare_BSSN_gridfunctions_if_not_declared_already()
    DIM = 3

    detgbarOverdetghat = sp.sympify(1)
    detgbarOverdetghat_dD = ixp.zerorank1()
    detgbarOverdetghat_dDD = ixp.zerorank2()

    if par.parval_from_str(thismodule + "::detgbarOverdetghat_equals_one") == "False":
        print("Error: detgbarOverdetghat_equals_one=\"False\" is not fully implemented yet.")
        exit(1)
    ## Approach for implementing detgbarOverdetghat_equals_one=False:
    #     detgbarOverdetghat = gri.register_gridfunctions("AUX", ["detgbarOverdetghat"])
    #     detgbarOverdetghatInitial = gri.register_gridfunctions("AUX", ["detgbarOverdetghatInitial"])
    #     detgbarOverdetghat_dD = ixp.declarerank1("detgbarOverdetghat_dD")
    #     detgbarOverdetghat_dDD = ixp.declarerank2("detgbarOverdetghat_dDD", "sym01")

    # Step 5.b: Define detgammabar, detgammabar_dD, and detgammabar_dDD (needed for \partial_t \bar{\Lambda}^i below)
    detgammabar = detgbarOverdetghat * rfm.detgammahat

    detgammabar_dD = ixp.zerorank1()
    for i in range(DIM):
        detgammabar_dD[i] = detgbarOverdetghat_dD[i] * rfm.detgammahat + detgbarOverdetghat * rfm.detgammahatdD[i]

    detgammabar_dDD = ixp.zerorank2()
    for i in range(DIM):
        for j in range(DIM):
            detgammabar_dDD[i][j] = detgbarOverdetghat_dDD[i][j] * rfm.detgammahat + \
                                    detgbarOverdetghat_dD[i] * rfm.detgammahatdD[j] + \
                                    detgbarOverdetghat_dD[j] * rfm.detgammahatdD[i] + \
                                    detgbarOverdetghat * rfm.detgammahatdDD[i][j]

# Step 6: Quantities related to conformal traceless
#         extrinsic curvature AbarDD:
#         AbarUU, AbarUD, and trAbar
def AbarUU_AbarUD_trAbar_AbarDD_dD():
    # Step 6.a: Declare as globals all expressions that may be used
    #           outside this function, declare BSSN gridfunctions
    #           if not defined already, and set DIM=3.
    global AbarUU,AbarUD,trAbar,AbarDD_dD,AbarDD_dupD
    hDD, aDD, lambdaU, vetU, betU, trK, cf, alpha = declare_BSSN_gridfunctions_if_not_declared_already()
    DIM = 3

    # Define AbarDD and gammabarDD in terms of BSSN gridfunctions
    BSSN_basic_tensors()
    # Define gammabarUU in terms of BSSN gridfunctions
    gammabar__inverse_and_derivs()


    # Step 6.a.i: Compute Abar^{ij} in terms of Abar_{ij} and gammabar^{ij}
    AbarUU = ixp.zerorank2()
    for i in range(DIM):
        for j in range(DIM):
            for k in range(DIM):
                for l in range(DIM):
                    # Abar^{ij} = gammabar^{ik} gammabar^{jl} Abar_{kl}
                    AbarUU[i][j] += gammabarUU[i][k] * gammabarUU[j][l] * AbarDD[k][l]

    # Step 6.a.ii: Compute Abar^i_j in terms of Abar_{ij} and gammabar^{ij}
    AbarUD = ixp.zerorank2()
    for i in range(DIM):
        for j in range(DIM):
            for k in range(DIM):
                # Abar^i_j = gammabar^{ik} Abar_{kj}
                AbarUD[i][j] += gammabarUU[i][k] * AbarDD[k][j]

    # Step 6.a.iii: Compute Abar^k_k = trace of Abar:
    trAbar = sp.sympify(0)
    for k in range(DIM):
        for j in range(DIM):
            # Abar^k_k = gammabar^{kj} Abar_{jk}
            trAbar += gammabarUU[k][j] * AbarDD[j][k]
            
    # Step 6.a.iv: Compute Abar_{ij,k}
    AbarDD_dD = ixp.zerorank3()
    AbarDD_dupD = ixp.zerorank3()
    aDD_dD   = ixp.declarerank3("aDD_dD"  ,"sym01")
    aDD_dupD = ixp.declarerank3("aDD_dupD","sym01")
    for i in range(DIM):
        for j in range(DIM):
            for k in range(DIM):
                AbarDD_dupD[i][j][k] = rfm.ReDDdD[i][j][k]*aDD[i][j] + rfm.ReDD[i][j]*aDD_dupD[i][j][k]
                AbarDD_dD[i][j][k]   = rfm.ReDDdD[i][j][k]*aDD[i][j] + rfm.ReDD[i][j]*aDD_dD[  i][j][k]

# Step 7: The conformal ("barred") Ricci tensor RbarDD
#         and associated quantities
def RicciBar__gammabarDD_dHatD__DGammaUDD__DGammaU():
    # Step 7.a: Declare as globals all expressions that may be used
    #           outside this function, declare BSSN gridfunctions
    #           if not defined already, and set DIM=3.
    global RbarDD,DGammaUDD,gammabarDD_dHatD,DGammaU
    hDD, aDD, lambdaU, vetU, betU, trK, cf, alpha = declare_BSSN_gridfunctions_if_not_declared_already()
    DIM = 3
    # GammabarUDD is used below, defined in
    #    gammabar__inverse_and_derivs()
    gammabar__inverse_and_derivs()

    # Step 7.a.i: Define \varepsilon_{ij} = epsDD[i][j]
    epsDD = ixp.zerorank3()
    for i in range(DIM):
        for j in range(DIM):
            epsDD[i][j] = hDD[i][j] * rfm.ReDD[i][j]

    # Step 7.a.ii: Define epsDD_dD[i][j][k]
    hDD_dD = ixp.declarerank3("hDD_dD", "sym01")
    epsDD_dD = ixp.zerorank3()
    for i in range(DIM):
        for j in range(DIM):
            for k in range(DIM):
                epsDD_dD[i][j][k] = hDD_dD[i][j][k] * rfm.ReDD[i][j] + hDD[i][j] * rfm.ReDDdD[i][j][k]

    # Step 7.a.iii: Define epsDD_dDD[i][j][k][l]
    hDD_dDD = ixp.declarerank4("hDD_dDD", "sym01_sym23")
    epsDD_dDD = ixp.zerorank4()
    for i in range(DIM):
        for j in range(DIM):
            for k in range(DIM):
                for l in range(DIM):
                    epsDD_dDD[i][j][k][l] = hDD_dDD[i][j][k][l] * rfm.ReDD[i][j] + \
                                            hDD_dD[i][j][k] * rfm.ReDDdD[i][j][l] + \
                                            hDD_dD[i][j][l] * rfm.ReDDdD[i][j][k] + \
                                            hDD[i][j] * rfm.ReDDdDD[i][j][k][l]

    # Step 7.a.iv: DhatgammabarDDdD[i][j][l] = \bar{\gamma}_{ij;\hat{l}}
    # \bar{\gamma}_{ij;\hat{l}} = \varepsilon_{i j,l}
    #                           - \hat{\Gamma}^m_{i l} \varepsilon_{m j}
    #                           - \hat{\Gamma}^m_{j l} \varepsilon_{i m}
    gammabarDD_dHatD = ixp.zerorank3()
    for i in range(DIM):
        for j in range(DIM):
            for l in range(DIM):
                gammabarDD_dHatD[i][j][l] = epsDD_dD[i][j][l]
                for m in range(DIM):
                    gammabarDD_dHatD[i][j][l] += - rfm.GammahatUDD[m][i][l] * epsDD[m][j] \
                                                 - rfm.GammahatUDD[m][j][l] * epsDD[i][m]

    # Step 7.a.v: \bar{\gamma}_{ij;\hat{l},k} = DhatgammabarDD_dHatD_dD[i][j][l][k]:
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
                        gammabarDD_dHatD_dD[i][j][l][k] += -rfm.GammahatUDDdD[m][i][l][k] * epsDD[m][j] \
                                                           - rfm.GammahatUDD[m][i][l] * epsDD_dD[m][j][k] \
                                                           - rfm.GammahatUDDdD[m][j][l][k] * epsDD[i][m] \
                                                           - rfm.GammahatUDD[m][j][l] * epsDD_dD[i][m][k]

    # Step 7.a.vi: \bar{\gamma}_{ij;\hat{l}\hat{k}} = DhatgammabarDD_dHatDD[i][j][l][k]
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
                        gammabarDD_dHatDD[i][j][l][k] += - rfm.GammahatUDD[m][l][k] * gammabarDD_dHatD[i][j][m] \
                                                         - rfm.GammahatUDD[m][i][k] * gammabarDD_dHatD[m][j][l] \
                                                         - rfm.GammahatUDD[m][j][k] * gammabarDD_dHatD[i][m][l]

    # Step 7.b: Second term of RhatDD: compute \hat{D}_{j} \bar{\Lambda}^{k} = LambarU_dHatD[k][j]
    lambdaU_dD = ixp.declarerank2("lambdaU_dD", "nosym")
    LambarU_dHatD = ixp.zerorank2()
    for j in range(DIM):
        for k in range(DIM):
            LambarU_dHatD[k][j] = lambdaU_dD[k][j] * rfm.ReU[k] + lambdaU[k] * rfm.ReUdD[k][j]
            for m in range(DIM):
                LambarU_dHatD[k][j] += rfm.GammahatUDD[k][m][j] * lambdaU[m] * rfm.ReU[m]
                
    # Step 7.c: Conformal Ricci tensor, part 3: The \Delta^{k} \Delta_{(i j) k}  
    #           + \bar{\gamma}^{k l}*(2 \Delta_{k(i}^{m} \Delta_{j) m l} 
    #           + \Delta_{i k}^{m} \Delta_{m j l}) terms

    # Step 7.c.i: Define \Delta^i_{jk} = \bar{\Gamma}^i_{jk} - \hat{\Gamma}^i_{jk} = DGammaUDD[i][j][k]
    DGammaUDD = ixp.zerorank3()
    for i in range(DIM):
        for j in range(DIM):
            for k in range(DIM):
                DGammaUDD[i][j][k] = GammabarUDD[i][j][k] - rfm.GammahatUDD[i][j][k]

    # Step 7.c.ii: Define \Delta^i = \bar{\gamma}^{jk} \Delta^i_{jk}
    DGammaU = ixp.zerorank1()
    for i in range(DIM):
        for j in range(DIM):
            for k in range(DIM):
                DGammaU[i] += gammabarUU[j][k] * DGammaUDD[i][j][k]

    # Step 7.c.iii: Define \Delta_{ijk} = \bar{\gamma}_{im} \Delta^m_{jk}
    DGammaDDD = ixp.zerorank3()
    for i in range(DIM):
        for j in range(DIM):
            for k in range(DIM):
                for m in range(DIM):
                    DGammaDDD[i][j][k] += gammabarDD[i][m] * DGammaUDD[m][j][k]
    
    # Step 7.d: Summing the terms and defining \bar{R}_{ij}
    # Step 7.d.i: Add the first term to RbarDD:
    #         Rbar_{ij} += - \frac{1}{2} \bar{\gamma}^{k l} \hat{D}_{k} \hat{D}_{l} \bar{\gamma}_{i j}
    RbarDD = ixp.zerorank2()
    RbarDDpiece = ixp.zerorank2()
    for i in range(DIM):
        for j in range(DIM):
            for k in range(DIM):
                for l in range(DIM):
                    RbarDD[i][j] += -sp.Rational(1, 2) * gammabarUU[k][l] * gammabarDD_dHatDD[i][j][l][k]
                    RbarDDpiece[i][j] += -sp.Rational(1, 2) * gammabarUU[k][l] * gammabarDD_dHatDD[i][j][l][k]

    # Step 7.d.ii: Add the second term to RbarDD:
    #         Rbar_{ij} += (1/2) * (gammabar_{ki} Lambar^k_{;\hat{j}} + gammabar_{kj} Lambar^k_{;\hat{i}})
    for i in range(DIM):
        for j in range(DIM):
            for k in range(DIM):
                RbarDD[i][j] += sp.Rational(1, 2) * (gammabarDD[k][i] * LambarU_dHatD[k][j] + \
                                                     gammabarDD[k][j] * LambarU_dHatD[k][i])

    # Step 7.d.iii: Add the remaining term to RbarDD:
    #      Rbar_{ij} += \Delta^{k} \Delta_{(i j) k} = 1/2 \Delta^{k} (\Delta_{i j k} + \Delta_{j i k})
    for i in range(DIM):
        for j in range(DIM):
            for k in range(DIM):
                RbarDD[i][j] += sp.Rational(1, 2) * DGammaU[k] * (DGammaDDD[i][j][k] + DGammaDDD[j][i][k])

    # Step 7.d.iv: Add the final term to RbarDD:
    #      Rbar_{ij} += \bar{\gamma}^{k l} (\Delta^{m}_{k i} \Delta_{j m l}
    #                   + \Delta^{m}_{k j} \Delta_{i m l}
    #                   + \Delta^{m}_{i k} \Delta_{m j l})
    for i in range(DIM):
        for j in range(DIM):
            for k in range(DIM):
                for l in range(DIM):
                    for m in range(DIM):
                        RbarDD[i][j] += gammabarUU[k][l] * (DGammaUDD[m][k][i] * DGammaDDD[j][m][l] +
                                                            DGammaUDD[m][k][j] * DGammaDDD[i][m][l] +
                                                            DGammaUDD[m][i][k] * DGammaDDD[m][j][l])

# Step 8: The unrescaled shift vector betaU spatial derivatives:
#         betaUdD & betaUdDD, written in terms of the
#         rescaled shift vector vetU
def betaU_derivs():
    # Step 8.i: Declare as globals all expressions that may be used
    #           outside this function, declare BSSN gridfunctions
    #           if not defined already, and set DIM=3.
    global betaU_dD,betaU_dupD,betaU_dDD
    hDD, aDD, lambdaU, vetU, betU, trK, cf, alpha = declare_BSSN_gridfunctions_if_not_declared_already()
    DIM = 3

    # Step 8.ii: Compute the unrescaled shift vector beta^i = ReU[i]*vet^i
    vetU_dD = ixp.declarerank2("vetU_dD", "nosym")
    vetU_dupD = ixp.declarerank2("vetU_dupD", "nosym")  # Needed for upwinded \beta^i_{,j}
    vetU_dDD = ixp.declarerank3("vetU_dDD", "sym12")  # Needed for \beta^i_{,j}
    betaU_dD = ixp.zerorank2()
    betaU_dupD = ixp.zerorank2()  # Needed for, e.g., \beta^i RHS
    betaU_dDD = ixp.zerorank3()  # Needed for, e.g., \bar{\Lambda}^i RHS
    for i in range(DIM):
        for j in range(DIM):
            betaU_dD[i][j] = vetU_dD[i][j] * rfm.ReU[i] + vetU[i] * rfm.ReUdD[i][j]
            betaU_dupD[i][j] = vetU_dupD[i][j] * rfm.ReU[i] + vetU[i] * rfm.ReUdD[i][j]  # Needed for \beta^i RHS
            for k in range(DIM):
                # Needed for, e.g., \bar{\Lambda}^i RHS:
                betaU_dDD[i][j][k] = vetU_dDD[i][j][k] * rfm.ReU[i] + vetU_dD[i][j] * rfm.ReUdD[i][k] + \
                                     vetU_dD[i][k] * rfm.ReUdD[i][j] + vetU[i] * rfm.ReUdDD[i][j][k]

# Step 9: Standard BSSN conformal factor phi,
#         and its partial and covariant derivatives,
#         all in terms of BSSN gridfunctions like cf
def phi_and_derivs():
    # Step 9.a: Declare as globals all expressions that may be used
    #           outside this function, declare BSSN gridfunctions
    #           if not defined already, and set DIM=3.
    global phi_dD,phi_dupD,phi_dDD,exp_m4phi,phi_dBarD,phi_dBarDD
    hDD, aDD, lambdaU, vetU, betU, trK, cf, alpha = declare_BSSN_gridfunctions_if_not_declared_already()
    DIM = 3

    # GammabarUDD is used below, defined in
    #    gammabar__inverse_and_derivs()
    gammabar__inverse_and_derivs()

    # Step 9.a.i: Define partial derivatives of \phi in terms of evolved quantity "cf":
    cf_dD = ixp.declarerank1("cf_dD")
    cf_dupD = ixp.declarerank1("cf_dupD")  # Needed for \partial_t \phi next.
    cf_dDD = ixp.declarerank2("cf_dDD", "sym01")
    phi_dD = ixp.zerorank1()
    phi_dupD = ixp.zerorank1()
    phi_dDD = ixp.zerorank2()
    exp_m4phi = sp.sympify(0)

    # Step 9.a.ii: Assuming cf=phi, define exp_m4phi, phi_dD,
    #              phi_dupD (upwind finite-difference version of phi_dD), and phi_DD
    if par.parval_from_str("BSSN.BSSN_quantities::EvolvedConformalFactor_cf") == "phi":
        for i in range(DIM):
            phi_dD[i] = cf_dD[i]
            phi_dupD[i] = cf_dupD[i]
            for j in range(DIM):
                phi_dDD[i][j] = cf_dDD[i][j]
        exp_m4phi = sp.exp(-4 * cf)

    # Step 9.a.iii: Assuming cf=W=e^{-2 phi}, define exp_m4phi, phi_dD,
    #               phi_dupD (upwind finite-difference version of phi_dD), and phi_DD
    if par.parval_from_str("BSSN.BSSN_quantities::EvolvedConformalFactor_cf") == "W":
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

    # Step 9.a.iv: Assuming cf=chi=e^{-4 phi}, define exp_m4phi, phi_dD,
    #              phi_dupD (upwind finite-difference version of phi_dD), and phi_DD
    if par.parval_from_str("BSSN.BSSN_quantities::EvolvedConformalFactor_cf") == "chi":
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

    # Step 9.a.v: Error out if unsupported EvolvedConformalFactor_cf choice is made:
    cf_choice = par.parval_from_str("BSSN.BSSN_quantities::EvolvedConformalFactor_cf")
    if not (cf_choice == "phi" or cf_choice == "W" or cf_choice == "chi"):
        print("Error: EvolvedConformalFactor_cf == " + par.parval_from_str("BSSN.BSSN_quantities::EvolvedConformalFactor_cf") + " unsupported!")
        exit(1)

    # Step 9.b: Define phi_dBarD = phi_dD (since phi is a scalar) and phi_dBarDD (covariant derivative)
    #          \bar{D}_i \bar{D}_j \phi = \phi_{;\bar{i}\bar{j}} = \bar{D}_i \phi_{,j}
    #                                   = \phi_{,ij} - \bar{\Gamma}^k_{ij} \phi_{,k}
    phi_dBarD = phi_dD
    phi_dBarDD = ixp.zerorank2()
    for i in range(DIM):
        for j in range(DIM):
            phi_dBarDD[i][j] = phi_dDD[i][j]
            for k in range(DIM):
                phi_dBarDD[i][j] += - GammabarUDD[k][i][j] * phi_dD[k]
