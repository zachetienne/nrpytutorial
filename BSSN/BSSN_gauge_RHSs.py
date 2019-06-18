# As documented in the NRPy+ tutorial module
#   Tutorial-BSSN_time_evolution-BSSN_gauge_RHSs.ipynb
#   this module will construct the right-hand sides (RHSs)
#   expressions for the time evolution equations of the
#   BSSN gauge quantities alpha and beta^i (i.e., the
#   lapse and shift)
#
# Non-gauge BSSN time-evolution equations are documented in the
#   NRPy+ tutorial module
#   Tutorial-BSSN_time_evolution-BSSN_RHSs.ipynb,
#   and separately constructed in the BSSN.BSSN_RHSs
#   Python module.

# Author: Zachariah B. Etienne
#         zachetie **at** gmail **dot* com

# Step 1: Import all needed modules from NRPy+:
import sympy as sp
import NRPy_param_funcs as par
import indexedexp as ixp
import grid as gri
import finite_difference as fin
import reference_metric as rfm

# Step 1.a: Load BSSN quantities module
import BSSN.BSSN_quantities as Bq
# Step 1.b: Load BSSN RHSs module
import BSSN.BSSN_RHSs as Brhs

# Step 1.c: Declare/initialize parameters for this module
thismodule = __name__
par.initialize_param(par.glb_param("char", thismodule, "LapseEvolutionOption", "OnePlusLog"))
par.initialize_param(par.glb_param("char", thismodule, "ShiftEvolutionOption", "GammaDriving2ndOrder_Covariant"))

def BSSN_gauge_RHSs():
    # Step 1.d: Set spatial dimension (must be 3 for BSSN, as BSSN is
    #           a 3+1-dimensional decomposition of the general
    #           relativistic field equations)
    DIM = 3

    # Step 1.e: Given the chosen coordinate system, set up
    #           corresponding reference metric and needed
    #           reference metric quantities
    # The following function call sets up the reference metric
    #    and related quantities, including rescaling matrices ReDD,
    #    ReU, and hatted quantities.
    rfm.reference_metric()

    # Step 1.f: Define needed BSSN quantities:
    # Declare scalars & tensors (in terms of rescaled BSSN quantities)
    Bq.BSSN_basic_tensors()
    Bq.betaU_derivs()
    # Declare BSSN_RHSs (excluding the time evolution equations for the gauge conditions)
    Brhs.BSSN_RHSs()

    # Step 2.a: The 1+log lapse condition:
    #   \partial_t \alpha = \beta^i \alpha_{,i} - 2*\alpha*K
    # First import expressions from BSSN_quantities
    cf = Bq.cf
    trK = Bq.trK
    alpha = Bq.alpha
    betaU = Bq.betaU

    # Implement the 1+log lapse condition
    global alpha_rhs
    alpha_rhs = sp.sympify(0)
    if par.parval_from_str(thismodule + "::LapseEvolutionOption") == "OnePlusLog":
        alpha_rhs = -2 * alpha * trK
        alpha_dupD = ixp.declarerank1("alpha_dupD")
        for i in range(DIM):
            alpha_rhs += betaU[i] * alpha_dupD[i]

    # Step 2.b: Implement the harmonic slicing lapse condition
    elif par.parval_from_str(thismodule + "::LapseEvolutionOption") == "HarmonicSlicing":
        if par.parval_from_str("BSSN.BSSN_quantities::EvolvedConformalFactor_cf") == "W":
            alpha_rhs = -3 * cf ** (-4) * Brhs.cf_rhs
        elif par.parval_from_str("BSSN.BSSN_quantities::EvolvedConformalFactor_cf") == "phi":
            alpha_rhs = 6 * sp.exp(6 * cf) * Brhs.cf_rhs
        else:
            print("Error LapseEvolutionOption==HarmonicSlicing unsupported for EvolvedConformalFactor_cf!=(W or phi)")
            exit(1)

    # Step 2.c: Frozen lapse
    #    \partial_t \alpha = 0
    elif par.parval_from_str(thismodule + "::LapseEvolutionOption") == "Frozen":
        alpha_rhs = sp.sympify(0)

    else:
        print("Error: "+thismodule + "::LapseEvolutionOption == "+
              par.parval_from_str(thismodule + "::LapseEvolutionOption")+" not supported!")
        exit(1)

    # Step 3.a: Set \partial_t \beta^i
    # First import expressions from BSSN_quantities
    BU = Bq.BU
    betU = Bq.betU
    betaU_dupD = Bq.betaU_dupD
    # Define needed quantities
    beta_rhsU = ixp.zerorank1()
    B_rhsU = ixp.zerorank1()
    if par.parval_from_str(thismodule + "::ShiftEvolutionOption") == "GammaDriving2ndOrder_NoCovariant":
        # Step 3.a.i: Compute right-hand side of beta^i
        # *  \partial_t \beta^i = \beta^j \beta^i_{,j} + B^i
        for i in range(DIM):
            beta_rhsU[i] += BU[i]
            for j in range(DIM):
                beta_rhsU[i] += betaU[j] * betaU_dupD[i][j]
        # Compute right-hand side of B^i:
        eta = par.Cparameters("REAL", thismodule, ["eta"])

        # Step 3.a.ii: Compute right-hand side of B^i
        # *  \partial_t B^i     = \beta^j B^i_{,j} + 3/4 * \partial_0 \Lambda^i - eta B^i
        # Step 3.a.iii: Define BU_dupD, in terms of derivative of rescaled variable \bet^i
        BU_dupD = ixp.zerorank2()
        betU_dupD = ixp.declarerank2("betU_dupD", "nosym")
        for i in range(DIM):
            for j in range(DIM):
                BU_dupD[i][j] = betU_dupD[i][j] * rfm.ReU[i] + betU[i] * rfm.ReUdD[i][j]

        # Step 3.a.iv: Compute \partial_0 \bar{\Lambda}^i = (\partial_t - \beta^i \partial_i) \bar{\Lambda}^j
        Lambdabar_partial0 = ixp.zerorank1()
        for i in range(DIM):
            Lambdabar_partial0[i] = Brhs.Lambdabar_rhsU[i]
        for i in range(DIM):
            for j in range(DIM):
                Lambdabar_partial0[j] += -betaU[i] * Brhs.LambdabarU_dupD[j][i]

        # Step 3.a.v: Evaluate RHS of B^i:
        for i in range(DIM):
            B_rhsU[i] += sp.Rational(3, 4) * Lambdabar_partial0[i] - eta * BU[i]
            for j in range(DIM):
                B_rhsU[i] += betaU[j] * BU_dupD[i][j]
                
    # Step 3.b: The right-hand side of the \partial_t \beta^i equation
    if par.parval_from_str(thismodule + "::ShiftEvolutionOption") == "GammaDriving2ndOrder_Covariant":
        # Step 3.b Option 2: \partial_t \beta^i = \left[\beta^j \bar{D}_j \beta^i\right] + B^{i}
        # First we need GammabarUDD, defined in Bq.gammabar__inverse_and_derivs()
        Bq.gammabar__inverse_and_derivs()
        GammabarUDD = Bq.GammabarUDD
        # Then compute right-hand side:
        # Term 1: \beta^j \beta^i_{,j}
        for i in range(DIM):
            for j in range(DIM):
                beta_rhsU[i] += betaU[j] * betaU_dupD[i][j]

        # Term 2: \beta^j \bar{\Gamma}^i_{mj} \beta^m
        for i in range(DIM):
            for j in range(DIM):
                for m in range(DIM):
                    beta_rhsU[i] += betaU[j] * GammabarUDD[i][m][j] * betaU[m]
        # Term 3: B^i
        for i in range(DIM):
            beta_rhsU[i] += BU[i]

    if par.parval_from_str(thismodule + "::ShiftEvolutionOption") == "GammaDriving2ndOrder_Covariant":
        # Step 3.c: Covariant option:
        #  \partial_t B^i = \beta^j \bar{D}_j B^i
        #               + \frac{3}{4} ( \partial_t \bar{\Lambda}^{i} - \beta^j \bar{D}_j \bar{\Lambda}^{i} )
        #               - \eta B^{i}
        #                 = \beta^j B^i_{,j} + \beta^j \bar{\Gamma}^i_{mj} B^m
        #               + \frac{3}{4}[ \partial_t \bar{\Lambda}^{i}
        #                            - \beta^j (\bar{\Lambda}^i_{,j} + \bar{\Gamma}^i_{mj} \bar{\Lambda}^m)]
        #               - \eta B^{i}
        # Term 1, part a: First compute B^i_{,j} using upwinded derivative
        BU_dupD = ixp.zerorank2()
        betU_dupD = ixp.declarerank2("betU_dupD", "nosym")
        for i in range(DIM):
            for j in range(DIM):
                BU_dupD[i][j] = betU_dupD[i][j] * rfm.ReU[i] + betU[i] * rfm.ReUdD[i][j]
        # Term 1: \beta^j B^i_{,j}
        for i in range(DIM):
            for j in range(DIM):
                B_rhsU[i] += betaU[j] * BU_dupD[i][j]
        # Term 2: \beta^j \bar{\Gamma}^i_{mj} B^m
        for i in range(DIM):
            for j in range(DIM):
                for m in range(DIM):
                    B_rhsU[i] += betaU[j] * GammabarUDD[i][m][j] * BU[m]
        # Term 3: \frac{3}{4}\partial_t \bar{\Lambda}^{i}
        for i in range(DIM):
            B_rhsU[i] += sp.Rational(3, 4) * Brhs.Lambdabar_rhsU[i]
        # Term 4: -\frac{3}{4}\beta^j \bar{\Lambda}^i_{,j}
        for i in range(DIM):
            for j in range(DIM):
                B_rhsU[i] += -sp.Rational(3, 4) * betaU[j] * Brhs.LambdabarU_dupD[i][j]
        # Term 5: -\frac{3}{4}\beta^j \bar{\Gamma}^i_{mj} \bar{\Lambda}^m
        for i in range(DIM):
            for j in range(DIM):
                for m in range(DIM):
                    B_rhsU[i] += -sp.Rational(3, 4) * betaU[j] * GammabarUDD[i][m][j] * Bq.LambdabarU[m]
        # Term 6: - \eta B^i
        # eta is a free parameter; we declare it here:
        eta = par.Cparameters("REAL", thismodule, ["eta"])
        for i in range(DIM):
            B_rhsU[i] += -eta * BU[i]

    # Step 4: Rescale the BSSN gauge RHS quantities so that the evolved
    #         variables may remain smooth across coord singularities
    global vet_rhsU,bet_rhsU
    vet_rhsU = ixp.zerorank1()
    bet_rhsU = ixp.zerorank1()
    for i in range(DIM):
        vet_rhsU[i] = beta_rhsU[i] / rfm.ReU[i]
        bet_rhsU[i] = B_rhsU[i] / rfm.ReU[i]
    # print(str(Abar_rhsDD[2][2]).replace("**","^").replace("_","").replace("xx","x").replace("sin(x2)","Sin[x2]").replace("sin(2*x2)","Sin[2*x2]").replace("cos(x2)","Cos[x2]").replace("detgbaroverdetghat","detg"))
    # print(str(Dbarbetacontraction).replace("**","^").replace("_","").replace("xx","x").replace("sin(x2)","Sin[x2]").replace("detgbaroverdetghat","detg"))
    # print(betaU_dD)
    # print(str(trK_rhs).replace("xx2","xx3").replace("xx1","xx2").replace("xx0","xx1").replace("**","^").replace("_","").replace("sin(xx2)","Sinx2").replace("xx","x").replace("sin(2*x2)","Sin2x2").replace("cos(x2)","Cosx2").replace("detgbaroverdetghat","detg"))
    # print(str(bet_rhsU[0]).replace("xx2","xx3").replace("xx1","xx2").replace("xx0","xx1").replace("**","^").replace("_","").replace("sin(xx2)","Sinx2").replace("xx","x").replace("sin(2*x2)","Sin2x2").replace("cos(x2)","Cosx2").replace("detgbaroverdetghat","detg"))
