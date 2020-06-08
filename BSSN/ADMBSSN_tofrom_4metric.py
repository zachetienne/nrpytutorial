# As documented in the NRPy+ tutorial module
#   Tutorial-ADMBSSN_tofrom_4metric.ipynb,
#   this module will construct expressions for
#   ADM or BSSN quantities in terms of the
#   4-metric g4DD, and g4DD/g4UU in terms of
#   ADM/BSSN quantities.

# Author: Zachariah B. Etienne
#         zachetie **at** gmail **dot* com

import sympy as sp                # SymPy: The Python computer algebra package upon which NRPy+ depends
import indexedexp as ixp          # NRPy+: Symbolic indexed expression (e.g., tensors, vectors, etc.) support
import sys                        # Standard Python modules for multiplatform OS-level functions

def setup_ADM_quantities(inputvars):
    if inputvars == "ADM":
        gammaDD = ixp.declarerank2("gammaDD", "sym01")
        betaU = ixp.declarerank1("betaU")
        alpha = sp.symbols("alpha", real=True)
    elif inputvars == "BSSN":
        import BSSN.ADM_in_terms_of_BSSN as AitoB

        # Construct gamma_{ij} in terms of cf & gammabar_{ij}
        AitoB.ADM_in_terms_of_BSSN()
        gammaDD = AitoB.gammaDD
        # Next construct beta^i in terms of vet^i and reference metric quantities
        import BSSN.BSSN_quantities as Bq

        Bq.BSSN_basic_tensors()
        betaU = Bq.betaU
        alpha = sp.symbols("alpha", real=True)
    else:
        print("inputvars = " + str(inputvars) + " not supported. Please choose ADM or BSSN.")
        sys.exit(1)
    return gammaDD,betaU,alpha

# g_{mu nu} in terms of BSSN (if inputvars=="BSSN") or ADM (if inputvars=="ADM") variables.
def g4DD_ito_BSSN_or_ADM(inputvars,gammaDD=None,betaU=None,alpha=None):
    # Step 0: Declare g4DD as globals, to make interfacing with other modules/functions easier
    global g4DD

    # Step 1: Set gammaDD, betaU, and alpha if not already input.
    if gammaDD is None and betaU is None and alpha is None: # Use "is None" instead of "==None", as the former is more correct.
        gammaDD,betaU,alpha = setup_ADM_quantities(inputvars)

    # Step 2: Compute g4DD = g_{mu nu}:
    # To get \gamma_{\mu \nu} = gamma4DD[mu][nu], we'll need to construct the 4-metric, using Eq. 2.122 in B&S:
    g4DD = ixp.zerorank2(DIM=4)

    # Step 2.a: Compute beta_i via Eq. 2.121 in B&S
    betaD = ixp.zerorank1()
    for i in range(3):
        for j in range(3):
            betaD[i] += gammaDD[i][j] * betaU[j]

    # Step 2.b: Compute beta_i beta^i, the beta contraction.
    beta2 = sp.sympify(0)
    for i in range(3):
        beta2 += betaU[i] * betaD[i]

    # Step 2.c: Construct g4DD via Eq. 2.122 in B&S
    g4DD[0][0] = -alpha ** 2 + beta2
    for mu in range(1, 4):
        g4DD[mu][0] = g4DD[0][mu] = betaD[mu - 1]
    for mu in range(1, 4):
        for nu in range(1, 4):
            g4DD[mu][nu] = gammaDD[mu - 1][nu - 1]

# g^{mu nu} in terms of BSSN (if inputvars=="BSSN") or ADM (if inputvars=="ADM") variables.
def g4UU_ito_BSSN_or_ADM(inputvars,gammaDD=None,betaU=None,alpha=None, gammaUU=None):
    # Step 0: Declare g4UU as globals, to make interfacing with other modules/functions easier
    global g4UU

    # Step 1: Set gammaDD, betaU, and alpha if not already input.
    if gammaDD is None and betaU is None and alpha is None:  # Use "is None" instead of "==None", as the former is more correct.
        gammaDD,betaU,alpha = setup_ADM_quantities(inputvars)

    # Step 2: Compute g4UU = g_{mu nu}:
    # To get \gamma^{\mu \nu} = gamma4UU[mu][nu], we'll need to use Eq. 2.119 in B&S.
    g4UU = ixp.zerorank2(DIM=4)

    # Step 3: Construct g4UU = g^{mu nu}
    # Step 3.a: Compute gammaUU based on provided gammaDD:
    if gammaUU is None: # Use "is None" instead of "==None", as the former is more correct.
        gammaUU, _gammaDET = ixp.symm_matrix_inverter3x3(gammaDD) # _gammaDET is unused.

    # Then evaluate g4UU:
    g4UU = ixp.zerorank2(DIM=4)

    g4UU[0][0] = -1 / alpha ** 2
    for mu in range(1, 4):
        g4UU[0][mu] = g4UU[mu][0] = betaU[mu - 1] / alpha ** 2
    for mu in range(1, 4):
        for nu in range(1, 4):
            g4UU[mu][nu] = gammaUU[mu - 1][nu - 1] - betaU[mu - 1] * betaU[nu - 1] / alpha ** 2

# BSSN (if inputvars=="BSSN") or ADM (if inputvars=="ADM") metric variables in terms of g_{mu nu}
def BSSN_or_ADM_ito_g4DD(inputvars,g4DD=None):
    # Step 0: Declare output variables as globals, to make interfacing with other modules/functions easier
    if inputvars == "ADM":
        global gammaDD, betaU, alpha
    elif inputvars == "BSSN":
        global hDD, cf, vetU, alpha
    else:
        print("inputvars = " + str(inputvars) + " not supported. Please choose ADM or BSSN.")
        sys.exit(1)

    # Step 1: declare g4DD as symmetric rank-4 tensor:
    g4DD_is_input_into_this_function = True
    if g4DD is None: # Use "is None" instead of "==None", as the former is more correct.
        g4DD = ixp.declarerank2("g4DD", "sym01", DIM=4)
        g4DD_is_input_into_this_function = False

    # Step 2: Compute gammaDD & betaD
    betaD = ixp.zerorank1()
    gammaDD = ixp.zerorank2()
    for i in range(3):
        betaD[i] = g4DD[0][i]
        for j in range(3):
            gammaDD[i][j] = g4DD[i + 1][j + 1]

    # Step 3: Compute betaU
    # Step 3.a: Compute gammaUU based on provided gammaDD
    gammaUU, _gammaDET = ixp.symm_matrix_inverter3x3(gammaDD) # _gammaDET is unused.

    # Step 3.b: Use gammaUU to raise betaU
    betaU = ixp.zerorank1()
    for i in range(3):
        for j in range(3):
            betaU[i] += gammaUU[i][j] * betaD[j]

    # Step 4:   Compute alpha = sqrt(beta^2 - g_{00}):
    # Step 4.a: Compute beta^2 = beta^k beta_k:
    beta_squared = sp.sympify(0)
    for k in range(3):
        beta_squared += betaU[k] * betaD[k]

    # Step 4.b: alpha = sqrt(beta^2 - g_{00}):
    if g4DD_is_input_into_this_function == False:
        alpha = sp.sqrt(sp.simplify(beta_squared) - g4DD[0][0])
    else:
        alpha = sp.sqrt(beta_squared - g4DD[0][0])

    # Step 5: If inputvars == "ADM", we are finished. Return.
    if inputvars == "ADM":
        return

    # Step 6: If inputvars == "BSSN", convert ADM to BSSN
    import BSSN.BSSN_in_terms_of_ADM as BitoA
    dummyBU = ixp.zerorank1()
    BitoA.gammabarDD_hDD(gammaDD)
    BitoA.cf_from_gammaDD(gammaDD)
    BitoA.betU_vetU(betaU, dummyBU)
    hDD  = BitoA.hDD
    cf   = BitoA.cf
    vetU = BitoA.vetU
