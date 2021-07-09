# This module provides functions for setting up the energy-momentum
#     tensor of a massless Scalar Field as documented in
#     Tutorial-ScalarField_Tmunu.ipynb

# Authors: Leonardo R. Werneck
#          wernecklr **at** gmail **dot* com
#          Zachariah B. Etienne

# First we import needed core NRPy+ modules
import sympy as sp                         # SymPy: The Python computer algebra package upon which NRPy+ depends
import indexedexp as ixp                   # NRPy+: Symbolic indexed expression (e.g., tensors, vectors, etc.) support
import reference_metric as rfm             # NRPy+: Reference metric support
import BSSN.BSSN_quantities as Bq          # NRPy+: BSSN quantities
import BSSN.ADM_in_terms_of_BSSN as BtoA   # NRPy+: ADM quantities in terms of BSSN quantities
import BSSN.ADMBSSN_tofrom_4metric as ADMg # NRPy+: ADM 4-metric to/from ADM or BSSN quantities

def ScalarField_Tmunu():

    global T4UU

    # Step 1.c: Set spatial dimension (must be 3 for BSSN, as BSSN is
    #           a 3+1-dimensional decomposition of the general
    #           relativistic field equations)
    DIM = 3

    # Step 1.d: Given the chosen coordinate system, set up
    #           corresponding reference metric and needed
    #           reference metric quantities
    #    The following function call sets up the reference metric
    #    and related quantities, including rescaling matrices ReDD,
    #    ReU, and hatted quantities.
    rfm.reference_metric()

    # Step 1.e: Import all basic (unrescaled) BSSN scalars & tensors
    Bq.BSSN_basic_tensors()
    alpha = Bq.alpha
    betaU = Bq.betaU

    # Step 1.g: Define ADM quantities in terms of BSSN quantities
    BtoA.ADM_in_terms_of_BSSN()
    gammaDD = BtoA.gammaDD
    gammaUU = BtoA.gammaUU

    # Step 1.h: Define scalar field quantitites
    sf_dD   = ixp.declarerank1("sf_dD")
    Pi      = sp.Symbol("sfM",real=True)

    # Step 2a: Set up \partial^{t}\varphi = Pi/alpha
    sf4dU = ixp.zerorank1(DIM=4)
    sf4dU[0] = Pi/alpha

    # Step 2b: Set up \partial^{i}\varphi = -Pi*beta^{i}/alpha + gamma^{ij}\partial_{j}\varphi
    for i in range(DIM):
        sf4dU[i+1] = -Pi*betaU[i]/alpha
        for j in range(DIM):
            sf4dU[i+1] += gammaUU[i][j]*sf_dD[j]

    # Step 2c: Set up \partial^{i}\varphi\partial_{i}\varphi = -Pi**2 + gamma^{ij}\partial_{i}\varphi\partial_{j}\varphi
    sf4d2 = -Pi**2
    for i in range(DIM):
        for j in range(DIM):
            sf4d2 += gammaUU[i][j]*sf_dD[i]*sf_dD[j]

    # Step 3a: Setting up g^{\mu\nu}
    ADMg.g4UU_ito_BSSN_or_ADM("ADM",gammaDD=gammaDD,betaU=betaU,alpha=alpha, gammaUU=gammaUU)
    g4UU = ADMg.g4UU

    # Step 3b: Setting up T^{\mu\nu} for a massless scalar field
    T4UU = ixp.zerorank2(DIM=4)
    for mu in range(4):
        for nu in range(4):
            T4UU[mu][nu] = sf4dU[mu]*sf4dU[nu] - g4UU[mu][nu]*sf4d2/2
