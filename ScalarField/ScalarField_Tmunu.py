# This module provides functions for setting up the energy-momentum
#     tensor of a massless Scalar Field as documented in
#     Tutorial-ScalarField_Tmunu.ipynb

# Authors: Leonardo R. Werneck
#          wernecklr **at** gmail **dot* com
#          Zachariah B. Etienne

# First we import needed core NRPy+ modules
from indexedexp import zerorank1, zerorank2      # NRPy+: Symbolic indexed expression support
from indexedexp import symm_matrix_inverter3x3   #        (e.g., tensors, vectors, etc.) 
from reference_metric import reference_metric    # NRPy+: Reference metric support
import BSSN.BSSN_quantities as Bq                # NRPy+: BSSN quantities
import BSSN.ADM_in_terms_of_BSSN as BtoA         # NRPy+: ADM quantities in terms of BSSN quantities
import BSSN.ADMBSSN_tofrom_4metric as ADMg       # NRPy+: ADM 4-metric to/from ADM or BSSN quantities
import ScalarField.ScalarField_quantities as SFq # NRPy+: Scalar field quantities

def ScalarField_T4UU(alpha=None, betaU=None, gammaDD=None):

    # Step 1: Set needed ADM quantities
    if alpha is None or betaU is None or gammaDD is None:
        # Step 1.a: Set ADM metric based on BSSN quantities
        # Step 1.a.1: Set reference metric
        reference_metric()

        # Step 1.a.2: Import all basic (unrescaled) BSSN scalars & tensors
        Bq.BSSN_basic_tensors()
        alpha = Bq.alpha
        betaU = Bq.betaU

        # Step 1.a.3: Define ADM quantities in terms of BSSN quantities
        BtoA.ADM_in_terms_of_BSSN()
        gammaDD = BtoA.gammaDD
        gammaUU = BtoA.gammaUU
    else:
        # Step 1.b: Compute gammaUU from input
        gammaUU, _ = symm_matrix_inverter3x3(gammaDD)
        
    # Step 2: Define scalar field quantitites
    SFq.ScalarField_quantities()
    sf_dD = SFq.sf_dD
    Pi    = SFq.sfM

    # Step 3: Set up \partial^{\mu}\varphi
    # Step 3.a: Set up \partial^{t}\varphi = Pi/alpha
    sf4dU = zerorank1(DIM=4)
    sf4dU[0] = Pi/alpha

    # Step 3.b: Set up \partial^{i}\varphi = -Pi*beta^{i}/alpha + gamma^{ij}\partial_{j}\varphi
    for i in range(3):
        sf4dU[i+1] = -Pi*betaU[i]/alpha
        for j in range(3):
            sf4dU[i+1] += gammaUU[i][j]*sf_dD[j]

    # Step 4: Set up \partial^{i}\varphi\partial_{i}\varphi = -Pi**2 + gamma^{ij}\partial_{i}\varphi\partial_{j}\varphi
    sf4d2 = -Pi**2
    for i in range(3):
        for j in range(3):
            sf4d2 += gammaUU[i][j]*sf_dD[i]*sf_dD[j]

    # Step 5: Set up g^{\mu\nu}
    ADMg.g4UU_ito_BSSN_or_ADM("ADM", alpha=alpha, betaU=betaU, gammaDD=gammaDD, gammaUU=gammaUU)
    g4UU = ADMg.g4UU

    # Step 6: Set up T^{\mu\nu} for a massless scalar field
    T4UU = zerorank2(DIM=4)
    for mu in range(4):
        for nu in range(4):
            T4UU[mu][nu] = sf4dU[mu]*sf4dU[nu] - g4UU[mu][nu]*sf4d2/2

    return T4UU
