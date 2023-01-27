# This module provides functions for setting up the right-hand sides
#     of the evolution equations a massless Scalar Field, i.e. the
#     Klein-Gordon equations, as documented in
#     Tutorial-ScalarField_RHSs.ipynb

# Authors: Leonardo R. Werneck
#          wernecklr **at** gmail **dot* com
#          Zachariah B. Etienne

# First we import needed core NRPy+ modules
import sympy as sp                # SymPy: The Python computer algebra package upon which NRPy+ depends
import indexedexp as ixp          # NRPy+: Symbolic indexed expression (e.g., tensors, vectors, etc.) support
import reference_metric as rfm    # NRPy+: Reference metric support
import BSSN.BSSN_quantities as Bq # NRPy+: BSSN quantities

# Checking Python version for correct import syntax
import sys
if sys.version_info[0] == 3:
    import ScalarField.ScalarField_declare_gridfunctions as sfgfs
elif sys.version_info[0] == 2:
    import ScalarField_declare_gridfunctions as sfgfs

###########################################
# Part B: The scalarfield_RHSs() function #
###########################################
def ScalarField_RHSs():

    # Step B.4: Set spatial dimension (must be 3 for BSSN, as BSSN is
    #           a 3+1-dimensional decomposition of the general
    #           relativistic field equations)
    DIM = 3

    # Step B.5: Import all basic (unrescaled) BSSN scalars & tensors
    Bq.BSSN_basic_tensors()
    trK   = Bq.trK
    alpha = Bq.alpha
    betaU = Bq.betaU
    Bq.gammabar__inverse_and_derivs()
    gammabarUU = Bq.gammabarUU

    global sf_rhs, sfM_rhs

    # Step B.5.a: Declare grid functions for varphi and Pi
    sf, sfM = sfgfs.declare_scalar_field_gridfunctions_if_not_declared_already()

    # Step 2.a: Add Term 1 to sf_rhs: -alpha*Pi
    sf_rhs = - alpha * sfM

    # Step 2.b: Add Term 2 to sf_rhs: beta^{i}\partial_{i}\varphi
    sf_dupD = ixp.declarerank1("sf_dupD")
    for i in range(DIM):
        sf_rhs += betaU[i] * sf_dupD[i]

    # Step 3a: Add Term 1 to sfM_rhs: alpha * K * Pi
    sfM_rhs = alpha * trK * sfM

    # Step 3b: Add Term 2 to sfM_rhs: beta^{i}\partial_{i}Pi
    sfM_dupD = ixp.declarerank1("sfM_dupD")
    for i in range(DIM):
        sfM_rhs += betaU[i] * sfM_dupD[i]

    # Step 3c: Adding Term 3 to sfM_rhs
    # Step 3c.i: Term 3a: gammabar^{ij}\partial_{i}\alpha\partial_{j}\varphi
    alpha_dD    = ixp.declarerank1("alpha_dD")
    sf_dD       = ixp.declarerank1("sf_dD")
    sfMrhsTerm3 = sp.sympify(0)
    for i in range(DIM):
        for j in range(DIM):
            sfMrhsTerm3 += - gammabarUU[i][j] * alpha_dD[i] * sf_dD[j]

    # Step 3c.ii: Term 3b: \alpha*gammabar^{ij}\partial_{i}\partial_{j}\varphi
    sf_dDD = ixp.declarerank2("sf_dDD","sym01")
    for i in range(DIM):
        for j in range(DIM):
            sfMrhsTerm3 += - alpha * gammabarUU[i][j] * sf_dDD[i][j]

    # Step 3c.iii: Term 3c: 2*alpha*gammabar^{ij}\partial_{j}\varphi\partial_{i}\phi
    Bq.phi_and_derivs() # sets exp^{-4phi} = exp_m4phi and \partial_{i}phi = phi_dD[i]
    for i in range(DIM):
        for j in range(DIM):
            sfMrhsTerm3 += - 2 * alpha * gammabarUU[i][j] * sf_dD[j] * Bq.phi_dD[i]

    # Step 3c.iv: Multiplying Term 3 by e^{-4phi} and adding it to sfM_rhs
    sfMrhsTerm3 *= Bq.exp_m4phi
    sfM_rhs     += sfMrhsTerm3

    # Step 3d: Adding Term 4 to sfM_rhs
    # Step 3d.i: Term 4a: \alpha \bar\Lambda^{i}\partial_{i}\varphi
    LambdabarU  = Bq.LambdabarU
    sfMrhsTerm4 = sp.sympify(0)
    for i in range(DIM):
        sfMrhsTerm4 += alpha * LambdabarU[i] * sf_dD[i]

    # Step 3d.ii: Evaluating \bar\gamma^{ij}\hat\Gamma^{k}_{ij}
    GammahatUDD = rfm.GammahatUDD
    gammabarGammahatContractionU = ixp.zerorank1()
    for k in range(DIM):
        for i in range(DIM):
            for j in range(DIM):
                gammabarGammahatContractionU[k] += gammabarUU[i][j] * GammahatUDD[k][i][j]

    # Step 3d.iii: Term 4b: \alpha \bar\gamma^{ij}\hat\Gamma^{k}_{ij}\partial_{k}\varphi
    for i in range(DIM):
        sfMrhsTerm4 += alpha * gammabarGammahatContractionU[i] * sf_dD[i]

    # Step 3d.iii: Multplying Term 4 by e^{-4phi} and adding it to sfM_rhs
    sfMrhsTerm4 *= Bq.exp_m4phi
    sfM_rhs     += sfMrhsTerm4

    return
