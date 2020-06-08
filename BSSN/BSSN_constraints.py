# As documented in the NRPy+ tutorial module
#   Tutorial-BSSN_constraints.ipynb,
#   this module will construct expressions for the
#   BSSN Hamiltonian and momentum constraint equations.

# Author: Zachariah B. Etienne
#         zachetie **at** gmail **dot* com

# Step 1: Initialize needed Python/NRPy+ modules
import sympy as sp                # SymPy: The Python computer algebra package upon which NRPy+ depends
import NRPy_param_funcs as par    # NRPy+: Parameter interface
import indexedexp as ixp          # NRPy+: Symbolic indexed expression (e.g., tensors, vectors, etc.) support
import reference_metric as rfm    # NRPy+: Reference metric support
import BSSN.BSSN_quantities as Bq # NRPy+: Computes useful BSSN quantities
import BSSN.BSSN_T4UUmunu_vars as BTmunu

thismodule = __name__

def BSSN_constraints(add_T4UUmunu_source_terms=False):
    # Step 1.a: Set spatial dimension (must be 3 for BSSN, as BSSN is
    #           a 3+1-dimensional decomposition of the general
    #           relativistic field equations)
    DIM = 3

    # Step 1.b: Given the chosen coordinate system, set up
    #           corresponding reference metric and needed
    #           reference metric quantities
    # The following function call sets up the reference metric
    #    and related quantities, including rescaling matrices ReDD,
    #    ReU, and hatted quantities.
    rfm.reference_metric()


    # Step 2: Hamiltonian constraint.
    # First declare all needed variables
    Bq.declare_BSSN_gridfunctions_if_not_declared_already()  # Sets trK
    Bq.BSSN_basic_tensors()  # Sets AbarDD
    Bq.gammabar__inverse_and_derivs()  # Sets gammabarUU
    Bq.AbarUU_AbarUD_trAbar_AbarDD_dD()  # Sets AbarUU and AbarDD_dD
    Bq.RicciBar__gammabarDD_dHatD__DGammaUDD__DGammaU()  # Sets RbarDD
    Bq.phi_and_derivs()  # Sets phi_dBarD & phi_dBarDD

    ###############################
    ###############################
    #  HAMILTONIAN CONSTRAINT
    ###############################
    ###############################

    # Term 1: 2/3 K^2
    global H
    H = sp.Rational(2, 3) * Bq.trK ** 2

    # Term 2: -A_{ij} A^{ij}
    for i in range(DIM):
        for j in range(DIM):
            H += -Bq.AbarDD[i][j] * Bq.AbarUU[i][j]

    # Term 3a: trace(Rbar)
    Rbartrace = sp.sympify(0)
    for i in range(DIM):
        for j in range(DIM):
            Rbartrace += Bq.gammabarUU[i][j] * Bq.RbarDD[i][j]

    # Term 3b: -8 \bar{\gamma}^{ij} \bar{D}_i \phi \bar{D}_j \phi = -8*phi_dBar_times_phi_dBar
    # Term 3c: -8 \bar{\gamma}^{ij} \bar{D}_i \bar{D}_j \phi      = -8*phi_dBarDD_contraction
    phi_dBar_times_phi_dBar = sp.sympify(0)  # Term 3b
    phi_dBarDD_contraction = sp.sympify(0)  # Term 3c
    for i in range(DIM):
        for j in range(DIM):
            phi_dBar_times_phi_dBar += Bq.gammabarUU[i][j] * Bq.phi_dBarD[i] * Bq.phi_dBarD[j]
            phi_dBarDD_contraction += Bq.gammabarUU[i][j] * Bq.phi_dBarDD[i][j]

    # Add Term 3:
    H += Bq.exp_m4phi * (Rbartrace - 8 * (phi_dBar_times_phi_dBar + phi_dBarDD_contraction))

    if add_T4UUmunu_source_terms:
        M_PI = par.Cparameters("#define", thismodule, "M_PI","")  # M_PI is pi as defined in C
        BTmunu.define_BSSN_T4UUmunu_rescaled_source_terms()
        rho = BTmunu.rho
        H += -16*M_PI*rho

    # FIXME: ADD T4UUmunu SOURCE TERMS TO MOMENTUM CONSTRAINT!

    # Step 3: M^i, the momentum constraint

    ###############################
    ###############################
    #  MOMENTUM CONSTRAINT
    ###############################
    ###############################

    # SEE Tutorial-BSSN_constraints.ipynb for full documentation.
    global MU
    MU = ixp.zerorank1()

    # Term 2: 6 A^{ij} \partial_j \phi:
    for i in range(DIM):
        for j in range(DIM):
            MU[i] += 6*Bq.AbarUU[i][j]*Bq.phi_dD[j]

    # Term 3: -2/3 \bar{\gamma}^{ij} K_{,j}
    trK_dD = ixp.declarerank1("trK_dD") # Not defined in BSSN_RHSs; only trK_dupD is defined there.
    for i in range(DIM):
        for j in range(DIM):
            MU[i] += -sp.Rational(2,3)*Bq.gammabarUU[i][j]*trK_dD[j]

    # Next evaluate the conformal covariant derivative \bar{D}_j \bar{A}_{lm}
    AbarDD_dBarD = ixp.zerorank3()
    for i in range(DIM):
        for j in range(DIM):
            for k in range(DIM):
                AbarDD_dBarD[i][j][k] = Bq.AbarDD_dD[i][j][k]
                for l in range(DIM):
                    AbarDD_dBarD[i][j][k] += -Bq.GammabarUDD[l][k][i] * Bq.AbarDD[l][j]
                    AbarDD_dBarD[i][j][k] += -Bq.GammabarUDD[l][k][j] * Bq.AbarDD[i][l]

    # Term 1: Contract twice with the metric to make \bar{D}_{j} \bar{A}^{ij}
    for i in range(DIM):
        for j in range(DIM):
            for k in range(DIM):
                for l in range(DIM):
                    MU[i] += Bq.gammabarUU[i][k] * Bq.gammabarUU[j][l] * AbarDD_dBarD[k][l][j]

    # Finally, we multiply by e^{-4 phi} and rescale the momentum constraint:
    for i in range(DIM):
        MU[i] *= Bq.exp_m4phi / rfm.ReU[i]
