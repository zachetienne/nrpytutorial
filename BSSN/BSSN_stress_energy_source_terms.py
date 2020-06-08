# As documented in the NRPy+ tutorial module
#   Tutorial-BSSN_stress_energy_source_terms.ipynb,
#   this module will construct expressions for
#   BSSN stress-energy source terms, in terms of
#   elements of T^{mu nu}.

# Author: Zachariah B. Etienne
#         zachetie **at** gmail **dot* com

import sympy as sp                         # SymPy: The Python computer algebra package upon which NRPy+ depends
import NRPy_param_funcs as par             # NRPy+: Parameter interface
import indexedexp as ixp                   # NRPy+: Symbolic indexed expression (e.g., tensors, vectors, etc.) support
import reference_metric as rfm             # NRPy+: Reference metric support
import BSSN.ADMBSSN_tofrom_4metric as AB4m # NRPy+: ADM/BSSN <-> 4-metric conversions
import sys                                 # Standard Python modules for multiplatform OS-level functions

thismodule = __name__

# Define BSSN source terms in terms of T^{mu nu} or T_{mu nu}
def stress_energy_source_terms_ito_T4UU_and_ADM_or_BSSN_metricvars(inputvars,custom_T4UU=None):
    # Step 1: Check if rfm.reference_metric() already called. If not, BSSN
    #         quantities are not yet defined, so cannot proceed!
    if rfm.have_already_called_reference_metric_function == False:
        print("BSSN_source_terms_ito_T4UU(): Must call reference_metric() first!")
        sys.exit(1)

    # Step 2.a: Define gamma4DD[mu][nu] = g_{mu nu} + n_{mu} n_{nu}
    alpha = sp.symbols("alpha", real=True)
    zero = sp.sympify(0)
    n4D = [sp.sympify(-1)*alpha, zero, zero, zero]
    AB4m.g4DD_ito_BSSN_or_ADM(inputvars)

    gamma4DD = ixp.zerorank2(DIM=4)
    for mu in range(4):
        for nu in range(4):
            gamma4DD[mu][nu] = AB4m.g4DD[mu][nu] + n4D[mu] * n4D[nu]

    # Step 2.b: If expression for components of T4UU not given, declare T4UU here
    if custom_T4UU is None: # Use "is None" instead of "==None", as the former is more correct.
        T4UU = ixp.declarerank2("T4UU","sym01",DIM=4)
    else:
        T4UU = custom_T4UU

    # Step 2.c: Define BSSN source terms
    global SDD,SD,S,rho
    # Step 2.c.i: S_{ij} = gamma_{i mu} gamma_{j nu} T^{mu nu}
    SDD = ixp.zerorank2()
    for i in range(3):
        for j in range(3):
            for mu in range(4):
                for nu in range(4):
                    SDD[i][j] += gamma4DD[i+1][mu] * gamma4DD[j+1][nu] * T4UU[mu][nu]
    # Step 2.c.ii: S_{i} = -gamma_{i mu} n_{nu} T^{mu nu}
    SD = ixp.zerorank1()
    for i in range(3):
        for mu in range(4):
            for nu in range(4):
                SD[i] += - gamma4DD[i+1][mu] * n4D[nu] * T4UU[mu][nu]
    # Step 2.c.iii: S = gamma^{ij} S_{ij}
    if inputvars == "ADM":
        gammaDD = ixp.declarerank2("gammaDD", "sym01")
        gammaUU, dummydet = ixp.symm_matrix_inverter3x3(gammaDD)  # Set gammaUU
    elif inputvars == "BSSN":
        import BSSN.ADM_in_terms_of_BSSN as AitoB  # NRPy+: ADM quantities in terms of BSSN quantities
        AitoB.ADM_in_terms_of_BSSN()
        gammaUU = AitoB.gammaUU

    S = zero
    for i in range(3):
        for j in range(3):
            S += gammaUU[i][j] * SDD[i][j]
    # Step 2.c.iv: rho = n_{mu} n_{nu} T^{mu nu}
    rho = zero
    for mu in range(4):
        for nu in range(4):
            rho += n4D[mu] * n4D[nu] * T4UU[mu][nu]
    return SDD,SD,S,rho

# Step 3: Add BSSN stress-energy source terms to BSSN RHSs
def BSSN_source_terms_for_BSSN_RHSs(custom_T4UU=None):
    global sourceterm_trK_rhs, sourceterm_a_rhsDD, sourceterm_lambda_rhsU, sourceterm_Lambdabar_rhsU

    # Step 3.a: Call BSSN_source_terms_ito_T4UU to get SDD, SD, S, & rho

    if custom_T4UU == "unrescaled BSSN source terms already given":
        SDD = ixp.declarerank2("SDD", "sym01")
        SD = ixp.declarerank1("SD")
        S = sp.symbols("S", real=True)
        rho = sp.symbols("rho", real=True)
    else:
        SDD,SD,S,rho = stress_energy_source_terms_ito_T4UU_and_ADM_or_BSSN_metricvars("BSSN", custom_T4UU)
    PI = par.Cparameters("REAL", thismodule, ["PI"], "3.14159265358979323846264338327950288")
    alpha = sp.symbols("alpha", real=True)

    # Step 3.b: trK_rhs
    sourceterm_trK_rhs = 4 * PI * alpha * (rho + S)

    # Step 3.c: Abar_rhsDD:
    # Step 3.c.i: Compute trace-free part of S_{ij}:
    import BSSN.BSSN_quantities as Bq
    Bq.BSSN_basic_tensors()  # Sets gammabarDD
    gammabarUU, dummydet = ixp.symm_matrix_inverter3x3(Bq.gammabarDD)  # Set gammabarUU
    tracefree_SDD = ixp.zerorank2()
    for i in range(3):
        for j in range(3):
            tracefree_SDD[i][j] = SDD[i][j]
    for i in range(3):
        for j in range(3):
            for k in range(3):
                for m in range(3):
                    tracefree_SDD[i][j] += -sp.Rational(1, 3) * Bq.gammabarDD[i][j] * gammabarUU[k][m] * SDD[k][m]
    # Step 3.c.ii: Define exp_m4phi = e^{-4 phi}
    Bq.phi_and_derivs()
    # Step 3.c.iii: Evaluate stress-energy part of AbarDD's RHS
    sourceterm_a_rhsDD = ixp.zerorank2()
    for i in range(3):
        for j in range(3):
            Abar_rhsDDij = -8 * PI * alpha * Bq.exp_m4phi * tracefree_SDD[i][j]
            sourceterm_a_rhsDD[i][j] = Abar_rhsDDij / rfm.ReDD[i][j]

    # Step 3.d: Stress-energy part of Lambdabar_rhsU = stressenergy_Lambdabar_rhsU
    sourceterm_Lambdabar_rhsU = ixp.zerorank1()
    for i in range(3):
        for j in range(3):
            sourceterm_Lambdabar_rhsU[i] += -16 * PI * alpha * gammabarUU[i][j] * SD[j]
    sourceterm_lambda_rhsU = ixp.zerorank1()
    for i in range(3):
        sourceterm_lambda_rhsU[i] = sourceterm_Lambdabar_rhsU[i] / rfm.ReU[i]


# Step 4: Add BSSN stress-energy source terms to BSSN constraints
def BSSN_source_terms_for_BSSN_constraints(custom_T4UU=None):
    global sourceterm_H, sourceterm_MU

    # Step 4.a: Call BSSN_source_terms_ito_T4UU to get SDD, SD, S, & rho
    if custom_T4UU == "unrescaled BSSN source terms already given":
        # SDD and S unused, so we ignore their return values from ixp.declarerankN() below
        ixp.declarerank2("SDD", "sym01")
        SD = ixp.declarerank1("SD")
        sp.symbols("S", real=True)
        rho = sp.symbols("rho", real=True)
    else:
        _SDD,SD,_S,rho = stress_energy_source_terms_ito_T4UU_and_ADM_or_BSSN_metricvars("BSSN", custom_T4UU) #_SDD,_S unused.
    PI = par.Cparameters("REAL", thismodule, ["PI"], "3.14159265358979323846264338327950288")

    # Step 4.b: Add source term to the Hamiltonian constraint H
    sourceterm_H = -16 * PI * rho

    # Step 4.c: Add source term to the momentum constraint M^i
    # Step 4.c.i: Compute gammaUU in terms of BSSN quantities
    import BSSN.ADM_in_terms_of_BSSN as AitoB
    AitoB.ADM_in_terms_of_BSSN()  # Provides gammaUU
    # Step 4.c.ii: Raise S_i
    SU = ixp.zerorank1()
    for i in range(3):
        for j in range(3):
            SU[i] += AitoB.gammaUU[i][j] * SD[j]
    # Step 4.c.iii: Add source term to momentum constraint & rescale:
    sourceterm_MU = ixp.zerorank1()
    for i in range(3):
        sourceterm_MU[i] = -8 * PI * SU[i] / rfm.ReU[i]
