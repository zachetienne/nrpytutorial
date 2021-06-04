# As documented in the NRPy+ tutorial module
#   Tutorial-GRFFE_Equations-Cartesian-new.ipynb
#   this module will construct useful quantities
#   for the GiRaFFE implementation of
#   general relativistic force-free
#   electrodynamics (GRFFE)
#
# Authors: Zachariah B. Etienne
#          zachetie **at** gmail **dot* com
#          Patrick Nelson

# Step 1: Import needed core NRPy+ modules
import sympy as sp               # SymPy: The Python computer algebra package upon which NRPy+ depends
import indexedexp as ixp         # NRPy+: Symbolic indexed expression (e.g., tensors, vectors, etc.) support

# Step 2: Define needed quantities for T_{EM}^{mu nu}, the EM part of the stress-energy tensor
# Step 2.a: Define B^i = Btilde^i / sqrt(gamma)
def compute_B_notildeU(sqrtgammaDET, B_tildeU):
    global B_notildeU
    B_notildeU = ixp.zerorank1(DIM=3)
    for i in range(3):
        B_notildeU[i] = B_tildeU[i]/sqrtgammaDET

# Step 2.b.i: Define b^mu.
def compute_smallb4U(gammaDD, betaU, alpha, u4U, B_notildeU, sqrt4pi):
    global smallb4U
    import BSSN.ADMBSSN_tofrom_4metric as AB4m
    AB4m.g4DD_ito_BSSN_or_ADM("ADM", gammaDD, betaU, alpha)

    u4D = ixp.zerorank1(DIM=4)
    for mu in range(4):
        for nu in range(4):
            u4D[mu] += AB4m.g4DD[mu][nu] * u4U[nu]
    smallb4U = ixp.zerorank1(DIM=4)
    u4_dot_B_notilde = sp.sympify(0)
    for i in range(3):
        u4_dot_B_notilde += u4D[i + 1] * B_notildeU[i]

    # b^0 = (u_j B^j)/[alpha * sqrt(4 pi)]
    smallb4U[0] = u4_dot_B_notilde / (alpha * sqrt4pi)
    # b^i = [B^i + (u_j B^j) u^i]/[alpha * u^0 * sqrt(4 pi)]
    for i in range(3):
        smallb4U[i + 1] = (B_notildeU[i] + u4_dot_B_notilde * u4U[i + 1]) / (alpha * u4U[0] * sqrt4pi)

# Step 2.b.ii: Define b^mu when u4 and B are orthogonal
def compute_smallb4U_with_driftvU_for_FFE(gammaDD,betaU,alpha, u4U,B_notildeU, sqrt4pi):
    global smallb4_with_driftv_for_FFE_U
    import BSSN.ADMBSSN_tofrom_4metric as AB4m
    AB4m.g4DD_ito_BSSN_or_ADM("ADM",gammaDD,betaU,alpha)

    u4D = ixp.zerorank1(DIM=4)
    for mu in range(4):
        for nu in range(4):
            u4D[mu] += AB4m.g4DD[mu][nu]*u4U[nu]
    smallb4_with_driftv_for_FFE_U = ixp.zerorank1(DIM=4)

    # b^0 = 0
    smallb4_with_driftv_for_FFE_U[0] = 0
    # b^i = B^i / [alpha * u^0 * sqrt(4 pi)]
    for i in range(3):
        smallb4_with_driftv_for_FFE_U[i+1] = B_notildeU[i] / (alpha*u4U[0]*sqrt4pi)

# Step 2.c: Define b^2.
def compute_smallbsquared(gammaDD, betaU, alpha, smallb4U):
    global smallbsquared
    import BSSN.ADMBSSN_tofrom_4metric as AB4m
    AB4m.g4DD_ito_BSSN_or_ADM("ADM", gammaDD, betaU, alpha)

    smallbsquared = sp.sympify(0)
    for mu in range(4):
        for nu in range(4):
            smallbsquared += AB4m.g4DD[mu][nu] * smallb4U[mu] * smallb4U[nu]

# Step 3: Define the electromagnetic stress-energy tensor T_{EM}^{mu nu}
# Step 3.a: Define T_{EM}^{mu nu} (a 4-dimensional tensor)
def compute_TEM4UU(gammaDD,betaU,alpha, smallb4U, smallbsquared,u4U):
    global TEM4UU

    # Then define g^{mu nu} in terms of the ADM quantities:
    import BSSN.ADMBSSN_tofrom_4metric as AB4m
    AB4m.g4UU_ito_BSSN_or_ADM("ADM",gammaDD,betaU,alpha)
    # Finally compute T^{mu nu}
    TEM4UU = ixp.zerorank2(DIM=4)
    for mu in range(4):
        for nu in range(4):
            TEM4UU[mu][nu] = smallbsquared*u4U[mu]*u4U[nu] \
                             + sp.Rational(1,2)*smallbsquared*AB4m.g4UU[mu][nu] \
                             - smallb4U[mu]*smallb4U[nu]


# Step 3.b: Define T^{mu}_{nu} (a 4-dimensional tensor)
def compute_TEM4UD(gammaDD,betaU,alpha, TEM4UU):
    global TEM4UD
    # Next compute T^mu_nu = T^{mu delta} g_{delta nu}, needed for S_tilde flux.
    # First we'll need g_{alpha nu} in terms of ADM quantities:
    import BSSN.ADMBSSN_tofrom_4metric as AB4m
    AB4m.g4DD_ito_BSSN_or_ADM("ADM",gammaDD,betaU,alpha)
    TEM4UD = ixp.zerorank2(DIM=4)
    for mu in range(4):
        for nu in range(4):
            for delta in range(4):
                TEM4UD[mu][nu] += TEM4UU[mu][delta]*AB4m.g4DD[delta][nu]

def compute_AD_flux_term(sqrtgammaDET,driftvU,BU):
    # Levi-Civita tensor for cross products
    LeviCivitaDDD = ixp.LeviCivitaTensorDDD_dim3_rank3(sqrtgammaDET)
    global A_fluxD
    A_fluxD = ixp.zerorank1()
    for i in range(3):
        for j in range(3):
            for k in range(3):
                # \epsilon_{ijk} v^j B^k
                A_fluxD[i] += LeviCivitaDDD[i][j][k]*driftvU[j]*BU[k]

def compute_AD_source_term_operand_for_FD(sqrtgammaDET,betaU,alpha,psi6Phi,AD):
    Phi = psi6Phi/sqrtgammaDET
    global AevolParen
    # \alpha \Phi - \beta^j A_j
    AevolParen = alpha * Phi
    for j in range(3):
        AevolParen += -betaU[j] * AD[j]

def compute_psi6Phi_rhs_flux_term_operand(gammaDD,sqrtgammaDET,betaU,alpha,AD,psi6Phi):
    gammaUU,_gammaDET = ixp.symm_matrix_inverter3x3(gammaDD) # _gammaDET unused.
    AU = ixp.zerorank1()
    # Raise the index on A in the usual way:
    for i in range(3):
        for j in range(3):
            AU[i] += gammaUU[i][j] * AD[j]

    global PhievolParenU
    PhievolParenU = ixp.zerorank1(DIM=3)

    for j in range(3):
        # \alpha\sqrt{\gamma}A^j - \beta^j [\sqrt{\gamma} \Phi]
        PhievolParenU[j] += alpha*sqrtgammaDET*AU[j] - betaU[j]*psi6Phi

def compute_psi6Phi_rhs_damping_term(alpha,psi6Phi,xi_damping):
    # - \xi \alpha [\sqrt{\gamma} \Phi]
    # Combine the divergence and the damping term
    global psi6Phi_damping
    psi6Phi_damping = - xi_damping * alpha * psi6Phi

def generate_everything_for_UnitTesting():
    # First define hydrodynamical quantities
    u4U = ixp.declarerank1("u4U", DIM=4)
    B_tildeU = ixp.declarerank1("B_tildeU", DIM=3)

    # Then ADM quantities
    gammaDD = ixp.declarerank2("gammaDD", "sym01", DIM=3)
    betaU = ixp.declarerank1("betaU", DIM=3)
    alpha = sp.symbols('alpha', real=True)

    # Then numerical constant
    sqrt4pi = sp.symbols('sqrt4pi', real=True)

    # First compute stress-energy tensor T4UU and T4UD:
    import GRHD.equations as GHeq
    GHeq.compute_sqrtgammaDET(gammaDD)
    compute_B_notildeU(GHeq.sqrtgammaDET, B_tildeU)
    compute_smallb4U(gammaDD, betaU, alpha, u4U, B_notildeU, sqrt4pi)
    compute_smallb4U_with_driftvU_for_FFE(gammaDD, betaU, alpha, u4U, B_notildeU, sqrt4pi)
    compute_smallbsquared(gammaDD, betaU, alpha, smallb4U)

    compute_TEM4UU(gammaDD, betaU, alpha, smallb4U, smallbsquared, u4U)
    compute_TEM4UD(gammaDD, betaU, alpha, TEM4UU)

    # Compute conservative variables in terms of primitive variables
    GHeq.compute_S_tildeD(alpha, GHeq.sqrtgammaDET, TEM4UD)
    global S_tildeD
    S_tildeD = GHeq.S_tildeD

    # Next compute fluxes of conservative variables
    GHeq.compute_S_tilde_fluxUD(alpha, GHeq.sqrtgammaDET, TEM4UD)
    global S_tilde_fluxUD
    S_tilde_fluxUD = GHeq.S_tilde_fluxUD

    # Then declare derivatives & compute g4DDdD
    gammaDD_dD = ixp.declarerank3("gammaDD_dD", "sym01", DIM=3)
    betaU_dD = ixp.declarerank2("betaU_dD", "nosym", DIM=3)
    alpha_dD = ixp.declarerank1("alpha_dD", DIM=3)
    GHeq.compute_g4DD_zerotimederiv_dD(gammaDD, betaU, alpha, gammaDD_dD, betaU_dD, alpha_dD)

    # Finally compute source terms on tau_tilde and S_tilde equations
    GHeq.compute_S_tilde_source_termD(alpha, GHeq.sqrtgammaDET, GHeq.g4DD_zerotimederiv_dD, TEM4UU)
    global S_tilde_source_termD
    S_tilde_source_termD = GHeq.S_tilde_source_termD
