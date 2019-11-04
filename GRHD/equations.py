# As documented in the NRPy+ tutorial module
#   Tutorial-GRHD_Equations-Cartesian-new.ipynb
#   this module will construct useful quantities
#   for the IllinoisGRMHD implementation of
#   general relativistic hydrodynamics (GRHD)
#
# Authors: Zachariah B. Etienne
#          zachetie **at** gmail **dot* com
#          Patrick Nelson

# Step 1: import all needed modules from NRPy+/Python:
from outputC import *            # NRPy+: Core C code output module
import indexedexp as ixp         # NRPy+: Symbolic indexed expression (e.g., tensors, vectors, etc.) support

# Step 2: Define the stress-energy tensor

# Step 2.a: First define h, the enthalpy:
def compute_enthalpy(rho_b,P,epsilon):
    global h
    h = 1 + epsilon + P/rho_b

# Step 2.b: Define T^{mu nu} (a 4-dimensional tensor)
def compute_T4UU(gammaDD,betaU,alpha, rho_b,P,epsilon,u4U):
    global T4UU

    compute_enthalpy(rho_b,P,epsilon)
    # Then define g^{mu nu} in terms of the ADM quantities:
    import BSSN.ADMBSSN_tofrom_4metric as AB4m
    AB4m.g4UU_ito_BSSN_or_ADM("ADM",gammaDD,betaU,alpha)

    # Finally compute T^{mu nu}
    T4UU = ixp.zerorank2(DIM=4)
    for mu in range(4):
        for nu in range(4):
            T4UU[mu][nu] = rho_b * h * u4U[mu]*u4U[nu] + P*AB4m.g4UU[mu][nu]

# Step 2.c: Define T^{mu}_{nu} (a 4-dimensional tensor)
def compute_T4UD(gammaDD,betaU,alpha, T4UU):
    global T4UD
    # Next compute T^mu_nu = T^{mu delta} g_{delta nu}, needed for S_tilde flux.
    # First we'll need g_{alpha nu} in terms of ADM quantities:
    import BSSN.ADMBSSN_tofrom_4metric as AB4m
    AB4m.g4DD_ito_BSSN_or_ADM("ADM",gammaDD,betaU,alpha)
    T4UD = ixp.zerorank2(DIM=4)
    for mu in range(4):
        for nu in range(4):
            for delta in range(4):
                T4UD[mu][nu] += T4UU[mu][delta]*AB4m.g4DD[delta][nu]

# Step 3: Writing the conservative variables in terms of the primitive variables
def compute_sqrtgammaDET(gammaDD):
    global sqrtgammaDET
    gammaUU, gammaDET = ixp.symm_matrix_inverter3x3(gammaDD)
    sqrtgammaDET = sp.sqrt(gammaDET)

def compute_rho_star(alpha, sqrtgammaDET, rho_b,u4U):
    global rho_star
    # Compute rho_star:
    rho_star = alpha*sqrtgammaDET*rho_b*u4U[0]

def compute_tau_tilde(alpha, sqrtgammaDET, T4UU,rho_star):
    global tau_tilde
    tau_tilde = alpha**2*sqrtgammaDET*T4UU[0][0] - rho_star

def compute_S_tildeD(alpha, sqrtgammaDET, T4UD):
    global S_tildeD
    S_tildeD = ixp.zerorank1(DIM=3)
    for i in range(3):
        S_tildeD[i] = alpha*sqrtgammaDET*T4UD[0][i+1]

# Step 4: Define the fluxes for the GRHD equations
# Step 4.a: vU from u4U may be needed for computing rho_star_flux from u4U
def compute_vU_from_u4U__no_speed_limit(u4U):
    global vU
    # Now compute v^i = u^i/u^0:
    vU = ixp.zerorank1(DIM=3)
    for j in range(3):
        vU[j] = u4U[j+1]/u4U[0]

# Step 4.b: rho_star flux
def compute_rho_star_fluxU(vU, rho_star):
    global rho_star_fluxU
    rho_star_fluxU = ixp.zerorank1(DIM=3)
    for j in range(3):
        rho_star_fluxU[j] = rho_star*vU[j]

# Step 4.c: tau_tilde flux
def compute_tau_tilde_fluxU(alpha, sqrtgammaDET, vU,T4UU):
    global tau_tilde_fluxU
    tau_tilde_fluxU = ixp.zerorank1(DIM=3)
    for j in range(3):
        tau_tilde_fluxU[j] = alpha**2*sqrtgammaDET*T4UU[0][j+1] - rho_star*vU[j]

# Step 4.d: S_tilde flux
def compute_S_tilde_fluxUD(alpha, sqrtgammaDET, T4UD):
    global S_tilde_fluxUD
    S_tilde_fluxUD = ixp.zerorank2(DIM=3)
    for j in range(3):
        for i in range(3):
            S_tilde_fluxUD[j][i] = alpha*sqrtgammaDET*T4UD[j+1][i+1]

# Step 5: Define source terms on RHSs of GRHD equations
# Step 5.a: tau_tilde RHS source term s
def compute_s_source_term(KDD,betaU,alpha, sqrtgammaDET,alpha_dD, T4UU):
    global s_source_term
    s_source_term = sp.sympify(0)
    for i in range(3):
        for j in range(3):
            s_source_term += (T4UU[0][0]*betaU[i]*betaU[j] + 2*T4UU[0][i+1]*betaU[j] + T4UU[i+1][j+1])*KDD[i][j]

    for i in range(3):
        s_source_term += -(T4UU[0][0]*betaU[i] + T4UU[0][i+1])*alpha_dD[i]

    s_source_term *= alpha*sqrtgammaDET

# Step 5.b: Define source term on RHS of $\tilde{S}_i$ equation
# Step 5.b.i: Compute g_{mu nu, i}, needed for the S tilde source term
def compute_g4DD_zerotimederiv_dD(gammaDD,betaU,alpha, gammaDD_dD,betaU_dD,alpha_dD):
    global g4DD_zerotimederiv_dD
    # Eq. 2.121 in B&S
    betaD = ixp.zerorank1()
    for i in range(3):
        for j in range(3):
            betaD[i] += gammaDD[i][j]*betaU[j]

    # gammaDD_dD = ixp.declarerank3("gammaDD_dDD","sym12",DIM=3)
    # betaU_dD   = ixp.declarerank2("betaU_dD"   ,"nosym",DIM=3)
    betaDdD = ixp.zerorank2()
    for i in range(3):
        for j in range(3):
            for k in range(3):
                # Recall that betaD[i] = gammaDD[i][j]*betaU[j] (Eq. 2.121 in B&S)
                betaDdD[i][k] += gammaDD_dD[i][j][k]*betaU[j] + gammaDD[i][j]*betaU_dD[j][k]

    # Eq. 2.122 in B&S
    g4DD_zerotimederiv_dD = ixp.zerorank3(DIM=4)
    for k in range(3):
        # Recall that g4DD[0][0] = -alpha^2 + betaU[j]*betaD[j]
        g4DD_zerotimederiv_dD[0][0][k+1] += -2*alpha*alpha_dD[k]
        for j in range(3):
            g4DD_zerotimederiv_dD[0][0][k+1] += betaU_dD[j][k]*betaD[j] + betaU[j]*betaDdD[j][k]

    for i in range(3):
        for k in range(3):
            # Recall that g4DD[i][0] = g4DD[0][i] = betaD[i]
            g4DD_zerotimederiv_dD[i+1][0][k+1] = g4DD_zerotimederiv_dD[0][i+1][k+1] = betaDdD[i][k]
    for i in range(3):
        for j in range(3):
            for k in range(3):
                # Recall that g4DD[i][j] = gammaDD[i][j]
                g4DD_zerotimederiv_dD[i+1][j+1][k+1] = gammaDD_dD[i][j][k]

# Step 5.b.ii: S_tilde source terms
def compute_S_tilde_source_termD(alpha, sqrtgammaDET,g4DD_zerotimederiv_dD, T4UU):
    global S_tilde_source_termD
    S_tilde_source_termD = ixp.zerorank1(DIM=3)
    for i in range(3):
        for mu in range(4):
            for nu in range(4):
                S_tilde_source_termD[i] += sp.Rational(1,2)*alpha*sqrtgammaDET*T4UU[mu][nu]*g4DD_zerotimederiv_dD[mu][nu][i+1]

# Step 6: Convert v^i into u^\mu, applying a speed limiter, courtesy Patrick Nelson
def u4U_in_terms_of_vU_apply_speed_limit(alpha, betaU, gammaDD, vU):
    global u4_ito_3velsU

    import GiRaFFE_HO.Stilde_flux as GSf
    ValenciavU = ixp.zerorank1()
    for i in range(3):
        ValenciavU[i] = (vU[i] + betaU[i]) / alpha

    GSf.compute_u0_noif(gammaDD, alpha, ValenciavU)
    u4_ito_3velsU = ixp.zerorank1(DIM=4)
    u4_ito_3velsU[0] = GSf.rescaledu0
    for i in range(3):
        u4_ito_3velsU[i+1] = GSf.rescaledu0 * (alpha * GSf.rescaledValenciavU[i] - betaU[i])

def generate_everything_for_UnitTesting():
    # First define hydrodynamical quantities
    u4U = ixp.declarerank1("u4U", DIM=4)
    rho_b, P, epsilon = sp.symbols('rho_b P epsilon', real=True)

    # Then ADM quantities
    gammaDD = ixp.declarerank2("gammaDD", "sym01", DIM=3)
    KDD = ixp.declarerank2("KDD", "sym01", DIM=3)
    betaU = ixp.declarerank1("betaU", DIM=3)
    alpha = sp.symbols('alpha', real=True)

    # First compute stress-energy tensor T4UU and T4UD:
    compute_T4UU(gammaDD, betaU, alpha, rho_b, P, epsilon, u4U)
    compute_T4UD(gammaDD, betaU, alpha, T4UU)

    # Next sqrt(gamma)
    compute_sqrtgammaDET(gammaDD)

    # Compute conservative variables in terms of primitive variables
    compute_rho_star(alpha, sqrtgammaDET, rho_b, u4U)
    compute_tau_tilde(alpha, sqrtgammaDET, T4UU, rho_star)
    compute_S_tildeD(alpha, sqrtgammaDET, T4UD)

    # Then compute v^i from u^mu
    compute_vU_from_u4U__no_speed_limit(u4U)

    # Next compute fluxes of conservative variables
    compute_rho_star_fluxU(vU, rho_star)
    compute_tau_tilde_fluxU(alpha, sqrtgammaDET, vU, T4UU)
    compute_S_tilde_fluxUD( alpha, sqrtgammaDET,     T4UD)

    # Then declare derivatives & compute g4DD_zerotimederiv_dD
    gammaDD_dD = ixp.declarerank3("gammaDD_dD", "sym01", DIM=3)
    betaU_dD = ixp.declarerank2("betaU_dD", "nosym", DIM=3)
    alpha_dD = ixp.declarerank1("alpha_dD", DIM=3)
    compute_g4DD_zerotimederiv_dD(gammaDD, betaU, alpha, gammaDD_dD, betaU_dD, alpha_dD)

    # Then compute source terms on tau_tilde and S_tilde equations
    compute_s_source_term(KDD, betaU, alpha, sqrtgammaDET, alpha_dD, T4UU)
    compute_S_tilde_source_termD(alpha, sqrtgammaDET, g4DD_zerotimederiv_dD, T4UU)

    # Finally compute the 4-velocities in terms of an input 3-velocity testvU
    testvU = ixp.declarerank1("testvU", DIM=3)
    u4U_in_terms_of_vU_apply_speed_limit(alpha, betaU, gammaDD, testvU)