# As documented in the NRPy+ tutorial module
#   Tutorial-GRMHD_Equations-Cartesian.ipynb
#   this module will construct useful quantities
#   for the IllinoisGRMHD implementation of
#   general relativistic hydrodynamics (GRHD)
#
# Author: Zachariah B. Etienne
#         zachetie **at** gmail **dot* com

# Step 1: import all needed modules from NRPy+/Python:
import sympy as sp               # SymPy: The Python computer algebra package upon which NRPy+ depends
import indexedexp as ixp         # NRPy+: Symbolic indexed expression (e.g., tensors, vectors, etc.) support
import GRHD.equations as GRHD    # NRPy+: Useful functions for general relativistic hydrodynamics
import GRFFE.equations as GRFFE  # NRPy+: Useful functions for general relativistic force-free electrodynamics

# Step 2: Define the stress-energy tensor
# Step 2.a: Define the GRMHD T^{mu nu} (a 4-dimensional tensor)
def compute_GRMHD_T4UU(gammaDD, betaU, alpha, rho_b, P, epsilon, u4U, smallb4U, smallbsquared):
    global GRHDT4UU
    global GRFFET4UU
    global T4UU

    GRHD.compute_T4UU(gammaDD, betaU, alpha, rho_b, P, epsilon, u4U)
    GRFFE.compute_TEM4UU(gammaDD, betaU, alpha, smallb4U, smallbsquared, u4U)

    GRHDT4UU = ixp.zerorank2(DIM=4)
    GRFFET4UU = ixp.zerorank2(DIM=4)
    T4UU = ixp.zerorank2(DIM=4)
    for mu in range(4):
        for nu in range(4):
            GRHDT4UU[mu][nu] = GRHD.T4UU[mu][nu]
            GRFFET4UU[mu][nu] = GRFFE.TEM4UU[mu][nu]
            T4UU[mu][nu] = GRHD.T4UU[mu][nu] + GRFFE.TEM4UU[mu][nu]


# Step 2.b: Define T^{mu}_{nu} (a 4-dimensional tensor)
def compute_GRMHD_T4UD(gammaDD, betaU, alpha, GRHDT4UU, GRFFET4UU):
    global T4UD

    GRHD.compute_T4UD(gammaDD, betaU, alpha, GRHDT4UU)
    GRFFE.compute_TEM4UD(gammaDD, betaU, alpha, GRFFET4UU)

    T4UD = ixp.zerorank2(DIM=4)
    for mu in range(4):
        for nu in range(4):
            T4UD[mu][nu] = GRHD.T4UD[mu][nu] + GRFFE.TEM4UD[mu][nu]

# Step 3: Declare ADM and hydrodynamical input variables, and
#         construct all terms in GRMHD equations
def generate_everything_for_UnitTesting():
    # First define hydrodynamical quantities
    u4U = ixp.declarerank1("u4U", DIM=4)
    rho_b,P,epsilon = sp.symbols('rho_b P epsilon',real=True)
    B_tildeU = ixp.declarerank1("B_tildeU", DIM=3)

    # Then ADM quantities
    gammaDD = ixp.declarerank2("gammaDD","sym01",DIM=3)
    KDD     = ixp.declarerank2("KDD"    ,"sym01",DIM=3)
    betaU   = ixp.declarerank1("betaU", DIM=3)
    alpha   = sp.symbols('alpha', real=True)

    # Then numerical constant
    sqrt4pi = sp.symbols('sqrt4pi', real=True)

    # First compute smallb4U & smallbsquared from BtildeU, which are needed
    #      for GRMHD stress-energy tensor T4UU and T4UD:
    GRHD.compute_sqrtgammaDET(gammaDD)
    GRFFE.compute_B_notildeU(GRHD.sqrtgammaDET,      B_tildeU)
    GRFFE.compute_smallb4U(     gammaDD,betaU,alpha, u4U,GRFFE.B_notildeU, sqrt4pi)
    GRFFE.compute_smallbsquared(gammaDD,betaU,alpha, GRFFE.smallb4U)

    # Then compute the GRMHD stress-energy tensor:
    compute_GRMHD_T4UU(gammaDD,betaU,alpha, rho_b,P,epsilon,u4U, GRFFE.smallb4U, GRFFE.smallbsquared)
    compute_GRMHD_T4UD(gammaDD,betaU,alpha, GRHDT4UU,GRFFET4UU)

    # Compute conservative variables in terms of primitive variables
    global rho_star,tau_tilde,S_tildeD
    GRHD.compute_rho_star( alpha, GRHD.sqrtgammaDET, rho_b,u4U)
    GRHD.compute_tau_tilde(alpha, GRHD.sqrtgammaDET, T4UU,GRHD.rho_star)
    GRHD.compute_S_tildeD( alpha, GRHD.sqrtgammaDET, T4UD)
    rho_star = GRHD.rho_star
    tau_tilde = GRHD.tau_tilde
    S_tildeD = GRHD.S_tildeD

    # Then compute v^i from u^mu
    GRHD.compute_vU_from_u4U__no_speed_limit(u4U)

    # Next compute fluxes of conservative variables
    global rho_star_fluxU,tau_tilde_fluxU,S_tilde_fluxUD
    GRHD.compute_rho_star_fluxU(                           GRHD.vU,     GRHD.rho_star)
    GRHD.compute_tau_tilde_fluxU(alpha, GRHD.sqrtgammaDET, GRHD.vU,T4UU,GRHD.rho_star)
    GRHD.compute_S_tilde_fluxUD( alpha, GRHD.sqrtgammaDET,         T4UD)
    rho_star_fluxU = GRHD.rho_star_fluxU
    tau_tilde_fluxU = GRHD.tau_tilde_fluxU
    S_tilde_fluxUD = GRHD.S_tilde_fluxUD

    # Then declare derivatives & compute g4DD_zerotimederiv_dD
    gammaDD_dD = ixp.declarerank3("gammaDD_dD","sym01",DIM=3)
    betaU_dD   = ixp.declarerank2("betaU_dD"  ,"nosym",DIM=3)
    alpha_dD   = ixp.declarerank1("alpha_dD"          ,DIM=3)
    GRHD.compute_g4DD_zerotimederiv_dD(gammaDD,betaU,alpha, gammaDD_dD,betaU_dD,alpha_dD)

    # Then compute source terms on tau_tilde and S_tilde equations
    global s_source_term,S_tilde_source_termD
    GRHD.compute_s_source_term(KDD,betaU,alpha, GRHD.sqrtgammaDET,alpha_dD,                   T4UU)
    GRHD.compute_S_tilde_source_termD(   alpha, GRHD.sqrtgammaDET,GRHD.g4DD_zerotimederiv_dD, T4UU)
    s_source_term = GRHD.s_source_term
    S_tilde_source_termD = GRHD.S_tilde_source_termD
