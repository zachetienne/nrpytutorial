# As documented in the NRPy+ tutorial module
#   Tutorial-GRHD_Equations-Cartesian-new.ipynb
#   this module will construct useful quantities
#   for the IllinoisGRMHD implementation of
#   general relativistic force-free
#   electrodynamics (GRFFE)
#
# Authors: Zachariah B. Etienne
#          zachetie **at** gmail **dot* com
#          Patrick Nelson

# Step 1: Import needed core NRPy+ modules
from outputC import *            # NRPy+: Core C code output module
import indexedexp as ixp         # NRPy+: Symbolic indexed expression (e.g., tensors, vectors, etc.) support

# Step 2: Define needed quantities for T_{EM}^{mu nu}, the EM part of the stress-energy tensor
# Step 2.a: Define B^i = Btilde^i / sqrt(gamma)
def compute_B_notildeU(sqrtgammaDET, B_tildeU):
    global B_notildeU
    B_notildeU = ixp.zerorank1(DIM=3)
    for i in range(3):
        B_notildeU[i] = B_tildeU[i]/sqrtgammaDET

# Step 2.b: Define b^mu.
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
    # b^i = [B^i + (u_j B^j)]/[alpha * u^0 * sqrt(4 pi)]
    for i in range(3):
        smallb4U[i + 1] = (B_notildeU[i] + u4_dot_B_notilde * u4U[i + 1]) / (alpha * u4U[0] * sqrt4pi)


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
                             + smallb4U[mu]*smallb4U[nu]

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