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

# Step 2.a: Define T^{mu nu} (a 4-dimensional tensor)
def compute_TEM4UU(gammaDD,betaU,alpha, smallb4U, smallb2,u4U):
    global TEM4UU

    # Then define g^{mu nu} in terms of the ADM quantities:
    import BSSN.ADMBSSN_tofrom_4metric as AB4m
    AB4m.g4UU_ito_BSSN_or_ADM("ADM",gammaDD,betaU,alpha)

    # Finally compute T^{mu nu}
    TEM4UU = ixp.zerorank2(DIM=4)
    for mu in range(4):
        for nu in range(4):
            TEM4UU[mu][nu] = smallb2*u4U[mu]*u4U[nu] \
                             + sp.Rational(1,2)*smallb2*AB4m.g4UU[mu][nu] \
                             + smallb4U[mu]*smallb4U[nu]

# Step 2.b: Define T^{mu}_{nu} (a 4-dimensional tensor)
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
                
# Step 3: Writing the conservative variables in terms of the primitive variables
def compute_sqrtgammaDET(gammaDD):
    global sqrtgammaDET
    gammaUU, gammaDET = ixp.symm_matrix_inverter3x3(gammaDD)
    sqrtgammaDET = sp.sqrt(gammaDET)

def compute_S_tildeD(alpha, sqrtgammaDET, T4UD):
    global S_tildeD
    S_tildeD = ixp.zerorank1(DIM=3)
    for i in range(3):
        S_tildeD[i] = alpha*sqrtgammaDET*T4UD[0][i+1]
        
# Step 4.d: S_tilde flux
def compute_S_tilde_fluxUD(gammaDD,betaU,alpha, sqrtgammaDET, T4UD):
    global S_tilde_fluxUD
    S_tilde_fluxUD = ixp.zerorank2(DIM=3)
    for j in range(3):
        for i in range(3):
            S_tilde_fluxUD[j][i] = alpha*sqrtgammaDET*T4UD[j+1][i+1]
            
def compute_g4DDdD(gammaDD,betaU,alpha, gammaDD_dD,betaU_dD,alpha_dD):
    global g4DDdD
    # Eq. 2.121 in B&S
    betaD = ixp.zerorank1()
    for i in range(3):
        for j in range(3):
            betaD[i] += gammaDD[i][j]*betaU[j]

    betaDdD = ixp.zerorank2()
    for i in range(3):
        for j in range(3):
            for k in range(3):
                # Recall that betaD[i] = gammaDD[i][j]*betaU[j] (Eq. 2.121 in B&S)
                betaDdD[i][k] += gammaDD_dD[i][j][k]*betaU[j] + gammaDD[i][j]*betaU_dD[j][k]

    # Eq. 2.122 in B&S
    g4DDdD = ixp.zerorank3(DIM=4)
    for k in range(3):
        # Recall that g4DD[0][0] = -alpha^2 + betaU[j]*betaD[j]
        g4DDdD[0][0][k+1] += -2*alpha*alpha_dD[k]
        for j in range(3):
            g4DDdD[0][0][k+1] += betaU_dD[j][k]*betaD[j] + betaU[j]*betaDdD[j][k]

    for i in range(3):
        for k in range(3):
            # Recall that g4DD[i][0] = g4DD[0][i] = betaD[i]
            g4DDdD[i+1][0][k+1] = g4DDdD[0][i+1][k+1] = betaDdD[i][k]
    for i in range(3):
        for j in range(3):
            for k in range(3):
                # Recall that g4DD[i][j] = gammaDD[i][j]
                g4DDdD[i+1][j+1][k+1] = gammaDD_dD[i][j][k]
                
# Step 5.b.ii: Compute S_tilde source term
def compute_S_tilde_source_termD(alpha, sqrtgammaDET,g4DDdD, T4UU):
    global S_tilde_source_termD
    S_tilde_source_termD = ixp.zerorank1(DIM=3)
    for i in range(3):
        for mu in range(4):
            for nu in range(4):
                S_tilde_source_termD[i] += sp.Rational(1,2)*alpha*sqrtgammaDET*T4UU[mu][nu]*g4DDdD[mu][nu][i+1]
                
# Step 6: Convert v^i into u^\mu, applying a speed limiter
def u4U_in_terms_of_vU_apply_speed_limit(alpha,betaU,gammaDD, vU):
    global u4_ito_3velsU
    
    import GiRaFFE_HO.Stilde_flux as GSf
    ValenciavU = ixp.zerorank1()
    for i in range(3):
        ValenciavU[i] = (vU[i] + betaU[i])/alpha

    GSf.compute_u0_noif(gammaDD,alpha,ValenciavU)
    u4_ito_3velsU = ixp.zerorank1(DIM=4)
    u4_ito_3velsU[0] = GSf.rescaledu0
    for i in range(3):
        u4_ito_3velsU[i+1] = GSf.rescaledu0 * (alpha * GSf.rescaledValenciavU[i] - betaU[i])