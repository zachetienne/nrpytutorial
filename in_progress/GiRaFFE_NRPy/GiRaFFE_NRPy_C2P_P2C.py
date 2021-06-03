from outputC import nrpyAbs      # NRPy+: Core C code output module
import sympy as sp               # SymPy: The Python computer algebra package upon which NRPy+ depends
import NRPy_param_funcs as par   # NRPy+: Parameter interface
import grid as gri               # NRPy+: Functions having to do with numerical grids
import indexedexp as ixp         # NRPy+: Symbolic indexed expression (e.g., tensors, vectors, etc.) support
import reference_metric as rfm   # NRPy+: Reference metric support
import GRHD.equations as GRHD    # NRPy+: Generate general relativistic hydrodynamics equations
import GRFFE.equations as GRFFE  # NRPy+: Generate general relativisitic force-free electrodynamics equations

thismodule = __name__

# There are several C parameters that we will need in this module:
M_PI  = par.Cparameters("#define",thismodule,["M_PI"], "")
GAMMA_SPEED_LIMIT = par.Cparameters("REAL",thismodule,"GAMMA_SPEED_LIMIT",10.0) # Default value based on
                                                                                # IllinoisGRMHD.
                                                                                # GiRaFFE default = 2000.0

# There are three fixes in this module; we might not always want to do all (or even any!) of them.
# So, we initialize some NRPy+ parameters to control this behavior
par.initialize_param(par.glb_param(type="bool", module=thismodule, parname="enforce_orthogonality_StildeD_BtildeU", defaultval=True))
par.initialize_param(par.glb_param(type="bool", module=thismodule, parname="enforce_speed_limit_StildeD", defaultval=True))
par.initialize_param(par.glb_param(type="bool", module=thismodule, parname="enforce_current_sheet_prescription", defaultval=True))

def GiRaFFE_NRPy_C2P(StildeD,BU,gammaDD,betaU,alpha):
    GRHD.compute_sqrtgammaDET(gammaDD)
    gammaUU,unusedgammadet = ixp.symm_matrix_inverter3x3(gammaDD)
    BtildeU = ixp.zerorank1()
    for i in range(3):
        # \tilde{B}^i = B^i \sqrt{\gamma}
        BtildeU[i] = GRHD.sqrtgammaDET*BU[i]

    BtildeD = ixp.zerorank1()
    for i in range(3):
        for j in range(3):
            BtildeD[j] += gammaDD[i][j]*BtildeU[i]

    Btilde2 = sp.sympify(0)
    for i in range(3):
        Btilde2 += BtildeU[i]*BtildeD[i]

    global outStildeD
    outStildeD = StildeD
    # Then, enforce the orthogonality:
    if par.parval_from_str("enforce_orthogonality_StildeD_BtildeU"):
        StimesB = sp.sympify(0)
        for i in range(3):
            StimesB += StildeD[i]*BtildeU[i]

        for i in range(3):
            # {\tilde S}_i = {\tilde S}_i - ({\tilde S}_j {\tilde B}^j) {\tilde B}_i/{\tilde B}^2
            outStildeD[i] -= StimesB*BtildeD[i]/Btilde2

    # Calculate \tilde{S}^2:
    Stilde2 = sp.sympify(0)
    for i in range(3):
        for j in range(3):
            Stilde2 += gammaUU[i][j]*outStildeD[i]*outStildeD[j]

    # First we need to compute the factor f:
    # f = \sqrt{(1-\Gamma_{\max}^{-2}){\tilde B}^4/(16 \pi^2 \gamma {\tilde S}^2)}
    speed_limit_factor = sp.sqrt((sp.sympify(1)-GAMMA_SPEED_LIMIT**(-2.0))*Btilde2*Btilde2*sp.Rational(1,16)/\
                                 (M_PI*M_PI*GRHD.sqrtgammaDET*GRHD.sqrtgammaDET*Stilde2))

    import Min_Max_and_Piecewise_Expressions as noif

    # Calculate B^2
    B2 = sp.sympify(0)
    for i in range(3):
        for j in range(3):
            B2 += gammaDD[i][j]*BU[i]*BU[j]

    # Enforce the speed limit on StildeD:
    if par.parval_from_str("enforce_speed_limit_StildeD"):
        for i in range(3):
            outStildeD[i] *= noif.min_noif(sp.sympify(1),speed_limit_factor)

    global ValenciavU
    ValenciavU = ixp.zerorank1()
    # Recompute 3-velocity:
    for i in range(3):
        for j in range(3):
            # \bar{v}^i = 4 \pi \gamma^{ij} {\tilde S}_j / (\sqrt{\gamma} B^2)
            ValenciavU[i] += sp.sympify(4)*M_PI*gammaUU[i][j]*outStildeD[j]/(GRHD.sqrtgammaDET*B2)

    # This number determines how far away (in grid points) we will apply the fix.
    grid_points_from_z_plane = par.Cparameters("REAL",thismodule,"grid_points_from_z_plane",4.0)

    if par.parval_from_str("enforce_current_sheet_prescription"):
        # Calculate the drift velocity
        driftvU = ixp.zerorank1()
        for i in range(3):
            driftvU[i] = alpha*ValenciavU[i] - betaU[i]

        # The direct approach, used by the original GiRaFFE:
        # v^z = -(\gamma_{xz} v^x + \gamma_{yz} v^y) / \gamma_{zz}
        newdriftvU2 = -(gammaDD[0][2]*driftvU[0] + gammaDD[1][2]*driftvU[1])/gammaDD[2][2]
        # Now that we have the z component, it's time to substitute its Valencia form in.
        # Remember, we only do this if abs(z) < (k+0.01)*dz. Note that we add 0.01; this helps
        # avoid floating point errors and division by zero. This is the same as abs(z) - (k+0.01)*dz<0
        coord = nrpyAbs(rfm.xx[2])
        bound =(grid_points_from_z_plane+sp.Rational(1,100))*gri.dxx[2]
        ValenciavU[2] = noif.coord_leq_bound(coord,bound)*(newdriftvU2+betaU[2])/alpha \
                      + noif.coord_greater_bound(coord,bound)*ValenciavU[2]

def GiRaFFE_NRPy_P2C(gammaDD,betaU,alpha,  ValenciavU,BU, sqrt4pi):
    # After recalculating the 3-velocity, we need to update the poynting flux:
    # We'll reset the Valencia velocity, since this will be part of a second call to outCfunction.

    # First compute stress-energy tensor T4UU and T4UD:

    GRHD.compute_sqrtgammaDET(gammaDD)
    GRHD.u4U_in_terms_of_ValenciavU__rescale_ValenciavU_by_applying_speed_limit(alpha, betaU, gammaDD, ValenciavU)

    GRFFE.compute_smallb4U_with_driftvU_for_FFE(gammaDD, betaU, alpha, GRHD.u4U_ito_ValenciavU, BU, sqrt4pi)
    GRFFE.compute_smallbsquared(gammaDD, betaU, alpha, GRFFE.smallb4_with_driftv_for_FFE_U)

    GRFFE.compute_TEM4UU(gammaDD, betaU, alpha, GRFFE.smallb4_with_driftv_for_FFE_U, GRFFE.smallbsquared, GRHD.u4U_ito_ValenciavU)
    GRFFE.compute_TEM4UD(gammaDD, betaU, alpha, GRFFE.TEM4UU)

    # Compute conservative variables in terms of primitive variables
    GRHD.compute_S_tildeD(alpha, GRHD.sqrtgammaDET, GRFFE.TEM4UD)

    global StildeD
    StildeD = GRHD.S_tildeD
