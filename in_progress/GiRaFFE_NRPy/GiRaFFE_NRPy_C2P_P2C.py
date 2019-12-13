from outputC import *            # NRPy+: Core C code output module
import finite_difference as fin  # NRPy+: Finite difference C code generation module
import NRPy_param_funcs as par   # NRPy+: Parameter interface
import grid as gri               # NRPy+: Functions having to do with numerical grids
import loop as lp                # NRPy+: Generate C code loops
import indexedexp as ixp         # NRPy+: Symbolic indexed expression (e.g., tensors, vectors, etc.) support
import reference_metric as rfm   # NRPy+: Reference metric support
import cmdline_helper as cmd     # NRPy+: Multi-platform Python command-line interface

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

def GiRaFFE_NRPy_C2P(StildeD,BU,gammaDD,gammaUU,gammadet,betaU,alpha):
    sqrtgammadet = sp.sqrt(gammadet)
    BtildeU = ixp.zerorank1()
    for i in range(3):
        # \tilde{B}^i = B^i \sqrt{\gamma}
        BtildeU[i] = sqrtgammadet*BU[i]

    BtildeD = ixp.zerorank1()
    for i in range(3):
        for j in range(3):
            BtildeD[j] += gammaDD[i][j]*BtildeU[i]

    Btilde2 = sp.sympify(0)
    for i in range(3):
        Btilde2 += BtildeU[i]*BtildeD[i]
        
    StimesB = sp.sympify(0)
    for i in range(3):
        StimesB += StildeD[i]*BtildeU[i]

    global outStildeD
    outStildeD = StildeD
    # Then, enforce the orthogonality:
    if par.parval_from_str("enforce_orthogonality_StildeD_BtildeU"):
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
    speed_limit_factor = sp.sqrt((1.0-GAMMA_SPEED_LIMIT**(-2.0))\
                                 *Btilde2*Btilde2*sp.Rational(1,16)/(M_PI*M_PI*gammadet*Stilde2))

    def min_noif(a,b):
        # This returns the minimum of a and b
        # If a>b, then we get 0.5*(a+b-a+b) = b
        # If b>a, then we get 0.5*(a+b+a-b) = a
        return sp.Rational(1,2)*(a+b-nrpyAbs(a-b))

    # Calculate B^2
    B2 = sp.sympify(0)
    for i in range(3):
        for j in range(3):
            B2 += gammaDD[i][j]*BU[i]*BU[j]

    # Enforce the speed limit on StildeD:
    if par.parval_from_str("enforce_speed_limit_StildeD"):
        for i in range(3):
            outStildeD[i] *= min_noif(1.0,speed_limit_factor)

    global ValenciavU
    ValenciavU = ixp.zerorank1()
    if par.parval_from_str("enforce_orthogonality_StildeD_BtildeU") or par.parval_from_str("enforce_speed_limit_StildeD"):
        # Recompute 3-velocity:
        for i in range(3):
            for j in range(3):
                # \bar{v}^i = 4 \pi \gamma^{ij} {\tilde S}_j / (\sqrt{\gamma} B^2)
                ValenciavU[i] += sp.sympify(4.0)*M_PI*gammaUU[i][j]*outStildeD[j]/(sqrtgammadet*B2)

    # We will use once more the trick from above with min and max without if. However, we we'll need a function
    # that returns either 0 or 1, so as to choose between two otherwise mathetmatically unrelated branches. 
    def max_normal0(a):
        # If a>0, return 1. Otherwise, return 0. This defines a 'greater than' branch.
        # WILL BREAK if a = 0. 
        return (a+nrpyAbs(a))/(a+a)

    def min_normal0(a):
        # If a<0, return 1. Otherwise, return 0. This defines a 'less than' branch.
        # WILL BREAK if a = 0. 
        return (a-nrpyAbs(a))/(a+a)

    # This number determines how far away (in grid points) we will apply the fix.
    grid_points_from_z_plane = par.Cparameters("REAL",thismodule,"grid_points_from_z_plane",4.0)
    # Set the Cartesian normal vector. This can be expanded later to arbitrary sheets and coordinate systems.
    nU = ixp.zerorank1()
    nU[2] = 1
    # Lower the index, as usual:
    nD = ixp.zerorank1()
    for i in range(3):
        for j in range(3):
            nD[i] = gammaDD[i][j]*nU[j]

    if par.parval_from_str("enforce_current_sheet_prescription"):
        # Calculate the drift velocity
        driftvU = ixp.declarerank1("driftvU")

        inner_product = sp.sympify(0)
        for i in range(3):
            inner_product += driftvU[i]*nD[i] # This is the portion of the drift velocity normal to the z plane
                                              # In flat space, this is just v^z
        # We'll use a sympy utility to solve for v^z. This should make it easier to generalize later
        newdriftvU2 = sp.solve(inner_product,driftvU[2])
        newdriftvU2 = newdriftvU2[0] # In flat space this reduces to v^z=0
        for i in range(3):
            # Now, we substitute drift velocity in terms of our preferred Valencia velocity
            newdriftvU2 = newdriftvU2.subs(driftvU[i],alpha*ValenciavU[i]-betaU[i])
        # Now that we have the z component, it's time to substitute its Valencia form in.
        # Remember, we only do this if abs(z) < (k+0.01)*dz. Note that we add 0.01; this helps
        # avoid floating point errors and division by zero. This is the same as abs(z) - (k+0.01)*dz<0
        boundary = nrpyAbs(rfm.xx[2]) - (grid_points_from_z_plane+sp.sympify(0.01))*gri.dxx[2]
        ValenciavU[2] = min_normal0(boundary)*(newdriftvU2+betaU[2])/alpha \
                         + max_normal0(boundary)*ValenciavU[2]

import GRFFE.equations as GRFFE
import GRHD.equations as GRHD

def GiRaFFE_NRPy_P2C(gammadet,gammaDD,betaU,alpha,  ValenciavU,BU, sqrt4pi):
    # After recalculating the 3-velocity, we need to update the poynting flux:
    # We'll reset the Valencia velocity, since this will be part of a second call to outCfunction.

    # First compute stress-energy tensor T4UU and T4UD:
    
    GRHD.u4U_in_terms_of_ValenciavU__rescale_ValenciavU_by_applying_speed_limit(alpha, betaU, gammaDD, ValenciavU)
    GRFFE.compute_smallb4U(gammaDD, betaU, alpha, GRHD.u4U_ito_ValenciavU, BU, sqrt4pi)
#     GRFFE.compute_smallbsquared(gammaDD, betaU, alpha, GRFFE.smallb4U)

#     GRFFE.compute_TEM4UU(gammaDD, betaU, alpha, GRFFE.smallb4U, GRFFE.smallbsquared, GRHD.u4U_ito_ValenciavU)
#     GRFFE.compute_TEM4UD(gammaDD, betaU, alpha, GRFFE.TEM4UU)

#     # Compute conservative variables in terms of primitive variables
#     GRHD.compute_S_tildeD(alpha, sp.sqrt(gammadet), GRFFE.TEM4UD)

    B2 = sp.sympify(0)
    for i in range(3):
        for j in range(3):
            B2 += gammaDD[i][j]*BU[i]*BU[j]

    global StildeD
#     StildeD = GRHD.S_tildeD
    StildeD = ixp.zerorank1()
    # \gamma_{ij} \frac{\bar{v}^j \sqrt{\gamma}B^2}{4 \pi}
    for i in range(3):
        for j in range(3):
            StildeD[i] += gammaDD[i][j]*ValenciavU[j]*sp.sqrt(gammadet)*B2/(sqrt4pi*sqrt4pi)
            
    # Missing term?
    # define u_k B^k = g{\mu k} u^\mu B^k
    import BSSN.ADMBSSN_tofrom_4metric as AB4m
    AB4m.g4DD_ito_BSSN_or_ADM("ADM",gammaDD,betaU,alpha)
    udotB = sp.sympify(0)
    for mu in range(4):
        for k in range(3):
            udotB += AB4m.g4DD[mu][k+1]*GRHD.u4U_ito_ValenciavU[mu]*BU[k]
            
    BD = ixp.zerorank1()
    for i in range(3):
        for j in range(3):
            BD[i] += gammaDD[i][j]*BU[j]
            
    for i in range(3):
        # - \frac{\sqrt{\gamma}B_i (u_k B^k)}{4 \pi \alpha} 
        StildeD[i] += -sp.sqrt(gammadet)*BD[i]*udotB/(sqrt4pi*sqrt4pi*alpha)
