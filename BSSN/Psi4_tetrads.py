# As documented in the NRPy+ tutorial module
#   Tutorial-Psi4_tetrads.ipynb,
#   this module will construct tetrads
#   needed to compute \psi_4 (as well as other
#   Weyl scalars and invariants in principle)

# Authors: Zachariah B. Etienne
#          (zachetie **at** gmail **dot* com),
#          and Patrick Nelson

# Step 1.a: import all needed modules from NRPy+:
import sympy as sp                # SymPy: The Python computer algebra package upon which NRPy+ depends
import NRPy_param_funcs as par    # NRPy+: Parameter interface
import indexedexp as ixp          # NRPy+: Symbolic indexed expression (e.g., tensors, vectors, etc.) support
import reference_metric as rfm    # NRPy+: Reference metric support
import sys                        # Standard Python modules for multiplatform OS-level functions

# Step 1.b: Initialize TetradChoice parameter
thismodule = __name__
# Current option: QuasiKinnersley = choice made in Baker, Campanelli, and Lousto. PRD 65, 044001 (2002)
par.initialize_param(par.glb_param("char", thismodule, "TetradChoice", "QuasiKinnersley"))
par.initialize_param(par.glb_param("char", thismodule, "UseCorrectUnitNormal", "False")) # False = consistent with WeylScal4 ETK thorn.

def Psi4_tetrads():
    global l4U, n4U, mre4U, mim4U

    # Step 1.c: Check if tetrad choice is implemented:
    if par.parval_from_str(thismodule+"::TetradChoice") != "QuasiKinnersley":
        print("ERROR: "+thismodule+"::TetradChoice = "+par.parval_from_str("TetradChoice")+" currently unsupported!")
        sys.exit(1)

    # Step 1.d: Given the chosen coordinate system, set up
    #           corresponding reference metric and needed
    #           reference metric quantities
    # The following function call sets up the reference metric
    #    and related quantities, including rescaling matrices ReDD,
    #    ReU, and hatted quantities.
    rfm.reference_metric()

    # Step 1.e: Set spatial dimension (must be 3 for BSSN, as BSSN is
    #           a 3+1-dimensional decomposition of the general
    #           relativistic field equations)
    DIM = 3

    # Step 1.f: Import all ADM quantities as written in terms of BSSN quantities
    import BSSN.ADM_in_terms_of_BSSN as AB
    AB.ADM_in_terms_of_BSSN()

    # Step 2.a: Declare the Cartesian x,y,z in terms of
    #           xx0,xx1,xx2.
    x = rfm.xx_to_Cart[0]
    y = rfm.xx_to_Cart[1]
    z = rfm.xx_to_Cart[2]

    # Step 2.b: Declare detgamma and gammaUU from
    #           BSSN.ADM_in_terms_of_BSSN;
    #           simplify detgamma & gammaUU expressions,
    #           which expedites Psi4 codegen.
    detgamma = sp.simplify(AB.detgamma)
    gammaUU = ixp.zerorank2()
    for i in range(DIM):
        for j in range(DIM):
            gammaUU[i][j] = sp.simplify(AB.gammaUU[i][j])

    # Step 2.c: Define v1U and v2U
    v1UCart = [-y, x, sp.sympify(0)]
    v2UCart = [x, y, z]

    # Step 2.d: Construct the Jacobian d x_Cart^i / d xx^j
    Jac_dUCart_dDrfmUD = ixp.zerorank2()
    for i in range(DIM):
        for j in range(DIM):
            Jac_dUCart_dDrfmUD[i][j] = sp.simplify(sp.diff(rfm.xx_to_Cart[i], rfm.xx[j]))

    # Step 2.e: Invert above Jacobian to get needed d xx^j / d x_Cart^i
    Jac_dUrfm_dDCartUD, dummyDET = ixp.generic_matrix_inverter3x3(Jac_dUCart_dDrfmUD)

    # Step 2.e.i: Simplify expressions for d xx^j / d x_Cart^i:
    for i in range(DIM):
        for j in range(DIM):
            Jac_dUrfm_dDCartUD[i][j] = sp.simplify(Jac_dUrfm_dDCartUD[i][j])

    # Step 2.f: Transform v1U and v2U from the Cartesian to the xx^i basis
    v1U = ixp.zerorank1()
    v2U = ixp.zerorank1()
    for i in range(DIM):
        for j in range(DIM):
            v1U[i] += Jac_dUrfm_dDCartUD[i][j] * v1UCart[j]
            v2U[i] += Jac_dUrfm_dDCartUD[i][j] * v2UCart[j]

    # Step 2.g: Define v3U
    v3U = ixp.zerorank1()
    LeviCivitaSymbolDDD = ixp.LeviCivitaSymbol_dim3_rank3()
    for a in range(DIM):
        for b in range(DIM):
            for c in range(DIM):
                for d in range(DIM):
                    v3U[a] += sp.sqrt(detgamma) * gammaUU[a][d] * LeviCivitaSymbolDDD[d][b][c] * v1U[b] * v2U[c]

    # Step 2.g.i: Simplify expressions for v1U,v2U,v3U. This greatly expedites the C code generation (~10x faster)
    #             Drat. Simplification with certain versions of SymPy & coord systems results in a hang. Let's just
    #             evaluate the expressions so the most trivial optimizations can be performed.
    for a in range(DIM):
        v1U[a] = v1U[a].doit()  # sp.simplify(v1U[a])
        v2U[a] = v2U[a].doit()  # sp.simplify(v2U[a])
        v3U[a] = v3U[a].doit()  # sp.simplify(v3U[a])

    # Step 2.h: Define omega_{ij}
    omegaDD = ixp.zerorank2()
    gammaDD = AB.gammaDD
    def v_vectorDU(v1U,v2U,v3U,  i,a):
        if i==0:
            return v1U[a]
        if i==1:
            return v2U[a]
        if i==2:
            return v3U[a]
        print("ERROR: unknown vector!")
        sys.exit(1)

    def update_omega(omegaDD, i,j, v1U,v2U,v3U,gammaDD):
        omegaDD[i][j] = sp.sympify(0)
        for a in range(DIM):
            for b in range(DIM):
                omegaDD[i][j] += v_vectorDU(v1U,v2U,v3U, i,a)*v_vectorDU(v1U,v2U,v3U, j,b)*gammaDD[a][b]

    # Step 2.i: Define e^a_i. Note that:
    #           omegaDD[0][0] = \omega_{11} above;
    #           omegaDD[1][1] = \omega_{22} above, etc.
    # First e_1^a: Orthogonalize & normalize:
    e1U = ixp.zerorank1()
    update_omega(omegaDD, 0,0, v1U,v2U,v3U,gammaDD)
    for a in range(DIM):
        e1U[a] = v1U[a]/sp.sqrt(omegaDD[0][0])

    # Next e_2^a: First orthogonalize:
    e2U = ixp.zerorank1()
    update_omega(omegaDD, 0,1, e1U,v2U,v3U,gammaDD)
    for a in range(DIM):
        e2U[a] = (v2U[a] - omegaDD[0][1]*e1U[a])
    # Then normalize:
    update_omega(omegaDD, 1,1, e1U,e2U,v3U,gammaDD)
    for a in range(DIM):
        e2U[a] /= sp.sqrt(omegaDD[1][1])

    # Next e_3^a: First orthogonalize:
    e3U = ixp.zerorank1()
    update_omega(omegaDD, 0,2, e1U,e2U,v3U,gammaDD)
    update_omega(omegaDD, 1,2, e1U,e2U,v3U,gammaDD)
    for a in range(DIM):
        e3U[a] = (v3U[a] - omegaDD[0][2]*e1U[a] - omegaDD[1][2]*e2U[a])
    # Then normalize:
    update_omega(omegaDD, 2,2, e1U,e2U,e3U,gammaDD)
    for a in range(DIM):
        e3U[a] /= sp.sqrt(omegaDD[2][2])

    # Step 2.j: Construct l^mu, n^mu, and m^mu, based on r^mu, theta^mu, phi^mu, and u^mu:
    r4U     = ixp.zerorank1(DIM=4)
    u4U     = ixp.zerorank1(DIM=4)
    theta4U = ixp.zerorank1(DIM=4)
    phi4U   = ixp.zerorank1(DIM=4)

    for a in range(DIM):
        r4U[    a+1] = e2U[a]
        theta4U[a+1] = e3U[a]
        phi4U[  a+1] = e1U[a]

    # FIXME? assumes alpha=1, beta^i = 0
    if par.parval_from_str(thismodule+"::UseCorrectUnitNormal") == "False":
        u4U[0] = 1
    else:
        # Eq. 2.116 in Baumgarte & Shapiro:
        #  n^mu = {1/alpha, -beta^i/alpha}. Note that n_mu = {alpha,0}, so n^mu n_mu = -1.
        import BSSN.BSSN_quantities as Bq
        Bq.declare_BSSN_gridfunctions_if_not_declared_already()
        Bq.BSSN_basic_tensors()
        u4U[0] = 1/Bq.alpha
        for i in range(DIM):
            u4U[i+1] = -Bq.betaU[i]/Bq.alpha

    l4U = ixp.zerorank1(DIM=4)
    n4U = ixp.zerorank1(DIM=4)
    mre4U  = ixp.zerorank1(DIM=4)
    mim4U  = ixp.zerorank1(DIM=4)

    # M_SQRT1_2 = 1 / sqrt(2) (defined in math.h on Linux)
    M_SQRT1_2 = par.Cparameters("#define",thismodule,"M_SQRT1_2","")
    isqrt2 = M_SQRT1_2 #1/sp.sqrt(2) <- SymPy drops precision to 15 sig. digits in unit tests
    for mu in range(4):
        l4U[mu]   = isqrt2*(u4U[mu] + r4U[mu])
        n4U[mu]   = isqrt2*(u4U[mu] - r4U[mu])
        mre4U[mu] = isqrt2*theta4U[mu]
        mim4U[mu] = isqrt2*  phi4U[mu]
