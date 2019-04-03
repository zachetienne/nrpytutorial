# As documented in the NRPy+ tutorial module
#   Tutorial-Psi4_tetrads.ipynb,
#   this module will construct a generic
#   expression for \psi_4

# Author: Zachariah B. Etienne
#         zachetie **at** gmail **dot* com

# Step 1.a: import all needed modules from NRPy+:
import sympy as sp
import NRPy_param_funcs as par
import indexedexp as ixp
import grid as gri
import finite_difference as fin
import reference_metric as rfm

# Step 1.b: Initialize TetradChoice parameter
thismodule = __name__
# Current option: QuasiKinnersley = choice made in Baker, Campanelli, and Lousto. PRD 65, 044001 (2002)
par.initialize_param(par.glb_param("char", thismodule, "TetradChoice", "QuasiKinnersley"))

def Psi4_tetrads():
    global l4U, n4U, mre4U, mim4U

    # Step 1.c: Check if tetrad choice is implemented:
    if par.parval_from_str(thismodule+"::TetradChoice") != "QuasiKinnersley":
        print("ERROR: "+thismodule+"::TetradChoice = "+par.parval_from_str("TetradChoice")+" currently unsupported!")
        exit(1)

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

    # Step 2.a: Declare the Cartesian x,y,z as input parameters
    #           and v_1^a, v_2^a, and v_3^a tetrads,
    #           as well as detgamma and gammaUU from
    #           BSSN.ADM_in_terms_of_BSSN
    x, y, z = par.Cparameters("REAL", thismodule, ["x", "y", "z"])

    v1U = ixp.zerorank1()
    v2U = ixp.zerorank1()
    v3U = ixp.zerorank1()

    detgamma = AB.detgamma
    gammaUU = AB.gammaUU

    # Step 2.b: Define the rank-3 version of the Levi-Civita symbol. Amongst
    #         other uses, this is needed for the construction of the approximate
    #         quasi-Kinnersley tetrad.
    def define_LeviCivitaSymbol_rank3(DIM=-1):
        if DIM == -1:
            DIM = par.parval_from_str("DIM")

        LeviCivitaSymbol = ixp.zerorank3()

        for i in range(DIM):
            for j in range(DIM):
                for k in range(DIM):
                    # From https://codegolf.stackexchange.com/questions/160359/levi-civita-symbol :
                    LeviCivitaSymbol[i][j][k] = (i - j) * (j - k) * (k - i) / 2
        return LeviCivitaSymbol

    # Step 2.c: define v1U and v2U
    v1U = [-y, x, sp.sympify(0)]
    v2U = [x, y, z]

    # Step 2.d: define v3U
    LeviCivitaSymbolDDD = define_LeviCivitaSymbol_rank3(DIM=3)
    for a in range(DIM):
        for b in range(DIM):
            for c in range(DIM):
                for d in range(DIM):
                    v3U[a] += sp.sqrt(detgamma) * gammaUU[a][d] * LeviCivitaSymbolDDD[d][b][c] * v1U[b] * v2U[c]

    # Step 2.e: Define omega_{ij}
    omegaDD = ixp.zerorank2()
    gammaDD = AB.gammaDD

    def v_vectorDU(v1U, v2U, v3U, i, a):
        if i == 0:
            return v1U[a]
        elif i == 1:
            return v2U[a]
        elif i == 2:
            return v3U[a]
        else:
            print("ERROR: unknown vector!")
            exit(1)

    for i in range(DIM):
        for j in range(DIM):
            for a in range(DIM):
                for b in range(DIM):
                    omegaDD[i][j] += v_vectorDU(v1U, v2U, v3U, i, a) * v_vectorDU(v1U, v2U, v3U, j, b) * gammaDD[a][b]

    # Step 2.f: Define e^a_i. Note that:
    #           omegaDD[0][0] = \omega_{11} above;
    #           omegaDD[1][1] = \omega_{22} above, etc.
    e1U = ixp.zerorank1()
    e2U = ixp.zerorank1()
    e3U = ixp.zerorank1()
    for a in range(DIM):
        e1U[a] = v1U[a] / omegaDD[0][0]
        e2U[a] = (v2U[a] - omegaDD[0][1] * e1U[a]) / omegaDD[1][1]
        e3U[a] = (v3U[a] - omegaDD[0][2] * e1U[a] - omegaDD[1][2] * e2U[a]) / omegaDD[2][2]

    # Step 2.g: Construct l^a, n^a, and m^a
    isqrt2 = 1 / sp.sqrt(2)
    l4U = ixp.zerorank1(DIM=4)
    n4U = ixp.zerorank1(DIM=4)
    mre4U = ixp.zerorank1(DIM=4)
    mim4U = ixp.zerorank1(DIM=4)
    l4U[0] = isqrt2
    n4U[0] = isqrt2
    mre4U[0] = sp.sympify(0)
    mim4U[0] = sp.sympify(0)
    for a in range(DIM):
        l4U[a + 1] = isqrt2 * e2U[a]
        n4U[a + 1] = -isqrt2 * e2U[a]
        mre4U[a + 1] = isqrt2 * e3U[a]
        mim4U[a + 1] = isqrt2 * e1U[a]
