# As documented in the NRPy+ tutorial module
#   Tutorial-Psi4.ipynb,
#   this module will construct a generic
#   expression for \psi_4

# Author: Zachariah B. Etienne
#         zachetie **at** gmail **dot* com

# Step 1.a: import all needed modules from NRPy+:
import sympy as sp                # SymPy: The Python computer algebra package upon which NRPy+ depends
import indexedexp as ixp          # NRPy+: Symbolic indexed expression (e.g., tensors, vectors, etc.) support
import reference_metric as rfm    # NRPy+: Reference metric support

def Psi4(specify_tetrad=True):

    global psi4_im_pt, psi4_re_pt

    # Step 1.b: Given the chosen coordinate system, set up
    #           corresponding reference metric and needed
    #           reference metric quantities
    # The following function call sets up the reference metric
    #    and related quantities, including rescaling matrices ReDD,
    #    ReU, and hatted quantities.
    rfm.reference_metric()

    # Step 1.c: Set spatial dimension (must be 3 for BSSN, as BSSN is
    #           a 3+1-dimensional decomposition of the general
    #           relativistic field equations)
    DIM = 3

    # Step 1.d: Import all ADM quantities as written in terms of BSSN quantities
    import BSSN.ADM_in_terms_of_BSSN as AB
    AB.ADM_in_terms_of_BSSN()

    # Step 1.e: Set up tetrad vectors
    if specify_tetrad==True:
        import BSSN.Psi4_tetrads as BP4t
        BP4t.Psi4_tetrads()
        mre4U = BP4t.mre4U
        mim4U = BP4t.mim4U
        n4U   = BP4t.n4U
    else:
        # For code validation against NRPy+ psi_4 tutorial module (Tutorial-Psi4.ipynb);
        #   ensures a more complete code validation.
        mre4U = ixp.declarerank1("mre4U", DIM=4)
        mim4U = ixp.declarerank1("mim4U", DIM=4)
        n4U = ixp.declarerank1("n4U", DIM=4)

    # Step 2: Construct the (rank-4) Riemann curvature tensor associated with the ADM 3-metric:
    RDDDD = ixp.zerorank4()
    gammaDDdDD = AB.gammaDDdDD

    for i in range(DIM):
        for k in range(DIM):
            for l in range(DIM):
                for m in range(DIM):
                    RDDDD[i][k][l][m] = sp.Rational(1, 2) * \
                                        (gammaDDdDD[i][m][k][l] + gammaDDdDD[k][l][i][m] - gammaDDdDD[i][l][k][m] -
                                         gammaDDdDD[k][m][i][l])

    # ... then we add the term on the right:
    gammaDD = AB.gammaDD
    GammaUDD = AB.GammaUDD

    for i in range(DIM):
        for k in range(DIM):
            for l in range(DIM):
                for m in range(DIM):
                    for n in range(DIM):
                        for p in range(DIM):
                            RDDDD[i][k][l][m] += gammaDD[n][p] * \
                                                 (GammaUDD[n][k][l] * GammaUDD[p][i][m] - GammaUDD[n][k][m] * GammaUDD[p][i][l])

    # Step 3: Construct the (rank-4) tensor in term 1 of psi_4 (referring to Eq 5.1 in
    #   Baker, Campanelli, Lousto (2001); https://arxiv.org/pdf/gr-qc/0104063.pdf
    rank4term1DDDD = ixp.zerorank4()
    KDD = AB.KDD

    for i in range(DIM):
        for j in range(DIM):
            for k in range(DIM):
                for l in range(DIM):
                    rank4term1DDDD[i][j][k][l] = RDDDD[i][j][k][l] + KDD[i][k] * KDD[l][j] - KDD[i][l] * KDD[k][j]

    # Step 4: Construct the (rank-3) tensor in term 2 of psi_4 (referring to Eq 5.1 in
    #   Baker, Campanelli, Lousto (2001); https://arxiv.org/pdf/gr-qc/0104063.pdf
    rank3term2DDD = ixp.zerorank3()
    KDDdD = AB.KDDdD

    for j in range(DIM):
        for k in range(DIM):
            for l in range(DIM):
                rank3term2DDD[j][k][l] = sp.Rational(1, 2) * (KDDdD[j][k][l] - KDDdD[j][l][k])

    # ... then we construct the second term in this sum:
    #  \Gamma^{p}_{j[k} K_{l]p} = \frac{1}{2} (\Gamma^{p}_{jk} K_{lp}-\Gamma^{p}_{jl} K_{kp}):
    for j in range(DIM):
        for k in range(DIM):
            for l in range(DIM):
                for p in range(DIM):
                    rank3term2DDD[j][k][l] += sp.Rational(1, 2) * (
                                GammaUDD[p][j][k] * KDD[l][p] - GammaUDD[p][j][l] * KDD[k][p])

    # Finally, we multiply the term by $-8$:
    for j in range(DIM):
        for k in range(DIM):
            for l in range(DIM):
                rank3term2DDD[j][k][l] *= sp.sympify(-8)

    # Step 5: Construct the (rank-2) tensor in term 3 of psi_4 (referring to Eq 5.1 in
    #   Baker, Campanelli, Lousto (2001); https://arxiv.org/pdf/gr-qc/0104063.pdf

    # Step 5.1: Construct 3-Ricci tensor R_{ij} = gamma^{im} R_{ijml}
    RDD = ixp.zerorank2()
    gammaUU = AB.gammaUU
    for j in range(DIM):
        for l in range(DIM):
            for i in range(DIM):
                for m in range(DIM):
                    RDD[j][l] += gammaUU[i][m]*RDDDD[i][j][m][l]

    # Step 5.2: Construct K^p_l = gamma^{pi} K_{il}
    KUD = ixp.zerorank2()
    for p in range(DIM):
        for l in range(DIM):
            for i in range(DIM):
                KUD[p][l] += gammaUU[p][i]*KDD[i][l]

    # Step 5.3: Construct trK = gamma^{ij} K_{ij}
    trK = sp.sympify(0)
    for i in range(DIM):
        for j in range(DIM):
            trK += gammaUU[i][j]*KDD[i][j]

    # Next we put these terms together to construct the entire term in parentheses:
    # +4 \left(R_{jl} - K_{jp} K^p_l + K K_{jl} \right),
    rank2term3DD = ixp.zerorank2()
    for j in range(DIM):
        for l in range(DIM):
            rank2term3DD[j][l] = RDD[j][l] + trK*KDD[j][l]
            for p in range(DIM):
                rank2term3DD[j][l] += - KDD[j][p]*KUD[p][l]
    # Finally we multiply by +4:
    for j in range(DIM):
        for l in range(DIM):
            rank2term3DD[j][l] *= sp.sympify(4)

    # Step 6: Construct real & imaginary parts of psi_4
    #         by contracting constituent rank 2, 3, and 4
    #         tensors with input tetrads mre4U, mim4U, & n4U.

    def tetrad_product__Real_psi4(n, Mre, Mim, mu, nu, eta, delta):
        return +n[mu] * Mre[nu] * n[eta] * Mre[delta] - n[mu] * Mim[nu] * n[eta] * Mim[delta]

    def tetrad_product__Imag_psi4(n, Mre, Mim, mu, nu, eta, delta):
        return -n[mu] * Mre[nu] * n[eta] * Mim[delta] - n[mu] * Mim[nu] * n[eta] * Mre[delta]

    # We split psi_4 into three pieces, to expedite & possibly parallelize C code generation.
    psi4_re_pt = [sp.sympify(0),sp.sympify(0),sp.sympify(0)]
    psi4_im_pt = [sp.sympify(0),sp.sympify(0),sp.sympify(0)]
    # First term:
    for i in range(DIM):
        for j in range(DIM):
            for k in range(DIM):
                for l in range(DIM):
                    psi4_re_pt[0] += rank4term1DDDD[i][j][k][l] * tetrad_product__Real_psi4(n4U, mre4U, mim4U,
                                                                                            i + 1, j + 1,k + 1, l + 1)
                    psi4_im_pt[0] += rank4term1DDDD[i][j][k][l] * tetrad_product__Imag_psi4(n4U, mre4U, mim4U,
                                                                                            i + 1, j + 1, k + 1, l + 1)

    # Second term:
    for j in range(DIM):
        for k in range(DIM):
            for l in range(DIM):
                psi4_re_pt[1] += rank3term2DDD[j][k][l] * \
                           sp.Rational(1, 2) * (+tetrad_product__Real_psi4(n4U, mre4U, mim4U, 0, j + 1, k + 1, l + 1)
                                                - tetrad_product__Real_psi4(n4U, mre4U, mim4U, j + 1, 0, k + 1, l + 1))
                psi4_im_pt[1] += rank3term2DDD[j][k][l] * \
                           sp.Rational(1, 2) * (+tetrad_product__Imag_psi4(n4U, mre4U, mim4U, 0, j + 1, k + 1, l + 1)
                                                - tetrad_product__Imag_psi4(n4U, mre4U, mim4U, j + 1, 0, k + 1, l + 1))
    # Third term:
    for j in range(DIM):
        for l in range(DIM):
            psi4_re_pt[2] += rank2term3DD[j][l] * \
                       (sp.Rational(1, 4) * (+tetrad_product__Real_psi4(n4U, mre4U, mim4U, 0, j + 1, 0, l + 1)
                                             - tetrad_product__Real_psi4(n4U, mre4U, mim4U, j + 1, 0, 0, l + 1)
                                             - tetrad_product__Real_psi4(n4U, mre4U, mim4U, 0, j + 1, l + 1, 0)
                                             + tetrad_product__Real_psi4(n4U, mre4U, mim4U, j + 1, 0, l + 1, 0)))
            psi4_im_pt[2] += rank2term3DD[j][l] * \
                       (sp.Rational(1, 4) * (+tetrad_product__Imag_psi4(n4U, mre4U, mim4U, 0, j + 1, 0, l + 1)
                                             - tetrad_product__Imag_psi4(n4U, mre4U, mim4U, j + 1, 0, 0, l + 1)
                                             - tetrad_product__Imag_psi4(n4U, mre4U, mim4U, 0, j + 1, l + 1, 0)
                                             + tetrad_product__Imag_psi4(n4U, mre4U, mim4U, j + 1, 0, l + 1, 0)))
