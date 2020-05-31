# This code calculates the Weyl scalars psi0, psi1, psi2, psi3, and psi4. It does so by following the paper
# Baker, Campanelli, and Lousto. PRD 65, 044001 (2002), gr-qc/0104063 and the example set by the Kranc-
# generated ETK thorn which can be found at https://bitbucket.org/einsteintoolkit/einsteinanalysis/src
# Step 1: import all needed modules from NRPy+:
import sympy as sp                # SymPy: The Python computer algebra package upon which NRPy+ depends
import indexedexp as ixp          # NRPy+: Symbolic indexed expression (e.g., tensors, vectors, etc.) support
import grid as gri                # NRPy+: Functions having to do with numerical grids
import NRPy_param_funcs as par    # NRPy+: Parameter interface

# Step 2: Initialize WeylScalars parameters
thismodule = __name__
# Current option: Approx_QuasiKinnersley = choice made in Baker, Campanelli, and Lousto. PRD 65, 044001 (2002)
par.initialize_param(par.glb_param("char *", thismodule, "TetradChoice", "Approx_QuasiKinnersley"))
# This controls what gets output. Acceptable values are "psi4_only", "all_psis", and "all_psis_and_invariants"
par.initialize_param(par.glb_param("char *", thismodule, "output_scalars", "all_psis_and_invariants"))

# Step 3: Define the rank-3 version of the Levi-Civita symbol. Amongst
#         other possible uses, this is needed for the construction of the approximate
#         quasi-Kinnersley tetrad.
def define_LeviCivitaSymbol_rank3(DIM=-1):
    if DIM == -1:
        DIM = par.parval_from_str("DIM")

    LeviCivitaSymbol = ixp.zerorank3()

    for i in range(DIM):
        for j in range(DIM):
            for k in range(DIM):
                # From https://codegolf.stackexchange.com/questions/160359/levi-civita-symbol :
                LeviCivitaSymbol[i][j][k] = (i - j) * (j - k) * (k - i) * sp.Rational(1,2)
    return LeviCivitaSymbol

# Step 1: Call BSSNs. This module computes many different quantities related to the metric,
#         many of which will be essential here. We must first change to our desired coordinate
#         system, however.
def WeylScalars_Cartesian():
    # We do not need the barred or hatted quantities calculated when using Cartesian coordinates.
    # Instead, we declare the PHYSICAL metric and extrinsic curvature as grid functions.
    gammaDD = ixp.register_gridfunctions_for_single_rank2("AUX","gammaDD", "sym01")
    kDD = ixp.register_gridfunctions_for_single_rank2("AUX","kDD", "sym01")
    gammaUU, detgamma = ixp.symm_matrix_inverter3x3(gammaDD)

    global psi4r,psi4i,psi3r,psi3i,psi2r,psi2i,psi1r,psi1i,psi0r,psi0i
    psi4r,psi4i,psi3r,psi3i,psi2r,psi2i,psi1r,psi1i,psi0r,psi0i = gri.register_gridfunctions("AUX",["psi4r","psi4i",
                                                                                                    "psi3r","psi3i",
                                                                                                    "psi2r","psi2i",
                                                                                                    "psi1r","psi1i",
                                                                                                    "psi0r","psi0i"])

    # Step 2a: Set spatial dimension (must be 3 for BSSN)
    DIM = 3
    par.set_parval_from_str("grid::DIM",DIM)

    # Step 2b: Set the coordinate system to Cartesian
    x,y,z = gri.register_gridfunctions("AUX",["x","y","z"])

    # Step 2c: Set which tetrad is used; at the moment, only one supported option
    if par.parval_from_str("WeylScal4NRPy.WeylScalars_Cartesian::TetradChoice") == "Approx_QuasiKinnersley":
        # Step 3a: Choose 3 orthogonal vectors. Here, we choose one in the azimuthal
        #          direction, one in the radial direction, and the cross product of the two.
        # Eqs 5.6, 5.7 in https://arxiv.org/pdf/gr-qc/0104063.pdf:
        # v_1^a &= [-y,x,0] \\
        # v_2^a &= [x,y,z] \\
        # v_3^a &= {\rm det}(g)^{1/2} g^{ad} \epsilon_{dbc} v_1^b v_2^c,
        v1U = ixp.zerorank1()
        v2U = ixp.zerorank1()
        v3U = ixp.zerorank1()
        v1U[0] = -y
        v1U[1] = x
        v1U[2] = sp.sympify(0)
        v2U[0] = x
        v2U[1] = y
        v2U[2] = z
        LeviCivitaSymbol_rank3 = define_LeviCivitaSymbol_rank3()
        for a in range(DIM):
            for b in range(DIM):
                for c in range(DIM):
                    for d in range(DIM):
                        v3U[a] += sp.sqrt(detgamma) * gammaUU[a][d] * LeviCivitaSymbol_rank3[d][b][c] * v1U[b] *v2U[c]

        # Step 3b: Gram-Schmidt orthonormalization of the vectors.
        # The w_i^a vectors here are used to temporarily hold values on the way to the final vectors e_i^a
        # e_1^a &= \frac{v_1^a}{\omega_{11}} \\
        # e_2^a &= \frac{v_2^a - \omega_{12} e_1^a}{\omega_{22}} \\
        # e_3^a &= \frac{v_3^a - \omega_{13} e_1^a - \omega_{23} e_2^a}{\omega_{33}}, \\

        # Normalize the first vector
        w1U = ixp.zerorank1()
        for a in range(DIM):
            w1U[a] = v1U[a]
        omega11 = sp.sympify(0)
        for a in range(DIM):
            for b in range(DIM):
                omega11 += w1U[a] * w1U[b] * gammaDD[a][b]
        e1U = ixp.zerorank1()
        for a in range(DIM):
            e1U[a] = w1U[a] / sp.sqrt(omega11)

        # Subtract off the portion of the first vector along the second, then normalize
        omega12 = sp.sympify(0)
        for a in range(DIM):
            for b in range(DIM):
                omega12 += e1U[a] * v2U[b] * gammaDD[a][b]
        w2U = ixp.zerorank1()
        for a in range(DIM):
            w2U[a] = v2U[a] - omega12*e1U[a]
        omega22 = sp.sympify(0)
        for a in range(DIM):
            for b in range(DIM):
                omega22 += w2U[a] * w2U[b] *gammaDD[a][b]
        e2U = ixp.zerorank1()
        for a in range(DIM):
            e2U[a] = w2U[a] / sp.sqrt(omega22)

        # Subtract off the portion of the first and second vectors along the third, then normalize
        omega13 = sp.sympify(0)
        for a in range(DIM):
            for b in range(DIM):
                omega13 += e1U[a] * v3U[b] * gammaDD[a][b]
        omega23 = sp.sympify(0)
        for a in range(DIM):
            for b in range(DIM):
                omega23 += e2U[a] * v3U[b] * gammaDD[a][b]
        w3U = ixp.zerorank1()
        for a in range(DIM):
            w3U[a] = v3U[a] - omega13*e1U[a] - omega23*e2U[a]
        omega33 = sp.sympify(0)
        for a in range(DIM):
            for b in range(DIM):
                omega33 += w3U[a] * w3U[b] * gammaDD[a][b]
        e3U = ixp.zerorank1()
        for a in range(DIM):
            e3U[a] = w3U[a] / sp.sqrt(omega33)

        # Step 3c: Construct the tetrad itself.
        # Eqs. 5.6:
        # l^a &= \frac{1}{\sqrt{2}} e_2^a \\
        # n^a &= -\frac{1}{\sqrt{2}} e_2^a \\
        # m^a &= \frac{1}{\sqrt{2}} (e_3^a + i e_1^a) \\
        # \overset{*}{m}{}^a &= \frac{1}{\sqrt{2}} (e_3^a - i e_1^a)
        isqrt2 = 1/sp.sqrt(2)
        ltetU = ixp.zerorank1()
        ntetU = ixp.zerorank1()
        #mtetU = ixp.zerorank1()
        #mtetccU = ixp.zerorank1()
        remtetU = ixp.zerorank1() # SymPy does not like trying to take the real/imaginary parts of such a
        immtetU = ixp.zerorank1() # complicated expression as the Weyl scalars, so we will do it ourselves.
        for i in range(DIM):
            ltetU[i] = isqrt2 * e2U[i]
            ntetU[i] = -isqrt2 * e2U[i]
            remtetU[i] = isqrt2 * e3U[i]
            immtetU[i] = isqrt2 * e1U[i]
        nn = isqrt2

    else:
        print("Error: TetradChoice == "+par.parval_from_str("TetradChoice")+" unsupported!")
        exit(1)

    gammaDD_dD = ixp.declarerank3("gammaDD_dD","sym01")

    # Define the Christoffel symbols
    GammaUDD = ixp.zerorank3(DIM)
    for i in range(DIM):
        for k in range(DIM):
            for l in range(DIM):
                for m in range(DIM):
                    GammaUDD[i][k][l] += (sp.Rational(1,2))*gammaUU[i][m]*\
                                         (gammaDD_dD[m][k][l] + gammaDD_dD[m][l][k] - gammaDD_dD[k][l][m])


    # Step 4b: Declare and construct the Riemann curvature tensor:
    # R_{abcd} = \frac{1}{2} (\gamma_{ad,cb}+\gamma_{bc,da}-\gamma_{ac,bd}-\gamma_{bd,ac})
    #            + \gamma_{je} \Gamma^{j}_{bc}\Gamma^{e}_{ad} - \gamma_{je} \Gamma^{j}_{bd} \Gamma^{e}_{ac}
    gammaDD_dDD = ixp.declarerank4("gammaDD_dDD","sym01_sym23")
    RiemannDDDD = ixp.zerorank4()
    for a in range(DIM):
        for b in range(DIM):
            for c in range(DIM):
                for d in range(DIM):
                    RiemannDDDD[a][b][c][d] = (gammaDD_dDD[a][d][c][b] +
                                               gammaDD_dDD[b][c][d][a] -
                                               gammaDD_dDD[a][c][b][d] -
                                               gammaDD_dDD[b][d][a][c]) / 2
                    for e in range(DIM):
                        for j in range(DIM):
                            RiemannDDDD[a][b][c][d] +=  gammaDD[j][e] * GammaUDD[j][b][c] * GammaUDD[e][a][d] - \
                                                        gammaDD[j][e] * GammaUDD[j][b][d] * GammaUDD[e][a][c]


    # Step 4c: We also need the extrinsic curvature tensor $K_{ij}$.
    # In Cartesian coordinates, we already made the components gridfunctions.
    # We will, however, need to calculate the trace of K seperately:
    trK = sp.sympify(0)
    for i in range(DIM):
        for j in range(DIM):
            trK += gammaUU[i][j] * kDD[i][j]

    # Step 5: Build the formula for \psi_4.
    # Gauss equation: involving the Riemann tensor and extrinsic curvature.
    # GaussDDDD[i][j][k][l] =& R_{ijkl} + 2K_{i[k}K_{l]j}
    GaussDDDD = ixp.zerorank4()
    for i in range(DIM):
        for j in range(DIM):
            for k in range(DIM):
                for l in range(DIM):
                    GaussDDDD[i][j][k][l] = RiemannDDDD[i][j][k][l] + kDD[i][k]*kDD[l][j] - kDD[i][l]*kDD[k][j]

    # Codazzi equation: involving partial derivatives of the extrinsic curvature.
    # We will first need to declare derivatives of kDD
    # CodazziDDD[j][k][l] =& -2 (K_{j[k,l]} + \Gamma^p_{j[k} K_{l]p})
    kDD_dD = ixp.declarerank3("kDD_dD","sym01")
    CodazziDDD = ixp.zerorank3()
    for j in range(DIM):
        for k in range(DIM):
            for l in range(DIM):
                CodazziDDD[j][k][l] = kDD_dD[j][l][k] - kDD_dD[j][k][l]
                for p in range(DIM):
                    CodazziDDD[j][k][l] += GammaUDD[p][j][l]*kDD[k][p] - GammaUDD[p][j][k]*kDD[l][p]

    # Another piece. While not associated with any particular equation,
    # this is still useful for organizational purposes.
    # RojoDD[j][l]} = & R_{jl} - K_{jp} K^p_l + KK_{jl} \\
    #               = & \gamma^{pd} R_{jpld} - K_{jp} K^p_l + KK_{jl}
    RojoDD = ixp.zerorank2()
    for j in range(DIM):
        for l in range(DIM):
            RojoDD[j][l] = trK*kDD[j][l]
            for p in range(DIM):
                for d in range(DIM):
                    RojoDD[j][l] += gammaUU[p][d]*RiemannDDDD[j][p][l][d] - kDD[j][p]*gammaUU[p][d]*kDD[d][l]

    # Now we can calculate $\psi_4$ itself! We assume l^0 = n^0 = \frac{1}{\sqrt{2}}
    # and m^0 = \overset{*}{m}{}^0 = 0 to simplify these equations.
    # We calculate the Weyl scalars as defined in https://arxiv.org/abs/gr-qc/0104063
    # In terms of the above-defined quantites, the psis are defined as:
    # \psi_4 =&\ (\text{GaussDDDD[i][j][k][l]}) n^i \overset{*}{m}{}^j n^k \overset{*}{m}{}^l \\
    #         &+2 (\text{CodazziDDD[j][k][l]}) n^{0} \overset{*}{m}{}^{j} n^k \overset{*}{m}{}^l \\
    #         &+ (\text{RojoDD[j][l]}) n^{0} \overset{*}{m}{}^{j} n^{0} \overset{*}{m}{}^{l}.
    # \psi_3 =&\ (\text{GaussDDDD[i][j][k][l]}) l^i n^j \overset{*}{m}{}^k n^l \\
    #         &+ (\text{CodazziDDD[j][k][l]}) (l^{0} n^{j} \overset{*}{m}{}^k n^l - l^{j} n^{0} \overset{*}{m}{}^k n^l - l^k n^j\overset{*}{m}{}^l n^0) \\
    #         &- (\text{RojoDD[j][l]}) l^{0} n^{j} \overset{*}{m}{}^l n^0 - l^{j} n^{0} \overset{*}{m}{}^l n^0 \\
    # \psi_2 =&\ (\text{GaussDDDD[i][j][k][l]}) l^i m^j \overset{*}{m}{}^k n^l \\
    #         &+ (\text{CodazziDDD[j][k][l]}) (l^{0} m^{j} \overset{*}{m}{}^k n^l - l^{j} m^{0} \overset{*}{m}{}^k n^l - l^k m^l \overset{*}{m}{}^l n^0) \\
    #         &- (\text{RojoDD[j][l]}) l^0 m^j \overset{*}{m}{}^l n^0 \\
    # \psi_1 =&\ (\text{GaussDDDD[i][j][k][l]}) n^i l^j m^k l^l \\
    #         &+ (\text{CodazziDDD[j][k][l]}) (n^{0} l^{j} m^k l^l - n^{j} l^{0} m^k l^l - n^k l^l m^j l^0) \\
    #         &- (\text{RojoDD[j][l]}) (n^{0} l^{j} m^l l^0 - n^{j} l^{0} m^l l^0) \\
    # \psi_0 =&\ (\text{GaussDDDD[i][j][k][l]}) l^i m^j l^k m^l \\
    #         &+2 (\text{CodazziDDD[j][k][l]}) (l^0 m^j l^k m^l + l^k m^l l^0 m^j) \\
    #         &+ (\text{RojoDD[j][l]}) l^0 m^j l^0 m^j. \\

    psi4r = sp.sympify(0)
    psi4i = sp.sympify(0)
    psi3r = sp.sympify(0)
    psi3i = sp.sympify(0)
    psi2r = sp.sympify(0)
    psi2i = sp.sympify(0)
    psi1r = sp.sympify(0)
    psi1i = sp.sympify(0)
    psi0r = sp.sympify(0)
    psi0i = sp.sympify(0)
    for l in range(DIM):
        for j in range(DIM):
            psi4r += RojoDD[j][l] * nn * nn * (remtetU[j]*remtetU[l]-immtetU[j]*immtetU[l])
            psi4i += RojoDD[j][l] * nn * nn * (-remtetU[j]*immtetU[l]-immtetU[j]*remtetU[l])
            psi3r +=-RojoDD[j][l] * nn * nn * (ntetU[j]-ltetU[j]) * remtetU[l]
            psi3i += RojoDD[j][l] * nn * nn * (ntetU[j]-ltetU[j]) * immtetU[l]
            psi2r +=-RojoDD[j][l] * nn * nn * (remtetU[l]*remtetU[j]+immtetU[j]*immtetU[l])
            psi2i +=-RojoDD[j][l] * nn * nn * (immtetU[l]*remtetU[j]-remtetU[j]*immtetU[l])
            psi1r += RojoDD[j][l] * nn * nn * (ntetU[j]*remtetU[l]-ltetU[j]*remtetU[l])
            psi1i += RojoDD[j][l] * nn * nn * (ntetU[j]*immtetU[l]-ltetU[j]*immtetU[l])
            psi0r += RojoDD[j][l] * nn * nn * (remtetU[j]*remtetU[l]-immtetU[j]*immtetU[l])
            psi0i += RojoDD[j][l] * nn * nn * (remtetU[j]*immtetU[l]+immtetU[j]*remtetU[l])

    for l in range(DIM):
        for j in range(DIM):
            for k in range(DIM):
                psi4r += 2 * CodazziDDD[j][k][l] * ntetU[k] * nn * (remtetU[j]*remtetU[l]-immtetU[j]*immtetU[l])
                psi4i += 2 * CodazziDDD[j][k][l] * ntetU[k] * nn * (-remtetU[j]*immtetU[l]-immtetU[j]*remtetU[l])
                psi3r += 1 * CodazziDDD[j][k][l] * nn * ((ntetU[j]-ltetU[j])*remtetU[k]*ntetU[l]-remtetU[j]*ltetU[k]*ntetU[l])
                psi3i +=-1 * CodazziDDD[j][k][l] * nn * ((ntetU[j]-ltetU[j])*immtetU[k]*ntetU[l]-immtetU[j]*ltetU[k]*ntetU[l])
                psi2r += 1 * CodazziDDD[j][k][l] * nn * (ntetU[l]*(remtetU[j]*remtetU[k]+immtetU[j]*immtetU[k])-ltetU[k]*(remtetU[j]*remtetU[l]+immtetU[j]*immtetU[l]))
                psi2i += 1 * CodazziDDD[j][k][l] * nn * (ntetU[l]*(immtetU[j]*remtetU[k]-remtetU[j]*immtetU[k])-ltetU[k]*(remtetU[j]*immtetU[l]-immtetU[j]*remtetU[l]))
                psi1r += 1 * CodazziDDD[j][k][l] * nn * (ltetU[j]*remtetU[k]*ltetU[l]-remtetU[j]*ntetU[k]*ltetU[l]-ntetU[j]*remtetU[k]*ltetU[l])
                psi1i += 1 * CodazziDDD[j][k][l] * nn * (ltetU[j]*immtetU[k]*ltetU[l]-immtetU[j]*ntetU[k]*ltetU[l]-ntetU[j]*immtetU[k]*ltetU[l])
                psi0r += 2 * CodazziDDD[j][k][l] * nn * ltetU[k]*(remtetU[j]*remtetU[l]-immtetU[j]*immtetU[l])
                psi0i += 2 * CodazziDDD[j][k][l] * nn * ltetU[k]*(remtetU[j]*immtetU[l]+immtetU[j]*remtetU[l])


    for l in range(DIM):
        for j in range(DIM):
            for k in range(DIM):
                for i in range(DIM):
                    psi4r += GaussDDDD[i][j][k][l] * ntetU[i] * ntetU[k] * (remtetU[j]*remtetU[l]-immtetU[j]*immtetU[l])
                    psi4i += GaussDDDD[i][j][k][l] * ntetU[i] * ntetU[k] * (-remtetU[j]*immtetU[l]-immtetU[j]*remtetU[l])
                    psi3r += GaussDDDD[i][j][k][l] * ltetU[i] * ntetU[j] * remtetU[k] * ntetU[l]
                    psi3i +=-GaussDDDD[i][j][k][l] * ltetU[i] * ntetU[j] * immtetU[k] * ntetU[l]
                    psi2r += GaussDDDD[i][j][k][l] * ltetU[i] * ntetU[l] * (remtetU[j]*remtetU[k]+immtetU[j]*immtetU[k])
                    psi2i += GaussDDDD[i][j][k][l] * ltetU[i] * ntetU[l] * (immtetU[j]*remtetU[k]-remtetU[j]*immtetU[k])
                    psi1r += GaussDDDD[i][j][k][l] * ntetU[i] * ltetU[j] * remtetU[k] * ltetU[l]
                    psi1i += GaussDDDD[i][j][k][l] * ntetU[i] * ltetU[j] * immtetU[k] * ltetU[l]
                    psi0r += GaussDDDD[i][j][k][l] * ltetU[i] * ltetU[k] * (remtetU[j]*remtetU[l]-immtetU[j]*immtetU[l])
                    psi0i += GaussDDDD[i][j][k][l] * ltetU[i] * ltetU[k] * (remtetU[j]*immtetU[l]+immtetU[j]*remtetU[l])

