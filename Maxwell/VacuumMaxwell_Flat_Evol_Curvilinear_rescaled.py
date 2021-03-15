# Import needed core NRPy+ modules
import grid as gri               # NRPy+: Functions having to do with numerical grids
import NRPy_param_funcs as par   # NRPy+: Parameter interface
import reference_metric as rfm   # NRPy+: Reference metric support
import sympy as sp               # SymPy: The Python computer algebra package upon which NRPy+ depends
import indexedexp as ixp         # NRPy+: Symbolic indexed expression (e.g., tensors, vectors, etc.) support


thismodule = __name__ # "VacuumMaxwell_Flat_Evol-Curvilinear"

def VacuumMaxwellRHSs_rescaled():
    global erhsU, arhsU, psi_rhs, Gamma_rhs, C, G, EU_Cart, AU_Cart

    #Step 0: Set the spatial dimension parameter to 3.
    par.set_parval_from_str("grid::DIM", 3)
    DIM = par.parval_from_str("grid::DIM")

    # Set reference metric related quantities
    rfm.reference_metric()

    # Register gridfunctions that are needed as input.

    #  Declare the rank-1 indexed expressions E^{i}, A^{i},
    #  and \partial^{i} \psi, that are to be evolved in time.
    #  Derivative variables like these must have an underscore
    #  in them, so the finite difference module can parse
    #  the variable name properly.

    # e^i
    eU = ixp.register_gridfunctions_for_single_rank1("EVOL", "eU")

    # \partial_k ( E^i ) --> rank two tensor
    eU_dD = ixp.declarerank2("eU_dD", "nosym")

    # a^i
    aU = ixp.register_gridfunctions_for_single_rank1("EVOL", "aU")

    # \partial_k ( a^i ) --> rank two tensor
    aU_dD = ixp.declarerank2("aU_dD", "nosym")

    # \partial_k partial_m ( a^i ) --> rank three tensor
    aU_dDD = ixp.declarerank3("aU_dDD", "sym12")

    # \psi is a scalar function that is time evolved
    # psi is unused here
    _psi = gri.register_gridfunctions("EVOL", ["psi"])

    # \Gamma is a scalar function that is time evolved
    Gamma = gri.register_gridfunctions("EVOL", ["Gamma"])

    # \partial_i \psi
    psi_dD = ixp.declarerank1("psi_dD")

    # \partial_i \Gamma
    Gamma_dD = ixp.declarerank1("Gamma_dD")

    # partial_i \partial_j \psi
    psi_dDD = ixp.declarerank2("psi_dDD", "sym01")

    ghatUU = rfm.ghatUU
    GammahatUDD = rfm.GammahatUDD
    GammahatUDDdD = rfm.GammahatUDDdD
    ReU = rfm.ReU
    ReUdD = rfm.ReUdD
    ReUdDD = rfm.ReUdDD

    # \partial_t a^i = -e^i - \frac{\hat{g}^{ij}\partial_j \varphi}{\text{ReU}[i]}
    arhsU = ixp.zerorank1()
    for i in range(DIM):
        arhsU[i] -= eU[i]
        for j in range(DIM):
            arhsU[i] -= (ghatUU[i][j]*psi_dD[j])/ReU[i]

    # A^i

    AU = ixp.zerorank1()

    # \partial_k ( A^i ) --> rank two tensor
    AU_dD = ixp.zerorank2()

    # \partial_k partial_m ( A^i ) --> rank three tensor
    AU_dDD = ixp.zerorank3()

    for i in range(DIM):
        AU[i] = aU[i]*ReU[i]
        for j in range(DIM):
            AU_dD[i][j] = aU_dD[i][j]*ReU[i] + aU[i]*ReUdD[i][j]
            for k in range(DIM):
                AU_dDD[i][j][k] = aU_dDD[i][j][k]*ReU[i] + aU_dD[i][j]*ReUdD[i][k] +\
                                  aU_dD[i][k]*ReUdD[i][j] + aU[i]*ReUdDD[i][j][k]


    # Term 1 = \hat{g}^{ij}\partial_j \Gamma
    Term1U = ixp.zerorank1()
    for i in range(DIM):
        for j in range(DIM):
            Term1U[i] += ghatUU[i][j]*Gamma_dD[j]

    # Term 2: A^i_{,kj}
    Term2UDD = ixp.zerorank3()
    for i in range(DIM):
        for j in range(DIM):
            for k in range(DIM):
                Term2UDD[i][j][k] += AU_dDD[i][k][j]

    # Term 3: \hat{\Gamma}^i_{mk,j} A^m + \hat{\Gamma}^i_{mk} A^m_{,j}
    #        + \hat{\Gamma}^i_{dj}A^d_{,k} - \hat{\Gamma}^d_{kj} A^i_{,d}
    Term3UDD = ixp.zerorank3()
    for i in range(DIM):
        for j in range(DIM):
            for k in range(DIM):
                for m in range(DIM):
                    Term3UDD[i][j][k] += GammahatUDDdD[i][m][k][j]*AU[m]    \
                                          + GammahatUDD[i][m][k]*AU_dD[m][j] \
                                          + GammahatUDD[i][m][j]*AU_dD[m][k] \
                                          - GammahatUDD[m][k][j]*AU_dD[i][m]
    # Term 4: \hat{\Gamma}^i_{dj}\hat{\Gamma}^d_{mk} A^m -
    #         \hat{\Gamma}^d_{kj} \hat{\Gamma}^i_{md} A^m
    Term4UDD = ixp.zerorank3()
    for i in range(DIM):
        for j in range(DIM):
            for k in range(DIM):
                for m in range(DIM):
                    for d in range(DIM):
                        Term4UDD[i][j][k] += ( GammahatUDD[i][d][j]*GammahatUDD[d][m][k] \
                                               -GammahatUDD[d][k][j]*GammahatUDD[i][m][d])*AU[m]

    # \partial_t E^i = \hat{g}^{ij}\partial_j \Gamma - \hat{\gamma}^{jk}*
    #    (A^i_{,kj}
    #   + \hat{\Gamma}^i_{mk,j} A^m + \hat{\Gamma}^i_{mk} A^m_{,j}
    #   + \hat{\Gamma}^i_{dj} A^d_{,k} - \hat{\Gamma}^d_{kj} A^i_{,d}
    #   + \hat{\Gamma}^i_{dj}\hat{\Gamma}^d_{mk} A^m
    #   - \hat{\Gamma}^d_{kj} \hat{\Gamma}^i_{md} A^m)

    ErhsU = ixp.zerorank1()
    for i in range(DIM):
        ErhsU[i] += Term1U[i]
        for j in range(DIM):
            for k in range(DIM):
                ErhsU[i] -= ghatUU[j][k]*(Term2UDD[i][j][k] + Term3UDD[i][j][k] + Term4UDD[i][j][k])

    erhsU = ixp.zerorank1()
    for i in range(DIM):
        erhsU[i] = ErhsU[i]/ReU[i]

    # \partial_t \Gamma = -\hat{g}^{ij} (\partial_i \partial_j \varphi -
    # \hat{\Gamma}^k_{ji} \partial_k \varphi)
    Gamma_rhs = sp.sympify(0)
    for i in range(DIM):
        for j in range(DIM):
            Gamma_rhs -= ghatUU[i][j]*psi_dDD[i][j]
            for k in range(DIM):
                Gamma_rhs += ghatUU[i][j]*GammahatUDD[k][j][i]*psi_dD[k]

    # \partial_t \varphi = -\Gamma
    psi_rhs = -Gamma

    # \mathcal{G} \equiv \Gamma - \partial_i A^i + \hat{\Gamma}^i_{ji} A^j
    G = Gamma
    for i in range(DIM):
        G -= AU_dD[i][i]
        for j in range(DIM):
            G += GammahatUDD[i][j][i]*AU[j]

    # E^i
    EU = ixp.zerorank1()

    EU_dD = ixp.zerorank2()
    for i in range(DIM):
        EU[i] = eU[i]*ReU[i]
        for j in range(DIM):
            EU_dD[i][j] = eU_dD[i][j]*ReU[i] + eU[i]*ReUdD[i][j]

    C = sp.sympify(0)
    for i in range(DIM):
        C += EU_dD[i][i]
        for j in range(DIM):
            C += GammahatUDD[i][j][i]*EU[j]

    def Convert_to_Cartesian_basis(VU):
        # Coordinate transformation from original basis to Cartesian
        rfm.reference_metric()

        VU_Cart = ixp.zerorank1()
        Jac_dxCartU_dxOrigD = ixp.zerorank2()
        for i in range(DIM):
            for j in range(DIM):
                Jac_dxCartU_dxOrigD[i][j] = sp.diff(rfm.xx_to_Cart[i], rfm.xx[j])

        for i in range(DIM):
            for j in range(DIM):
                VU_Cart[i] += Jac_dxCartU_dxOrigD[i][j]*VU[j]
        return VU_Cart

    AU_Cart = Convert_to_Cartesian_basis(AU)
    EU_Cart = Convert_to_Cartesian_basis(EU)
