# Import needed core NRPy+ modules
import grid as gri               # NRPy+: Functions having to do with numerical grids
import NRPy_param_funcs as par   # NRPy+: Parameter interface
import reference_metric as rfm   # NRPy+: Reference metric support
import sympy as sp               # SymPy: The Python computer algebra package upon which NRPy+ depends
import indexedexp as ixp         # NRPy+: Symbolic indexed expression (e.g., tensors, vectors, etc.) support


thismodule = __name__ # "VacuumMaxwell_Flat_Evol"

def VacuumMaxwellRHSs():
    system = par.parval_from_str("Maxwell.InitialData::System_to_use")
    global ArhsU, ErhsU, C, psi_rhs

    #Step 0: Set the spatial dimension parameter to 3.
    DIM = par.parval_from_str("grid::DIM")

    rfm.reference_metric()

        # Register gridfunctions that are needed as input.

    #          Declare the rank-1 indexed expressions E_{i}, A_{i},
    #          that are to be evolved in time.
    #          Derivative variables like these must have an underscore
    #          in them, so the finite difference module can parse
    #          the variable name properly.

    # E^i
    EU = ixp.register_gridfunctions_for_single_rank1("EVOL", "EU")

    # A^i, _AU is unused
    _AU = ixp.register_gridfunctions_for_single_rank1("EVOL", "AU")

    # \psi is a scalar function that is time evolved
    # _psi is unused
    _psi = gri.register_gridfunctions("EVOL", ["psi"])

    # \partial_i \psi
    psi_dD = ixp.declarerank1("psi_dD")

    # \partial_k ( A^i ) --> rank two tensor
    AU_dD = ixp.declarerank2("AU_dD", "nosym")

    # \partial_k partial_m ( A^i ) --> rank three tensor
    AU_dDD = ixp.declarerank3("AU_dDD", "sym12")

    EU_dD = ixp.declarerank2("EU_dD","nosym")

    C = gri.register_gridfunctions("AUX", "DivE")

    # Equation 12 of https://arxiv.org/abs/gr-qc/0201051
    C = EU_dD[0][0] + EU_dD[1][1] + EU_dD[2][2]

    if system == "System_I":
        print('Currently using ' + system + ' RHSs \n')

        # Define right-hand sides for the evolution.
        # Equations 10 and 11 from https://arxiv.org/abs/gr-qc/0201051

        # \partial_t A^i = E^i - \partial_i \psi
        ArhsU = ixp.zerorank1()

        # \partial_t E^i = -\partial_j^2 A^i + \partial_j \partial_i A^j
        ErhsU = ixp.zerorank1()

        # Lorenz gauge condition
        # \partial_t \psi = -\partial_i A^i
        psi_rhs = sp.sympify(0)

        for i in range(DIM):
            ArhsU[i] = -EU[i] - psi_dD[i]
            psi_rhs -= AU_dD[i][i]
            for j in range(DIM):
                ErhsU[i] += -AU_dDD[i][j][j] + AU_dDD[j][j][i]

    elif system == "System_II":
        global Gamma_rhs, G
        print('Currently using ' + system + ' RHSs \n')
        # We inherit here all of the definitions from System I, above

        # Register the scalar auxiliary variable \Gamma
        Gamma = gri.register_gridfunctions("EVOL", ["Gamma"])

        # Declare the ordinary gradient \partial_{i} \Gamma
        Gamma_dD = ixp.declarerank1("Gamma_dD")

        # partial_i \partial_j \psi
        psi_dDD = ixp.declarerank2("psi_dDD", "sym01")

        # Lorenz gauge condition
        psi_rhs = -Gamma

        # Define right-hand sides for the evolution.
        # Equations 10 and 14 https://arxiv.org/abs/gr-qc/0201051
        ArhsU = ixp.zerorank1()
        ErhsU = ixp.zerorank1()

        # Equation 13 of https://arxiv.org/abs/gr-qc/0201051
        # G = \Gamma - \partial_i A^i
        G = gri.register_gridfunctions("AUX", ["G"])
        G = Gamma - AU_dD[0][0] - AU_dD[1][1] - AU_dD[2][2]

        # Equation 15 https://arxiv.org/abs/gr-qc/0201051
        # Gamma_rhs = -DivE
        Gamma_rhs = sp.sympify(0)

        for i in range(DIM):
            ArhsU[i] = -EU[i] - psi_dD[i]
            ErhsU[i] = Gamma_dD[i]
            Gamma_rhs -= psi_dDD[i][i]
            for j in range(DIM):
                ErhsU[i] -= AU_dDD[i][j][j]

    else:
        print("Invalid choice of system: System_to_use must be either System_I or System_II")
