import NRPy_param_funcs as par
import indexedexp as ixp
import grid as gri
import finite_difference as fin
from outputC import *

par.initialize_param(par.glb_param("char", __name__, "System_to_use", "System_II"))

def MaxwellCartesian_Evol():
    #Step 0: Set the spatial dimension parameter to 3.
    par.set_parval_from_str("grid::DIM", 3)
    DIM = par.parval_from_str("grid::DIM")

    # Step 1: Set the finite differencing order to 4.
    par.set_parval_from_str("finite_difference::FD_CENTDERIVS_ORDER", 4)

    # Step 2: Register gridfunctions that are needed as input.
    psi = gri.register_gridfunctions("EVOL", ["psi"])

    # Step 3a: Declare the rank-1 indexed expressions E_{i}, A_{i},
    #          and \partial_{i} \psi. Derivative variables like these
    #          must have an underscore in them, so the finite
    #          difference module can parse the variable name properly.
    ED = ixp.register_gridfunctions_for_single_rank1("EVOL", "ED")
    AD = ixp.register_gridfunctions_for_single_rank1("EVOL", "AD")
    psi_dD = ixp.declarerank1("psi_dD")
    x,y,z = gri.register_gridfunctions("AUX",["x","y","z"])

    ## Step 3b: Declare the conformal metric tensor and its first 
    #           derivative. These are needed to find the Christoffel
    #           symbols, which we need for covariant derivatives.
    gammaDD = ixp.register_gridfunctions_for_single_rank2("AUX","gammaDD", "sym01") # The AUX or EVOL designation is *not*
                                                                                    # used in diagnostic modules.
    gammaDD_dD = ixp.declarerank3("gammaDD_dD","sym01")
    gammaDD_dDD = ixp.declarerank4("gammaDD_dDD","sym01_sym23")

    gammaUU = ixp.declarerank3("gammaUU","sym01")
    detgamma = gri.register_gridfunctions("AUX",["detgamma"])
    gammaUU, detgamma = ixp.symm_matrix_inverter3x3(gammaDD)
    gammaUU_dD = ixp.declarerank3("gammaDD_dD","sym01")

    # Define the Christoffel symbols
    GammaUDD = ixp.zerorank3(DIM)
    for i in range(DIM):
        for k in range(DIM):
            for l in range(DIM):
                for m in range(DIM):
                    GammaUDD[i][k][l] += (sp.Rational(1,2))*gammaUU[i][m]*\
                                         (gammaDD_dD[m][k][l] + gammaDD_dD[m][l][k] - gammaDD_dD[k][l][m])

    # Step 3b: Declare the rank-2 indexed expression \partial_{j} A_{i},
    #          which is not symmetric in its indices.
    #          Derivative variables like these must have an underscore
    #          in them, so the finite difference module can parse the
    #          variable name properly.
    AD_dD = ixp.declarerank2("AD_dD", "nosym")

    # Step 3c: Declare the rank-3 indexed expression \partial_{jk} A_{i},
    #          which is symmetric in the two {jk} indices.
    AD_dDD = ixp.declarerank3("AD_dDD", "sym12")

    # Step 4: Calculate first and second covariant derivatives, and the
    #         necessary contractions.
    # First covariant derivative
    # D_{j} A_{i} = A_{i,j} - \Gamma^{k}_{ij} A_{k}
    AD_dcovD = ixp.zerorank2()
    for i in range(DIM):
        for j in range(DIM):
            AD_dcovD[i][j] = AD_dD[i][j]
            for k in range(DIM):
                AD_dcovD[i][j] -= GammaUDD[k][i][j] * AD[k]

    # First, we must construct the lowered Christoffel symbols:
    # \Gamma_{ijk} = \gamma_{il} \Gamma^l_{jk}
    # And raise the index on A:
    # A^j = \gamma^{ij} A_i
    GammaDDD = ixp.zerorank3()
    AU = ixp.zerorank1()
    for i in range(DIM):
        for j in range(DIM):
            AU[j] = gammaUU[i][j] * AD[i]
            for k in range(DIM):
                for l in range(DIM):
                    GammaDDD[i][j][k] = gammaDD[i][l] * GammaUDD[l][j][k]

    # Covariant second derivative (the bracketed terms):
    # D_j D^j A_i = \gamma^{jk} [A_{i,jk} - A^l (\gamma_{li,kj} + \gamma_{kl,ij} - \gamma_{ik,lj})
    #               + \Gamma_{lik} (\gamma^{lm} A_{m,j} + A_m \gamma^{lm}{}_{,j})
    #               - (\Gamma^l_{ij} A_{l,k} + \Gamma^l_{jk} A_{i,l})
    #               + (\Gamma^m_{ij} \Gamma^l_{mk} A_l + \Gamma ^m_{jk} \Gamma^l_{im} A_l)
    AD_dcovDD = ixp.zerorank3()
    for i in range(DIM):
        for j in range(DIM):
            for k in range(DIM):
                AD_dcovDD[i][j][k] = AD_dDD[i][j][k]
                for l in range(DIM):
                    # Terms 1 and 3
                    AD_dcovDD[i][j][k] -= AU[l] * (gammaDD_dDD[l][i][k][j] + gammaDD_dDD[k][l][i][j] - \
                                                   gammaDD_dDD[i][k][l][j]) \
                                        + GammaUDD[l][i][j] * AD_dD[l][k] + GammaUDD[l][j][k] * AD_dD[i][l]
                    for m in range(DIM):
                        # Terms 2 and 4
                        AD_dcovDD[i][j][k] += GammaDDD[l][i][k] * (gammaUU[l][m] * AD_dD[m][j] + AD[m] * gammaUU_dD[l][m][j]) \
                                            + GammaUDD[m][i][j] * GammaUDD[l][m][k] * AD[l] \
                                            + GammaUDD[m][j][k] * GammaUDD[l][i][m] * AD[l]

    # Covariant divergence
    # D_{i} A^{i} = \gamma^{ij} D_{j} A_{i}
    DivA = 0
    # Gradient of covariant divergence
    # DivA_dD_{i} = \gamma^{jk} A_{k;\hat{j}\hat{i}}
    DivA_dD = ixp.zerorank1()
    # Covariant Laplacian
    # LapAD_{i} = \gamma^{jk} A_{i;\hat{j}\hat{k}}
    LapAD = ixp.zerorank1()
    for i in range(DIM):
        for j in range(DIM):
            DivA += gammaUU[i][j] * AD_dcovD[i][j]
            for k in range(DIM):
                DivA_dD[i] += gammaUU[j][k] * AD_dcovDD[k][j][i]
                LapAD[i]   += gammaUU[j][k] * AD_dcovDD[i][j][k]
    
    global ArhsD, ErhsD, psi_rhs
    system = par.parval_from_str("System_to_use")
    if system is "System_I":
        # Step 5: Define right-hand sides for the evolution.
        print("Warning: System I is less stable!")
        ArhsD = ixp.zerorank1()
        ErhsD = ixp.zerorank1()
        for i in range(DIM):
            ArhsD[i] = -ED[i] - psi_dD[i]
            ErhsD[i] = -LapAD[i] + DivA_dD[i]
        psi_rhs = -DivA

    elif system is "System_II":
        # We inherit here all of the definitions from System I, above

        # Step 7a: Register the scalar auxiliary variable \Gamma
        Gamma = gri.register_gridfunctions("EVOL", ["Gamma"])

        # Step 7b: Declare the ordinary gradient \partial_{i} \Gamma
        Gamma_dD = ixp.declarerank1("Gamma_dD")

        # Step 8a: Construct the second covariant derivative of the scalar \psi
        # \psi_{;\hat{i}\hat{j}} = \psi_{,i;\hat{j}}
        #                        = \psi_{,ij} - \Gamma^{k}_{ij} \psi_{,k}
        psi_dDD = ixp.declarerank2("psi_dDD", "sym01")
        psi_dcovDD = ixp.zerorank2()
        for i in range(DIM):
            for j in range(DIM):
                psi_dcovDD[i][j] = psi_dDD[i][j]
                for k in range(DIM):
                    psi_dcovDD[i][j] += - GammaUDD[k][i][j] * psi_dD[k]

        # Step 8b: Construct the covariant Laplacian of \psi
        # Lappsi = ghat^{ij} D_{j} D_{i} \psi
        Lappsi = 0
        for i in range(DIM):
            for j in range(DIM):
                Lappsi += gammaUU[i][j] * psi_dcovDD[i][j]

        # Step 9: Define right-hand sides for the evolution.
        global Gamma_rhs
        ArhsD = ixp.zerorank1()
        ErhsD = ixp.zerorank1()
        for i in range(DIM):
            ArhsD[i] = -ED[i] - psi_dD[i]
            ErhsD[i] = -LapAD[i] + Gamma_dD[i]
        psi_rhs = -Gamma
        Gamma_rhs = -Lappsi

    else:
        print("Invalid choice of system: System_to_use must be either System_I or System_II")
    
    ED_dD = ixp.declarerank2("ED_dD","nosym")
    global Cviola
    Cviola = gri.register_gridfunctions("AUX", ["Cviola"])
    Cviola = sp.sympify(0)
    for i in range(DIM):
        for j in range(DIM):
            for b in range(DIM):
                Cviola += gammaUU[i][j] * (ED_dD[j][i] - GammaUDD[b][i][j]*ED[b])
