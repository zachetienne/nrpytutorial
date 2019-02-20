import indexedexp as ixp
import sympy as sp

'''
Conversion FROM ADM variables gamma_{ij}, beta^i, and alpha TO 4-metric g_{mu nu} and/or inverse 4-metric g^{mu nu} 
'''
def ADM_to_four_metric(gammaDD,betaU,alpha, returng4DD=True, returng4UU=False):
    # The ADM formulation decomposes Einstein's 4D equations into 3+1 dimensional form, with
    #     the 3 spatial dimensions separated from the 1 temporal dimension. DIM here refers to
    #     the spatial dimension.
    DIM = 3

    # Eq 4.47 in [Gourgoulhon](https://arxiv.org/pdf/gr-qc/0703035.pdf):

    # g_{tt} = -\alpha^2 + \beta^k \beta_k
    # g_{ti} = \beta_i
    # g_{ij} = \gamma_{ij}
    
    # Eq. 2.121 in B&S
    betaD = ixp.zerorank1()
    for i in range(DIM):
        for j in range(DIM):
            betaD[i] += gammaDD[i][j] * betaU[j]
    
    # Now compute the beta contraction.
    beta2 = sp.sympify(0)
    for i in range(DIM):
        beta2 += betaU[i] * betaD[i]
    
    g4DD = ixp.zerorank2(DIM=4)
    g4DD[0][0] = -alpha ** 2 + beta2
    for i in range(DIM):
        g4DD[i + 1][0] = g4DD[0][i + 1] = betaD[i]
        for j in range(DIM):
            g4DD[i + 1][j + 1] = gammaDD[i][j]

    if returng4DD == True and returng4UU == False:
        return g4DD

    gammaUU, gammaDET = ixp.symm_matrix_inverter3x3(gammaDD)
    # Eq. 4.49 in [Gourgoulhon](https://arxiv.org/pdf/gr-qc/0703035.pdf):

    # g^{tt} = -1 / alpha^2
    # g^{ti} = beta^i / alpha^2
    # g^{ij} = gamma^{ij} - beta^i beta^j / alpha^2

    g4UU = ixp.zerorank2(DIM=4)
    g4UU[0][0] = - 1 / alpha ** 2
    for i in range(DIM):
        g4UU[i + 1][0] = g4UU[0][i + 1] = betaU[i] / alpha ** 2
        for j in range(DIM):
            g4UU[i + 1][j + 1] = gammaUU[i][j] - betaU[i]*betaU[j] / alpha ** 2

    if returng4DD == True and returng4UU == True:
        return g4DD,g4UU
    if returng4DD == False and returng4UU == True:
        return g4UU
    print("Error: ADM_to_four_metric() called without requesting anything being returned!")
    exit(1)

'''
Conversion FROM 4-metric g_{mu nu} TO ADM variables gamma_{ij}, beta^i, and alpha 
'''
def four_metric_to_ADM(g4DD):
    # The ADM formulation decomposes Einstein's 4D equations into 3+1 dimensional form, with
    #     the 3 spatial dimensions separated from the 1 temporal dimension. DIM here refers to
    #     the spatial dimension.
    DIM = 3

    # Eq 4.47 in [Gourgoulhon](https://arxiv.org/pdf/gr-qc/0703035.pdf):
    # g_{ij} = \gamma_{ij}

    gammaDD = ixp.zerorank2()
    for i in range(DIM):
        for j in range(DIM):
            gammaDD[i][j] = g4DD[i + 1][j + 1]

    # Eq 4.47 in [Gourgoulhon](https://arxiv.org/pdf/gr-qc/0703035.pdf):
    # g_{ti} = \beta_i
    betaD = ixp.zerorank1()
    for i in range(DIM):
        betaD[i] = g4DD[i+1][0]

    #Eq. 2.121 in B&S:
    # beta^i = gamma^{ij} beta_j
    betaU = ixp.zerorank1()
    gammaUU, gammaDET = ixp.symm_matrix_inverter3x3(gammaDD)
    for i in range(DIM):
        for j in range(DIM):
            betaU[i] += gammaUU[i][j]*betaD[j]

    # Eq 4.47 in [Gourgoulhon](https://arxiv.org/pdf/gr-qc/0703035.pdf):
    # g_{tt} = -\alpha^2 + \beta^k \beta_k
    # -> alpha = sqrt(beta^2 - g_tt)
    # Now compute the beta contraction.
    beta2 = sp.sympify(0)
    for i in range(DIM):
        beta2 += betaU[i] * betaD[i]

    alpha = sp.sqrt(beta2 - g4DD[0][0])

    return gammaDD,betaU,alpha