# This module computes the the Weyl scalars based on
# Baker, Campanelli, and Lousto. PRD 65, 044001 (2002);
#  https://arxiv.org/abs/gr-qc/0104063

# Step 1: import all needed modules from NRPy+:
import NRPy_param_funcs as par
import indexedexp as ixp
import grid as gri
import finite_difference as fin
import reference_metric as rfm
rfm.reference_metric()
from outputC import *
import BSSN_RHSs_new as bssn
import sympy as sp

# Step 1: Initialize WeylScalar parameters
thismodule = __name__
# Use proper names for Tetrad Choices. If no name given (hunt the literature), then use the literature reference as the name.
TetradChoice = par.initialize_param(par.glb_param("char", thismodule, "TetradChoice", "Approx_QuasiKinnersley"))
# Why are these needed?
# xorig = par.initialize_param(par.glb_param("REAL", thismodule, "xorig", "0.0"))
# yorig = par.initialize_param(par.glb_param("REAL", thismodule, "yorig", "0.0"))
# zorig = par.initialize_param(par.glb_param("REAL", thismodule, "zorig", "0.0"))
# offset = par.initialize_param(par.glb_param("REAL", thismodule, "offset", "1.0e-15"))

# Step 2: Define the Levi-Civita symbol, used with tetrads <- better description needed.
def define_LeviCivitaSymbol(DIM=-1):
    if DIM == -1:
        DIM = par.parval_from_str("DIM")

    LeviCivitaSymbol = ixp.zerorank3()

    for i in range(DIM):
        for j in range(DIM):
            for k in range(DIM):
                # From https://codegolf.stackexchange.com/questions/160359/levi-civita-symbol :
                LeviCivitaSymbol[i][j][k] = (i - j) * (j - k) * (k - i) / 2
    return LeviCivitaSymbol

# Step 3: Compute the Weyl scalars
def WeylScalars():
    # Step 1:
    bssn.BSSN_RHSs()

    # Step 2b: Set spatial dimension (must be 3 for BSSN)
    DIM = 3
    par.set_parval_from_str("grid::DIM",DIM)

    x = rfm.xxCart[0]
    y = rfm.xxCart[1]
    z = rfm.xxCart[2]

    if TetradChoice == "Approx_QuasiKinnersley":
        # Eqs 5.6 in https://arxiv.org/pdf/gr-qc/0104063.pdf
        xmoved = x - xorig  # Make sure I'm handling coordinates correctly
        ymoved = y - yorig
        zmoved = z - zorig

        # Eqs 5.7
        vbU = ixp.zerorank1("vbU")
        vaU = ixp.zerorank1("vaU")
        vcU = ixp.zerorank1("vcU")
        vaU[0] = -ymoved
        vaU[1] = xmoved + offset
        vaU[2] = 0
        vbU[0] = xmoved + offset
        vbU[1] = ymoved
        vbU[2] = zmoved
        LeviCivitaSymbol = define_LeviCivitaSymbol()
        for a in range(DIM):
            for b in range(DIM):
                for c in range(DIM):
                    for d in range(DIM):
                        vcU[a] += sp.sqrt(bssn.detgammabar) * bssn.gammabarUU[a][d] * LeviCivitaSymbol[d][b][c] * vaU[b] *vbU[c]

        # Graham-Schmidt orthonormalization of the tetrad
        waU = ixp.zerorank1("waU")
        wbU = ixp.zerorank1("wbU")
        wcU = ixp.zerorank1("wcU")
        eaU = ixp.zerorank1("eaU")
        ebU = ixp.zerorank1("ebU")
        ecU = ixp.zerorank1("ecU")

        waU[a] = vaU[a]
        omega11 = 0
        for a in range(DIM):
            for b in range(DIM):
                omega11 += waU[a] * waU[b] * bssn.gammabarDD[a][b]
        eaU = waU / sp.sqrt(omega11)

        omega12 = 0
        for a in range(DIM):
            for b in range(DIM):
                omega12 += eaU[a] * vaU[b] * bssn.gammabarDD[a][b]
        wbU = vbU - omega12*eaU
        omega22 = 0
        for a in range(DIM):
            for b in range(DIM):
                omega22 += wbU[a] * wbU[b] *bssn.gammabarDD[a][b]
        ebU = wbU / sqrt(omega22)

        omega13 = 0
        for a in range(DIM):
            for b in range(DIM):
                omega13 += eaU[a] * vcU[b] * bssn.gammabarDD[a][b]
        omega23 = 0
        for a in range(DIM):
            for b in range(DIM):
                omega23 += ebU[a] * vcU[b] * bssn.gammabarDD[a][b]
        wcU = vcU - omega13*eaU - omega23*ebU

        # Construct the tetrad
        isqrt2 = 1/sp.sqrt(2)
        ltetU = isqrt2 * ebU
        ntetU = -isqrt2 * ebU
        mtetU = isqrt2 * (ecU + I*eaU)
        mtetbarU = sp.conjugate(mtetU)
