
# This module computes the the Weyl scalar $\psi_4$ based on Baker, Campanelli, and Lousto. PRD 65,044001 (2002).
# Step 1: import all needed modules from NRPy+:
import NRPy_param_funcs as par
import indexedexp as ixp
import grid as gri
import finite_difference as fin
import reference_metric as rfm
rfm.reference_metric()
from outputC import *
import BSSN_rhs as bsn
bsn.bssn_metric()
from sympy import conjugate,I

# Step 1a: define a useful function
def levicivita(i,j,k):
    i = int(i)
    j = int(j)
    k = int(k)
    a = j-i
    b = k-j
    c = i-k
    if abs(a) == 2:
        a = a/2

# Step 2a: Initialize Psi4_RHS parameters
thismodule = __name__
TetradChoice = par.initialize_param(par.glb_param("int", thismodule, "TetradChoice", "0"))
xorig = par.initialize_param(par.glb_param("double", thismodule, "xorig", "0.0"))
yorig = par.initialize_param(par.glb_param("double", thismodule, "yorig", "0.0"))
zorig = par.initialize_param(par.glb_param("double", thismodule, "zorig", "0.0"))
offset = par.initialize_param(par.glb_param("double", thismodule, "offset", "1.0e-15"))

# Step 2b: Set spatial dimension (must be 3 for BSSN)
DIM = 3
par.set_parval_from_str("grid::DIM",DIM)

x,y,z = gri.register_gridfunctions("AUX",["x","y","z"])

if TetradChoice is 0:
    xmoved = x - xorig  # Make sure I'm handling coordinates correctly
    ymoved = y - yorig
    zmoved = z - zorig
    vaU = ixp.declarerank1("vaU")
    vbU = ixp.declarerank1("vbU")
    vcU = ixp.declarerank1("vcU")
    waU = ixp.declarerank1("waU")
    wbU = ixp.declarerank1("wbU")
    wcU = ixp.declarerank1("wcU")
    eaU = ixp.declarerank1("eaU")
    ebU = ixp.declarerank1("ebU")
    ecU = ixp.declarerank1("ecU")
    vaU[0] = -ymoved
    vaU[1] = xmoved + offset
    vaU[2] = 0
    vbU[0] = xmoved + offset
    vbU[1] = ymoved
    vbU[2] = zmoved
    for a in range(DIM):
        for b in range(DIM):
            for c in range(DIM)
                for d in range(DIM)
                    vcU[a] += sqrt(detgammabar) * gammabarUU[a][d] * levicivita(d,b,c) * vaU[b] *vbU[c]

    # TO DO: code up the Levi-Civita tensor
    # Graham-Schmidt orthonormalization of the tetrad
    waU[a] = vaU[a]
    omega11 = 0
    for a in range(DIM):
        for b in range(DIM):
            omega11 += waU[a] * waU[b] * gammabarDD[a][b]
    eaU = waU / sqrt(omega11)

    omega12 = 0
    for a in range(DIM):
        for b in range(DIM):
            omega12 += eaU[a] * vaU[b] * gammabarDD[a][b]
    wbU = vbU - omega12*eaU
    omega22 = 0
    for a in range(DIM):
        for b in range(DIM):
            omega22 += wbU[a] * wbU[b] *gammabarDD[a][b]
    ebU = wbU / sqrt(omega22)

    omega13 = 0
    for a in range(DIM):
        for b in range(DIM):
            omega13 += eaU[a] * vcU[b] * gammabarDD[a][b]
    omega23 = 0
    for a in range(DIM):
        for b in range(DIM):
            omega23 += ebU[a] * vcU[b] * gammabarDD[a][b]
    wcU = vcU - omega13*eaU - omega23*ebU

    # Construct the tetrad
    isqrt2 = 1/sqrt(2)
    ltetU = isqrt2 * ebU
    ntetU = -isqrt2 * ebU
    mtetU = isqrt2 * (ecU + I*eaU)
    mtetbarU = conjugate(mtetU)