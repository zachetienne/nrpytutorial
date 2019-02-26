# # [BSSN](http://www2.yukawa.kyoto-u.ac.jp/~yuichiro.sekiguchi/3+1.pdf) Hamiltonian and momentum constraint equations, in ***curvilinear*** coordinates, using a covariant reference metric approach: C code generation
# 
# Python module containing the final expressions: [BSSN/BSSN_Constraints.py](../edit/BSSN/BSSN_Constraints.py)
# 
# Citations: Generic curvilinear coordinate reference metric approach matches that of
# [Ruchlin, Etienne, and Baumgarte (2018)](https://arxiv.org/abs/1712.07658),
# which is an extension of the spherical coordinate reference metric approach of
# [Baumgarte, Montero, Cordero-Carri\'on, and M\"uller (2012)](https://arxiv.org/abs/1211.6632),
# which builds upon the covariant "Lagrangian" BSSN formalism of
# [Brown (2009)](https://arxiv.org/abs/0902.3652).
#
# *See also citations within each article.*
# 
# We start by loading the needed modules. Notably, this module depends on several quantities
# defined in the BSSN/BSSN_RHSs.py Python code, documented in the NRPy+ "BSSN in curvilinear
# coordinates" module. Thus in Step 2 below we call BSSN_RHSs() to set these quantities.

# Step 1: Load SymPy and other needed core NRPy+ modules
import sympy as sp
import NRPy_param_funcs as par
import indexedexp as ixp
import grid as gri
import reference_metric as rfm
from outputC import *
import BSSN.BSSN_rescaled_vars as Brv
import BSSN.BSSN_unrescaled_and_barred_vars as Bubv
import finite_difference as fin
import BSSN.BSSN_T4UUmunu_vars as BTmunu

thismodule = __name__

def BSSNConstraints(add_T4UUmunu_source_terms=False):
    # Step 0: Declare as globals all quantities computed in this function.
    global H, MU

    # Step 1: All barred quantities are defined in terms of BSSN rescaled gridfunctions,
    #         which we declare here in case they haven't yet been declared elsewhere.
    Brv.declare_BSSN_rescaled_gridfunctions_if_not_declared_already()
    trK = Brv.trK
    aDD = Brv.aDD

    # Step 2: Evaluate all barred quantities in terms of BSSN rescaled gridfunctions.
    Bubv.BSSN_barred_variables()
    GammabarUDD = Bubv.GammabarUDD
    gammabarUU  = Bubv.gammabarUU
    AbarDD      = Bubv.AbarDD
    AbarUU      = Bubv.AbarUU

    Bubv.RicciBar__gammabarDD_dHatD__DGammaUDD__DGammaU()
    RbarDD      = Bubv.RbarDD

    Bubv.phi_and_derivs()
    phi_dD     = Bubv.phi_dD
    exp_m4phi  = Bubv.exp_m4phi
    phi_dBarD  = Bubv.phi_dBarD
    phi_dBarDD = Bubv.phi_dBarDD

    # Step 3: Set spatial dimension (must be 3 for BSSN)
    DIM = 3
    par.set_parval_from_str("grid::DIM",DIM)

    ###############################
    ###############################
    #  HAMILTONIAN CONSTRAINT
    ###############################
    ###############################

    # Next we define the Hamiltonian constraint. Eq. 13 of [Baumgarte *et al.*](https://arxiv.org/pdf/1211.6632.pdf) yields:
    # $$
    # H = {\underbrace {\textstyle \frac{2}{3} K^2}_{\rm Term\ 1}} -
    # {\underbrace {\textstyle \bar{A}_{ij} \bar{A}^{ij}}_{\rm Term\ 2}} +
    # {\underbrace {\textstyle e^{-4\phi} \left(\bar{R} - 8 \bar{D}^i \phi \bar{D}_i \phi - 8 \bar{D}^2 \phi\right)}_{\rm Term\ 3}}
    # $$

    # Term 1: 2/3 K^2
    H = sp.Rational(2,3)*trK**2

    # Term 2: -A_{ij} A^{ij}
    for i in range(DIM):
        for j in range(DIM):
            H += -AbarDD[i][j]*AbarUU[i][j]

    # Term 3a: trace(Rbar)
    Rbartrace = sp.sympify(0)
    for i in range(DIM):
        for j in range(DIM):
            Rbartrace += gammabarUU[i][j]*RbarDD[i][j]

    # Term 3b: -8 \bar{\gamma}^{ij} \bar{D}_i \phi \bar{D}_j \phi = -8*phi_dBar_times_phi_dBar
    # Term 3c: -8 \bar{\gamma}^{ij} \bar{D}_i \bar{D}_j \phi      = -8*phi_dBarDD_contraction
    phi_dBar_times_phi_dBar = sp.sympify(0) # Term 3b
    phi_dBarDD_contraction  = sp.sympify(0) # Term 3c
    for i in range(DIM):
        for j in range(DIM):
            phi_dBar_times_phi_dBar += gammabarUU[i][j]*phi_dBarD[i]*phi_dBarD[j]
            phi_dBarDD_contraction  += gammabarUU[i][j]*phi_dBarDD[i][j]

    # Add Term 3:
    H += exp_m4phi*(Rbartrace - 8*(phi_dBar_times_phi_dBar + phi_dBarDD_contraction))

    if add_T4UUmunu_source_terms:
        M_PI = par.Cparameters("REAL", thismodule, "M_PI")  # M_PI is pi as defined in C
        BTmunu.define_BSSN_T4UUmunu_rescaled_source_terms()
        rho = BTmunu.rho
        H += -16*M_PI*rho

    # FIXME: ADD T4UUmunu SOURCE TERMS TO MOMENTUM CONSTRAINT!

    ###############################
    ###############################
    #  MOMENTUM CONSTRAINT
    ###############################
    ###############################
    # SEE Tutorial-BSSNConstraints.ipynb for full documentation.
    
    # Let's first implement Terms 2 & 3 of the Momentum constraint:
    MU = ixp.zerorank1()

    # Term 2: 6 A^{ij} \partial_j \phi:
    for i in range(DIM):
        for j in range(DIM):
            MU[i] += 6*AbarUU[i][j]*phi_dD[j]

    # Term 3: -2/3 \bar{\gamma}^{ij} K_{,j}
    trK_dD = ixp.declarerank1("trK_dD") # Not defined in BSSN_RHSs; only trK_dupD is defined there.
    for i in range(DIM):
        for j in range(DIM):
            MU[i] += -sp.Rational(2,3)*gammabarUU[i][j]*trK_dD[j]

    # Finally Term 1: Dbar_j Abar^{ij}

    # We first compute:
    # \bar{D}_{k} \bar{A}_{i j} = \partial_{k} \bar{A}_{i j} - \bar{\Gamma}^{l}_{k i} \bar{A}_{l j} - \bar{\Gamma}^{l}_{k j} \bar{A}_{i l}
    # First define aDD_dD:
    aDD_dD = ixp.declarerank3("aDD_dD","sym01")

    # Then define AbarDD_dD in terms of aDD_dD and other derivatives
    AbarDD_dD = ixp.zerorank3()
    for i in range(DIM):
        for j in range(DIM):
            for k in range(DIM):
                AbarDD_dD[i][j][k] += aDD_dD[i][j][k]*rfm.ReDD[i][j] + aDD[i][j]*rfm.ReDDdD[i][j][k]

    # Then evaluate the conformal covariant derivative \bar{D}_j \bar{A}_{lm}
    AbarDD_dBarD = ixp.zerorank3()
    for i in range(DIM):
        for j in range(DIM):
            for k in range(DIM):
                AbarDD_dBarD[i][j][k] = AbarDD_dD[i][j][k]
                for l in range(DIM):
                    AbarDD_dBarD[i][j][k] += -GammabarUDD[l][k][i]*AbarDD[l][j]
                    AbarDD_dBarD[i][j][k] += -GammabarUDD[l][k][j]*AbarDD[i][l]

    # We then apply two raising operators:
    # \bar{D}_{k} \bar{A}^{i j} = gammabar^{ij} gammabar^{kl} \bar{D}_{k} \bar{A}^{k l}
    # Contract twice with the metric to make \bar{D}_{j} \bar{A}^{ij}
    for i in range(DIM):
        for j in range(DIM):
            for k in range(DIM):
                for l in range(DIM):
                    MU[i] += gammabarUU[i][k]*gammabarUU[j][l]*AbarDD_dBarD[k][l][j]

    # Finally, we multiply by e^{-4 phi} and the appropriate scale factor.
    for i in range(DIM):
        MU[i] *= exp_m4phi / rfm.ReU[i]

def output_C__Hamiltonian_h(add_T4UUmunu_source_terms=False,enable_verbose=True):
    # Calling BSSNConstraints() (defined above) computes H and MU:
    BSSNConstraints(add_T4UUmunu_source_terms)

    import time
    start = time.time()
    if enable_verbose:
        print("Generating optimized C code for Hamiltonian constraint. May take a while, depending on CoordSystem.")
    Hamiltonianstring = fin.FD_outputC("returnstring", lhrh(lhs=gri.gfaccess("aux_gfs", "H"), rhs=H),
                                       params="outCverbose=False")
    end = time.time()
    if enable_verbose:
        print("Finished in " + str(end - start) + " seconds.")

    import loop as lp
    with open("BSSN/Hamiltonian.h", "w") as file:
        file.write(lp.loop(["i2", "i1", "i0"], ["NGHOSTS", "NGHOSTS", "NGHOSTS"],
                           ["NGHOSTS+Nxx[2]", "NGHOSTS+Nxx[1]", "NGHOSTS+Nxx[0]"],
                           ["1", "1", "1"], ["const REAL invdx0 = 1.0/dxx[0];\n" +
                                             "const REAL invdx1 = 1.0/dxx[1];\n" +
                                             "const REAL invdx2 = 1.0/dxx[2];\n" +
                                             "#pragma omp parallel for",
                                             "    const REAL xx2 = xx[2][i2];",
                                             "        const REAL xx1 = xx[1][i1];"], "",
                           "const REAL xx0 = xx[0][i0];\n" + Hamiltonianstring))
    print("Output C implementation of Hamiltonian constraint to BSSN/Hamiltonian.h")