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
import BSSN.BSSN_RHSs as bssnrhs

def BSSNConstraints():
    # Step 2: Call BSSN_RHSs() to load needed quantities, but only
    #         if it has not already been called; calling
    #         BSSNConstraints() after a BSSN_RHSs() call will result
    #         in a doubly-declared gridfunction error.
    BSSN_RHSs_has_been_called = False
    for i in range(len(gri.glb_gridfcs_list)):
        if "hDD00" in gri.glb_gridfcs_list[i].name:
            BSSN_RHSs_has_been_called = True
    if BSSN_RHSs_has_been_called == False:
        bssnrhs.BSSN_RHSs()

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
    global H
    H = sp.Rational(2,3)*bssnrhs.trK**2

    # Term 2: -A_{ij} A^{ij}
    for i in range(DIM):
        for j in range(DIM):
            H += -bssnrhs.AbarDD[i][j]*bssnrhs.AbarUU[i][j]

    # Term 3a: trace(Rbar)
    Rbartrace = sp.sympify(0)
    for i in range(DIM):
        for j in range(DIM):
            Rbartrace += bssnrhs.gammabarUU[i][j]*bssnrhs.RbarDD[i][j]

    # Term 3b: -8 \bar{\gamma}^{ij} \bar{D}_i \phi \bar{D}_j \phi = -8*phi_dBar_times_phi_dBar
    # Term 3c: -8 \bar{\gamma}^{ij} \bar{D}_i \bar{D}_j \phi      = -8*phi_dBarDD_contraction
    phi_dBar_times_phi_dBar = sp.sympify(0) # Term 3b
    phi_dBarDD_contraction  = sp.sympify(0) # Term 3c
    for i in range(DIM):
        for j in range(DIM):
            phi_dBar_times_phi_dBar += bssnrhs.gammabarUU[i][j]*bssnrhs.phi_dBarD[i]*bssnrhs.phi_dBarD[j]
            phi_dBarDD_contraction  += bssnrhs.gammabarUU[i][j]*bssnrhs.phi_dBarDD[i][j]

    # Add Term 3:
    H += bssnrhs.exp_m4phi*(Rbartrace - 8*(phi_dBar_times_phi_dBar + phi_dBarDD_contraction))

    ###############################
    ###############################
    #  MOMENTUM CONSTRAINT
    ###############################
    ###############################
    # Eq. 14 of [Baumgarte *et al.*](https://arxiv.org/pdf/1211.6632.pdf) for the momentum constraint is *allegedly* missing the
    # term $e^{-4\phi} \bar{A}^{ik} \Delta\Gamma^j_{jk}$, which we include below:

    # \mathcal{M}^i = e^{-4\phi} \left(
    # {\underbrace {\textstyle \frac{1}{\sqrt{\bar{\gamma}}} \hat{D}_j\left(\sqrt{\bar{\gamma}}\bar{A}^{ij}\right)}_{\rm Term\ 1}} +
    # {\underbrace {\textstyle 6 \bar{A}^{ij}\partial_j \phi}_{\rm Term\ 2}} -
    # {\underbrace {\textstyle \frac{2}{3} \bar{\gamma}^{ij}\partial_j K}_{\rm Term\ 3}} +
    # {\underbrace {\textstyle \bar{A}^{jk} \Delta\Gamma^i_{jk} {\color{red}{+ \bar{A}^{ik} \Delta\Gamma^j_{jk}}}}_{\rm Term\ 4}}\right)

    # Let's first implement Terms 2-4 of the Momentum constraint:
    global MU
    MU = ixp.zerorank1()

    # Term 2: 6 A^{ij} \partial_j \phi:
    for i in range(DIM):
        for j in range(DIM):
            MU[i] += 6*bssnrhs.AbarUU[i][j]*bssnrhs.phi_dD[j]

    # Term 3: -2/3 \bar{\gamma}^{ij} K_{,j}
    trK_dD = ixp.declarerank1("trK_dD") # Not defined in BSSN_RHSs; only trK_dupD is defined there.
    for i in range(DIM):
        for j in range(DIM):
            MU[i] += -sp.Rational(2,3)*bssnrhs.gammabarUU[i][j]*trK_dD[j]

    # Term 4: \bar{A}^{jk} \Delta\Gamma^i_{jk} + \bar{A}^{ik} \Delta\Gamma^j_{jk}
    for i in range(DIM):
        for j in range(DIM):
            for k in range(DIM):
                MU[i] += bssnrhs.AbarUU[j][k]*bssnrhs.DGammaUDD[i][j][k]
                # Pending further investigation: + bssnrhs.AbarUU[i][k]*bssnrhs.DGammaUDD[j][j][k]

    # Term 1: \frac{1}{\bar{\gamma}} \hat{D}_j ( \sqrt{\bar{\gamma}}\bar{A^{ij}} )

    # \begin{align}
    # \text{Term 1} &= \frac{1}{\sqrt{\bar{\gamma}}} \hat{D}_j\left(\sqrt{\bar{\gamma}}\bar{A}^{ij}\right)\\
    # &= {\underbrace {\textstyle \hat{D}_j \bar{A}^{ij}}_{\rm Term\ 1a}} +
    # {\underbrace {\textstyle \frac{1}{2} \bar{A}^{ij} \frac{\bar{\gamma}_{,j}}{\bar{\gamma}}}_{\rm Term\ 1b}}
    # \end{align}
    #
    # Let's first implement Term 1a, $\hat{D}_j \bar{A}^{ij}$:
    #
    # Since the "up-up" tensor $\bar{A}^{ij}$ cannot be easily expressed in terms of the BSSN gridfunctions (i.e., the
    # functions we actually sample to take derivatives), we must rewrite Term 1a in terms of derivatives of
    # $\bar{A}_{ij}$:
    # \begin{align}
    # \hat{D}_j \bar{A}^{ij} &= \hat{D}_j \left(\bar{\gamma}^{i\ell}\bar{\gamma}^{jm} \bar{A}_{\ell m}\right)\\
    # &=
    # {\underbrace {\textstyle \bar{A}_{\ell m} \hat{D}_j \left(\bar{\gamma}^{i\ell}\bar{\gamma}^{jm}\right)}_{\rm Term\ 1a.i}} +
    # {\underbrace {\textstyle \bar{\gamma}^{i\ell}\bar{\gamma}^{jm}\hat{D}_j \bar{A}_{\ell m}}_{\rm Term\ 1a.ii}}.
    # \end{align}
    #
    # Similarly, the "up-up" tensor $\bar{\gamma}^{ij}$ cannot be easily expressed in terms of the BSSN gridfunctions,
    # so we must express terms like $\hat{D}_j\bar{\gamma}^{i\ell}$ in terms of $\hat{D}_j \bar{\gamma}_{i\ell}$,
    # as computed in the [BSSN RHS tutorial module (needed for $\bar{R}_{ij}$)](Tutorial-BSSNCurvilinear.ipynb).
    # This is straightforward given the following identity:
    # \begin{align}
    # 0 &= \hat{D}_{k} \delta_{i}^{j} \\
    # &= \hat{D}_{k} (\bar{\gamma}_{i l} \bar{\gamma}^{l j}) \\
    # &= \bar{\gamma}^{l j} \hat{D}_{k} \bar{\gamma}_{i l} + \bar{\gamma}_{i l} \hat{D}_{k} \bar{\gamma}^{l j} \\
    # \implies \bar{\gamma}_{i l} \hat{D}_{k} \bar{\gamma}^{l j} &= -\bar{\gamma}^{l j} \hat{D}_{k} \bar{\gamma}_{i l}\\
    # \implies \bar{\gamma}^{i m} \bar{\gamma}_{i l} \hat{D}_{k} \bar{\gamma}^{l j}
    # &= -\bar{\gamma}^{i m} \bar{\gamma}^{l j} \hat{D}_{k} \bar{\gamma}_{i l}\\
    # \implies \hat{D}_{k} \bar{\gamma}^{m j} &= -\bar{\gamma}^{i m} \bar{\gamma}^{l j} \hat{D}_{k} \bar{\gamma}_{i l}.
    # \end{align}
    #
    # Next, the covariant derivative $\hat{D}_j \bar{A}_{\ell m}$ is, by definition:

    # \hat{D}_j \bar{A}_{\ell m} = \partial_j \bar{A}_{\ell m}
    # - \hat{\Gamma}^k_{j\ell} \bar{A}_{km}
    # - \hat{\Gamma}^k_{jm}    \bar{A}_{\ell k}.

    # First we implement Term 1a.i:
    gammabarUU_dHatD = ixp.zerorank3()
    for i in range(DIM):
        for j in range(DIM):
            for k in range(DIM):
                for l in range(DIM):
                    for m in range(DIM):
                        gammabarUU_dHatD[m][j][k] += \
                            -bssnrhs.gammabarUU[i][m]*bssnrhs.gammabarUU[l][j]*bssnrhs.gammabarDD_dHatD[i][l][k]

    # Next we implement Term 1a.ii:

    # First define aDD_dD:
    AbarDD_dD = ixp.zerorank3()
    aDD_dD = ixp.declarerank3("aDD_dD","sym01")
    for i in range(DIM):
        for j in range(DIM):
            for k in range(DIM):
                AbarDD_dD[i][j][k] += aDD_dD[i][j][k]*rfm.ReDD[i][j] + bssnrhs.aDD[i][j]*rfm.ReDDdD[i][j][k]

    # Then evaluate \hat{D}_j \bar{A}_{lm}
    AbarDD_dHatD = ixp.zerorank3()
    for j in range(DIM):
        for l in range(DIM):
            for m in range(DIM):
                AbarDD_dHatD[l][m][j] = AbarDD_dD[l][m][j]
                for k in range(DIM):
                    AbarDD_dHatD[l][m][j] += -rfm.GammahatUDD[k][j][l]*bssnrhs.AbarDD[k][m]
                    AbarDD_dHatD[l][m][j] += -rfm.GammahatUDD[k][j][m]*bssnrhs.AbarDD[l][k]


    for i in range(DIM):
        for j in range(DIM):
            for l in range(DIM):
                for m in range(DIM):
                    MU[i] += bssnrhs.AbarDD[l][m]*(gammabarUU_dHatD[i][l][j]*bssnrhs.gammabarUU[j][m] +
                                                   bssnrhs.gammabarUU[i][l]*gammabarUU_dHatD[j][m][j]) # Term 1a.i
                    MU[i] += bssnrhs.gammabarUU[i][l]*bssnrhs.gammabarUU[j][m]*AbarDD_dHatD[l][m][j]   # Term 1a.ii

    # Next we implement Term 1b:
    for i in range(DIM):
        for j in range(DIM):
            MU[i] += sp.Rational(1,2)*bssnrhs.AbarUU[i][j]*bssnrhs.detgammabar_dD[j]/bssnrhs.detgammabar

    # Finally, we multiply by e^{-4 phi} and the appropriate scale factor.
    for i in range(DIM):
        MU[i] *= bssnrhs.exp_m4phi / rfm.ReU[i]