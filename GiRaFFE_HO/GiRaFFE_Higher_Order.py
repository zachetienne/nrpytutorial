#!/usr/bin/env python
# coding: utf-8

# $\newcommand{\giraffe}{\texttt{GiRaFFE}}$
# # $\giraffe$: General Relativistic Force-Free Electrodynamics
# ## Porting the original $\giraffe$ code to NRPy+
# 
# Porting the original $\giraffe$ code, as presented in [the original paper](https://arxiv.org/pdf/1704.00599.pdf), to NRPy+ will generally make it easier to maintain, as well as to make changes. Specifically, it will make it nearly trivial to increase the finite-differencing order.

# Our ultimate goal here will be to code the evolution equation (from eq. 13 of the [original paper](https://arxiv.org/pdf/1704.00599.pdf)):
# \begin{equation}
# \partial_t \tilde{S}_i = - \partial_j \left( \alpha \sqrt{\gamma} T^j_{{\rm EM} i} \right) + \frac{1}{2} \alpha \sqrt{\gamma} T^{\mu \nu}_{\rm EM} \partial_i g_{\mu \nu},
# \end{equation}
# where the densitized spatial Poynting flux one-form $\tilde{S}_i = \sqrt{\gamma} S_i$ (and $S_i$ comes from $S_{\mu} -n_{\nu} T^{\nu}_{{\rm EM} \mu}$) and
# \begin{align}
# T^{\mu \nu}_{\rm EM} &= b^2 u^{\mu} u^{\nu} + \frac{b^2}{2} g^{\mu \nu} - b^{\mu} b^{\nu}, \\
# \sqrt{4\pi} b^0 &= B^0_{\rm (u)} = \frac{u_j B^j}{\alpha}, \\
# \sqrt{4\pi} b^i &= B^i_{\rm (u)} = \frac{B^i + (u_j B^j) u^i}{\alpha u^0}, \\
# \end{align}
# and 
# $$B^i = \frac{\tilde{B}^i}{\sqrt{\gamma}}.$$
# Furthermore, the four-metric $g_{\mu\nu}$ is related to the three-metric $\gamma_{ij}$, lapse $\beta_i$, and shift $\alpha$ by
# $$
# g_{\mu\nu} = \begin{pmatrix} 
# -\alpha^2 + \beta^k \beta_k & \beta_j \\
# \beta_i & \gamma_{ij}
# \end{pmatrix}.
# $$
# Most of these are computed in the module u0_smallb_Poynting__Cartesian.py, and we will import that module to save effort.
# Note that as usual, Greek indices refer to four-dimensional quantities where the zeroth component indicates $t$ components, while Latin indices refer to three-dimensional quantities. Since Python always indexes its lists from 0, however, the zeroth component will indicate a spatial direction, and any expressions involving mixed Greek and Latin indices will need to offset one set of indices by one.

# ## Preliminaries
# First, we will import the core modules of NRPy that we will need and specify the main gridfunctions we will need. The metric quantities $\alpha$, $\beta^i$, and $\gamma_{ij}$ will be initially set by the $\text{ShiftedKerrSchild}$ thorn, while the Valencia 3-velocity $v^i_{(n)}$ and vector potential $A_i$ will initially be set in the separate thorn $\text{GiRaFFEFood_HO}$.

# In[1]:


import NRPy_param_funcs as par
import indexedexp as ixp
import grid as gri
import finite_difference as fin
from outputC import *

def GiRaFFE_Higher_Order():
    #Step 0: Set the spatial dimension parameter to 3.
    par.set_parval_from_str("grid::DIM", 3)
    DIM = par.parval_from_str("grid::DIM")

    # Step 1: Set the finite differencing order to 4.
    par.set_parval_from_str("finite_difference::FD_CENTDERIVS_ORDER", 4)

    thismodule = "GiRaFFE_NRPy"

    # M_PI will allow the C code to substitute the correct value
    M_PI = par.Cparameters("REAL",thismodule,"M_PI")
    # Here, we declare the 3-metric, shift, and lapse as usual
    gammaDD = ixp.register_gridfunctions_for_single_rank2("AUX","gammaDD", "sym01",DIM=3)
    betaU   = ixp.register_gridfunctions_for_single_rank1("AUX","betaU",DIM=3)
    alpha   = gri.register_gridfunctions("AUX","alpha")
    # These two gridfunctions are the basic inputs into StildeD_rhs
    ValenciavU = ixp.register_gridfunctions_for_single_rank1("AUX","ValenciavU",DIM=3)
    AD = ixp.register_gridfunctions_for_single_rank1("AUX","AD",DIM=3)

    # Step 2: Import the four metric
    global gammaUU,gammadet
    gammaUU = ixp.register_gridfunctions_for_single_rank2("AUX","gammaUU","sym01")
    gammadet = gri.register_gridfunctions("AUX","gammadet")
    gammaUU, gammadet = ixp.symm_matrix_inverter3x3(gammaDD)


    # Then, we will declare the gridfunctions related to the metric and build the four metric using code from [Tutorial-smallb2_Poynting_vector-Cartesian.ipynb](Tutorial-smallb2_Poynting_vector-Cartesian.ipynb), which requires the magnetic field $B^i$. So, we will also define the magnetic field as in [eq. 18](https://arxiv.org/pdf/1704.00599.pdf):
    # $$B^i = \epsilon^{ijk} \partial_j A_k,$$
    # where $\epsilon^{ijk}$ the rank-3 Levi-Civita tensor, related to the rank-3 Levi-Civita symbol $[ijk]$ and determinant of the three-metric $\gamma$ by $$\epsilon^{ijk} = [ijk]/\sqrt{\gamma}.$$ 


    # We already have a handy function to define the Levi-Civita symbol in WeylScalars
    import WeylScal4NRPy.WeylScalars_Cartesian as weyl
    LeviCivitaDDD = weyl.define_LeviCivitaSymbol_rank3()
    LeviCivitaUUU = ixp.zerorank3()
    for i in range(DIM):
        for j in range(DIM):
            for k in range(DIM):
                LCijk = LeviCivitaDDD[i][j][k]
                LeviCivitaDDD[i][j][k] = LCijk * sp.sqrt(gammadet)
                LeviCivitaUUU[i][j][k] = LCijk / sp.sqrt(gammadet)

    AD_dD = ixp.declarerank2("AD_dD","nosym")
    # With Levi-Civita and the derivative of A_i, we set B^i to the curl of A^i
    BU = ixp.zerorank1()
    for i in range(DIM):
        for j in range(DIM):
            for k in range(DIM):
                BU[i] += LeviCivitaUUU[i][j][k] * AD_dD[k][j]

    # We use previous work to set u0 and smallb along with the four metric and its inverse
    import u0_smallb_Poynting__Cartesian.u0_smallb_Poynting__Cartesian as u0b
    u0b.compute_u0_smallb_Poynting__Cartesian(gammaDD,betaU,alpha,ValenciavU,BU)


    # Recall that the four-metric $g_{\mu\nu}$ is related to the three-metric $\gamma_{ij}$, lapse $\beta_i$, and shift $\alpha$ by  
    # $$
    # g_{\mu\nu} = \begin{pmatrix} 
    # -\alpha^2 + \beta^k \beta_k & \beta_j \\
    # \beta_i & \gamma_{ij}
    # \end{pmatrix}.
    # $$
    # This tensor and its inverse have already been built by the u0_smallb_Poynting__Cartesian.py module, so we can simply import the variables.


    betaD = ixp.zerorank1()
    for i in range(DIM):
        for j in range(DIM):
            betaD[i] = gammaDD[i][j] * betaU[j]

    # We will now pull in the four metric and its inverse.
    g4DD = ixp.zerorank2(DIM=4)
    g4UU = ixp.zerorank2(DIM=4)
    for mu in range(4):
        for nu in range(4):
            g4DD[mu][nu] = u0b.g4DD[mu][nu]
            g4UU[mu][nu] = u0b.g4UU[mu][nu]


    # We will also need spatial derivatives of the metric, $\partial_i g_{\mu\nu} = g_{\mu\nu,i}$. To do this, we will need derivatives of the shift with its indexed lowered, $\beta_{i,j} = \partial_j \beta_i$. This becomes 
    # \begin{align}
    # \beta_{i,j} &= \partial_j \beta_i \\
    #             &= \partial_j (\gamma_{ik} \beta^k) \\
    #             &= \gamma_{ik} \partial_j\beta^k + \beta^k \partial_j \gamma_{ik} \\
    # \beta_{i,j} &= \gamma_{ik} \beta^k_{\ ,j} + \beta^k \gamma_{ik,j} \\
    # \end{align}
    # 
    # In terms of the three-metric, lapse, and shift, we find
    # $$
    # g_{\mu\nu,l} = \begin{pmatrix} 
    # -2\alpha \alpha_{,l} + \beta^k_{\ ,l} \beta_k + \beta^k \beta_{k,l} & \beta_{i,l} \\
    # \beta_{j,l} & \gamma_{ij,l}
    # \end{pmatrix}.
    # $$
    # 
    # Since this expression mixes Greek and Latin indices, we will need to store the expressions for each of the three spatial derivatives as separate variables. 
    # Also, consider the term $\beta_{i,j} = \partial_j \beta_i = \partial_j (\gamma_{ik} \beta^k) =  \gamma_{ik} \partial_j\beta^k + \beta^k \partial_j \gamma_{ik} = \gamma_{ik} \beta^k_{\ ,j} + \beta^k \gamma_{ik,j}$, that is, 
    # \begin{align}
    # \beta_{i,j} &= \gamma_{ik} \beta^k_{\ ,j} + \beta^k \gamma_{ik,j}. \\
    # \end{align}
    # 
    # So, we will first set 
    # $$ g_{00,l} = \underbrace{-2\alpha \alpha_{,l}}_{\rm Term\ 1} + \underbrace{\beta^k_{\ ,l} \beta_k}_{\rm Term\ 2} + \underbrace{\beta^k \beta_{k,l}}_{\rm Term\ 3} $$


    # Step 3: Build the spatial derivative of the four metric
    # Step 3a: Declare derivatives of grid functions. These will be handled by FD_outputC
    alpha_dD   = ixp.declarerank1("alpha_dD")
    betaU_dD   = ixp.declarerank2("betaU_dD","nosym")
    gammaDD_dD = ixp.declarerank3("gammaDD_dD","sym01")

    # Step 3b: These derivatives will be constructed analytically.
    betaDdD    = ixp.zerorank2()

    g4DDdD     = ixp.zerorank3(DIM=4)

    for i in range(DIM):
        for j in range(DIM):
            for k in range(DIM):
                betaDdD[i][j] = gammaDD[i][k] * betaU_dD[k][j] + betaU[k] * gammaDD_dD[i][k][j]

    # Step 3c: Set the 00 components
    # Step 3c.i: Term 1: -2\alpha \alpha_{,l}
    for l in range(DIM):
        g4DDdD[0][0][l+1] = -2*alpha*alpha_dD[l]

    # Step 3c.ii: Term 2: \beta^k_{\ ,l} \beta_k
    for l in range(DIM):
        for k in range(DIM):
            g4DDdD[0][0][l+1] += betaU_dD[k][l] * betaD[k]

    # Step 3c.iii: Term 3: \beta^k \beta_{k,l}
    for l in range(DIM):
        for k in range(DIM):
            g4DDdD[0][0][l+1] += betaU[k] * betaDdD[k][l]



    # Now we will contruct the other components of $g_{\mu\nu,l}$. We will first construct 
    # $$ g_{i0,l} = g_{0i,l} = \beta_{i,l}, $$
    # then
    # $$ g_{ij,l} = \gamma_{ij,l} $$


    # Step 3d: Set the i0 and 0j components
    for l in range(DIM):
        for i in range(DIM):
            g4DDdD[i+1][0][l+1] = g4DDdD[0][i+1][l+1] = betaDdD[i][l]

    #Step 3e: Set the ij components
    for l in range(DIM):
        for i in range(DIM):
            for j in range(DIM):
                g4DDdD[i+1][j+1][l+1] = gammaDD_dD[i][j][l]


    # ## $T^{\mu\nu}_{\rm EM}$ and its derivatives
    # 
    # Now that the metric and its derivatives are out of the way, we will return our attention to the electromagnetic stress-energy tensor, drawn from eq. 27 of [this paper](https://arxiv.org/pdf/1310.3274.pdf):
    # $$T^{\mu \nu}_{\rm EM} = b^2 u^{\mu} u^{\nu} + \frac{b^2}{2} g^{\mu \nu} - b^{\mu} b^{\nu}.$$
    # We will need the four-velocity $u^\mu$, where 
    # \begin{align}
    # u^i &= u^0 (\alpha v^i_{(n)} - \beta^i), \\
    # u_j &= \alpha u^0 \gamma_{ij} v^i_{(n)}, \\
    # \end{align}
    # and $v^i_{(n)}$ is the Valencia three-velocity, as shown in [Duez, et al, eqs. 53 and 56](https://arxiv.org/pdf/astro-ph/0503420.pdf). These have already been built by the u0_smallb_Poynting__Cartesian.py module, so we can simply import the variables.
    # 


    # Step 4a: import the four-velocity terms
    #u0 = par.Cparameters("REAL",thismodule,"u0")
    #u0 = gri.register_gridfunctions("AUX","u0")
    global uD,uU
    uD = ixp.register_gridfunctions_for_single_rank1("AUX","uD")
    uU = ixp.register_gridfunctions_for_single_rank1("AUX","uU")

    u0 = u0b.u0
    for i in range(DIM):
        uD[i] = u0b.uD[i]
        uU[i] = u0b.uU[i]


    # We also need the vector $b^{\mu}$ before we can compute this, which is related to the magnetic field by 
    # \begin{align}
    # b^0 &= \frac{1}{\sqrt{4\pi}} B^0_{\rm (u)} = \frac{u_j B^j}{\sqrt{4\pi}\alpha}, \\
    # b^i &= \frac{1}{\sqrt{4\pi}} B^i_{\rm (u)} = \frac{B^i + (u_j B^j) u^i}{\sqrt{4\pi}\alpha u^0}, \\
    # \end{align} and \begin{align}
    # B^i &= \frac{\tilde{B}^i}{\sqrt{\gamma}},
    # \end{align}
    # where $B^i$ is the variable tracked by the HydroBase thorn in the Einstein Toolkit. These have already been built by the u0_smallb_Poynting__Cartesian.py module, so we can simply import the variables.


    # Step 4b: import the small b terms
    smallb4U = ixp.zerorank1(DIM=4)
    smallb4D = ixp.zerorank1(DIM=4)
    for mu in range(4):
        smallb4U[mu] = u0b.smallb4U[mu]
        smallb4D[mu] = u0b.smallb4D[mu]

    smallb2 = u0b.smallb2


    # We now have all the pieces to calculate the stress-energy tensor,
    # $$T^{\mu \nu}_{\rm EM} = \underbrace{b^2 u^{\mu} u^{\nu}}_{\rm Term\ 1} + 
    # \underbrace{\frac{b^2}{2} g^{\mu \nu}}_{\rm Term\ 2}
    # - \underbrace{b^{\mu} b^{\nu}}_{\rm Term\ 3}.$$
    # Because $u^0$ is a separate variable, we could build the $00$ component separately, then the $\mu0$ and $0\nu$ components, and finally the $\mu\nu$ components. Alternatively, for clarity, we could create a temporary variable $u^\mu=\left( u^0, u^i \right)$


    # Step 5: Construct the electromagnetic stress-energy tensor
    # Step 5a: Set up the four-velocity vector
    u4U = ixp.zerorank1(DIM=4)
    u4U[0] = u0
    for i in range(DIM):
        u4U[i+1] = uU[i]

    # Step 5b: Build T4EMUU itself
    T4EMUU = ixp.zerorank2(DIM=4)
    for mu in range(4):
        for nu in range(4):
            # Term 1: \underbrace{b^2 u^{\mu} u^{\nu}}
            T4EMUU[mu][nu] = smallb2*u4U[mu]*u4U[nu]

    for mu in range(4):
        for nu in range(4):
            # Term 2: \frac{b^2}{2} g^{\mu \nu}
            T4EMUU[mu][nu] += smallb2*g4UU[mu][nu]/2

    for mu in range(4):
        for nu in range(4):
            # Term 3: b^{\mu} b^{\nu}
            T4EMUU[mu][nu] += -smallb4U[mu]*smallb4U[nu]


    # If we look at the evolution equation, we see that we will need spatial  derivatives of $T^{\mu\nu}_{\rm EM}$. In previous tutorials, when confronted with derivatives of complicated expressions, we have declared those expressions as gridfunctions themselves, thus allowing NRPy+ to take finite-difference derivatives of the expressions. This generally reduces the error, because the alternative is to use a function of several finite-difference derivatives, allowing more error to accumulate than the extra gridfunction will introduce. While we will use that technique for some of the subexpressions of $T^{\mu\nu}_{\rm EM}$, we don't want to rely on it for the whole expression; doing so would require us to take the derivative of the magnetic field $B^i$, which is itself found by finite-differencing the vector potential $A_i$. It requires some finesse, then, to find $B^i$ inside the ghost zones, and those values of $B^i$ are necessarily less accurate; taking another derivative from them would only compound the problem. Instead, we will take analytic derivatives of $T^{\mu\nu}_{\rm EM}$.
    # 
    # We will now now take these spatial  derivatives of $T^{\mu\nu}_{\rm EM}$, applying the chain rule until it is only in terms of basic gridfunctions: $\alpha$, $\beta^i$, $\gamma_{ij}$, $A_i$, and the Valencia 3-velocity, $v^i_{(n)}$. We will need the definitions of $\tilde{B}^i$ and $B^i$ in terms of $B^i$ and $A_i$:
    # \begin{align}
    # \tilde{B}^i &= \sqrt{\gamma} B^i \\
    # B^i &= \epsilon^{ijk} \partial_j A_k \\
    # \end{align}
    # So then, 
    # \begin{align}
    # \partial_j T^{j}_{{\rm EM} i} &= \partial_j (\gamma_{ki} T^{kj}_{\rm EM}) \\
    # &= \partial_j [\gamma_{ki} (b^2 u^j u^k + \frac{b^2}{2} g^{jk} - b^j b^k)] \\
    # &= \underbrace{\gamma_{ki,j} T^{kj}_{\rm EM}}_{\rm Term\ 1} + \gamma_{ki} \left( \underbrace{\partial_j \left(b^2 u^j u^k\right)}_{\rm Term\ 2} + \underbrace{\partial_j \left(\frac{b^2}{2} g^{jk}\right)}_{\rm Term\ 3} - \underbrace{\partial_j \left(b^j b^k\right)}_{\rm Term\ 4} \right) \\
    # \end{align}
    # Following the product and chain rules for each term (and recalling that $b^2 = \gamma_{ij} b^i b^j$, so $\partial_k b^2 = \gamma_{ij,k} b^i b^j + \gamma_{ij} b^i_{,k} b^j + \gamma_{ij} b^i b^j_{,k} = \gamma_{ij,k} b^i b^j + 2 b_j b^j_{,k}$), we find that 
    # \begin{align}
    # {\rm Term\ 2} &= \partial_j (b^2 u^j u^k) \\
    #               &= \partial_j b^2 u^j u^k + b^2 \partial_j u^j u^k + b^2 u^j \partial_j u_k \\
    #               &= \left(\gamma_{lm,j} b^l b^m + 2 b_l b^l_{,j}\right) u^k u^k + b^2 u^j_{,j} u^k + b^2 u^j u^k_{,j} \\
    # {\rm Term\ 3} &= \partial_j \left(\frac{b^2}{2} g^{jk}\right) \\
    #               &= \frac{1}{2} \left( \partial_j b^2 g^{jk} + b^2 \partial_j g^{jk} \right) \\
    #               &= \frac{1}{2}\left(\gamma_{lm,j} b^l b^m + 2 b_l b^l_{,j}\right)g^{jk} + \frac{b^2}{2} g^{jk}_{\ ,j} \\
    # {\rm Term\ 4} &= \partial_j (b^j b^k) \\
    #               &= b^j_{,j} b^k + b^j b^k_{,j}\\
    # \end{align}
    # 
    # So, 
    # \begin{align}
    # \partial_j T^{j}_{{\rm EM} i} &= \gamma_{ki,j} T^{kj}_{\rm EM} \\
    # &+ \gamma_{ki} \left(\left(\gamma_{lm,j} b^l b^m + 2 b_l b^l_{,j}\right) u^j u^k +b^2 u^j_{,j} u^k + b^2 u^j u^k_{,j} + \frac{1}{2}\left(\gamma_{lm,j} b^l b^m + 2 b_l b^l_{,j}\right)g^{jk} + \frac{b^2}{2} g^{jk}_{\ ,j} + b^j_{,j} b^k + b^j b^k_{,j}\right),
    # \end{align}
    # where 
    # \begin{align}
    # %u^i_{,j} &= u^0_{,j} (\alpha v^i_{(n)} - \beta^i) + u^0 (\alpha_{,j} v^i_{(n)} + \alpha v^i_{(n),j} - \beta^i_{,j}) \\
    # %u_{j,k} &= \alpha_{,k} u^0 \gamma_{ij} v^i_{(n)} + \alpha u^0_{,k} \gamma_{ij} v^i_{(n)} + \alpha u^0 \gamma_{ij,k} v^i_{(n)} + \alpha u^0 \gamma_{ij} v^i_{(n),k} \\
    # b^i_{,k} &= \frac{1}{\sqrt{4 \pi}} \frac{\left(\alpha u^0\right)  \left(B^i_{,k} + u_{j,k} B^j u^i + u_j B^j_{,k} u^i + u_j B^j u^i_{,k}\right) - \left(B^i + (u_j B^j) u^i\right) \partial_k \left(\alpha u^0\right)}{\left(\alpha u^0\right)^2} \\
    # B^i_{,i} &= -\frac{\gamma_{,i} B^i}{2\gamma} \\
    # B^i_{,l} &= -\frac{\gamma_{,l}}{2\gamma} B^i + \epsilon^{ijk} A_{k,jl}, i \neq l, \\
    # \end{align}
    # 
    # First, we will handle the derivatives of the velocity $u^i$ and its lowered form. We will simply use finite-difference derivatives here by setting $u^i$ and $u_i$ as gridfunctions (done above). This will reduce the chances for finite differencing errors to accumulate, because the alternative is to use a function of several finite-difference derivatives.


    # Step 6: Derivatives of the electromagnetic stress-energy tensor
    # Step 6a: Declare gridfunctions and their derivatives that will be useful for TEMUD_dD
    # We already handled the ADMBase variables' derivatives when we built g4DDdD.
    # That leaves the valencia 3 velocity and tilde B field.
    ValenciavU_dD = ixp.declarerank2("ValenciavU_dD","nosym")
    #BtildeU_dD    = ixp.declarerank2("BtildeU_dD",   "nosym")
    global alphau0
    alphau0 = gri.register_gridfunctions("AUX","alphau0")
    alphau0 = alpha * u0
    alphau0_dD = ixp.declarerank1("alphau0_dD")

    # We'll let these be handled by finite differencing.
    uU_dD = ixp.declarerank2("uU_dD","nosym")
    #for i in range(DIM):
    #     for j in range(DIM):
    #            uU_dD[i][j] = u0_dD[j]*(alpha*ValenciavU[i]-betaU[i]) + u0*(alpha_dD[j]*ValenciavU[i]\
    #                                                                        +alpha*ValenciavU_dD[i][j]\
    #                                                                        +betaU_dD[i][j])

    uD_dD = ixp.declarerank2("uD_dD","nosym")
    #for i in range(DIM):
    #    for j in range(DIM):
    #        for k in range(DIM):
    #            uD_dD[j][k] =  alpha_dD[k]*u0*gammaDD[i][j]*ValenciavU[i]\
    #                         + alpha*u0_dD[k]*gammaDD[i][j]*ValenciavU[i]\
    #                         + alpha*u0*gammaDD_dD[i][j][k]*ValenciavU[i]\
    #                         + alpha*u0*gammaDD[i][j]*ValenciavU_dD[i][k]\


    # Now, we will build the derivatives of the magnetic field. We will build one expression for the divergence of $b^i$, $b^i_{,i}$ (reducing to functions of $\tilde{B}^i$ (which is itself a function of $B^i$) since $\tilde{B}^i$ is divergenceless) and another for $b^i_{,l}$, reducing to functions of second derivatives of the vector potential $A_i$. That is, since 
    # $$
    # b^i_{,k} = \frac{1}{\sqrt{4 \pi}} \frac{\left(\alpha u^0\right)  \left(B^i_{,k} + u_{j,k} B^j u^i + u_j B^j_{,k} u^i + u_j B^j u^i_{,k}\right) - \left(B^i + (u_j B^j) u^i\right) \partial_k \left(\alpha u^0\right)}{\left(\alpha u^0\right)^2},
    # $$
    # 
    # we will need to define $B^i_{,i}$ and $B^i_{,l}$. For the divergence $B^i_{,i}$, since $B^i = \tilde{B}^i/\sqrt{\gamma}$,
    # the quotient rule tells us that 
    # \begin{align}
    # B^i_{,i} &= \frac{\sqrt{\gamma}\tilde{B}^i_{,i} - \tilde{B}^i \partial_i \sqrt{\gamma}}{\sqrt{\gamma}^2} \\
    #          &= \frac{\sqrt{\gamma}\tilde{B}^i_{,i} - \tilde{B}^i \frac{1}{2\sqrt{\gamma}} \gamma_{,i}}{\gamma}. \\
    # \end{align}
    # Since we have defined the divergence $\tilde{B}^i_{,i} = 0$, the first term in the numerator vanishes, and we are left with
    # \begin{align}
    # B^i_{,i} &= \frac{-\tilde{B}^i \frac{1}{2\sqrt{\gamma}} \gamma_{,i}}{\gamma} \\
    #          &= \frac{-\gamma_{,i} \tilde{B}^i}{2\gamma^{3/2}} \\
    #          &= \frac{-\gamma_{,i} B^i}{2\gamma}, \\
    # \end{align}
    # because we have defined $\tilde{B}^i = \sqrt{\gamma} B^i$. Let us now give a similar treatment to the derivative $B^i_{,l}$; we will start from the definition of $B^i$ in terms of $A_i$, $B^i = \frac{[ijk]}{\sqrt{\gamma}} \partial_j A_k$. We will first apply the product rule, noting that the symbol $[ijk]$ consists purely of the integers $-1, 0, 1$ and thus can be treated as a constant in this process.
    # \begin{align}
    # B^i_{,l} &= \partial_l \left( \frac{[ijk]}{\sqrt{\gamma}} \partial_j A_k \right)  \\
    #          &= [ijk] \partial_l \left( \frac{1}{\sqrt{\gamma}}\right) \partial_j A_k + \frac{[ijk]}{\sqrt{\gamma}} \partial_l \partial_j A_k \\
    #          &= [ijk]\left(-\frac{\gamma_{,l}}{2\gamma^{3/2}}\right) \partial_j A_k + \frac{[ijk]}{\sqrt{\gamma}} \partial_l \partial_j A_k \\
    # \end{align}
    # Now, we will substitute back in for the definition of the Levi-Civita tensor: $\epsilon^{ijk} = [ijk] / \sqrt{\gamma}$. Then will substitute the magnetic field $B^i$ back in.
    # \begin{align}
    # B^i_{,l} &= -\frac{\gamma_{,l}}{2\gamma} \epsilon^{ijk} \partial_j A_k + \epsilon^{ijk} \partial_l \partial_j A_k \\
    #          &= -\frac{\gamma_{,l}}{2\gamma} B^i + \epsilon^{ijk} A_{k,jl}, i \neq l, \\
    # \end{align}
    # Note that the first term here is identical to the divergence-like terms $B^i_{,i}$, so these expressions are consistent given that the divergence of a curl must be zero. Thus, we will not need to code up those terms separately.
    # 
    # Thus, the expression we are left with for the derivatives of the magnetic field is:
    # \begin{align}
    # B^i_{,l} &= \underbrace{-\frac{\gamma_{,l}}{2\gamma} B^i}_{\rm Term\ 1} + \underbrace{\epsilon^{ijk} A_{k,jl}}_{\rm Term\ 2}, \\
    # \end{align}
    # where $\epsilon^{ijk} = [ijk] / \sqrt{\gamma}$ is the antisymmetric Levi-Civita tensor and $\gamma$ is the determinant of the three-metric.
    # 


    #Step 6b: Construct the derivatives of the magnetic field.
    gammadet_dD = ixp.declarerank1("gammadet_dD")

    AD_dDD = ixp.declarerank3("AD_dDD","sym12")           
    # The other partial derivatives of B^i
    BU_dD = ixp.zerorank2()
    for i in range(DIM):
        for l in range(DIM):
            # Term 1: -\frac{\gamma_{,l}}{2\gamma} B^i
            BU_dD[i][l] = -gammadet_dD[l]*BU[i]/(2*gammadet)

    for i in range(DIM):
        for l in range(DIM):
            for j in range(DIM):
                for k in range(DIM):
                    # Term 2: \epsilon^{ijk} A_{k,jl}
                    BU_dD[i][l] += LeviCivitaUUU[i][j][k] * AD_dDD[k][j][l]


    # Now, we will code the derivatives of the spatial components of $b^{\mu}$, $b^i$.
    # $$
    # b^i_{,k} = \frac{1}{\sqrt{4 \pi}} \frac{\left(\alpha u^0\right)  \left(B^i_{,k} + u_{j,k} B^j u^i + u_j B^j_{,k} u^i + u_j B^j u^i_{,k}\right) - \left(B^i + (u_j B^j) u^i\right) \partial_k \left(\alpha u^0\right)}{\left(\alpha u^0\right)^2} 
    # $$
    # 
    # Let's go into a little more detail on where this comes from. We start from the definition $$b^i = \frac{B^i + (u_j B^j) u^i}{\sqrt{4\pi}\alpha u^0};$$ We then apply the quotient rule: 
    # \begin{align}
    # b^i_{,k} &= \frac{\left(\sqrt{4\pi}\alpha u^0\right) \partial_k \left(B^i + (u_j B^j) u^i\right) - \left(B^i + (u_j B^j) u^i\right) \partial_k \left(\sqrt{4\pi}\alpha u^0\right)}{\left(\sqrt{4\pi}\alpha u^0\right)^2} \\
    # &= \frac{1}{\sqrt{4 \pi}} \frac{\left(\alpha u^0\right) \partial_k \left(B^i + (u_j B^j) u^i\right) - \left(B^i + (u_j B^j) u^i\right) \partial_k \left(\alpha u^0\right)}{\left(\alpha u^0\right)^2} \\
    # \end{align}
    # Note that $\left( \alpha u^0 \right)$ is being used as its own gridfunction, so we can be done with that term. We will now apply the product rule to the term $\partial_k \left(B^i + (u_j B^j) u^i\right) = B^i_{,k} + u_{j,k} B^j u^i + u_j B^j_{,k} u^i + u_j B^j u^i_{,k}$. So, 
    # $$ b^i_{,k} = \frac{1}{\sqrt{4 \pi}} \frac{\left(\alpha u^0\right)  \left(B^i_{,k} + u_{j,k} B^j u^i + u_j B^j_{,k} u^i + u_j B^j u^i_{,k}\right) - \left(B^i + (u_j B^j) u^i\right) \partial_k \left(\alpha u^0\right)}{\left(\alpha u^0\right)^2}. $$
    # 
    # It will be easier to code this up if we rearrange these terms to group together the terms that involve contractions over $j$. Doing that, we find
    # $$
    # b^i_{,k} = \frac{\overbrace{\alpha u^0 B^i_{,k} - B^i \partial_k (\alpha u^0)}^{\rm Term\ Num1} + \overbrace{\left( \alpha u^0 \right) \left( u_{j,k} B^j u^i + u_j B^j_{,k} u^i + u_j B^j u^i_{,k} \right)}^{\rm Term\ Num2.a} - \overbrace{\left( u_j B^j u^i \right) \partial_k \left( \alpha u^0 \right) }^{\rm Term\ Num2.b}}{\underbrace{\sqrt{4 \pi} \left( \alpha u^0 \right)^2}_{\rm Term\ Denom}}.
    # $$


    # Step 6c: Construct derivatives of the small b vector
    smallbU_dD = ixp.zerorank2()
    for i in range(DIM):
        for k in range(DIM):
            # Term Num1: \alpha u^0 B^i_{,k} - B^i \partial_k (\alpha u^0)
            smallbU_dD[i][k] += alphau0*BU_dD[i][k]-BU[i]*alphau0_dD[k]

    for i in range(DIM):
        for k in range(DIM):
            for j in range(DIM):
                # Term Num2.a: terms that require contractions over k, and thus an extra loop.
                # \left( \alpha u^0 \right) \left( u_{j,k} B^j u^i + u_j B^j_{,k} u^i + u_j B^j u^i_{,k} \right)
                smallbU_dD[i][k] += alphau0*(uD_dD[j][k]*BU[j]*uU[i]+uD[j]*BU_dD[j][k]*uU[i]+uD[j]*BU[j]*uU_dD[i][k])

    for i in range(DIM):
        for k in range(DIM):
            for j in range(DIM):
                #Term 2.b (More contractions over k): \left( u_j B^j u^i \right) \partial_k \left( \alpha u^0 \right)
                smallbU_dD[i][k] += -(uD[j]*BU[j]*uU[i])*alphau0_dD[k]

    for i in range(DIM):
        for k in range(DIM):
            # Term Denom requires us to divide the whole expressions through by sqrt(4 pi) * (alpha u^0)^2
            smallbU_dD[i][k] /= sp.sqrt(4*M_PI) * alphau0 * alphau0


    # We will also need derivatives of the spatial part of the inverse four-metric: since $g^{ij} = \gamma^{ij} - \frac{\beta^i \beta^j}{\alpha^2}$ ([Gourgoulhon, eq. 4.49](https://arxiv.org/pdf/gr-qc/0703035.pdf)), 
    # \begin{align}
    # g^{ij}_{\ ,k} &= \gamma^{ij}_{\ ,k} - \frac{\alpha^2 \partial_k (\beta^i \beta^j) - \beta^i \beta^j \partial_k \alpha^2}{(\alpha^2)^2} \\
    # &= \gamma^{ij}_{\ ,k} - \frac{\alpha^2\beta^i \beta^j_{,k}+\alpha^2\beta^i_{,k} \beta^j-2\beta^i \beta^j \alpha \alpha_{,k}}{\alpha^4}. \\
    # &= \gamma^{ij}_{\ ,k} - \frac{\alpha\beta^i \beta^j_{,k}+\alpha\beta^i_{,k} \beta^j-2\beta^i \beta^j \alpha_{,k}}{\alpha^3} \\
    # g^{ij}_{\ ,k} &= \underbrace{\gamma^{ij}_{\ ,k}}_{\rm Term\ 1} - \underbrace{\frac{\beta^i \beta^j_{,k}}{\alpha^2}}_{\rm Term\ 2} - \underbrace{\frac{\beta^i_{,k} \beta^j}{\alpha^2}}_{\rm Term\ 3} + \underbrace{\frac{2\beta^i \beta^j \alpha_{,k}}{\alpha^3}}_{\rm Term\ 4}. \\
    # \end{align}
    # 


    # Step 6d: Construct derivatives of the spatial components of g^{ij}
    gammaUU_dD = ixp.declarerank3("gammaUU_dD","sym01")

    # The spatial derivatives of the spatial components of the four metric:
    # Term 1: \gamma^{ij}_{\ ,k}
    gSpatialUU_dD = ixp.zerorank3()
    for i in range(DIM):
        for j in range(DIM):
            for k in range(DIM):
                gSpatialUU_dD[i][j][k] = gammaUU_dD[i][j][k]

    # Term 2: - \frac{\beta^i \beta^j_{,k}}{\alpha^2}
    for i in range(DIM):
        for j in range(DIM):
            for k in range(DIM):
                gSpatialUU_dD[i][j][k] += -betaU[i]*betaU_dD[j][k]/alpha**2

    # Term 3: - \frac{\beta^i_{,k} \beta^j}{\alpha^2}
    for i in range(DIM):
        for j in range(DIM):
            for k in range(DIM):
                gSpatialUU_dD[i][j][k] += -betaU_dD[i][k]*betaU[j]/alpha**2

    # Term 4: \frac{2\beta^i \beta^j \alpha_{,k}}{\alpha^3}
    for i in range(DIM):
        for j in range(DIM):
            for k in range(DIM):
                gSpatialUU_dD[i][j][k] += 2*betaU[i]*betaU[j]*alpha_dD[k]/alpha**3


    # So, we can now put it all together:
    # \begin{align}
    # \partial_j  T^{j}_{{\rm EM} i} &= \gamma_{ki,j} T^{kj}_{\rm EM} + \gamma_{ki} \left(\left(\gamma_{lm,j} b^l b^m + 2 b_l b^l_{,j}\right) u^j u^k +b^2 u^j_{,j} u^k + b^2 u^j u^k_{,j} + \frac{1}{2}\left(\gamma_{lm,j} b^l b^m + 2 b_l b^l_{,j}\right)g^{jk} + \frac{b^2}{2} g^{jk}_{\ ,j} + b^j_{,j} b^k + b^j b^k_{,j}\right).
    # \end{align}
    # It should be noted that due to the way our indexing conventions have have fallen, the Python indices for $T^{ij}_{\rm EM}$, $g^{ij}$, $b^i$ and $b_i$ will need to be incremented to correctly use the spatial components. We will also quickly rearrange the terms of the expression to better mimic the loop structure we will need to create:
    # \begin{align}
    # \partial_j  T^{j}_{{\rm EM} i} =& \ \gamma_{ki,j} T^{kj}_{\rm EM} \\
    # & + \gamma_{ki} \left( b^2 u^j_{,j} u^k + b^2 u^j u^k_{,j} + \frac{b^2}{2} g^{jk}_{\ ,j} + b^j_{,j} b^k + b^j b^k_{,j} \right) \\
    # & + \gamma_{ki} \left( 2b_l b^l_{,j} u^j u^k + 2 b_l b^l_{,j} g^{jk} \right) \\
    # & + \gamma_{ki} \left( \gamma_{lm,j} b^l b^m u^j u^k + \frac{1}{2} \gamma_{lm,j} b^l b^m g^{jk} \right). \\
    # \end{align}
    # 
    # However, we can simplify this a bit more and pull out some common factors, leaving us with 
    # \begin{align}
    # \partial_j  T^{j}_{{\rm EM} i} =& \ 
    # \underbrace{\gamma_{ki,j} T^{kj}_{\rm EM}}_{\rm Term\ 1} \\
    # & + \underbrace{\gamma_{ki} \left( b^2 u^j_{,j} u^k + b^2 u^j u^k_{,j} + \frac{b^2}{2} g^{jk}_{\ ,j} + b^j_{,j} b^k + b^j b^k_{,j} \right)}_{\rm Term\ 2} \\
    # & + \underbrace{2 b_l b^l_{,j} \gamma_{ki} \left( u^j u^k +  g^{jk} \right)}_{\rm Term\ 3} \\
    # & + \underbrace{\gamma_{ki} \gamma_{lm,j} b^l b^m \left( u^j u^k + \frac{1}{2} g^{jk} \right)}_{\rm Term\ 4}. \\
    # \end{align}
    # We will now construct this term by term. Term 1 is straightforward: $${\rm Term\ 1} = \gamma_{ki,j} T^{kj}_{\rm EM}.$$


    # Step 6e: Construct TEMUD_dD itself
    # We will only set the divergence-like components that we need.
    TEMUD_dD_contracted = ixp.zerorank1()
    for i in range(DIM):
        for j in range(DIM):
            for k in range(DIM):
                TEMUD_dD_contracted[i] += gammaDD_dD[k][i][j] * T4EMUU[k+1][j+1]


    # We will now add $${\rm Term\ 2} = \gamma_{ki} \left( \underbrace{b^2 u^j_{,j} u^k}_{\rm Term\ 2a} + \underbrace{b^2 u^j u^k_{,j}}_{\rm Term\ 2b} + \underbrace{\frac{b^2}{2} g^{jk}_{\ ,j}}_{\rm Term\ 2c} + \underbrace{b^j_{,j} b^k}_{\rm Term\ 2d} + \underbrace{b^j b^k_{,j}}_{\rm Term\ 2e} \right)$$ to $\partial_j  T^{j}_{{\rm EM} i}$. These are the terms that involve contractions over $k$ (but no metric derivatives like Term 1 had). 
    # 


    for i in range(DIM):
        for j in range(DIM):
            for k in range(DIM):
                # Term 2a: \gamma_{ki} b^2 u^j_{,j} u^k
                TEMUD_dD_contracted[i] += gammaDD[k][i]*smallb2*uU_dD[j][j]*uU[k]

    for i in range(DIM):
        for j in range(DIM):
            for k in range(DIM):
                # Term 2b: \gamma_{ki} b^2 u^j u^k_{,j}
                TEMUD_dD_contracted[i] += gammaDD[k][i]*smallb2*uU[j]*uU_dD[k][j]

    for i in range(DIM):
        for j in range(DIM):
            for k in range(DIM):
                # Term 2c: \gamma_{ki} \frac{b^2}{2} g^{jk}_{\ ,j}
                TEMUD_dD_contracted[i] += gammaDD[k][i]*smallb2*gSpatialUU_dD[j][k][j]/2

    for i in range(DIM):
        for j in range(DIM):
            for k in range(DIM):
                # Term 2d: \gamma_{ki} b^j_{,j} b^k
                TEMUD_dD_contracted[i] += gammaDD[k][i]*smallbU_dD[j][j]*smallb4U[k+1]

    for i in range(DIM):
        for j in range(DIM):
            for k in range(DIM):
                # Term 2e: \gamma_{ki} b^j b^k_{,j}
                TEMUD_dD_contracted[i] += gammaDD[k][i]*smallb4U[j+1]*smallbU_dD[k][j]


    # Now, we will add $${\rm Term\ 3} = 2 b_l b^l_{,j} \gamma_{ki} \left( \underbrace{u^j u^k}_{\rm Term\ 3a} + \underbrace{g^{jk}}_{\rm Term\ 3b} \right),$$ which involves contractions over $k$ and $l$.


    for i in range(DIM):
        for j in range(DIM):
            for k in range(DIM):
                for l in range(DIM):
                    # Term 3a: 2 b_l b^l_{,j} \gamma_{ki} u^j u^k
                    TEMUD_dD_contracted[i] += 2*smallb4D[l+1]*smallbU_dD[l][j]*gammaDD[k][i]*uU[j]*uU[k]

    for i in range(DIM):
        for j in range(DIM):
            for k in range(DIM):
                for l in range(DIM):
                    # Term 3b: b_l b^l_{,j} \gamma_{ki} g^{jk}
                    TEMUD_dD_contracted[i] += smallb4D[l+1]*smallbU_dD[l][j]*gammaDD[k][i]*g4UU[j+1][k+1]


    # Finally, we will add contractions over $k$, $l$, and $m$: $${\rm Term\ 4} = \gamma_{ki} \gamma_{lm,j} b^l b^m \left( \underbrace{u^j u^k}_{\rm Term\ 4a} + \underbrace{\frac{1}{2} g^{jk}}_{\rm Term\ 4b} \right)$$


    for i in range(DIM):
        for j in range(DIM):
            for k in range(DIM):
                for l in range(DIM):
                    for m in range(DIM):
                        # Term 4a: \gamma_{ki} \gamma_{lm,j} b^l b^m u^j u^k
                        TEMUD_dD_contracted[i] += gammaDD[k][i]*gammaDD_dD[l][m][j]*smallb4U[l+1]*smallb4U[m+1]*uU[j]*uU[k]

    for i in range(DIM):
        for j in range(DIM):
            for k in range(DIM):
                for l in range(DIM):
                    for m in range(DIM):
                        # Term 4b: \gamma_{ki} \gamma_{lm,j} b^l b^m \frac{1}{2} g^{jk}
                        TEMUD_dD_contracted[i] += gammaDD[k][i]*gammaDD_dD[l][m][j]*smallb4U[l+1]*smallb4U[m+1]*g4UU[j+1][k+1]


    # ## Evolution equation for $\tilde{S}_i$
    # Finally, we will return our attention to the time evolution equation (from eq. 13 of the [original paper](https://arxiv.org/pdf/1704.00599.pdf)),
    # \begin{align}
    # \partial_t \tilde{S}_i &= - \partial_j \left( \alpha \sqrt{\gamma} T^j_{{\rm EM} i} \right) + \frac{1}{2} \alpha \sqrt{\gamma} T^{\mu \nu}_{\rm EM} \partial_i g_{\mu \nu} \\
    #                        &= -T^j_{{\rm EM} i} \partial_j (\alpha \sqrt{\gamma}) - \alpha \sqrt{\gamma} \partial_j T^j_{{\rm EM} i} + \frac{1}{2} \alpha \sqrt{\gamma} T^{\mu \nu}_{\rm EM} \partial_i g_{\mu \nu} \\
    #                        &= \underbrace{\frac{1}{2} \alpha \sqrt{\gamma} T^{\mu \nu}_{\rm EM} \partial_i g_{\mu \nu}}_{\rm Term\ 1} - \underbrace{\alpha \sqrt{\gamma} \partial_j T^j_{{\rm EM} i}}_{\rm Term\ 2} - \underbrace{T^j_{{\rm EM} i} \partial_j (\alpha \sqrt{\gamma})}_{\rm Term\ 3} .
    # \end{align}
    # We construct the first term separately at first, to reduce the complication of dealing with mixed Greek and Latin indices.
    # Then we will take derivatives of $\alpha \sqrt{\gamma}$.


    # Step 7: Construct the evolution equation for \tilde{S}_i
    # Here, we set up the necessary machinery to take FD derivatives of alpha * sqrt(gamma)
    global alpsqrtgam
    alpsqrtgam = gri.register_gridfunctions("AUX","alpsqrtgam")
    alpsqrtgam = alpha*sp.sqrt(gammadet)
    alpsqrtgam_dD = ixp.declarerank1("alpsqrtgam_dD")

    global Stilde_rhsD
    Stilde_rhsD = ixp.zerorank1()
    # The first term: \frac{1}{2} \alpha \sqrt{\gamma} T^{\mu \nu}_{\rm EM} \partial_i g_{\mu \nu}
    for i in range(DIM):
        for mu in range(4):
            for nu in range(4):
                Stilde_rhsD[i] += alpsqrtgam * T4EMUU[mu][nu] * g4DDdD[mu][nu][i+1] / 2

    # The second term: \alpha \sqrt{\gamma} \partial_j T^j_{{\rm EM} i}
    for i in range(DIM):
        Stilde_rhsD[i] += -alpsqrtgam * TEMUD_dD_contracted[i]

    # The third term: T^j_{{\rm EM} i} \partial_j (\alpha \sqrt{\gamma})
    for i in range(DIM):
        for j in range(DIM):
            for k in range(DIM):
                Stilde_rhsD[i] += -gammaDD[i][k]*T4EMUU[k][j]*alpsqrtgam_dD[j]


    # ## Evolution equations for $A_i$ and $\Phi$
    # 
    # We will also need to evolve the vector potential $A_i$. This evolution is given as eq. 17 in the [$\giraffe$](https://arxiv.org/pdf/1704.00599.pdf) paper:
    # $$\partial_t A_i = \epsilon_{ijk} v^j B^k - \partial_i (\alpha \Phi - \beta^j A_j),$$
    # where $\epsilon_{ijk} = [ijk] \sqrt{\gamma}$ is the antisymmetric Levi-Civita tensor, the drift velocity $v^i = u^i/u^0$, and $\gamma$ is the determinant of the three metric. 
    # The scalar electric potential $\Phi$ is also evolved by eq. 19:
    # $$\partial_t [\sqrt{\gamma} \Phi] = -\partial_j (\alpha\sqrt{\gamma}A^j - \beta^j [\sqrt{\gamma} \Phi]) - \xi \alpha [\sqrt{\gamma} \Phi],$$
    # with $\xi$ chosen as a damping factor. 
    # 
    # After declaring a  some needed quantities, we will also define the parenthetical terms that we need to take derivatives of. For $A_i$, we will define $${\rm AevolParen} = \alpha \Phi - \beta^j A_j,$$ and for $\sqrt{\gamma} \Phi$, we will define $${\rm PevolParenU[j]} = \alpha\sqrt{\gamma}A^j - \beta^j [\sqrt{\gamma} \Phi].$$


    #Step 8: Construct some useful auxiliary gridfunctions for the other evolution equations
    xi = par.Cparameters("REAL",thismodule,"xi") # The damping factor

    # Call sqrt(gamma)Phi psi6Phi
    psi6Phi = gri.register_gridfunctions("AUX","psi6Phi")
    Phi = psi6Phi / sp.sqrt(gammadet)

    # We'll define a few extra gridfunctions to avoid complicated derivatives
    global AevolParen,PevolParenU
    AevolParen = gri.register_gridfunctions("AUX","AevolParen")
    PevolParenU = ixp.register_gridfunctions_for_single_rank1("AUX","PevolParenU")

    # {\rm AevolParen} = \alpha \Phi - \beta^j A_j
    AevolParen = alpha*Phi
    for j in range(DIM):
        AevolParen     = -betaU[j] * AD[j]

    # {\rm PevolParenU[j]} = \alpha\sqrt{\gamma}A^j - \beta^j [\sqrt{\gamma} \Phi]
    for j in range(DIM):
        PevolParenU[j] = -betaU[j] * psi6Phi
        for i in range(DIM):
            PevolParenU[j] += alpha * sp.sqrt(gammadet) * gammaUU[i][j] * AD[i]

    AevolParen_dD = ixp.declarerank1("AevolParen_dD")
    PevolParenU_dD = ixp.declarerank2("PevolParenU_dD","nosym")


    # Now to set the evolution equations ([eqs. 17 and 19](https://arxiv.org/pdf/1704.00599.pdf)):
    # \begin{align}
    # \partial_t A_i &= \epsilon_{ijk} v^j B^k - \partial_i (\alpha \Phi - \beta^j A_j) \\
    #                &= \epsilon_{ijk} v^j B^k - {\rm AevolParen\_dD[i]} \\
    # \partial_t [\sqrt{\gamma} \Phi] &= -\partial_j (\alpha\sqrt{\gamma}A^j - \beta^j [\sqrt{\gamma} \Phi]) - \xi \alpha [\sqrt{\gamma} \Phi] \\
    #                                 &= {\rm PevolParenU\_dD[j][j]} - \xi \alpha [\sqrt{\gamma} \Phi]. \\
    # \end{align}


    # Step 8b: Construct the evolution equations for A_i and sqrt(gamma)Phi
    global A_rhsD,psi6Phi_rhs
    A_rhsD = ixp.zerorank1()
    psi6Phi_rhs = sp.sympify(0)

    for i in range(DIM):
        A_rhsD[i] = -AevolParen_dD[i]
        for j in range(DIM):
            for k in range(DIM):
                A_rhsD[i] += LeviCivitaDDD[i][j][k]*(uU[j]/u0)*BU[k]

    psi6Phi_rhs = -xi*alpha*psi6Phi
    for j in range(DIM):
        psi6Phi_rhs += -PevolParenU_dD[j][j]


    # ## Setting the initial $\tilde{S}_i$ from initial data
    # We will now find the densitized Poynting flux given by equation 21 in [the original paper](https://arxiv.org/pdf/1704.00599.pdf), $$\tilde{S}_i = \gamma_{ij} \frac{(v^j+\beta^j)\sqrt{\gamma}B^2}{4 \pi \alpha}.$$ This is needed to set initial data for $\tilde{S}_i$ after $B^i$ is set from the initial $A_i$. Note, however, that this expression uses the drift velocity $v^i = \alpha v^i_{(n)} - \beta^i$; substituting this into the definition of $\tilde{S}_i$ yields an expression in terms of the Valencia velocity: $$\tilde{S}_i = \gamma_{ij} \frac{v^i_{(n)} \sqrt{\gamma}B^2}{4 \pi}.$$


    # Step 9: Build the expression for \tilde{S}_i
    global StildeD
    StildeD = ixp.zerorank1()
    BU = ixp.declarerank1("BU") # Reset the values in BU so that the C code accesses the gridfunctions directly
    B2 = sp.sympify(0)
    for i in range(DIM):
        for j in range(DIM):
            B2 += gammaDD[i][j] * BU[i] * BU[j]
    for i in range(DIM):
        StildeD[i] = 0
        for j in range(DIM):
            StildeD[i] += gammaDD[i][j] * (ValenciavU[j])*sp.sqrt(gammadet)*B2/4/M_PI
