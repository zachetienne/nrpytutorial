#!/usr/bin/env python
# coding: utf-8

# $\newcommand{\giraffe}{\texttt{GiRaFFE}}$
# # $\giraffe$: General Relativistic Force-Free Electrodynamics
# ## Porting the original $\giraffe$ code to NRPy+
# 
# Porting the original $\giraffe$ code as presented in [the original paper](https://arxiv.org/pdf/1704.00599.pdf) to NRPy+ will generally make it easier to maintain, as well as to make changes. Specifically, it will make it nearly trivial to increase the finite-differencing order.

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
# $$B^i = \frac{\tilde{B}^i}{\gamma}.$$
# Furthermore, the four-metric $g_{\mu\nu}$ is related to the three-metric $\gamma_{ij}$, lapse $\beta_i$, and shift $\alpha$ by
# $$
# g_{\mu\nu} = \begin{pmatrix} 
# -\alpha^2 + \beta^k \beta_k & \beta_i \\
# \beta_j & \gamma_{ij}
# \end{pmatrix}.
# $$
# Most of these are computed in the module u0_smallb_Poynting__Cartesian.py, and we will import that module to save effort.
# Note that as usual, Greek indices refer to four-dimensional quantities where the zeroth component indicates $t$ components, while Latin indices refer to three-dimensional quantities. Since Python always indexes its lists from 0, however, the zeroth component will indicate a spatial direction, and any expressions involving mixed Greek and Latin indices will need to offset one set of indices by one.

# ## Preliminaries
# First, we will import the core modules of NRPy that we will need. Then, we will declare the gridfunctions related to the metric and build the four metric using code from [Tutorial-smallb2_Poynting_vector-Cartesian.ipynb](Tutorial-smallb2_Poynting_vector-Cartesian.ipynb). We will also define the magnetic field $B^k$ as in [eq. 18](https://arxiv.org/pdf/1704.00599.pdf):
# $$B^i = \epsilon^{ijk} \partial_j A_k.$$


import NRPy_param_funcs as par
import indexedexp as ixp
import grid as gri
import finite_difference as fin
from outputC import *

#Step 0: Set the spatial dimension parameter to 3.
par.set_parval_from_str("grid::DIM", 3)
DIM = par.parval_from_str("grid::DIM")

# Step 1: Set the finite differencing order to 4.
par.set_parval_from_str("finite_difference::FD_CENTDERIVS_ORDER", 4)

thismodule = __name__

def GiRaFFE_Higher_Order():
    M_PI = par.Cparameters("REAL",thismodule,"M_PI")
    gammaDD = ixp.register_gridfunctions_for_single_rank2("AUX","gammaDD", "sym01",DIM=3)
    betaU   = ixp.register_gridfunctions_for_single_rank1("AUX","betaU",DIM=3)
    alpha   = gri.register_gridfunctions("AUX","alpha")
    ValenciavU = ixp.register_gridfunctions_for_single_rank1("AUX","ValenciavU",DIM=3)
    AD = ixp.register_gridfunctions_for_single_rank1("AUX","AD",DIM=3)

    # Step 2: Import the four metric
    gammaUU = ixp.register_gridfunctions_for_single_rank2("AUX","gammaUU","sym01")
    global gammadet # Needed for the A to B routine
    gammadet = gri.register_gridfunctions("AUX","gammadet")
    gammaUU, gammadet = ixp.symm_matrix_inverter3x3(gammaDD)

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
    BU = ixp.zerorank1()
    for i in range(DIM):
        for j in range(DIM):
            for k in range(DIM):
                BU[i] += LeviCivitaUUU[i][j][k] * AD_dD[k][j]

    import u0_smallb_Poynting__Cartesian.u0_smallb_Poynting__Cartesian as u0b
    u0b.compute_u0_smallb_Poynting__Cartesian(gammaDD,betaU,alpha,ValenciavU,BU)


    # Recall that the four-metric $g_{\mu\nu}$ is related to the three-metric $\gamma_{ij}$, lapse $\beta_i$, and shift $\alpha$ by  
    # $$
    # g_{\mu\nu} = \begin{pmatrix} 
    # -\alpha^2 + \beta^k \beta_k & \beta_i \\
    # \beta_j & \gamma_{ij}
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


    # We will also need spatial derivatives of the metric, $\partial_i g_{\mu\nu} = g_{\mu\nu,i}$. In terms of the three-metric, lapse, and shift, we find
    # $$
    # g_{\mu\nu,i} = \begin{pmatrix} 
    # -2\alpha \alpha_{,i} + \beta^k_{\ ,i} \beta_k + \beta^k \beta_{k,i} & \beta_{i,i} \\
    # \beta_{j,i} & \gamma_{ij,i}
    # \end{pmatrix}.
    # $$
    # 
    # Since this expression mixes Greek and Latin indices, we will need to store the expressions for each of the three spatial derivatives as separate variables. 
    # Also, consider the term $\beta_{i,j} = \partial_j \beta_i = \partial_j (\gamma_{ik} \beta^k) =  \gamma_{ik} \partial_j\beta^k + \beta^k \partial_j \gamma_{ik} = \gamma_{ik} \beta^k_{\ ,j} + \beta^k \gamma_{ik,j}$


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
    for chi in range(1,4):
        g4DDdD[0][0][i] = -2*alpha*alpha_dD[chi-1]
        for k in range(DIM):
            g4DDdD[0][0][chi] += betaU_dD[k][chi-1] * betaD[k] + betaU[k] * betaDdD[k][chi-1]
        for mu in range(1,4):
            g4DDdD[mu][0][i] = g4DDdD[0][mu][chi-1] = betaDdD[mu-1][chi-1]
        for mu in range(1,4):
            for nu in range(1,4):
                g4DDdD[mu][nu][chi] = gammaDD_dD[mu-1][nu-1][chi-1]


    # Now that the metric and its derivatives are out of the way, we will return our attention to the electromagnetic stress-energy tensor, drawn from eq. 27 of [this paper](https://arxiv.org/pdf/1310.3274.pdf):
    # $$T^{\mu \nu}_{\rm EM} = b^2 u^{\mu} u^{\nu} + \frac{b^2}{2} g^{\mu \nu} - b^{\mu} b^{\nu}.$$
    # We will need the four-velocity $u^\mu$, where 
    # \begin{align}
    # u^i &= u^0 (\alpha v^i_{(n)} - \beta^i), \\
    # u_j &= \alpha u^0 \gamma_{ij} v^i_{(n)}, \\
    # \end{align}
    # and $v^i_{(n)}$ is the Valencia three-velocity, as shown in [Duez, et al, eqs. 53 and 56](https://arxiv.org/pdf/astro-ph/0503420.pdf). The values of $u^0$ and $v^i_{(n)}$ can be read from HydroBase and IllinoisGRMHD. These have already been built by the u0_smallb_Poynting__Cartesian.py module, so we can simply import the variables.
    # 


    #u0 = par.Cparameters("REAL",thismodule,"u0")
    #u0 = gri.register_gridfunctions("AUX","u0")
    uD = ixp.register_gridfunctions_for_single_rank1("AUX","uD")
    global uU
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
    # B^i &= \frac{\tilde{B}^i}{\gamma},
    # \end{align}
    # where $B^i$ is the variable tracked by the HydroBase thorn in the Einstein Toolkit. These have already been built by the u0_smallb_Poynting__Cartesian.py module, so we can simply import the variables.


    smallbU = ixp.zerorank1(DIM=4)
    smallbD = ixp.zerorank1(DIM=4)
    for mu in range(4):
        smallbU[mu] = u0b.smallb4U[mu]
        smallbD[mu] = u0b.smallb4D[mu]

    smallb2 = u0b.smallb2


    # We now have all the pieces to calculate the stress-energy tensor,
    # $$T^{\mu \nu}_{\rm EM} = b^2 u^{\mu} u^{\nu} + \frac{b^2}{2} g^{\mu \nu} - b^{\mu} b^{\nu}.$$
    # 


    TEMUU = ixp.register_gridfunctions_for_single_rank2("AUX","TEMUU","sym01",DIM=4)

    TEMUU[0][0] = smallb2*u0*u0 + smallb2*g4UU[0][0]/2 - smallbU[0]*smallbU[0]
    for mu in range(1,4):
        TEMUU[mu][0] = TEMUU[0][mu] = smallb2*uU[mu-1]*u0 + smallb2*g4UU[mu][0]/2 - smallbU[mu]*smallbU[0]
    for mu in range(1,4):
        for nu in range(1,4):
            TEMUU[mu][nu] = smallb2*uU[mu-1]*uU[nu-1] + smallb2*g4UU[mu][nu]/2 - smallbU[mu]*smallbU[nu]


    # If we look at the evolution equation, we see that we will need spatial  derivatives of $T^{\mu\nu}_{\rm EM}$; we will now now take these derivatives, applying the chain rule until it is only in terms of basic gridfunctions: $\alpha$, $\beta^i$, $\gamma_{ij}$, $A_i$, and the Valencia 3-velocity, $v^i_{(n)}$. We will need the definitions of $\tilde{B}^i$ and $B^i$ in terms of $B^i$ and $A_i$:
    # \begin{align}
    # \tilde{B}^i &= \gamma B^i \\
    # B^i &= \epsilon^{ijk} \partial_j A_k \\
    # \end{align}
    # So then, 
    # \begin{align}
    # \partial_j T^{j}_{{\rm EM} i} &= \partial_j (\gamma_{ki} T^{kj}_{\rm EM}) \\
    # &= \partial_j [\gamma_{ki} (b^2 u^j u^k + \frac{b^2}{2} g^{jk} - b^j b^k)] \\
    # &= \gamma_{ki,j} T^{kj}_{\rm EM} + \\
    # &\gamma_{ki} \left(\left(\gamma_{lm,j} b^l b^m + 2 b_l b^l_{,j}\right) u^j u^k +b^2 u^j_{,j} u^k + b^2 u^j u^k_{,j} + \frac{1}{2}\left(\gamma_{lm,j} b^l b^m + 2 b_l b^l_{,j}\right)g^{jk} + \frac{b^2}{2} g^{jk}_{\ ,j} + b^j_{,j} b^k + b^j b^k_{,j}\right),
    # \end{align}
    # where 
    # \begin{align}
    # u^i_{,j} &= u^0_{,j} (\alpha v^i_{(n)} - \beta^i) + u^0 (\alpha_{,j} v^i_{(n)} + \alpha v^i_{(n),j} - \beta^i_{,j}) \\
    # u_{j,k} &= \alpha_{,k} u^0 \gamma_{ij} v^i_{(n)} + \alpha u^0_{,k} \gamma_{ij} v^i_{(n)} + \alpha u^0 \gamma_{ij,k} v^i_{(n)} + \alpha u^0 \gamma_{ij} v^i_{(n),k} \\
    # b^i_{,k} &= \frac{1}{\sqrt{4 \pi}} \frac{\alpha u^0 (B^i_{,k} + u_{j,k} B^j u^i + u_j B^j_{,k} u^i + u_j B^j u^i_{,k}) - (B^i + u_j B^j u^i)(\alpha_{,k} u^0 + \alpha u^0_{,k})}{(\alpha u^0)^2} \\
    # B^i_{,i} &= \frac{\gamma_{,i} \tilde{B}^i}{\gamma^2} \\
    # B^i_{,l} &= \partial_l \left( \frac{[ijk]}{\sqrt{\gamma}} \partial_j A_k \right) = [ijk] \partial_l \left( \frac{1}{\sqrt{\gamma}}\right) \partial_j A_k + \frac{[ijk]}{\sqrt{\gamma}} \partial_l \partial_j A_k \\
    # &= [ijk]\left(-\frac{\gamma_{,l}}{2\gamma^{3/2}}\right) \partial_j A_k + \frac{[ijk]}{\sqrt{\gamma}} \partial_l \partial_j A_k \\
    # &= -\frac{\gamma_{,l}}{2\gamma} \epsilon^{ijk} \partial_j A_k + \epsilon^{ijk} \partial_l \partial_j A_k \\
    # &= \frac{\gamma_{,l}}{2\gamma} B^i + \epsilon^{ijk} A_{k,jl}, i \neq l, \\
    # \end{align}
    # 
    # First, we will handle the derivatives of the velocity $u^i$ and its lowered form.


    # We already handled the ADMBase variables' derivatives when we built g4DDdD.
    # That leaves the valencia 3 velocity and tilde B field.
    ValenciavU_dD = ixp.declarerank2("ValenciavU_dD","nosym")
    #BtildeU_dD    = ixp.declarerank2("BtildeU_dD",   "nosym")
    global alphau0
    alphau0 = gri.register_gridfunctions("AUX","alphau0")
    alphau0 = alpha * u0
    alphau0_dD = ixp.declarerank1("alphau0_dD")

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

    #alpsqrtgamTEMUD = ixp.register_gridfunctions_for_single_rank2("AUX","alpsqrtgamTEMUD","nosym")
    #for i in range(DIM):
    #    for j in range(DIM):
    #        for k in range(DIM):
    #            # Since TEMUU is a 4D quantity, we increment its indices to match gamma
    #            alpsqrtgamTEMUD[i][j] = alpha * sp.sqrt(gammadet) * gammaDD[k][i] * TEMUU[j+1][k+1]
    #
    #alpsqrtgamTEMUD_dD = ixp.declarerank3("alpsqrtgamTEMUD","nosym")


    # Now, we will build the derivatives of the magnetic field. We will build one expression for the divergence of $b^i$, $b^i_{,i}$ (reducing to functions of $\tilde{B}^i$ (which is itself a function of $B^i$) since $\tilde{B}^i$ is divergenceless) and another for $b^i_{,l}$, reducing to functions of second derivatives of the vector potential $A_i$:
    # \begin{align}
    # b^i_{,k} &= \frac{1}{\sqrt{4 \pi}} \frac{\alpha u^0 (B^i_{,k} + u_{j,k} B^j u^i + u_j B^j_{,k} u^i + u_j B^j u^i_{,k}) - (B^i + u_j B^j u^i)(\alpha_{,k} u^0 + \alpha u^0_{,k})}{(\alpha u^0)^2} \\
    # B^i_{,i} &= -\frac{\gamma_{,i} \tilde{B}^i}{\gamma^2} = -\frac{\gamma_{,i} B^i}{\gamma} \\
    # B^i_{,l} &= \frac{\gamma_{,l}}{2\gamma} B^i + \epsilon^{ijk} A_{k,jl}, i \neq l, \\
    # \end{align}
    # where $\epsilon_{ijk} = [ijk] \sqrt{\gamma}$ is the antisymmetric Levi-Civita tensor and $\gamma$ is the determinant of the three metric.
    # 


    gammadet_dD = ixp.declarerank1("gammadet_dD")

    # The divergence of the magnetic field can be expressed in terms of the B^i itself
    # since Btilde is what we choose to be divergenceless.
    divB = sp.sympify(0)
    for i in range(DIM):
        divB += -gammadet_dD[i]*BU[i]/gammadet

    AD_dDD = ixp.declarerank3("AD_dDD","sym12")           
    # The other partial derivatives of B^i
    BU_dD = ixp.zerorank2()
    for i in range(DIM):
        for l in range(DIM):
            BU_dD[i][l] = gammadet_dD[l]*BU[i]/(2*gammadet)
            for j in range(DIM):
                for k in range(DIM):
                    BU_dD[i][l] += LeviCivitaUUU[i][j][k] * AD_dDD[k][j][l]


    # Now, we will code the derivatives of the spatial componenets of $b^{\mu}$, $b^i$.
    # $$
    # b^i_{,k} = \frac{1}{\sqrt{4 \pi}} \frac{\alpha u^0 (B^i_{,k} + u_{j,k} B^j u^i + u_j B^j_{,k} u^i + u_j B^j u^i_{,k}) - (B^i + u_j B^j u^i)(\alpha_{,k} u^0 + \alpha u^0_{,k})}{(\alpha u^0)^2} 
    # $$


    divb = alphau0*divB
    for i in range(DIM):
        divb += -BU[i]*alphau0_dD[i]
        for j in range(DIM):
            divb += alphau0*(uD_dD[j][i]*BU[j]*uU[i]+uD[j]*BU_dD[j][i]*uU[i]+uD[j]*BU[j]*uU_dD[i][i])                -(uD[j]*BU[j]*uU[i])*alphau0_dD[i]
    divb /= sp.sqrt(4*M_PI) * alpha * u0 * alpha * u0

    smallbU_dD = ixp.zerorank2()
    for i in range(DIM):
        for k in range(DIM):
            smallbU_dD[i][k] += alphau0*BU_dD[i][k]-BU[i]*alphau0_dD[k]
            for j in range(DIM):
                divb += alphau0*(uD_dD[j][k]*BU[j]*uU[i]+uD[j]*BU_dD[j][k]*uU[i]+uD[j]*BU[j]*uU_dD[i][k])                    -(uD[j]*BU[j]*uU[i])*alphau0_dD[k]
            smallbU_dD[i][k] /= sp.sqrt(4*M_PI) * alpha * u0 * alpha * u0


    # We will also need derivatives of the spatial part of the inverse four-metric: since $g^{ij} = \gamma^{ij} - \frac{\beta^i \beta^j}{\alpha^2}$ ([Gourgoulhon, eq. 4.49](https://arxiv.org/pdf/gr-qc/0703035.pdf)), $$g^{ij}_{\ ,k} = \gamma^{ij}_{\ ,k} - \frac{2\alpha^2\beta^i \beta^j_{,k}-2\beta^i \beta^j \alpha \alpha_{,k}}{\alpha^4}.$$ 
    # So, we can now put it all together:
    # \begin{align}
    # \partial_j  T^{j}_{{\rm EM} i} &= \gamma_{ki,j} T^{kj}_{\rm EM} + \gamma_{ki} \left(\left(\gamma_{lm,j} b^l b^m + 2 b_l b^l_{,j}\right) u^j u^k +b^2 u^j_{,j} u^k + b^2 u^j u^k_{,j} + \frac{1}{2}\left(\gamma_{lm,j} b^l b^m + 2 b_l b^l_{,j}\right)g^{jk} + \frac{b^2}{2} g^{jk}_{\ ,j} + b^j_{,j} b^k + b^j b^k_{,j}\right).
    # \end{align}
    # It should be noted that due to the way our indexing conventions have have fallen, the Python indices for $T^{ij}_{\rm EM}$, $g^{ij}$, $b^i$ and $b_i$ will need to be incremented to correctly use the spatial components. We will also quickly rearrange the terms of the expression to better mimic the loop structure we will need to create.
    # \begin{align}
    # \partial_j  T^{j}_{{\rm EM} i} =& \ \gamma_{ki,j} T^{kj}_{\rm EM} \\
    # & + \gamma_{ki} \left( b^2 u^j_{,j} u^k + b^2 u^j u^k_{,j} + \frac{b^2}{2} g^{jk}_{\ ,j} + b^j_{,j} b^k + b^j b^k_{,j} \right) \\
    # & + \gamma_{ki} \left( 2b_l b^l_{,j} u^j u^k + 2 b_l b^l_{,j} g^{jk} \right) \\
    # & + \gamma_{ki} \left( \gamma_{lm,j} b^l b^m u^j u^k + \frac{1}{2} \gamma_{lm,j} b^l b^m g^{jk} \right) \\
    # \end{align}


    gammaUU_dD = ixp.declarerank3("gammaUU_dD","sym01")

    gSpatUU_dD = ixp.zerorank3()
    for i in range(DIM):
        for j in range(DIM):
            for k in range(DIM):
                gSpatUU_dD[i][j][k] = gammaUU_dD[i][j][k] - 2*betaU[i]*(alpha*betaU_dD[j][k]-betaU[j]*alpha_dD[k])/alpha**3

    # We will only set the divergence-like components that we need.
    TEMUD_dD = ixp.zerorank3()
    for i in range(DIM):
        for j in range(DIM):
            for k in range(DIM):
                TEMUD_dD[j][i][j]  = gammaDD_dD[k][i][j] * TEMUU[k+1][j+1]
                TEMUD_dD[j][i][j] += gammaDD[k][i]*(smallb2*uU_dD[j][j]*uU[k]+smallb2*uU[j]*uU_dD[k][j]+                                                smallb2*gSpatUU_dD[j][k][j]/2+smallbU_dD[j][j]*smallbU[k+1]+                                                smallbU[j+1]*smallbU_dD[k][j])
                for l in range(DIM):
                    TEMUD_dD[j][i][j] += gammaDD[k][i]*2*smallbD[l]*smallbU_dD[l][j]*(uU[j]*uU[k]+g4UU[j+1][k+1])
                    for m in range(DIM):
                        TEMUD_dD[j][i][j] += gammaDD[k][i]*gammaDD_dD[l][m][j]*smallbU[l+1]*smallbU[m+1]*(uU[j]*uU[k]                                                                                                      +g4UU[j+1][k+1])


    # Finally, we will return our attention to the time evolution equation (from eq. 13 of the [original paper](https://arxiv.org/pdf/1704.00599.pdf)),
    # \begin{equation}
    # \partial_t \tilde{S}_i = - T^j_{{\rm EM} i} \partial_j (\alpha \sqrt{\gamma}) - \alpha \sqrt{\gamma} \partial_j T^j_{{\rm EM} i} + \frac{1}{2} \alpha \sqrt{\gamma} T^{\mu \nu}_{\rm EM} \partial_i g_{\mu \nu}.
    # \end{equation}
    # We first construct the third term, to reduce the complication of dealing with mixed Greek and Latin indices.
    # Then we will take derivatives of $\alpha \sqrt{\gamma}$.


    thirdterm = ixp.zerorank1()
    for i in range(DIM):
        for mu in range(DIM):
            for nu in range(DIM):
                thirdterm[i] += alpha * sp.sqrt(gammadet) * TEMUU[mu][nu] * g4DDdD[mu][nu][i] / 2

    global alpsqrtgam
    alpsqrtgam = gri.register_gridfunctions("AUX","alpsqrtgam")
    alpsqrtgam = alpha*sp.sqrt(gammadet)
    alpsqrtgam_dD = ixp.register_gridfunctions_for_single_rank1("AUX","alpsqrtgam_dD")

    global Stilde_rhsD
    Stilde_rhsD = ixp.register_gridfunctions_for_single_rank1("EVOL","Stilde_rhsD")
    for i in range(DIM):
        Stilde_rhsD[i] = thirdterm[i]
        for j in range(DIM):
            Stilde_rhsD[i] += -alpsqrtgam * TEMUD_dD[j][i][j]
            for k in range(DIM):
                Stilde_rhsD[i] += -gammaDD[i][k]*TEMUU[k][j]*alpsqrtgam_dD[j]


    # We will also need to evolve the vector potential $A_i$. This evolution is given as eq. 17 in the [$\giraffe$](https://arxiv.org/pdf/1704.00599.pdf) paper:
    # $$\partial_t A_i = \epsilon_{ijk} v^j B^k - \partial_i (\alpha \Phi - \beta^j A_j),$$
    # where $\epsilon_{ijk} = [ijk] \sqrt{\gamma}$ is the antisymmetric Levi-Civita tensor, the drift velocity $v^i = u^i/u^0$, and $\gamma$ is the determinant of the three metric. 
    # The scalar electric potential $\Phi$ is also evolved by eq. 19:
    # $$\partial_t [\sqrt{\gamma} \Phi] = -\partial_j (\alpha\sqrt{\gamma}A^j - \beta^j [\sqrt{\gamma} \Phi]) - \xi \alpha [\sqrt{\gamma} \Phi],$$
    # with $\xi$ chosen as a damping factor. 


    xi = par.Cparameters("REAL",thismodule,"xi") # The damping factor

    # Call sqrt(gamma)Phi psi6Phi
    psi6Phi = gri.register_gridfunctions("AUX","psi6Phi")
    Phi = psi6Phi / sp.sqrt(gammadet)

    # We'll define a few extra gridfunctions to avoid complicated derivatives
    global AevolParen,PevolParenU
    AevolParen = gri.register_gridfunctions("AUX","AevolParen")
    PevolParenU = ixp.register_gridfunctions_for_single_rank1("AUX","PevolParenU")

    AevolParen = alpha*Phi
    for j in range(DIM):
        AevolParen     = -betaU[j] * AD[j]
        PevolParenU[j] = -betaU[j] * psi6Phi
        for i in range(DIM):
            PevolParenU[j] += alpha * sp.sqrt(gammadet) * gammaUU[i][j] * AD[i]

    AevolParen_dD = ixp.register_gridfunctions_for_single_rank1("AUX","AevolParen_dD")
    PevolParenU_dD = ixp.register_gridfunctions_for_single_rank2("AUX","PevolParenU_dD","nosym")


    # Now to set the evolution equations ([eqs. 17 and 19](https://arxiv.org/pdf/1704.00599.pdf)):
    # \begin{align}
    # \partial_t A_i &= \epsilon_{ijk} v^j B^k - \partial_i (\alpha \Phi - \beta^j A_j) \\
    # \partial_t [\sqrt{\gamma} \Phi] &= -\partial_j (\alpha\sqrt{\gamma}A^j - \beta^j [\sqrt{\gamma} \Phi]) - \xi \alpha [\sqrt{\gamma} \Phi]. \\
    # \end{align}


    global A_rhsD,psi6Phi_rhs
    A_rhsD = ixp.register_gridfunctions_for_single_rank1("EVOL","A_rhsD")
    psi6Phi_rhs = gri.register_gridfunctions("EVOL","psi6Phi_rhs")

    for i in range(DIM):
        A_rhsD[i] = -AevolParen_dD[i]
        for j in range(DIM):
            for k in range(DIM):
                A_rhsD[i] += LeviCivitaDDD[i][j][k]*(uU[j]/u0)*BU[k]

    psi6Phi_rhs = -xi*alpha*psi6Phi
    for j in range(DIM):
        psi6Phi_rhs += -PevolParenU_dD[j][j]

