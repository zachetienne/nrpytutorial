#!/usr/bin/env python
# coding: utf-8

# <a id='top'></a>
#
#
# # $\texttt{GiRaFFE}$: General Relativistic Force-Free Electrodynamics
# $$\label{top}$$
#
# **Authors: Patrick Nelson, Zachariah Etienne, George Vopal, and Maria Babiuc-Hamilton**
#
# <font color='red'>**This module is under active development -- do *not* use the resulting C code output for doing science.**</font>
#
# <font color='red'>**GiRaFFE_HO_v2 is an experiment that omits the analytic derivatives of $\partial_j T^j_{{\rm EM} i}$**</font>
#
# ## A new version of $\texttt{GiRaFFE}$, using NRPy+
#
# The original $\texttt{GiRaFFE}$ code, as presented in [the original paper](https://arxiv.org/pdf/1704.00599.pdf), exists as a significant modification to $\texttt{IllinoisGRMHD}$. As such, it used a third-order reconstruction algorithm with a slope limiter (Colella et al's piecewise parabolic method, or PPM) to handle spatial derivatives in the general relativistic force-free electrodynamics (GRFFE) equations. However the GRFFE equations do not generally permit shocks, so a more optimal approach would involve finite differencing all derivatives in the GRFFE equations. As it happens, NRPy+ was designed to generate C codes involving complex tensorial expressions and finite difference spatial derivatives, with finite-differencing order a freely-specifiable parameter.
#
# The purpose of this notebook is to rewrite the equations of GRFFE as used in the original $\texttt{GiRaFFE}$ code so that all derivatives that appear are represented numerically as finite difference derivatives. As we will see, the largest complication stems from derivatives of magnetic fields--requiring judicious application of the chain and product rules.
#
# The GRFFE evolution equations (from eq. 13 of the [original paper](https://arxiv.org/pdf/1704.00599.pdf)) we wish to encode in the NRPy+ version of $\texttt{GiRaFFE}$ are as follows:
#
# * $\partial_t \tilde{S}_i = - \partial_j \left( \alpha \sqrt{\gamma} T^j_{{\rm EM} i} \right) + \frac{1}{2} \alpha \sqrt{\gamma} T^{\mu \nu}_{\rm EM} \partial_i g_{\mu \nu}$: [Link to Step 3.0](#step7)
# * $\partial_t A_i = \epsilon_{ijk} v^j B^k - \partial_i (\alpha \Phi - \beta^j A_j)$: [Link to Step 4.0](#step8)
# * $\partial_t [\sqrt{\gamma} \Phi] = -\partial_j (\alpha\sqrt{\gamma}A^j - \beta^j [\sqrt{\gamma} \Phi]) - \xi \alpha [\sqrt{\gamma} \Phi]$: [Link to Step 4.0](#step8)
#
# Here, the densitized spatial Poynting flux one-form $\tilde{S}_i = \sqrt{\gamma} S_i$ (and $S_i$ comes from $S_{\mu} -n_{\nu} T^{\nu}_{{\rm EM} \mu}$), and $(\Phi, A_i)$ is the vector potential. We will solve these PDEs using the method of lines, where the right-hand sides of these equations (involving no explicit time derivatives) will be constructed using NRPy+.
#
# ### A Note on Notation
#
# As is standard in NRPy+,
#
# * Greek indices refer to four-dimensional quantities where the zeroth component indicates temporal (time) component.
# * Latin indices refer to three-dimensional quantities. This is somewhat counterintuitive since Python always indexes its lists starting from 0. As a result, the zeroth component of three-dimensional quantities will necessarily indicate the first *spatial* direction.
#
# As a corollary, any expressions involving mixed Greek and Latin indices will need to offset one set of indices by one: A Latin index in a four-vector will be incremented and a Greek index in a three-vector will be decremented (however, the latter case does not occur in this tutorial module).
#
# ### Table of Contents:
#
# 1. Preliminaries
#     1. [Steps 1.0-1.2](#steps0to2): Set up the basic NRPy+ infrastructure we need
#         * This part imports the core NRPy+ modules that we need, sets the dimensionality of the grid with parameter $\text{grid::DIM}$, and declares some of the basic gridfunctions
#     1. [Step 1.3](#step3): Build the spatial derivatives of the four metric
# 1. $T^{\mu\nu}_{\rm EM}$ and its derivatives
#     1. [Step 2.0](#step4): $u^i$ and $b^i$ and related quantities
#     1. [Step 2.1](#step5): Construct the electromagnetic stress-energy tensor
#     1. [Step 2.2](#step6): Derivatives of the electromagnetic stress-energy tensor
# 1. Evolution equation for $\tilde{S}_i$ (depends on quantities defined in previous steps)
#     1. [Step 3.0](#step7): Construct the evolution equation for $\tilde{S}_i$
# 1. Evolution equations for $A_i$ and $\Phi$ (depends on quantities defined in previous steps)
#     1. [Step 4.0](#step8): Construct the evolution equations for $A_i$ and $\sqrt{\gamma}\Phi$
# 1. Code Validation
#     1. [Step 5.0](#step9): NRPy+ Module Code Validation

# # Preliminaries
#
# First, we will import the core modules of NRPy that we will need and specify the main gridfunctions we will need.
#
# <a id='steps0to2'></a>
#
# ## Steps 1.0-1.1: Set up the needed NRPy+ infrastructure and declare core gridfunctions used by $\texttt{GiRaFFE}$ \[Back to [top](#top)\]
#
# 1. Set some basic NRPy+ parameters. E.g., set the spatial dimension parameter to 3 and the finite differencing order to 4.
# 1. Next, declare some gridfunctions that are provided as input to the equations:
#     1. $\alpha$, $\beta^i$, and $\gamma_{ij}$: These ADM 3+1 metric quantities are declared in the ADMBase Einstein Toolkit thorn, and are assumed to be made available to $\texttt{GiRaFFE}$ at this stage.
#     1. The Valencia 3-velocity $v^i_{(n)}$ and vector potential $A_i$: Declared by $\texttt{GiRaFFE}$, and will have their initial values set in the separate thorn **GiRaFFEfood_HO**.
#     1. The magnetic field as measured by a normal observer $B^i$: The quantities evolved forward in time in $\texttt{GiRaFFE}$ do not include the Valencia 3-velocity, so this quantity is not automatically updated. Instead, we compute it based on the evolved quantity $\tilde{S}_i$ and $B^i = \epsilon^{ijk} \partial_j A_k$ (where $A_k$ is another evolved quantity and $\epsilon^{ijk}$ is the Levi-Civita tensor). $B^i$ is evaluated using finite differences of $A_k$ in a separate function, though it can only be evaluated consistently **(Error check: how does this sentence end?)**
# $$\label{steps0to2}$$


import NRPy_param_funcs as par
import indexedexp as ixp
import grid as gri
import sympy as sp

def GiRaFFE_Higher_Order_v2():
    #Step 1.0: Set the spatial dimension parameter to 3.
    par.set_parval_from_str("grid::DIM", 3)
    DIM = par.parval_from_str("grid::DIM")

    # Step 1.1: Set the finite differencing order to 4.
    #par.set_parval_from_str("finite_difference::FD_CENTDERIVS_ORDER", 4)

    thismodule = "GiRaFFE_NRPy"

    # M_PI will allow the C code to substitute the correct value
    # M_PI = par.Cparameters("#define",thismodule,"M_PI","")
    # ADMBase defines the 4-metric in terms of the 3+1 spacetime metric quantities gamma_{ij}, beta^i, and alpha
    gammaDD = ixp.register_gridfunctions_for_single_rank2("AUX","gammaDD", "sym01",DIM=3)
    betaU   = ixp.register_gridfunctions_for_single_rank1("AUX","betaU",DIM=3)
    alpha   = gri.register_gridfunctions("AUX","alpha")
    # GiRaFFE uses the Valencia 3-velocity and A_i, which are defined in the initial data module(GiRaFFEfood)
    ValenciavU = ixp.register_gridfunctions_for_single_rank1("AUX","ValenciavU",DIM=3)
    AD = ixp.register_gridfunctions_for_single_rank1("EVOL","AD",DIM=3)
    # B^i must be computed at each timestep within GiRaFFE so that the Valencia 3-velocity can be evaluated
    BU = ixp.register_gridfunctions_for_single_rank1("AUX","BU",DIM=3)


    # <a id='step3'></a>
    #
    # ## Step 1.2: Build the four metric $g_{\mu\nu}$, its inverse $g^{\mu\nu}$ and spatial derivatives $g_{\mu\nu,i}$ from ADM 3+1 quantities $\gamma_{ij}$, $\beta^i$, and $\alpha$ \[Back to [top](#top)\]
    #
    # Notice that the time evolution equation for $\tilde{S}_i$
    # $$
    # \partial_t \tilde{S}_i = - \partial_j \left( \alpha \sqrt{\gamma} T^j_{{\rm EM} i} \right) + \frac{1}{2} \alpha \sqrt{\gamma} T^{\mu \nu}_{\rm EM} \partial_i g_{\mu \nu}
    # $$
    # contains $\partial_i g_{\mu \nu} = g_{\mu\nu,i}$. We will now focus on evaluating this term.
    #
    # The four-metric $g_{\mu\nu}$ is related to the three-metric $\gamma_{ij}$, index-lowered shift $\beta_i$, and lapse $\alpha$ by
    # $$
    # g_{\mu\nu} = \begin{pmatrix}
    # -\alpha^2 + \beta^k \beta_k & \beta_j \\
    # \beta_i & \gamma_{ij}
    # \end{pmatrix}.
    # $$
    # This tensor and its inverse have already been built by the u0_smallb_Poynting__Cartesian.py module ([documented here](Tutorial-u0_smallb_Poynting-Cartesian.ipynb)), so we can simply load the module and import the variables.
    # $$\label{step3}$$


    # Step 1.2: import u0_smallb_Poynting__Cartesian.py to set
    #           the four metric and its inverse. This module also sets b^2 and u^0.
    import u0_smallb_Poynting__Cartesian.u0_smallb_Poynting__Cartesian as u0b
    u0b.compute_u0_smallb_Poynting__Cartesian(gammaDD,betaU,alpha,ValenciavU,BU)

    betaD = ixp.zerorank1()
    for i in range(DIM):
        for j in range(DIM):
            betaD[i] += gammaDD[i][j] * betaU[j]
    # Error check: fixed = to +=

    # We will now pull in the four metric and its inverse.
    import BSSN.ADMBSSN_tofrom_4metric as AB4m # NRPy+: ADM/BSSN <-> 4-metric conversions
    AB4m.g4DD_ito_BSSN_or_ADM("ADM")
    g4DD = AB4m.g4DD
    AB4m.g4UU_ito_BSSN_or_ADM("ADM")
    g4UU = AB4m.g4UU


    # Next we compute spatial derivatives of the metric, $\partial_i g_{\mu\nu} = g_{\mu\nu,i}$, written in terms of the three-metric, shift, and lapse. Simply taking the derivative of the expression for $g_{\mu\nu}$ above, we find
    # $$
    # g_{\mu\nu,l} = \begin{pmatrix}
    # -2\alpha \alpha_{,l} + \beta^k_{\ ,l} \beta_k + \beta^k \beta_{k,l} & \beta_{j,l} \\
    # \beta_{i,l} & \gamma_{ij,l}
    # \end{pmatrix}.
    # $$
    #
    # Notice the derivatives of the shift vector with its indexed lowered, $\beta_{i,j} = \partial_j \beta_i$. This can be easily computed in terms of the given ADMBase quantities $\beta^i$ and $\gamma_{ij}$ via:
    # \begin{align}
    # \beta_{i,j} &= \partial_j \beta_i \\
    #             &= \partial_j (\gamma_{ik} \beta^k) \\
    #             &= \gamma_{ik} \partial_j\beta^k + \beta^k \partial_j \gamma_{ik} \\
    # \beta_{i,j} &= \gamma_{ik} \beta^k_{\ ,j} + \beta^k \gamma_{ik,j}.
    # \end{align}
    #
    # Since this expression mixes Greek and Latin indices, we will declare this as a 4 dimensional quantity, but only set the three spatial components of its last index (that is, leaving $l=0$ unset).
    #
    # So, we will first set
    # $$ g_{00,l} = \underbrace{-2\alpha \alpha_{,l}}_{\rm Term\ 1} + \underbrace{\beta^k_{\ ,l} \beta_k}_{\rm Term\ 2} + \underbrace{\beta^k \beta_{k,l}}_{\rm Term\ 3} $$


    # Step 1.2, cont'd: Build spatial derivatives of the four metric
    # Step 1.2.a: Declare derivatives of grid functions. These will be handled by FD_outputC
    alpha_dD   = ixp.declarerank1("alpha_dD")
    betaU_dD   = ixp.declarerank2("betaU_dD","nosym")
    gammaDD_dD = ixp.declarerank3("gammaDD_dD","sym01")

    # Step 1.2.b: These derivatives will be constructed analytically.
    betaDdD    = ixp.zerorank2()
    g4DDdD     = ixp.zerorank3(DIM=4)

    for i in range(DIM):
        for j in range(DIM):
            for k in range(DIM):
                # \gamma_{ik} \beta^k_{,j} + \beta^k \gamma_{ik,j}
                betaDdD[i][j] += gammaDD[i][k] * betaU_dD[k][j] + betaU[k] * gammaDD_dD[i][k][j]
    #Error check: changed = to += above

    # Step 1.2.c: Set the 00 components
    # Step 1.2.c.i: Term 1: -2\alpha \alpha_{,l}
    for l in range(DIM):
        g4DDdD[0][0][l+1] = -2*alpha*alpha_dD[l]

    # Step 1.2.c.ii: Term 2: \beta^k_{\ ,l} \beta_k
    for l in range(DIM):
        for k in range(DIM):
            g4DDdD[0][0][l+1] += betaU_dD[k][l] * betaD[k]

    # Step 1.2.c.iii: Term 3: \beta^k \beta_{k,l}
    for l in range(DIM):
        for k in range(DIM):
            g4DDdD[0][0][l+1] += betaU[k] * betaDdD[k][l]


    # Now we will contruct the other components of $g_{\mu\nu,l}$. We will first construct
    # $$ g_{i0,l} = g_{0i,l} = \beta_{i,l}, $$
    # then
    # $$ g_{ij,l} = \gamma_{ij,l} $$


    # Step 1.2.d: Set the i0 and 0j components
    for l in range(DIM):
        for i in range(DIM):
            # \beta_{i,l}
            g4DDdD[i+1][0][l+1] = g4DDdD[0][i+1][l+1] = betaDdD[i][l]

    #Step 1.2.e: Set the ij components
    for l in range(DIM):
        for i in range(DIM):
            for j in range(DIM):
                # \gamma_{ij,l}
                g4DDdD[i+1][j+1][l+1] = gammaDD_dD[i][j][l]


    # <a id='step4'></a>
    #
    # # $T^{\mu\nu}_{\rm EM}$ and its derivatives
    #
    # Now that the metric and its derivatives are out of the way, let's return to the evolution equation for $\tilde{S}_i$,
    # $$
    # \partial_t \tilde{S}_i = - \partial_j \left( \alpha \sqrt{\gamma} T^j_{{\rm EM} i} \right) + \frac{1}{2} \alpha \sqrt{\gamma} T^{\mu \nu}_{\rm EM} \partial_i g_{\mu \nu}.
    # $$
    # We turn our focus now to $T^j_{{\rm EM} i}$ and its derivatives. To this end, we start by computing $T^{\mu \nu}_{\rm EM}$ (from eq. 27 of [Paschalidis & Shapiro's paper on their GRFFE code](https://arxiv.org/pdf/1310.3274.pdf)):
    #
    # $$\boxed{T^{\mu \nu}_{\rm EM} = b^2 u^{\mu} u^{\nu} + \frac{b^2}{2} g^{\mu \nu} - b^{\mu} b^{\nu}.}$$
    #
    # Notice that $T^{\mu\nu}_{\rm EM}$ is written in terms of
    #
    # * $b^\mu$, the 4-component magnetic field vector, related to the comoving magnetic field vector $B^i_{(u)}$
    # * $u^\mu$, the 4-velocity
    # * $g^{\mu \nu}$, the inverse 4-metric
    #
    # However, $\texttt{GiRaFFE}$ has access to only the following quantities, requiring in the following sections that we write the above quantities in terms of the following ones:
    #
    # * $\gamma_{ij}$, the 3-metric
    # * $\alpha$, the lapse
    # * $\beta^i$, the shift
    # * $A_i$, the vector potential
    # * $B^i$, the magnetic field (we assume only in the grid interior, not the ghost zones)
    # * $\left[\sqrt{\gamma}\Phi\right]$, the zero-component of the vector potential $A_\mu$, times the square root of the determinant of the 3-metric
    # * $v_{(n)}^i$, the Valencia 3-velocity
    # * $u^0$, the zero-component of the 4-velocity
    #
    # ## Step 2.0: $u^i$ and $b^i$ and related quantities \[Back to [top](#top)\]
    #
    # We begin by importing what we can from u0_smallb_Poynting__Cartesian.py. We will need the four-velocity $u^\mu$, which is related to the Valencia 3-velocity $v^i_{(n)}$ used directly by $\texttt{GiRaFFE}$ (see also [Duez, et al, eqs. 53 and 56](https://arxiv.org/pdf/astro-ph/0503420.pdf))
    # \begin{align}
    # u^i &= u^0 (\alpha v^i_{(n)} - \beta^i), \\
    # u_j &= \alpha u^0 \gamma_{ij} v^i_{(n)},
    # \end{align}
    # where $v^i_{(n)}$ is the Valencia three-velocity. These have already been constructed in terms of the Valencia 3-velocity and other 3+1 ADM quantities by the u0_smallb_Poynting__Cartesian.py module, so we can simply import these variables:
    # $$\label{step4}$$


    # Step 2.0: u^i, b^i, and related quantities
    # Step 2.0.a: import the four-velocity, as written in terms of the Valencia 3-velocity
    uD = ixp.zerorank1()
    uU = ixp.zerorank1()
    u4upperZero = gri.register_gridfunctions("AUX","u4upperZero")

    for i in range(DIM):
        uD[i] = u0b.uD[i].subs(u0b.u0,u4upperZero)
        uU[i] = u0b.uU[i].subs(u0b.u0,u4upperZero)


    # We also need the magnetic field 4-vector $b^{\mu}$, which is related to the magnetic field by [eqs. 23, 24, and 31 in Duez, et al](https://arxiv.org/pdf/astro-ph/0503420.pdf):
    # \begin{align}
    # b^0 &= \frac{1}{\sqrt{4\pi}} B^0_{\rm (u)} = \frac{u_j B^j}{\sqrt{4\pi}\alpha}, \\
    # b^i &= \frac{1}{\sqrt{4\pi}} B^i_{\rm (u)} = \frac{B^i + (u_j B^j) u^i}{\sqrt{4\pi}\alpha u^0}, \\
    # \end{align}
    # where $B^i$ is the variable tracked by the HydroBase thorn in the Einstein Toolkit. Again, these have already been built by the u0_smallb_Poynting__Cartesian.py module, so we can simply import the variables.


    # Step 2.0.b: import the small b terms
    smallb4U = ixp.zerorank1(DIM=4)
    smallb4D = ixp.zerorank1(DIM=4)
    for mu in range(4):
        smallb4U[mu] = u0b.smallb4U[mu].subs(u0b.u0,u4upperZero)
        smallb4D[mu] = u0b.smallb4D[mu].subs(u0b.u0,u4upperZero)

    smallb2 = u0b.smallb2etk.subs(u0b.u0,u4upperZero)


    # <a id='step5'></a>
    #
    # ## Step 2.1: Construct all components of the electromagnetic stress-energy tensor $T^{\mu \nu}_{\rm EM}$ \[Back to [top](#top)\]
    #
    # We now have all the pieces to calculate the stress-energy tensor,
    # $$T^{\mu \nu}_{\rm EM} = \underbrace{b^2 u^{\mu} u^{\nu}}_{\rm Term\ 1} +
    # \underbrace{\frac{b^2}{2} g^{\mu \nu}}_{\rm Term\ 2}
    # - \underbrace{b^{\mu} b^{\nu}}_{\rm Term\ 3}.$$
    # Because $u^0$ is a separate variable, we could build the $00$ component separately, then the $\mu0$ and $0\nu$ components, and finally the $\mu\nu$ components. Alternatively, for clarity, we could create a temporary variable $u^\mu=\left( u^0, u^i \right)$
    # $$\label{step5}$$


    # Step 2.1: Construct the electromagnetic stress-energy tensor
    # Step 2.1.a: Set up the four-velocity vector
    u4U = ixp.zerorank1(DIM=4)
    u4U[0] = u4upperZero
    for i in range(DIM):
        u4U[i+1] = uU[i]

    # Step 2.1.b: Build T4EMUU itself
    T4EMUU = ixp.zerorank2(DIM=4)
    for mu in range(4):
        for nu in range(4):
            # Term 1: b^2 u^{\mu} u^{\nu}
            T4EMUU[mu][nu] = smallb2*u4U[mu]*u4U[nu]

    for mu in range(4):
        for nu in range(4):
            # Term 2: b^2 / 2 g^{\mu \nu}
            T4EMUU[mu][nu] += smallb2*g4UU[mu][nu]/2

    for mu in range(4):
        for nu in range(4):
            # Term 3: -b^{\mu} b^{\nu}
            T4EMUU[mu][nu] += -smallb4U[mu]*smallb4U[nu]


    #
    # # Evolution equation for $\tilde{S}_i$
    # <a id='step7'></a>
    #
    # ## Step 3.0: Construct the evolution equation for $\tilde{S}_i$ \[Back to [top](#top)\]
    #
    # Finally, we will return our attention to the time evolution equation (from eq. 13 of the [original paper](https://arxiv.org/pdf/1704.00599.pdf)),
    # \begin{align}
    # \partial_t \tilde{S}_i &= - \partial_j \underbrace{\left( \alpha \sqrt{\gamma} T^j_{{\rm EM} i} \right)}_{\rm SevolParenUD} + \frac{1}{2} \alpha \sqrt{\gamma} T^{\mu \nu}_{\rm EM} \partial_i g_{\mu \nu} \\
    #                        &= - \partial_j{\rm SevolParenUD[i][j]} + \frac{1}{2} \alpha \sqrt{\gamma} T^{\mu \nu}_{\rm EM} \partial_i g_{\mu \nu} .
    # \end{align}
    # We will first take construct ${\rm SevolParenUD}$, then use its derivatives to construct the evolution equation. Note that
    # \begin{align}
    # {\rm SevolParenUD[i][j]} &= \alpha \sqrt{\gamma} T^j_{{\rm EM} i} \\
    #                              &= \alpha \sqrt{\gamma} g_{\mu i} T^{\mu j}_{\rm EM}.
    # \end{align}
    # $$\label{step7}$$


    # Step 3.0: Construct the evolution equation for \tilde{S}_i
    # Here, we set up the necessary machinery to take FD derivatives of alpha * sqrt(gamma)
    global gammaUU,gammadet
    gammaUU = ixp.register_gridfunctions_for_single_rank2("AUX","gammaUU","sym01")
    gammadet = gri.register_gridfunctions("AUX","gammadet")
    gammaUU, gammadet = ixp.symm_matrix_inverter3x3(gammaDD)

    alpsqrtgam = alpha*sp.sqrt(gammadet)

    global SevolParenUD
    SevolParenUD = ixp.register_gridfunctions_for_single_rank2("AUX","SevolParenUD","nosym")
    SevolParenUD = ixp.zerorank2()

    for i in range(DIM):
        for j in range (DIM):
            for mu in range(4):
                SevolParenUD[j][i] += alpsqrtgam * g4DD[mu][i+1] * T4EMUU[mu][j+1]

    SevolParenUD_dD = ixp.declarerank3("SevolParenUD_dD","nosym")

    global Stilde_rhsD
    Stilde_rhsD = ixp.zerorank1()
    # The first term: \alpha \sqrt{\gamma} \partial_j T^j_{{\rm EM} i}
    for i in range(DIM):
        for j in range(DIM):
            Stilde_rhsD[i] += -SevolParenUD_dD[j][i][j]

    # The second term: \alpha \sqrt{\gamma} T^{\mu \nu}_{\rm EM} \partial_i g_{\mu \nu} / 2
    for i in range(DIM):
        for mu in range(4):
            for nu in range(4):
                Stilde_rhsD[i] += alpsqrtgam * T4EMUU[mu][nu] * g4DDdD[mu][nu][i+1] / 2


    # # Evolution equations for $A_i$ and $\Phi$
    # <a id='step8'></a>
    #
    # ## Step 4.0: Construct the evolution equations for $A_i$ and $[\sqrt{\gamma}\Phi]$ \[Back to [top](#top)\]
    #
    # We will also need to evolve the vector potential $A_i$. This evolution is given as eq. 17 in the [$\texttt{GiRaFFE}$](https://arxiv.org/pdf/1704.00599.pdf) paper:
    # $$\boxed{\partial_t A_i = \epsilon_{ijk} v^j B^k - \partial_i (\underbrace{\alpha \Phi - \beta^j A_j}_{\rm AevolParen}),}$$
    # where $\epsilon_{ijk} = [ijk] \sqrt{\gamma}$ is the antisymmetric Levi-Civita tensor, the drift velocity $v^i = u^i/u^0$, $\gamma$ is the determinant of the three metric, $B^k$ is the magnetic field, $\alpha$ is the lapse, and $\beta$ is the shift.
    # The scalar electric potential $\Phi$ is also evolved by eq. 19:
    # $$\boxed{\partial_t [\sqrt{\gamma} \Phi] = -\partial_j (\underbrace{\alpha\sqrt{\gamma} \gamma^{ij} A_i - \beta^j [\sqrt{\gamma} \Phi]}_{\rm PevolParenU[j]}) - \xi \alpha [\sqrt{\gamma} \Phi],}$$
    # with $\xi$ chosen as a damping factor.
    #
    # ### Step 4.0.a: Construct some useful auxiliary gridfunctions for the other evolution equations
    #
    # After declaring a  some needed quantities, we will also define the parenthetical terms (underbrace above) that we need to take derivatives of. That way, we can take finite-difference derivatives easily. Note that the above equations incorporate the fact that $\gamma^{ij}$ is the appropriate raising operator for $A_i$: $A^j = \gamma^{ij} A_i$. This is correct since $n_\mu A^\mu = 0$, where $n_\mu$ is a normal to the hypersurface, so $A^0=0$ (according to Sec. II, subsection C of [the "Improved EM gauge condition" paper of Etienne *et al*](https://arxiv.org/pdf/1110.4633.pdf)).
    # $$\label{step8}$$


    # Step 4.0: Construct the evolution equations for A_i and sqrt(gamma)Phi
    # Step 4.0.a: Construct some useful auxiliary gridfunctions for the other evolution equations
    xi = par.Cparameters("REAL",thismodule,"xi", 0.1) # The (dimensionful) Lorenz damping factor. Recommendation: set to ~1.5/max(delta t).

    # Define sqrt(gamma)Phi as psi6Phi
    psi6Phi = gri.register_gridfunctions("EVOL","psi6Phi")
    Phi = psi6Phi / sp.sqrt(gammadet)

    # We'll define a few extra gridfunctions to avoid complicated derivatives
    global AevolParen,PevolParenU
    AevolParen  = gri.register_gridfunctions("AUX","AevolParen")
    PevolParenU = ixp.register_gridfunctions_for_single_rank1("AUX","PevolParenU")

    # {\rm AevolParen} = \alpha \Phi - \beta^j A_j
    AevolParen = alpha*Phi
    for j in range(DIM):
        AevolParen    += -betaU[j] * AD[j]

    # {\rm PevolParenU[j]} = \alpha\sqrt{\gamma} \gamma^{ij} A_i - \beta^j [\sqrt{\gamma} \Phi]
    for j in range(DIM):
        PevolParenU[j] = -betaU[j] * psi6Phi
        for i in range(DIM):
            PevolParenU[j] += alpsqrtgam * gammaUU[i][j] * AD[i]

    AevolParen_dD  = ixp.declarerank1("AevolParen_dD")
    PevolParenU_dD = ixp.declarerank2("PevolParenU_dD","nosym")


    # ### Step 4.0.b: Complete the construction of the evolution equations for $A_i$ and $[\sqrt{\gamma}\Phi]$
    #
    # Now to set the evolution equations ([eqs. 17 and 19](https://arxiv.org/pdf/1704.00599.pdf)), recalling that the drift velocity $v^i = u^i/u^0$:
    # \begin{align}
    # \partial_t A_i &= \epsilon_{ijk} v^j B^k - \partial_i (\alpha \Phi - \beta^j A_j) \\
    #                &= \epsilon_{ijk} \frac{u^j}{u^0} B^k - {\rm AevolParen\_dD[i]} \\
    # \partial_t [\sqrt{\gamma} \Phi] &= -\partial_j \left(\left(\alpha\sqrt{\gamma}\right)A^j - \beta^j [\sqrt{\gamma} \Phi]\right) - \xi \alpha [\sqrt{\gamma} \Phi] \\
    #                                 &= -{\rm PevolParenU\_dD[j][j]} - \xi \alpha [\sqrt{\gamma} \Phi]. \\
    # \end{align}


    # Step 4.0.b: Construct the actual evolution equations for A_i and sqrt(gamma)Phi
    global A_rhsD,psi6Phi_rhs
    A_rhsD = ixp.zerorank1()
    psi6Phi_rhs = sp.sympify(0)

    # We already have a handy function to define the Levi-Civita symbol in WeylScalars
    import WeylScal4NRPy.WeylScalars_Cartesian as weyl
    # Initialize the Levi-Civita tensor by setting it equal to the Levi-Civita symbol
    LeviCivitaSymbolDDD = weyl.define_LeviCivitaSymbol_rank3()
    LeviCivitaTensorDDD = ixp.zerorank3()
    #LeviCivitaTensorUUU = ixp.zerorank3()

    for i in range(DIM):
        for j in range(DIM):
            for k in range(DIM):
                LeviCivitaTensorDDD[i][j][k] = LeviCivitaSymbolDDD[i][j][k] * sp.sqrt(gammadet)
                #LeviCivitaTensorUUU[i][j][k] = LeviCivitaSymbolDDD[i][j][k] / sp.sqrt(gammadet)

    for i in range(DIM):
        A_rhsD[i] = -AevolParen_dD[i]
        for j in range(DIM):
            for k in range(DIM):
                A_rhsD[i] += LeviCivitaTensorDDD[i][j][k]*(uU[j]/u4upperZero)*BU[k]

    psi6Phi_rhs = -xi*alpha*psi6Phi
    for j in range(DIM):
        psi6Phi_rhs += -PevolParenU_dD[j][j]
