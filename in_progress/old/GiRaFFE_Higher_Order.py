#!/usr/bin/env python
# coding: utf-8

# <a id='top'></a>
#
#
# # $\texttt{GiRaFFE}$: General Relativistic Force-Free Electrodynamics
# $$\label{top}$$
#
# **Authors: Patrick Nelson, Zachariah Etienne, and Maria Babiuc-Hamilton**
#
# <font color='red'>**This module is under active development -- do *not* use the resulting C code output for doing science.**</font>
#
# ## A new version of $\texttt{GiRaFFE}$, using NRPy+
#
# The original $\texttt{GiRaFFE}$ code, as presented in [the original paper](https://arxiv.org/pdf/1704.00599.pdf), exists as a significant modification to $\texttt{IllinoisGRMHD}$. As such, it used a third-order reconstruction algorithm with a slope limiter (Colella et al's piecewise parabolic method, or PPM) to handle derivatives. However the equations of general relativistic force-free electrodynamics (GRFFE) do not generally permit shocks, so a more optimal approach would involve finite differencing all derivatives in the GRFFE equations. As it happens, NRPy+ was designed to generate C codes involving complex tensorial expressions and finite difference derivatives, with finite-differencing order a freely-specifiable parameter.
#
# The purpose of this notebook is to rewrite the equations of GRFFE as used in the original $\texttt{GiRaFFE}$ code so that all derivatives that appear are represented numerically as finite difference derivatives. As we will see, the largest complication stems from derivatives of magnetic fields--requiring judicious application of the chain and product rules.
#
# The GRFFE evolution equations (from eq. 13 of the [original paper](https://arxiv.org/pdf/1704.00599.pdf)) we wish to encode in the NRPy+ version of $\texttt{GiRaFFE}$ are as follows:
#
# * $\partial_t \tilde{S}_i = - \partial_j \left( \alpha \sqrt{\gamma} T^j_{{\rm EM} i} \right) + \frac{1}{2} \alpha \sqrt{\gamma} T^{\mu \nu}_{\rm EM} \partial_i g_{\mu \nu}$: [Link to Step 3.0](#step7)
# * $\partial_t A_i = \epsilon_{ijk} v^j B^k - \partial_i (\alpha \Phi - \beta^j A_j)$: [Link to Step 4.0](#step8)
# * $\partial_t [\sqrt{\gamma} \Phi] = -\partial_j (\alpha\sqrt{\gamma}A^j - \beta^j [\sqrt{\gamma} \Phi]) - \xi \alpha [\sqrt{\gamma} \Phi]$: [Link to Step 4.0](#step8)
#
# Here, the densitized spatial Poynting flux one-form $\tilde{S}_i = \sqrt{\gamma} S_i$ (and $S_i$ comes from $S_{\mu} -n_{\nu} T^{\nu}_{{\rm EM} \mu}$), and $(\Phi, A_i)$ is the potential.
#
# ### A Note on Notation
#
# As is standard in NRPy+,
#
# * Greek indices refer to four-dimensional quantities where the zeroth component indicates $t$ components.
# * Latin indices refer to three-dimensional quantities. This is somewhat counterintuitive since Python always indexes its lists starting from 0. As a result, the zeroth component of three-dimensional quantities will necessarily indicate the first *spatial* direction.
#
# As a corollary, any expressions involving mixed Greek and Latin indices will need to offset one set of indices by one: A Latin index in a four-vector will be incremented and a Greek index in a three-vector will be decremented (however, the latter case does not occur in this tutorial module).
#
# ### Table of Contents:
#
# 1. Preliminaries
#     1. [Steps 1.0-1.2](#steps0to2): Set up the basic NRPy+ infrastructure we need
#         * This part imports the core NRPy+ modules that we need, sets the dimensionality of the grid with parameter $\text{grid::DIM}$, and declares the basic gridfunctions
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
# ## Steps 1.0-1.1: Set up the needed NRPy+ infrastructure and declare core gridfunctions used by $\texttt{GiRaFFE}$
#
# $$\label{steps0to2}$$
#
# \[Back to [top](#top)\]
#
# 1. Set some basic NRPy+ parameters. E.g., set the spatial dimension parameter to 3 and the finite differencing order to 4.
# 1. Next, declare some gridfunctions that are provided as input to the equations:
#     1. $\alpha$, $\beta^i$, and $\gamma_{ij}$: These ADM 3+1 metric quantities are declared in the ADMBase Einstein Toolkit thorn, and are assumed to be made available to $\texttt{GiRaFFE}$ at this stage.
#     1. The Valencia 3-velocity $v^i_{(n)}$ and vector potential $A_i$: Declared by $\texttt{GiRaFFE}$, and will have their initial values set in the separate thorn **GiRaFFEfood_HO**.
#     1. The magnetic field as measured by a normal observer $B^i$: The quantities evolved forward in time in $\texttt{GiRaFFE}$ do not include the Valencia 3-velocity, so this quantity is not automatically updated. Instead, we compute it based on the evolved quantity $\tilde{S}_i$ and $B^i = \epsilon^{ijk} \partial_j A_k$ (where $A_k$ is another evolved quantity and $\epsilon^{ijk}$ is the Levi-Civita tensor). $B^i$ is evaluated using finite differences of $A_k$ in a separate function, though it can only be evaluated consistently


import NRPy_param_funcs as par
import indexedexp as ixp
import grid as gri
import sympy as sp

def GiRaFFE_Higher_Order():
    #Step 1.0: Set the spatial dimension parameter to 3.
    par.set_parval_from_str("grid::DIM", 3)
    DIM = par.parval_from_str("grid::DIM")

    # Step 1.1: Set the finite differencing order to 4.
    #par.set_parval_from_str("finite_difference::FD_CENTDERIVS_ORDER", 4)

    thismodule = "GiRaFFE_NRPy"

    # M_PI will allow the C code to substitute the correct value
    M_PI = par.Cparameters("#define",thismodule,"M_PI","")
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
    # ## Step 1.2: Build the four metric $g_{\mu\nu}$, its inverse $g^{\mu\nu}$ and spatial derivatives $g_{\mu\nu,i}$ from ADM 3+1 quantities $\gamma_{ij}$, $\beta^i$, and $\alpha$
    #
    # $$\label{step3}$$
    # \[Back to [top](#top)\]
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


    # Step 1.2: import u0_smallb_Poynting__Cartesian.py to set
    #           the four metric and its inverse. This module also sets b^2 and u^0.
    import u0_smallb_Poynting__Cartesian.u0_smallb_Poynting__Cartesian as u0b
    u0b.compute_u0_smallb_Poynting__Cartesian(gammaDD,betaU,alpha,ValenciavU,BU)

    betaD = ixp.zerorank1()
    for i in range(DIM):
        for j in range(DIM):
            betaD[i] += gammaDD[i][j] * betaU[j]

    # We will now pull in the four metric and its inverse.
    import BSSN.ADMBSSN_tofrom_4metric as AB4m # NRPy+: ADM/BSSN <-> 4-metric conversions
    AB4m.g4DD_ito_BSSN_or_ADM("ADM")
    g4DD = AB4m.g4DD
    AB4m.g4UU_ito_BSSN_or_ADM("ADM")
    g4UU = AB4m.g4UU

    # Next we compute spatial derivatives of the metric, $\partial_i g_{\mu\nu} = g_{\mu\nu,i}$, written in terms of the three-metric, shift, and lapse. Simply taking the derivative of the expression for $g_{\mu\nu}$ above, we find
    # $$
    # g_{\mu\nu,l} = \begin{pmatrix}
    # -2\alpha \alpha_{,l} + \beta^k_{\ ,l} \beta_k + \beta^k \beta_{k,l} & \beta_{i,l} \\
    # \beta_{j,l} & \gamma_{ij,l}
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
    # Since this expression mixes Greek and Latin indices, we will need to store the expressions for each of the three spatial derivatives as separate variables.
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
    # ## Step 2.0: $u^i$ and $b^i$ and related quantities
    # $$\label{step4}$$
    # \[Back to [top](#top)\]
    #
    # We begin by importing what we can from u0_smallb_Poynting__Cartesian.py. We will need the four-velocity $u^\mu$, which is related to the Valencia 3-velocity $v^i_{(n)}$ used directly by $\texttt{GiRaFFE}$ (see also [Duez, et al, eqs. 53 and 56](https://arxiv.org/pdf/astro-ph/0503420.pdf))
    # \begin{align}
    # u^i &= u^0 (\alpha v^i_{(n)} - \beta^i), \\
    # u_j &= \alpha u^0 \gamma_{ij} v^i_{(n)},
    # \end{align}
    # and $v^i_{(n)}$ is the Valencia three-velocity. These have already been constructed in terms of the Valencia 3-velocity and other 3+1 ADM quantities by the u0_smallb_Poynting__Cartesian.py module, so we can simply import these variables:


    # Step 2.0: u^i, b^i, and related quantities
    # Step 2.0.a: import the four-velocity, as written in terms of the Valencia 3-velocity
    global uD,uU
    uD = ixp.register_gridfunctions_for_single_rank1("AUX","uD")
    uU = ixp.register_gridfunctions_for_single_rank1("AUX","uU")
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
    # ## Step 2.1: Construct all components of the electromagnetic stress-energy tensor $T^{\mu \nu}_{\rm EM}$
    # $$\label{step5}$$
    #
    # \[Back to [top](#top)\]
    #
    # We now have all the pieces to calculate the stress-energy tensor,
    # $$T^{\mu \nu}_{\rm EM} = \underbrace{b^2 u^{\mu} u^{\nu}}_{\rm Term\ 1} +
    # \underbrace{\frac{b^2}{2} g^{\mu \nu}}_{\rm Term\ 2}
    # - \underbrace{b^{\mu} b^{\nu}}_{\rm Term\ 3}.$$
    # Because $u^0$ is a separate variable, we could build the $00$ component separately, then the $\mu0$ and $0\nu$ components, and finally the $\mu\nu$ components. Alternatively, for clarity, we could create a temporary variable $u^\mu=\left( u^0, u^i \right)$


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


    # <a id='step6'></a>
    #
    # # Step 2.2: Derivatives of the electromagnetic stress-energy tensor
    # $$\label{step6}$$
    #
    # \[Back to [top](#top)\]
    #
    # If we look at the evolution equation, we see that we will need spatial  derivatives of $T^{\mu\nu}_{\rm EM}$. When confronted with derivatives of complicated expressions, it is generally convenient to declare those expressions as gridfunctions themselves, allowing NRPy+ to take finite-difference derivatives of the expressions. This can even reduce the truncation error associated with the finite differences, because the alternative is to use a function of several finite-difference derivatives, allowing more error to accumulate than the extra gridfunction will introduce. While we will use that technique for some of the subexpressions of $T^{\mu\nu}_{\rm EM}$, we don't want to rely on it for the whole expression; doing so would require us to take the derivative of the magnetic field $B^i$, which is itself found by finite-differencing the vector potential $A_i$. Thus $B^i$ cannot be *consistently* defined in ghost zones. To potentially reduce numerical errors induced by inconsistent finite differencing, we will differentiate $T^{\mu\nu}_{\rm EM}$ term-by-term so that finite-difference derivatives of $A_i$ appear.
    #
    # We will now now take these spatial derivatives of $T^{\mu\nu}_{\rm EM}$, applying the chain rule until it is only in terms of basic gridfunctions and their derivatives: $\alpha$, $\beta^i$, $\gamma_{ij}$, $A_i$, and the four-velocity $u^i$. Along the way, we will also set up useful temporary variables representing the steps of the chain rule. (Notably, *all* of these quantities will be written in terms of $A_i$ and its derivatives):
    #
    # * $B^i$ (already computed in terms of $A_k$, via $B^i = \epsilon^{ijk} \partial_j A_k$),
    # * $B^i_{,l}$,
    # * $b^i$ and $b_i$ (already computed),
    # * $b^i_{,k}$,
    # * $b^2$ (already computed),
    # * and $\left(b^2\right)_{,j}$.
    #
    # (The variables not already computed will not be seen by the ETK, as they are written in terms of $A_i$ and its derivatives; they simply help to organize the NRPy+ code.)
    #
    # So then,
    # \begin{align}
    # \partial_j T^{j}_{{\rm EM} i} &= \partial_j (g_{\mu i} T^{\mu j}_{\rm EM}) \\
    # &= \partial_j \left[g_{\mu i} \left(b^2 u^j u^\mu + \frac{b^2}{2} g^{j\mu} - b^j b^\mu\right)\right] \\
    # &= \underbrace{g_{\mu i,j} T^{\mu j}_{\rm EM}}_{\rm Term\ A} + g_{\mu i} \left( \underbrace{\partial_j \left(b^2 u^j u^\mu \right)}_{\rm Term\ B} + \underbrace{\partial_j \left(\frac{b^2}{2} g^{j\mu}\right)}_{\rm Term\ C} - \underbrace{\partial_j \left(b^j b^k\right)}_{\rm Term\ D} \right) \\
    # \end{align}
    # Following the product and chain rules for each term, we find that
    # \begin{align}
    # {\rm Term\ B} &= \partial_j (b^2 u^j u^\mu) \\
    #               &= \partial_j b^2 u^j u^\mu + b^2 \partial_j u^j u^\mu + b^2 u^j \partial_j u^\mu \\
    #               &= \underbrace{\left(b^2\right)_{,j} u^j u^\mu + b^2 u^j_{,j} u^\mu + b^2 u^j u^{\mu}_{,j}}_{\rm To\ Term\ 2\ below} \\
    # {\rm Term\ C} &= \partial_j \left(\frac{b^2}{2} g^{j\mu}\right) \\
    #               &= \frac{1}{2} \left( \partial_j b^2 g^{j\mu} + b^2 \partial_j g^{j\mu} \right) \\
    #               &= \underbrace{\frac{1}{2} \left(b^2\right)_{,j} g^{j\mu} + \frac{b^2}{2} g^{j\mu}_{\ ,j}}_{\rm To\ Term\ 3\ below} \\
    # {\rm Term\ D} &= \partial_j (b^j b^\mu) \\
    #               &= \underbrace{b^j_{,j} b^\mu + b^j b^\mu_{,j}}_{\rm To\ Term\ 2\ below}\\
    # \end{align}
    #
    # So,
    # \begin{align}
    # \partial_j T^{j}_{{\rm EM} i} &= g_{\mu i,j} T^{\mu j}_{\rm EM} \\
    # &+ g_{\mu i} \left(\left(b^2\right)_{,j} u^j u^\mu +b^2 u^j_{,j} u^\mu + b^2 u^j u^{\mu}_{,j} + \frac{1}{2}\left(b^2\right)_{,j} g^{j\mu} + \frac{b^2}{2} g^{j\mu}_{\ ,j} + b^j_{,j} b^\mu + b^j b^\mu_{,j}\right);
    # \end{align}
    # We will rearrange this once more, collecting the $b^2$ terms together, noting that Term A will become Term 1:
    # \begin{align}
    # \partial_j  T^{j}_{{\rm EM} i} =& \
    # \underbrace{g_{\mu i,j} T^{\mu j}_{\rm EM}}_{\rm Term\ 1} \\
    # & + \underbrace{g_{\mu i} \left( b^2 u^j_{,j} u^\mu + b^2 u^j u^\mu_{,j} + \frac{b^2}{2} g^{j\mu}_{\ ,j} + b^j_{,j} b^\mu + b^j b^\mu_{,j} \right)}_{\rm Term\ 2} \\
    # & + \underbrace{g_{\mu i} \left( \left(b^2\right)_{,j} u^j u^\mu + \frac{1}{2} \left(b^2\right)_{,j} g^{j\mu} \right).}_{\rm Term\ 3} \\
    # \end{align}
    #
    # <a id='table2'></a>
    #
    # **List of Derivatives**
    # $$\label{table2}$$
    #
    # Note that this is in terms of the derivatives of several other quantities:
    #
    # * [Step 2.2.a](#capitalBideriv): $B^i_{,l}$: Since $b^i$ is itself a function of $B^i$, we will first need the derivatives $B^i_{,l}$ in terms of the evolved quantity $A_i$ (the vector potential).
    # * [Step 2.2.b](#bideriv): $b^i_{,k}$: Once we have $B^i_{,l}$ we can evaluate derivatives of $b^i$, $b^i_{,k}$
    # * [Step 2.2.c](#b2deriv): The derivative of $b^2 = g_{\mu\nu} b^\mu b^\nu$, $\left(b^2\right)_{,j}$
    # * [Step 2.2.d](#gupijderiv): Derivatives of $g^{\mu\nu}$, $g^{\mu\nu}_{\ ,k}$
    # * [Step 2.2.e](#alltogether): Putting it together: $\partial_j T^{j}_{{\rm EM} i}$
    #     * [Step 2.2.e.i](#alltogether1): Putting it together: Term 1
    #     * [Step 2.2.e.ii](#alltogether2): Putting it together: Term 2
    #     * [Step 2.2.e.iii](#alltogether3): Putting it together: Term 3

    # <a id='capitalBideriv'></a>
    #
    # ## Step 2.2.a: Derivatives of $B^i$
    #
    # $$\label{capitalbideriv}$$
    #
    # \[Back to [List of Derivatives](#table2)\]
    #
    # First, we will build the derivatives of the magnetic field. Since $b^i$ is a function of $B^i$, we will start from the definition of $B^i$ in terms of $A_i$, $B^i = \frac{[ijk]}{\sqrt{\gamma}} \partial_j A_k$. We will first apply the product rule, noting that the symbol $[ijk]$ consists purely of the integers $-1, 0, 1$ and thus can be treated as a constant in this process.
    # \begin{align}
    # B^i_{,l} &= \partial_l \left( \frac{[ijk]}{\sqrt{\gamma}} \partial_j A_k \right)  \\
    #          &= [ijk] \partial_l \left( \frac{1}{\sqrt{\gamma}}\right) \partial_j A_k + \frac{[ijk]}{\sqrt{\gamma}} \partial_l \partial_j A_k \\
    #          &= [ijk]\left(-\frac{\gamma_{,l}}{2\gamma^{3/2}}\right) \partial_j A_k + \frac{[ijk]}{\sqrt{\gamma}} \partial_l \partial_j A_k \\
    # \end{align}
    # Now, we will substitute back in for the definition of the Levi-Civita tensor: $\epsilon^{ijk} = [ijk] / \sqrt{\gamma}$. Then we will substitute the magnetic field $B^i$ back in.
    # \begin{align}
    # B^i_{,l} &= -\frac{\gamma_{,l}}{2\gamma} \epsilon^{ijk} \partial_j A_k + \epsilon^{ijk} \partial_l \partial_j A_k \\
    #          &= -\frac{\gamma_{,l}}{2\gamma} B^i + \epsilon^{ijk} A_{k,jl}, \\
    # \end{align}
    #
    # Thus, the expression we are left with for the derivatives of the magnetic field is:
    # \begin{align}
    # B^i_{,l} &= \underbrace{-\frac{\gamma_{,l}}{2\gamma} B^i}_{\rm Term\ 1} + \underbrace{\epsilon^{ijk} A_{k,jl}}_{\rm Term\ 2}, \\
    # \end{align}
    # where $\epsilon^{ijk} = [ijk] / \sqrt{\gamma}$ is the antisymmetric Levi-Civita tensor and $\gamma$ is the determinant of the three-metric.
    #


    # Step 2.2: Derivatives of the electromagnetic stress-energy tensor
    ixp.register_gridfunctions_for_single_rank2("AUX","gammaUU","sym01")
    gri.register_gridfunctions("AUX","gammadet")
    gammaUU, gammadet = ixp.symm_matrix_inverter3x3(gammaDD)

    # We already have a handy function to define the Levi-Civita symbol in indexedexp.py
    # Initialize the Levi-Civita tensor by setting it equal to the Levi-Civita symbol
    LeviCivitaSymbolDDD = ixp.LeviCivitaSymbol_dim3_rank3()
    LeviCivitaTensorDDD = ixp.LeviCivitaTensorDDD_dim3_rank3(sp.sqrt(gammadet))
    LeviCivitaTensorUUU = ixp.LeviCivitaTensorUUU_dim3_rank3(sp.sqrt(gammadet))

    # AD_dD = ixp.declarerank2("AD_dD","nosym")

    # Step 2.2.a: Construct the derivatives of the magnetic field.
    gammadet_dD = ixp.declarerank1("gammadet_dD")

    AD_dDD = ixp.declarerank3("AD_dDD","sym12")
    # The other partial derivatives of B^i
    BUdD = ixp.zerorank2()
    for i in range(DIM):
        for l in range(DIM):
            # Term 1: -\gamma_{,l} / (2\gamma) B^i
            BUdD[i][l] = -gammadet_dD[l]*BU[i]/(2*gammadet)

    for i in range(DIM):
        for l in range(DIM):
            for j in range(DIM):
                for k in range(DIM):
                    # Term 2: \epsilon^{ijk} A_{k,jl}
                    BUdD[i][l] += LeviCivitaTensorUUU[i][j][k] * AD_dDD[k][j][l]


    # <a id='bideriv'></a>
    #
    # ## Step 2.2.b: Derivatives of $b^i$
    # $$\label{bideriv}$$
    #
    # \[Back to [List of Derivatives](#table2)\]
    #
    # Now, we will code the derivatives of the spatial components of $b^{\mu}$, $b^i$:
    # $$
    # b^i_{,k} = \frac{1}{\sqrt{4 \pi}} \frac{\left(\alpha u^0\right)  \left(B^i_{,k} + u_{j,k} B^j u^i + u_j B^j_{,k} u^i + u_j B^j u^i_{,k}\right) - \left(B^i + (u_j B^j) u^i\right) \partial_k \left(\alpha u^0\right)}{\left(\alpha u^0\right)^2}.
    # $$
    #
    # We should note that while $b^\mu$ is a four-vector (and the code reflects this: $\text{smallb4U}$ and $\text{smallb4U}$ have $\text{DIM=4}$), we only need the spatial components. We will only focus on the spatial components for the moment.
    #
    #
    # Let's go into a little more detail on where this comes from. We start from the definition $$b^i = \frac{B^i + (u_j B^j) u^i}{\sqrt{4\pi}\alpha u^0};$$ we then apply the quotient rule:
    # \begin{align}
    # b^i_{,k} &= \frac{\left(\sqrt{4\pi}\alpha u^0\right) \partial_k \left(B^i + (u_j B^j) u^i\right) - \left(B^i + (u_j B^j) u^i\right) \partial_k \left(\sqrt{4\pi}\alpha u^0\right)}{\left(\sqrt{4\pi}\alpha u^0\right)^2} \\
    # &= \frac{1}{\sqrt{4 \pi}} \frac{\left(\alpha u^0\right) \partial_k \left(B^i + (u_j B^j) u^i\right) - \left(B^i + (u_j B^j) u^i\right) \partial_k \left(\alpha u^0\right)}{\left(\alpha u^0\right)^2} \\
    # \end{align}
    # Note that $\left( \alpha u^0 \right)$ is being used as its own gridfunction, so $\partial_k \left(a u^0\right)$ will be finite-differenced by NRPy+ directly. We will also apply the product rule to the term $\partial_k \left(B^i + (u_j B^j) u^i\right) = B^i_{,k} + u_{j,k} B^j u^i + u_j B^j_{,k} u^i + u_j B^j u^i_{,k}$. So,
    # $$ b^i_{,k} = \frac{1}{\sqrt{4 \pi}} \frac{\left(\alpha u^0\right)  \left(B^i_{,k} + u_{j,k} B^j u^i + u_j B^j_{,k} u^i + u_j B^j u^i_{,k}\right) - \left(B^i + (u_j B^j) u^i\right) \partial_k \left(\alpha u^0\right)}{\left(\alpha u^0\right)^2}. $$
    #
    # It will be easier to code this up if we rearrange these terms to group together the terms that involve contractions over $j$. Doing that, we find
    # $$
    # b^i_{,k} = \frac{\overbrace{\alpha u^0 B^i_{,k} - B^i \partial_k (\alpha u^0)}^{\rm Term\ Num1} + \overbrace{\left( \alpha u^0 \right) \left( u_{j,k} B^j u^i + u_j B^j_{,k} u^i + u_j B^j u^i_{,k} \right)}^{\rm Term\ Num2.a} - \overbrace{\left( u_j B^j u^i \right) \partial_k \left( \alpha u^0 \right) }^{\rm Term\ Num2.b}}{\underbrace{\sqrt{4 \pi} \left( \alpha u^0 \right)^2}_{\rm Term\ Denom}}.
    # $$


    global u0alpha
    u0alpha = gri.register_gridfunctions("AUX","u0alpha")
    u0alpha = alpha * u4upperZero
    u0alpha_dD = ixp.declarerank1("u0alpha_dD")
    uU_dD = ixp.declarerank2("uU_dD","nosym")
    uD_dD = ixp.declarerank2("uD_dD","nosym")

    # Step 2.2.b: Construct derivatives of the small b vector
    # smallbUdD represents the derivative of smallb4U
    smallbUdD = ixp.zerorank2()
    for i in range(DIM):
        for k in range(DIM):
            # Term Num1: \alpha u^0 B^i_{,k} - B^i \partial_k (\alpha u^0)
            smallbUdD[i][k] += u0alpha*BUdD[i][k]-BU[i]*u0alpha_dD[k]

    for i in range(DIM):
        for k in range(DIM):
            for j in range(DIM):
                # Term Num2.a: terms that require contractions over k, and thus an extra loop.
                # ( \alpha u^0 ) (  u_{j,k} B^j u^i
                #                 + u_j B^j_{,k} u^i
                #                 + u_j B^j u^i_{,k} )
                smallbUdD[i][k] += u0alpha*(uD_dD[j][k]*BU[j]*uU[i]                                         +uD[j]*BUdD[j][k]*uU[i]                                         +uD[j]*BU[j]*uU_dD[i][k])

    for i in range(DIM):
        for k in range(DIM):
            for j in range(DIM):
                #Term 2.b (More contractions over k): ( u_j B^j u^i ) ( \alpha u^0 ),k
                smallbUdD[i][k] += -(uD[j]*BU[j]*uU[i])*u0alpha_dD[k]

    for i in range(DIM):
        for k in range(DIM):
            # Term Denom: Divide the numerator by sqrt(4 pi) * (alpha u^0)^2
            smallbUdD[i][k] /= sp.sqrt(4*M_PI) * u0alpha * u0alpha


    # <a id='b2deriv'></a>
    #
    # ## Step 2.2.c: Derivative of $b^2$
    # $$\label{b2deriv}$$
    #
    # \[Back to [List of Derivatives](#table2)\]
    #
    # Here, we will take the derivative of $b^2 = g_{\mu\nu} b^\mu b^\nu$. Using the product rule,
    # \begin{align}
    # \left(b^2\right)_{,j} &= \partial_j \left( g_{\mu\nu} b^\mu b^\nu \right) \\
    #                       &= g_{\mu\nu,j} b^\mu b^\nu + g_{\mu\nu} b^\mu_{,j} b^\nu + g_{\mu\nu} b^\mu b^\nu_{,j} \\
    #                       &= g_{\mu\nu,j} b^\mu b^\nu + 2 g_{\mu\nu} b^\mu_{,j} b^\nu.
    # \end{align}
    # We have already defined the spatial derivatives of the four-metric $g_{\mu\nu,j}$ in [this section](#step3); we have also defined the spatial derivatives of spatial components of $b^\mu$, $b^i_{,k}$ in [this section](#bideriv). Notice the above expression requires spatial derivatives of the *zeroth* component of $b^\mu$ as well, $b^0_{,j}$, which we will now compute. Starting with the definition, and applying the quotient rule:
    # \begin{align}
    # b^0 &= \frac{u_k B^k}{\sqrt{4\pi}\alpha}, \\
    # \rightarrow b^0_{,j} &= \frac{1}{\sqrt{4\pi}} \frac{\alpha \left( u_{k,j} B^k + u_k B^k_{,j} \right) - u_k B^k \alpha_{,j}}{\alpha^2} \\
    #     &= \frac{\alpha u_{k,j} B^k + \alpha u_k B^k_{,j} - \alpha_{,j} u_k B^k}{\sqrt{4\pi} \alpha^2}.
    # \end{align}
    # We will first code the numerator, and then divide through by the denominator.


    # Step 2.2.c: Construct the derivative of b^2
    # First construct the derivative b^0_{,j}
    # This four-vector will make b^2 simpler:
    smallb4UdD = ixp.zerorank2(DIM=4)
    # Fill in the zeroth component
    for j in range(DIM):
        for k in range(DIM):
            # The numerator:  \alpha u_{k,j} B^k
            #               + \alpha u_k B^k_{,j}
            #               - \alpha_{,j} u_k B^k
            smallb4UdD[0][j+1] +=   alpha*uD_dD[k][j]*BU[k]                                + alpha*uD[k]*BUdD[k][j]                                - alpha_dD[j]*uD[k]*BU[k]
    for j in range(DIM):
        # Divide through by the denominator: \sqrt{4\pi} \alpha^2
        smallb4UdD[0][j+1] /= sp.sqrt(4*M_PI)*alpha*alpha


    # At this point, both $b^0_{\ ,j}$ and $b^i_{\ ,j}$ have been computed, but one exists inconveniently in the $4\times 4$ component $\verb|smallb4UdD[][]|$ and the other in the $3\times 3$ component $\verb|smallbUdD[][]|$. So that we can perform full implied sums over $g_{\mu\nu} b^\mu_{,j} b^\nu$ more conveniently, we will now store all information from $\verb|smallbUdD[i][j]|$ into $\verb|smallb4UdD[i+1][j+1]|$:


    # Now, we'll fill out the rest of the four-vector with b^i_{,j} that we derived above.
    for i in range(DIM):
        for j in range(DIM):
            smallb4UdD[i+1][j+1] = smallbUdD[i][j]


    # Using 4-component (Greek-indexed) quantities, we can now complete our construction of
    # $$\left(b^2\right)_{,j} = g_{\mu\nu,j} b^\mu b^\nu + 2 g_{\mu\nu} b^\mu_{,j} b^\nu:$$


    smallb2_dD = ixp.zerorank1()
    for j in range(DIM):
        for mu in range(4):
            for nu in range(4):
                #   g_{\mu\nu,j} b^\mu b^\nu
                # + 2 g_{\mu\nu} b^\mu_{,j} b^\nu
                smallb2_dD[j] +=   g4DDdD[mu][nu][j+1]*smallb4U[mu]*smallb4U[nu]                              + 2*g4DD[mu][nu]*smallb4UdD[mu][j+1]*smallb4U[nu]


    # <a id='gupijderiv'></a>
    #
    # ## Step 2.2.d: Derivatives of $g^{\mu\nu}$
    # $$\label{gupijderiv}$$
    #
    # \[Back to [List of Derivatives](#table2)\]
    #
    # We will need derivatives of the inverse four-metric, as well. Let us begin with $g^{00}$: since $g^{00} = -1/\alpha^2$ ([Gourgoulhon, eq. 4.49](https://arxiv.org/pdf/gr-qc/0703035.pdf)), $$g^{00}_{\ ,k} = \frac{2 \alpha_{,k}}{\alpha^3}$$
    #


    # Step 2.2.d: Construct derivatives of the components of g^{\mu\nu}
    g4UUdD = ixp.zerorank3(DIM=4)

    for k in range(DIM):
        # 2 \alpha_{,k} / \alpha^3
        g4UUdD[0][0][k+1] = 2*alpha_dD[k]/alpha**3


    # Now, we will code the $g^{i0}_{\ ,k}$ and $g^{0j}_{\ ,k}$ components. According to [Gourgoulhon, eq. 4.49](https://arxiv.org/pdf/gr-qc/0703035.pdf), $g^{i0} = g^{0i} = \beta^i/\alpha^2$, so $$g^{i0}_{\ ,k} = g^{0i}_{\ ,k} = \frac{\alpha^2 \beta^i_{,k} - 2 \beta^i \alpha \alpha_{,k}}{\alpha^4}$$ by the quotient rule. So, we'll code
    # $$
    # g^{i0} = g^{0i} =
    # \underbrace{\frac{\beta^i_{,k}}{\alpha^2}}_{\rm Term\ 1}
    # - \underbrace{\frac{2 \beta^i \alpha_{,k}}{\alpha^3}}_{\rm Term\ 2}
    # $$


    for k in range(DIM):
        for i in range(DIM):
            # Term 1:                                    \beta^i_{,k} / \alpha^2
            g4UUdD[i+1][0][k+1] = g4UUdD[0][i+1][k+1] = betaU_dD[i][k] / alpha**2

    for k in range(DIM):
        for i in range(DIM):
            # Term 2:              -2 \beta^i \alpha_{,k} / \alpha^3
            g4UUdD[i+1][0][k+1] += -2 * betaU[i] * alpha_dD[k] / alpha**3
            g4UUdD[0][i+1][k+1] += -2 * betaU[i] * alpha_dD[k] / alpha**3


    # We will also need derivatives of the spatial part of the inverse four-metric: since $g^{ij} = \gamma^{ij} - \frac{\beta^i \beta^j}{\alpha^2}$ ([Gourgoulhon, eq. 4.49](https://arxiv.org/pdf/gr-qc/0703035.pdf)),
    # \begin{align}
    # g^{ij}_{\ ,k} &= \gamma^{ij}_{\ ,k} - \frac{\alpha^2 \partial_k (\beta^i \beta^j) - \beta^i \beta^j \partial_k \alpha^2}{(\alpha^2)^2} \\
    # &= \gamma^{ij}_{\ ,k} - \frac{\alpha^2\beta^i \beta^j_{,k}+\alpha^2\beta^i_{,k} \beta^j-2\beta^i \beta^j \alpha \alpha_{,k}}{\alpha^4}. \\
    # &= \gamma^{ij}_{\ ,k} - \frac{\alpha\beta^i \beta^j_{,k}+\alpha\beta^i_{,k} \beta^j-2\beta^i \beta^j \alpha_{,k}}{\alpha^3} \\
    # g^{ij}_{\ ,k} &= \underbrace{\gamma^{ij}_{\ ,k}}_{\rm Term\ 1} - \underbrace{\frac{\beta^i \beta^j_{,k}}{\alpha^2}}_{\rm Term\ 2} - \underbrace{\frac{\beta^i_{,k} \beta^j}{\alpha^2}}_{\rm Term\ 3} + \underbrace{\frac{2\beta^i \beta^j \alpha_{,k}}{\alpha^3}}_{\rm Term\ 4}. \\
    # \end{align}
    #


    gammaUU_dD = ixp.declarerank3("gammaUU_dD","sym01")

    # The spatial derivatives of the spatial components of the four metric:
    # Term 1: \gamma^{ij}_{\ ,k}
    for i in range(DIM):
        for j in range(DIM):
            for k in range(DIM):
                g4UUdD[i+1][j+1][k+1] = gammaUU_dD[i][j][k]

    # Term 2: - \beta^i \beta^j_{,k} / \alpha^2
    for i in range(DIM):
        for j in range(DIM):
            for k in range(DIM):
                g4UUdD[i+1][j+1][k+1] += -betaU[i]*betaU_dD[j][k]/alpha**2

    # Term 3: - \beta^i_{,k} \beta^j / \alpha^2
    for i in range(DIM):
        for j in range(DIM):
            for k in range(DIM):
                g4UUdD[i+1][j+1][k+1] += -betaU_dD[i][k]*betaU[j]/alpha**2

    # Term 4: 2\beta^i \beta^j \alpha_{,k}\alpha^3
    for i in range(DIM):
        for j in range(DIM):
            for k in range(DIM):
                g4UUdD[i+1][j+1][k+1] += 2*betaU[i]*betaU[j]*alpha_dD[k]/alpha**3


    # <a id='alltogether'></a>
    #
    # ## Step 2.2.e: Putting it together:
    # $$\label{alltogether}$$
    #
    # \[Back to [List of Derivatives](#table2)\]
    #
    # So, we can now put it all together, starting from the expression we derived above in [Step 2.2](#step6):
    # \begin{align}
    # \partial_j  T^{j}_{{\rm EM} i} =& \
    # \underbrace{g_{\mu i,j} T^{\mu j}_{\rm EM}}_{\rm Term\ 1} \\
    # & + \underbrace{g_{\mu i} \left( b^2 u^j_{,j} u^\mu + b^2 u^j u^\mu_{,j} + \frac{b^2}{2} g^{j\mu}_{\ ,j} + b^j_{,j} b^\mu + b^j b^\mu_{,j} \right)}_{\rm Term\ 2} \\
    # & + \underbrace{g_{\mu i} \left( \left(b^2\right)_{,j} u^j u^\mu + \frac{1}{2} \left(b^2\right)_{,j} g^{j\mu} \right).}_{\rm Term\ 3} \\
    # \end{align}
    #
    # <a id='alltogether1'></a>
    #
    # ### Step 2.2.e.i: Putting it together: Term 1
    # $$\label{alltogether1}$$
    #
    # \[Back to [List of Derivatives](#table2)\]
    #
    # We will now construct this term by term. Term 1 is straightforward: $${\rm Term\ 1} = \gamma_{\mu i,j} T^{\mu j}_{\rm EM}.$$


    # Step 2.2.e: Construct TEMUDdD_contracted itself
    # Step 2.2.e.i
    TEMUDdD_contracted = ixp.zerorank1()
    for i in range(DIM):
        for j in range(DIM):
            for mu in range(4):
                # Term 1:                g_{\mu i,j} T^{\mu j}_{\rm EM}
                TEMUDdD_contracted[i] += g4DDdD[mu][i+1][j+1] * T4EMUU[mu][j+1]

    # We'll need derivatives of u4U for the next part:
    u4UdD = ixp.zerorank2(DIM=4)
    u4upperZero_dD = ixp.declarerank1("u4upperZero_dD") # Note that derivatives can't be done in 4-D with the current version of NRPy
    for i in range(DIM):
        u4UdD[0][i+1] = u4upperZero_dD[i]
    for i in range(DIM):
        for j in range(DIM):
            u4UdD[i+1][j+1] = uU_dD[i][j]


    # <a id='alltogether2'></a>
    #
    # ### Step 2.2.e.ii: Putting it together: Term 2
    # $$\label{alltogether2}$$
    #
    # \[Back to [List of Derivatives](#table2)\]
    #
    # We will now add $${\rm Term\ 2} = g_{\mu i} \left( \underbrace{b^2 u^j_{,j} u^\mu}_{\rm Term\ 2a} + \underbrace{b^2 u^j u^\mu_{,j}}_{\rm Term\ 2b} + \underbrace{\frac{b^2}{2} g^{j\mu}_{\ ,j}}_{\rm Term\ 2c} + \underbrace{b^j_{,j} b^\mu}_{\rm Term\ 2d} + \underbrace{b^j b^\mu_{,j}}_{\rm Term\ 2e} \right)$$ to $\partial_j  T^{j}_{{\rm EM} i}$. These are the terms that involve contractions over $k$ (but no metric derivatives like Term 1 had).
    #


    # Step 2.2.e.ii
    for i in range(DIM):
        for j in range(DIM):
            for mu in range(4):
                # Term 2a: g_{\mu i} b^2 u^j_{,j} u^\mu
                TEMUDdD_contracted[i] += g4DD[mu][i+1]*smallb2*uU_dD[j][j]*u4U[mu]

    for i in range(DIM):
        for j in range(DIM):
            for mu in range(4):
                # Term 2b: g_{\mu i} b^2 u^j u^\mu_{,j}
                TEMUDdD_contracted[i] += g4DD[mu][i+1]*smallb2*uU[j]*u4UdD[mu][j+1]

    for i in range(DIM):
        for j in range(DIM):
            for mu in range(4):
                # Term 2c: g_{\mu i} b^2 g^{j \mu}_{,j} / 2
                TEMUDdD_contracted[i] += g4DD[mu][i+1]*smallb2*g4UUdD[j+1][mu][j+1]/2

    for i in range(DIM):
        for j in range(DIM):
            for mu in range(4):
                # Term 2d: g_{\mu i} b^j_{,j} b^\mu
                TEMUDdD_contracted[i] += g4DD[mu][i+1]*smallbUdD[j][j]*smallb4U[mu]

    for i in range(DIM):
        for j in range(DIM):
            for mu in range(4):
                # Term 2e: g_{\mu i} b^j b^\mu_{,j}
                TEMUDdD_contracted[i] += g4DD[mu][i+1]*smallb4U[j+1]*smallb4UdD[mu][j+1]


    # <a id='alltogether3'></a>
    #
    # ### Step 2.2.e.iii: Putting it together: Term 3
    # $$\label{alltogether3}$$
    #
    # \[Back to [List of Derivatives](#table2)\]
    #
    # Now, we will add $${\rm Term\ 3} = g_{\mu i} \left( \underbrace{\left(b^2\right)_{,j} u^j u^\mu}_{\rm Term\ 3a} + \underbrace{\frac{1}{2} \left(b^2\right)_{,j} g^{j\mu}}_{\rm Term\ 3b} \right).$$


    # Step 2.2.e.iii
    for i in range(DIM):
        for j in range(DIM):
            for mu in range(4):
                # Term 3a: g_{\mu i} ( b^2 )_{,j} u^j u^\mu
                TEMUDdD_contracted[i] += g4DD[mu][i+1]*smallb2_dD[j]*uU[j]*u4U[mu]

    for i in range(DIM):
        for j in range(DIM):
            for mu in range(4):
                # Term 3b: g_{mu i} ( b^2 )_{,j} g^{j\mu} / 2
                TEMUDdD_contracted[i] += g4DD[mu][i+1]*smallb2_dD[j]*g4UU[j+1][mu]/2


    #
    # # Evolution equation for $\tilde{S}_i$
    # <a id='step7'></a>
    #
    # ## Step 3.0: Construct the evolution equation for $\tilde{S}_i$
    # $$\label{step7}$$
    #
    # \[Back to [top](#top)\]
    #
    # Finally, we will return our attention to the time evolution equation (from eq. 13 of the [original paper](https://arxiv.org/pdf/1704.00599.pdf)),
    # \begin{align}
    # \partial_t \tilde{S}_i &= - \partial_j \left( \alpha \sqrt{\gamma} T^j_{{\rm EM} i} \right) + \frac{1}{2} \alpha \sqrt{\gamma} T^{\mu \nu}_{\rm EM} \partial_i g_{\mu \nu} \\
    #                        &= -T^j_{{\rm EM} i} \partial_j (\alpha \sqrt{\gamma}) - \alpha \sqrt{\gamma} \partial_j T^j_{{\rm EM} i} + \frac{1}{2} \alpha \sqrt{\gamma} T^{\mu \nu}_{\rm EM} \partial_i g_{\mu \nu} \\
    #                        &= \underbrace{-g_{i\mu} T^{\mu j}_{\rm EM} \partial_j (\alpha \sqrt{\gamma})
    # }_{\rm Term\ 1} - \underbrace{\alpha \sqrt{\gamma} \partial_j T^j_{{\rm EM} i}}_{\rm Term\ 2} + \underbrace{\frac{1}{2} \alpha \sqrt{\gamma} T^{\mu \nu}_{\rm EM} \partial_i g_{\mu \nu}}_{\rm Term\ 3} .
    # \end{align}
    # We will first take derivatives of $\alpha \sqrt{\gamma}$, then construct each term in turn.


    # Step 3.0: Construct the evolution equation for \tilde{S}_i
    # Here, we set up the necessary machinery to take FD derivatives of alpha * sqrt(gamma)
    global alpsqrtgam
    alpsqrtgam = gri.register_gridfunctions("AUX","alpsqrtgam")
    alpsqrtgam = alpha*sp.sqrt(gammadet)
    alpsqrtgam_dD = ixp.declarerank1("alpsqrtgam_dD")

    global Stilde_rhsD
    Stilde_rhsD = ixp.zerorank1()
    # The first term: g_{i\mu} T^{\mu j}_{\rm EM} \partial_j (\alpha \sqrt{\gamma})
    for i in range(DIM):
        for j in range(DIM):
            for mu in range(4):
                Stilde_rhsD[i] += -g4DD[i+1][mu]*T4EMUU[mu][j+1]*alpsqrtgam_dD[j]

    # The second term: \alpha \sqrt{\gamma} \partial_j T^j_{{\rm EM} i}
    for i in range(DIM):
        Stilde_rhsD[i] += -alpsqrtgam * TEMUDdD_contracted[i]

    # The third term: \alpha \sqrt{\gamma} T^{\mu \nu}_{\rm EM} \partial_i g_{\mu \nu} / 2
    for i in range(DIM):
        for mu in range(4):
            for nu in range(4):
                Stilde_rhsD[i] += alpsqrtgam * T4EMUU[mu][nu] * g4DDdD[mu][nu][i+1] / 2


    # # Evolution equations for $A_i$ and $\Phi$
    # <a id='step8'></a>
    #
    # ## Step 4.0: Construct the evolution equations for $A_i$ and $[\sqrt{\gamma}\Phi]$
    # $$\label{step8}$$
    #
    # \[Back to [top](#top)\]
    #
    # We will also need to evolve the vector potential $A_i$. This evolution is given as eq. 17 in the [$\texttt{GiRaFFE}$](https://arxiv.org/pdf/1704.00599.pdf) paper:
    # $$\boxed{\partial_t A_i = \epsilon_{ijk} v^j B^k - \partial_i (\underbrace{\alpha \Phi - \beta^j A_j}_{\rm AevolParen}),}$$
    # where $\epsilon_{ijk} = [ijk] \sqrt{\gamma}$ is the antisymmetric Levi-Civita tensor, the drift velocity $v^i = u^i/u^0$, $\gamma$ is the determinant of the three metric, $B^k$ is the magnetic field, $\alpha$ is the lapse, and $\beta$ is the shift.
    # The scalar electric potential $\Phi$ is also evolved by eq. 19:
    # $$\boxed{\partial_t [\sqrt{\gamma} \Phi] = -\partial_j (\underbrace{\alpha\sqrt{\gamma}A^j - \beta^j [\sqrt{\gamma} \Phi]}_{\rm PevolParenU[j]}) - \xi \alpha [\sqrt{\gamma} \Phi],}$$
    # with $\xi$ chosen as a damping factor.
    #
    # ### Step 4.0.a: Construct some useful auxiliary gridfunctions for the other evolution equations
    #
    # After declaring a  some needed quantities, we will also define the parenthetical terms (underbrace above) that we need to take derivatives of. That way, we can take finite-difference derivatives easily. Note that we use $A^j = \gamma^{ij} A_i$, while $A_i$ (with $\Phi$) is technically a four-vector; this is justified, however, since $n_\mu A^\mu = 0$, where $n_\mu$ is a normal to the hypersurface, $A^0=0$ (according to Sec. II, subsection C of [this paper](https://arxiv.org/pdf/1110.4633.pdf)).


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
            PevolParenU[j] += alpha * sp.sqrt(gammadet) * gammaUU[i][j] * AD[i]

    AevolParen_dD  = ixp.declarerank1("AevolParen_dD")
    PevolParenU_dD = ixp.declarerank2("PevolParenU_dD","nosym")


    # ### Step 4.0.b: Construct the actual evolution equations for $A_i$ and $[\sqrt{\gamma}\Phi]$
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

    for i in range(DIM):
        A_rhsD[i] = -AevolParen_dD[i]
        for j in range(DIM):
            for k in range(DIM):
                A_rhsD[i] += LeviCivitaTensorDDD[i][j][k]*(uU[j]/u4upperZero)*BU[k]

    psi6Phi_rhs = -xi*alpha*psi6Phi
    for j in range(DIM):
        psi6Phi_rhs += -PevolParenU_dD[j][j]
