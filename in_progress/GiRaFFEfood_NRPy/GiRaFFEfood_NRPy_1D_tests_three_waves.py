#!/usr/bin/env python
# coding: utf-8

# <a id='top'></a>
#
#
# # $\texttt{GiRaFFEfood}$: Initial data for $\texttt{GiRaFFE}$
#
# ## Alfv&eacute;n Wave
#
# $$\label{top}$$
#
# This module provides another initial data option for $\texttt{GiRaFFE}$, drawn from [this paper](https://arxiv.org/abs/1310.3274) . This is a flat-spacetime test with initial data
# \begin{align}
# A_x &= 0 \\
# A_y &= \left \{ \begin{array}{lll}\gamma_\mu x - 0.015 & \mbox{if} & x \leq -0.1/\gamma_\mu \\
# 1.15 \gamma_\mu x - 0.03g(x) & \mbox{if} & -0.1/\gamma_\mu \leq x \leq 0.1/\gamma_\mu \\
# 1.3 \gamma_\mu x - 0.015 & \mbox{if} & x \geq 0.1/\gamma_\mu \end{array} \right. , \\
#  A_z = &\ y - \gamma_\mu (1-\mu)x ,
# \end{align}
# which generates the magnetic field in the wave frame,
# \begin{align}
# B'^{x'}(x') = &\ 1.0,\ B'^y(x') = 1.0, \\
# B'^z(x') = &\ \left \{ \begin{array}{lll} 1.0 & \mbox{if} & x' \leq -0.1 \\
# 				1.0+0.15 f(x') & \mbox{if} & -0.1 \leq x' \leq 0.1 \\
# 				1.3 & \mbox{if} & x' \geq 0.1 \end{array} \right. ,
# \end{align}
# and the electric field in the wave frame,
# $$E'^{x'}(x') = -B'^z(0,x') \ \ , \ \ E'^y(x') = 0.0 \ \ , \ \ E'^z(x') = 1.0  .$$
#
# These are converted to the grid frame by
# \begin{align}
#   B^x(0,x) = &\ B'^{x'}(\gamma_\mu x) , \\
#   B^y(0,x) = &\ \gamma_\mu [ B'^y(\gamma_\mu x) - \mu E'^z(\gamma_\mu x) ] , \\
#   B^z(0,x) = &\ \gamma_\mu [ B'^z(\gamma_\mu x) + \mu E'^y(\gamma_\mu x) ] ,
# \end{align}
# and
# \begin{align}
#   E^x(0,x) = &\ E'^{x'}(\gamma_\mu x) , \\
#   E^y(0,x) = &\ \gamma_\mu [ E'^y(\gamma_\mu x) + \mu B'^z(\gamma_\mu x) ] ,\\
#   E^z(0,x) = &\ \gamma_\mu [ E'^z(\gamma_\mu x) - \mu B'^y(\gamma_\mu x) ],
# \end{align}
# and the velocity is given by $$\mathbf{v} = \frac{\mathbf{E} \times \mathbf{B}}{B^2}$$ in flat spacetime. Additionally, $f(x)=1+\sin (5\pi x)$, $-1<\mu<1$ is the wave speed relative to the grid frame and $\gamma_\mu = (1-\mu^2)^{-1/2}$, and $g(x) = \cos (5\pi \gamma_\mu x)/\pi$.
#
# For the eventual purpose of testing convergence, any quantity $Q$ evolves as $Q(t,x) = Q(0,x-\mu t)$
#
# See [previous NRPy+ tutorial module](Tutorial-GiRaFFEfood_NRPy.ipynb) for more general detail on how this is used.
#
# #### Table of Contents:
# 1. [Steps 0-1:](#preliminaries) Preliminaries
# 1. [Step 2:](#step2) Set the vector $A_k$
# 1. [Step 3:](#step3) Set the vectors $B^i$ and $E^i$ for the velocity
# 1. [Step 4:](#step4) Calculate $v^i$
# 1. [Step 5:](#step6) NRPy+ Module Code Validation
#

# <a id='preliminaries'></a>
#
# ### Steps 0-1: Preliminaries \[Back to [top](#top)\]
#
# Here, we will import the NRPy+ core modules and set the reference metric to Cartesian, set commonly used NRPy+ parameters, and set C parameters that will be set from outside the code eventually generated from these expressions. We will also set up a parameter to determine what initial data is set up, although it won't do much yet.
# $$\label{preliminaries}$$


# Step 0: Import the NRPy+ core modules and set the reference metric to Cartesian
import NRPy_param_funcs as par
import grid as gri               # NRPy+: Functions having to do with numerical grids
import indexedexp as ixp
import sympy as sp               # SymPy: The Python computer algebra package upon which NRPy+ depends
import reference_metric as rfm
par.set_parval_from_str("reference_metric::CoordSystem","Cartesian")
rfm.reference_metric()

# Step 1a: Set commonly used parameters.
thismodule = __name__


# <a id='step2'></a>
# ### Set the vector $A_k$
# The vector potential is given as
# \begin{align}
# A_x &= 0 \\
# A_y &= \left \{ \begin{array}{lll}\gamma_\mu x - 0.015 & \mbox{if} & x \leq -0.1/\gamma_\mu \\
# 1.15 \gamma_\mu x - 0.03g(x) & \mbox{if} & -0.1/\gamma_\mu \leq x \leq 0.1/\gamma_\mu \\
# 1.3 \gamma_\mu x - 0.015 & \mbox{if} & x \geq 0.1/\gamma_\mu \end{array} \right. , \\
# A_z &= y - \gamma_\mu (1-\mu)x .
# \end{align}
# First, however, we must set $$\gamma_\mu = (1-\mu^2)^{-1/2}$$ and $$g(x) = \cos (5\pi \gamma_\mu x)/\pi$$.
# $$\label{step2}$$

def GiRaFFEfood_NRPy_1D_tests_three_waves(stagger = False):

    # We'll use reference_metric.py to define x and y
    x = rfm.xx_to_Cart[0]
    y = rfm.xx_to_Cart[1]
    if stagger:
        x_p_half = x + sp.Rational(1,2)*gri.dxx[0]
        y_p_half = y + sp.Rational(1,2)*gri.dxx[1]

    # Now, we can define the vector potential. We will create three copies of this variable, because the potential is uniquely defined in three zones. Data for $x \leq -0.1/\gamma_\mu$ shall be referred to as "left", data for $-0.1/\gamma_\mu \leq x \leq 0.1/\gamma_\mu$ as "center", and data for $x \geq 0.1/\gamma_\mu$ as "right".

    global AD
    AD = ixp.zerorank1(DIM=3)

    import Min_Max_and_Piecewise_Expressions as noif

    AD[0] = sp.sympify(0)
    if stagger:
        AD[1] = sp.Rational(7,2)*x_p_half*noif.coord_greater_bound(-x_p_half,0) + sp.sympify(3)*x_p_half*noif.coord_greater_bound(x_p_half,0)
        AD[2] = y_p_half-sp.Rational(3,2)*x_p_half*noif.coord_greater_bound(-x_p_half,0) - sp.sympify(3)*x_p_half*noif.coord_greater_bound(x_p_half,0)
    else:
        AD[1] = sp.Rational(7,2)*x*noif.coord_greater_bound(-x,0) + sp.sympify(3)*x*noif.coord_greater_bound(x,0)
        AD[2] = y-sp.Rational(3,2)*x*noif.coord_greater_bound(-x,0) - sp.sympify(3)*x*noif.coord_greater_bound(x,0)

    # ### Set the vectors $B^i$ and $E^i$ for the velocity
    #
    # Now, we will set the magnetic and electric fields that we will need to define the initial velocities. First, we need to define $$f(x)=1+\sin (5\pi x);$$ note that in the definition of $B^i$, we need $f(x')$ where $x'=\gamma_\mu x$.
    # $$\label{step2}$$


    # We will now set the magnetic field in the wave frame:
    # \begin{align}
    # B'^{x'}(x') = &\ 1.0,\ B'^y(x') = 1.0, \\
    # B'^z(x') = &\ \left \{ \begin{array}{lll} 1.0 & \mbox{if} & x' \leq -0.1 \\
    # 				1.0+0.15 f(x') & \mbox{if} & -0.1 \leq x' \leq 0.1 \\
    # 				1.3 & \mbox{if} & x' \geq 0.1 \end{array} \right. .
    # \end{align}
    #

    B_aU = ixp.zerorank1(DIM=3)
    E_aU = ixp.zerorank1(DIM=3)
    B_pU = ixp.zerorank1(DIM=3)
    E_pU = ixp.zerorank1(DIM=3)
    B_mU = ixp.zerorank1(DIM=3)
    E_mU = ixp.zerorank1(DIM=3)

#     if stagger:
#         B_aU[0] = sp.sympify(1)
#         B_aU[1] = noif.coord_leq_bound(x_p_half,0) * sp.sympify(1) + noif.coord_greater_bound(x_p_half,0) * sp.Rational(3,2)
#         B_aU[2] = sp.sympify(2)

#         E_aU[0] = noif.coord_leq_bound(x_p_half,0) * sp.sympify(-1) + noif.coord_greater_bound(x_p_half,0) * sp.Rational(-3,2)
#         E_aU[1] = sp.sympify(1)
#         E_aU[2] = sp.sympify(0)

#         B_pU[0] = sp.sympify(0)
#         B_pU[1] = noif.coord_leq_bound(x_p_half,0) * sp.sympify(0) + noif.coord_greater_bound(x_p_half,0) * sp.Rational(3,2)
#         B_pU[2] = noif.coord_leq_bound(x_p_half,0) * sp.sympify(0) + noif.coord_greater_bound(x_p_half,0) * sp.sympify(1)

#         E_pU[0] = sp.sympify(0)
#         E_pU[1] = noif.coord_leq_bound(x_p_half,0) * sp.sympify(0) + noif.coord_greater_bound(x_p_half,0) * sp.sympify(1)
#         E_pU[2] = noif.coord_leq_bound(x_p_half,0) * sp.sympify(0) + noif.coord_greater_bound(x_p_half,0) * sp.Rational(-3,2)

#         B_mU[0] = sp.sympify(0)
#         B_mU[1] = noif.coord_leq_bound(x_p_half,0) * sp.Rational(1,2)  + noif.coord_greater_bound(x_p_half,0) * sp.sympify(0)
#         B_mU[2] = noif.coord_leq_bound(x_p_half,0) * sp.Rational(3,2) + noif.coord_greater_bound(x_p_half,0) * sp.sympify(0)

#         E_mU[0] = sp.sympify(0)
#         E_mU[1] = noif.coord_leq_bound(x_p_half,0) * sp.Rational(-3,2) + noif.coord_greater_bound(x_p_half,0) * sp.sympify(0)
#         E_mU[2] = noif.coord_leq_bound(x_p_half,0) * sp.Rational(1,2)  + noif.coord_greater_bound(x_p_half,0) * sp.sympify(0)
#     else:
    B_aU[0] = sp.sympify(1)
    B_aU[1] = noif.coord_leq_bound(x,0) * sp.sympify(1) + noif.coord_greater_bound(x,0) * sp.Rational(3,2)
    B_aU[2] = sp.sympify(2)

    E_aU[0] = noif.coord_leq_bound(x,0) * sp.sympify(-1) + noif.coord_greater_bound(x,0) * sp.Rational(-3,2)
    E_aU[1] = sp.sympify(1)
    E_aU[2] = sp.sympify(0)

    B_pU[0] = sp.sympify(0)
    B_pU[1] = noif.coord_leq_bound(x,0) * sp.sympify(0) + noif.coord_greater_bound(x,0) * sp.Rational(3,2)
    B_pU[2] = noif.coord_leq_bound(x,0) * sp.sympify(0) + noif.coord_greater_bound(x,0) * sp.sympify(1)

    E_pU[0] = sp.sympify(0)
    E_pU[1] = noif.coord_leq_bound(x,0) * sp.sympify(0) + noif.coord_greater_bound(x,0) * sp.sympify(1)
    E_pU[2] = noif.coord_leq_bound(x,0) * sp.sympify(0) + noif.coord_greater_bound(x,0) * sp.Rational(-3,2)

    B_mU[0] = sp.sympify(0)
    B_mU[1] = noif.coord_leq_bound(x,0) * sp.Rational(1,2)  + noif.coord_greater_bound(x,0) * sp.sympify(0)
    B_mU[2] = noif.coord_leq_bound(x,0) * sp.Rational(3,2) + noif.coord_greater_bound(x,0) * sp.sympify(0)

    E_mU[0] = sp.sympify(0)
    E_mU[1] = noif.coord_leq_bound(x,0) * sp.Rational(-3,2) + noif.coord_greater_bound(x,0) * sp.sympify(0)
    E_mU[2] = noif.coord_leq_bound(x,0) * sp.Rational(1,2)  + noif.coord_greater_bound(x,0) * sp.sympify(0)

    global BU
    BU = ixp.zerorank1(DIM=3)
    EU = ixp.zerorank1(DIM=3)
    for i in range(3):
        BU[i] = B_aU[i] + B_pU[i] + B_mU[i]
        EU[i] = E_aU[i] + E_pU[i] + E_mU[i]
# <a id='step3'></a>
    # ### Calculate $v^i$
    #
    # Now, we calculate $$\mathbf{v} = \frac{\mathbf{E} \times \mathbf{B}}{B^2},$$ which is equivalent to $$v^i = [ijk] \frac{E^j B^k}{B^2},$$ where $[ijk]$ is the Levi-Civita symbol and $B^2 = \gamma_{ij} B^i B^j$ is a trivial dot product in flat space.
    # $$\label{step3}$$


    LeviCivitaSymbolDDD = ixp.LeviCivitaSymbol_dim3_rank3()

    B2 = sp.sympify(0)
    for i in range(3):
        # In flat spacetime, gamma_{ij} is just a Kronecker delta
        B2 += BU[i]**2 # This is trivial to extend to curved spacetime

    global ValenciavU
    ValenciavU = ixp.zerorank1()
    for i in range(3):
        for j in range(3):
            for k in range(3):
                ValenciavU[i] += LeviCivitaSymbolDDD[i][j][k] * EU[j] * BU[k] / B2
