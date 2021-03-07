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
import Min_Max_and_Piecewise_Expressions as noif
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

M_PI  = par.Cparameters("#define",thismodule,["M_PI"], "")

#################
# Generic function for all 1D tests: Compute Ax,Ay,Az
def Axyz_func(Ax_func,Ay_func,Az_func, stagger_enable, **params):
    x = rfm.xx_to_Cart[0]
    y = rfm.xx_to_Cart[1]
    z = rfm.xx_to_Cart[2]
    # First Ax
    if stagger_enable:
        y += sp.Rational(1,2)*gri.dxx[1]
        z += sp.Rational(1,2)*gri.dxx[2]
    Ax = Ax_func(x,y,z, **params)
    # Then Ay
    if stagger_enable:
        x += sp.Rational(1,2)*gri.dxx[0]
        z += sp.Rational(1,2)*gri.dxx[2]
    Ay = Ay_func(x,y,z, **params)
    # Finally Az
    if stagger_enable:
        x += sp.Rational(1,2)*gri.dxx[0]
        y += sp.Rational(1,2)*gri.dxx[1]
    Az = Az_func(x,y,z, **params)

    return Ax,Ay,Az

#################
# Generic function for all 1D tests: Valencia 3-velocity from EU and BU
def compute_ValenciavU_from_EU_and_BU(EU, BU):
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

    ValenciavU = ixp.zerorank1()
    for i in range(3):
        for j in range(3):
            for k in range(3):
                ValenciavU[i] += LeviCivitaSymbolDDD[i][j][k] * EU[j] * BU[k] / B2

    return ValenciavU

########################################
########################################
## Degenerate Alfven wave:
def Ax_DAW(x,y,z, **params):
    return sp.sympify(0)

def Ay_DAW(x,y,z, **params):
    gammamu = params["gammamu"]
    bound = sp.Rational(1,10)/gammamu

    Ayleft = -sp.Rational(4,5)/M_PI
    Aycenter = -sp.Rational(4,5)/M_PI * h1_AW
    Ayright = sp.sympify(2)*(gammamu*x-sp.Rational(1,10))
    return noif.coord_leq_bound(x,-bound)*Ayleft\
           +noif.coord_greater_bound(x,-bound)*noif.coord_leq_bound(x,bound)*Aycenter\
           +noif.coord_greater_bound(x,bound)*Ayright

def Az_DAW(x,y,z, **params):
    gammamu = params["gammamu"]
    bound = sp.Rational(1,10)/gammamu

    Azleft = -sp.sympify(2)*(gammamu*x+sp.Rational(1,10))
    Azcenter = -sp.Rational(4,5)/M_PI * h2_AW
    Azright = -sp.Rational(4,5)/M_PI
    return noif.coord_leq_bound(x,-bound)*Azleft\
           +noif.coord_greater_bound(x,-bound)*noif.coord_leq_bound(x,bound)*Azcenter\
           +noif.coord_greater_bound(x,bound)*Azright

def ValenciavU_DAW(**params):
    gammamu = params["gammamu"]
    mu_AW   = params["mu_DAW"]

    x = rfm.xx_to_Cart[0]
    y = rfm.xx_to_Cart[1]
    z = rfm.xx_to_Cart[2]

    # ### Set the vectors $B^i$ and $E^i$ for the velocity
    #
    # Now, we will set the magnetic and electric fields that we will need to define the initial velocities. First, we need to define $$f(x)=1+\sin (5\pi x);$$ note that in the definition of $B^i$, we need $f(x')$ where $x'=\gamma_\mu x$.
    # $$\label{step2}$$
    xprime = gammamu*x
    bound = sp.Rational(1,10)

    phileft = sp.sympify(0)
    phicenter = sp.Rational(5,2)*M_PI*(xprime+sp.Rational(1,10))
    phiright = sp.Rational(1,2)*M_PI

    phi = noif.coord_leq_bound(xprime,-bound)*phileft\
         +noif.coord_greater_bound(xprime,-bound)*noif.coord_leq_bound(x,bound)*phicenter\
         +noif.coord_greater_bound(xprime,bound)*phiright

    # We will now set the magnetic field in the wave frame:
    # \begin{align}
    # B'^{x'}(x') = &\ 1.0,\ B'^y(x') = 1.0, \\
    # B'^z(x') = &\ \left \{ \begin{array}{lll} 1.0 & \mbox{if} & x' \leq -0.1 \\
    # 				1.0+0.15 f(x') & \mbox{if} & -0.1 \leq x' \leq 0.1 \\
    # 				1.3 & \mbox{if} & x' \geq 0.1 \end{array} \right. .
    # \end{align}
    #

    BpU = ixp.zerorank1()
    BpU[0] = sp.sympify(0)
    BpU[1] = sp.sympify(2)*sp.cos(phi)
    BpU[2] = sp.sympify(2)*sp.sin(phi)

    # Now, we will set the electric field in the wave frame:
    # \begin{align}
    # E'^{x'}(x') &= -B'^z(0,x'), \\
    # E'^y(x') &= 0.0, \\
    # E'^z(x') &= 1.0  .
    # \end{align}


    EpU = ixp.zerorank1()

    # Next, we must transform the fields into the grid frame. We'll do the magnetic fields first.
    # \begin{align}
    #   B^x(0,x) = &\ B'^{x'}(\gamma_\mu x) , \\
    #   B^y(0,x) = &\ \gamma_\mu [ B'^y(\gamma_\mu x) - \mu E'^z(\gamma_\mu x) ] , \\
    #   B^z(0,x) = &\ \gamma_\mu [ B'^z(\gamma_\mu x) + \mu E'^y(\gamma_\mu x) ] ,
    # \end{align}
    #
    global BU
    BU = ixp.zerorank1()
    BU[0] = BpU[0]
    BU[1] = gammamu*(BpU[1]-mu_AW*EpU[2])
    BU[2] = gammamu*(BpU[2]+mu_AW*EpU[1])


    # And now the electric fields:
    # \begin{align}
    #   E^x(0,x) = &\ E'^{x'}(\gamma_\mu x) , \\
    #   E^y(0,x) = &\ \gamma_\mu [ E'^y(\gamma_\mu x) + \mu B'^z(\gamma_\mu x) ] ,\\
    #   E^z(0,x) = &\ \gamma_\mu [ E'^z(\gamma_\mu x) - \mu B'^y(\gamma_\mu x) ],
    # \end{align}
    #
    EU = ixp.zerorank1()
    EU[0] = EpU[0]
    EU[1] = gammamu*(EpU[1]+mu_AW*BpU[2])
    EU[2] = gammamu*(EpU[2]-mu_AW*BpU[1])

    return compute_ValenciavU_from_EU_and_BU(EU, BU)
########################################
########################################
## Fast wave:
def Ax_FW(x,y,z, **params):
    return sp.sympify(0)

def Ay_FW(x,y,z, **params):
    return sp.sympify(0)

def Az_FW(x,y,z, **params):
    bound = sp.Rational(1,10)

    # A_x = 0, A_y = 0
    # A_z = y+ (-x-0.0075) if x <= -0.1
    #          (0.75x^2 - 0.85x) if -0.1 < x <= 0.1
    #          (-0.7x-0.0075) if x > 0.1
    Azleft = y - x - sp.Rational(75,10000)
    Azcenter = y + sp.Rational(75,100)*x*x - sp.Rational(85,100)*x
    Azright = y - sp.Rational(7,10)*x - sp.Rational(75,10000)

    return noif.coord_leq_bound(x,-bound)*Azleft\
           +noif.coord_greater_bound(x,-bound)*noif.coord_leq_bound(x,bound)*Azcenter\
           +noif.coord_greater_bound(x,bound)*Azright

def ValenciavU_FW(**params):
    # B^x(0,x) = 1.0
    # B^y(0,x) = 1.0 if x <= -0.1
    #            1.0-1.5(x+0.1) if -0.1 < x <= 0.1
    #            0.7 if x > 0.1
    # B^z(0,x) = 0

    Byleft = sp.sympify(1)
    Bycenter = sp.sympify(1) - sp.Rational(15,10)*(x+sp.Rational(1,10))
    Byright = sp.Rational(7,10)

    global BU
    BU = ixp.zerorank1()
    BU[0] = sp.sympify(1)
    BU[1] = noif.coord_leq_bound(x,-bound)*Byleft\
            +noif.coord_greater_bound(x,-bound)*noif.coord_leq_bound(x,bound)*Bycenter\
            +noif.coord_greater_bound(x,bound)*Byright
    BU[2] = sp.sympify(0)

    # E^x(0,x) = 0.0 , E^y(x) = 0.0 , E^z(x) = -B^y(0,x)
    EU = ixp.zerorank1()
    EU[0] = sp.sympify(0)
    EU[1] = sp.sympify(0)
    EU[2] = -BU[1]

    return compute_ValenciavU_from_EU_and_BU(EU, BU)
########################################
########################################

def GiRaFFEfood_NRPy_1D_tests(ID_type = "DegenAlfvenWave", stagger_enable = False):
    global AD, ValenciaVU
    AD = ixp.zerorank1()
    if ID_type == "DegenAlfvenWave":
        mu_DAW = par.Cparameters("REAL",thismodule,["mu_DAW"], -0.5) # The wave speed
        AD = Axyz_func(Ax_DAW, Ay_DAW, Az_DAW, stagger_enable,
                       gammamu=sp.sympify(1)/sp.sqrt(sp.sympify(1)-mu_AW**2))
        ValenciaVU = ValenciavU_DAW(mu_DAW=mu_DAW, gammamu=gammamu)
    elif ID_type == "FastWave":
        AD = Axyz_func(Ax_FW, Ay_FW, Az_FW, stagger_enable)
        ValenciaVU = ValenciavU_FW()
