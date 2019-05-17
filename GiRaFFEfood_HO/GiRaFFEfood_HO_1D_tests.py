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
# See [previous NRPy+ tutorial module](Tutorial-GiRaFFEfood_HO.ipynb) for more general detail on how this is used.
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
import indexedexp as ixp
import grid as gri
import finite_difference as fin
from outputC import *
import loop

import reference_metric as rfm
par.set_parval_from_str("reference_metric::CoordSystem","Cartesian")
rfm.reference_metric()

# Step 1a: Set commonly used parameters.
thismodule = "GiRaFFEfood_HO_1D"
# Set the spatial dimension parameter to 3.
par.set_parval_from_str("grid::DIM", 3)
DIM = par.parval_from_str("grid::DIM")


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


mu_AW,M_PI = par.Cparameters("REAL",thismodule,["mu_AW","M_PI"]) # The wave speed and pi in C

def GiRaFFEfood_HO_1D_tests():
    gammamu = (1-mu_AW**2)**(-1/2)

    # We'll use reference_metric.py to define x and y
    x = rfm.xxCart[0]
    y = rfm.xxCart[1]

    g_AW = sp.cos(5*M_PI*gammamu)/M_PI


    # Now, we can define the vector potential. We will create three copies of this variable, because the potential is uniquely defined in three zones. Data for $x \leq -0.1/\gamma_\mu$ shall be referred to as "left", data for $-0.1/\gamma_\mu \leq x \leq 0.1/\gamma_\mu$ as "center", and data for $x \geq 0.1/\gamma_\mu$ as "right".
    # 
    # Starting on the left, 
    # \begin{align}
    # A_x &= 0 \\
    # A_y &= \gamma_\mu x - 0.015 \\
    # A_z &= y - \gamma_\mu (1-\mu)x .
    # \end{align}


    AD = ixp.register_gridfunctions_for_single_rank1("AUX","AD")
    global AleftD,AcenterD,ArightD
    AleftD = ixp.zerorank1()

    AleftD[0] = sp.sympify(0)
    AleftD[1] = gammamu*x-0.015
    AleftD[2] = y-gammamu*(1-mu_AW)*x


    # In the the center,
    # \begin{align}
    # A_x &= 0 \\
    # A_y &= 1.15 \gamma_\mu x - 0.03g(x) \\
    # A_z &= y - \gamma_\mu (1-\mu)x .
    # \end{align}


    AcenterD = ixp.zerorank1()

    AcenterD[0] = sp.sympify(0)
    AcenterD[1] = 1.15*gammamu*x-0.03*g_AW
    AcenterD[2] = y-gammamu*(1-mu_AW)*x


    # And on the right,
    # \begin{align}
    # A_x &= 0 \\
    # A_y &= 1.3 \gamma_\mu x - 0.015 \\
    # A_z &= y - \gamma_\mu (1-\mu)x .
    # \end{align}


    ArightD = ixp.zerorank1()

    ArightD[0] = sp.sympify(0)
    ArightD[1] = 1.3*gammamu*x-0.015
    ArightD[2] = y-gammamu*(1-mu_AW)*x


    # <a id='step2'></a>
    # ### Set the vectors $B^i$ and $E^i$ for the velocity
    # 
    # Now, we will set the magnetic and electric fields that we will need to define the initial velocities. First, we need to define $$f(x)=1+\sin (5\pi x);$$ note that in the definition of $B^i$, we need $f(x')$ where $x'=\gamma_\mu x$.
    # $$\label{step2}$$


    xprime = gammamu*x
    f_AW = 1 + sp.sin(5*M_PI*xprime)


    # We will now set the magnetic field in the wave frame:
    # \begin{align}
    # B'^{x'}(x') = &\ 1.0,\ B'^y(x') = 1.0, \\
    # B'^z(x') = &\ \left \{ \begin{array}{lll} 1.0 & \mbox{if} & x' \leq -0.1 \\
    # 				1.0+0.15 f(x') & \mbox{if} & -0.1 \leq x' \leq 0.1 \\
    # 				1.3 & \mbox{if} & x' \geq 0.1 \end{array} \right. .
    # \end{align}
    # 


    BleftpU = ixp.zerorank1()
    BleftpU[0] = sp.sympify(1.0)
    BleftpU[1] = sp.sympify(1.0)
    BleftpU[2] = sp.sympify(1.0)

    BcenterpU = ixp.zerorank1()
    BcenterpU[0] = sp.sympify(1.0)
    BcenterpU[1] = sp.sympify(1.0)
    BcenterpU[2] = 1.0 + 0.15*f_AW

    BrightpU = ixp.zerorank1()
    BrightpU[0] = sp.sympify(1.0)
    BrightpU[1] = sp.sympify(1.0)
    BrightpU[2] = sp.sympify(1.3)


    # Now, we will set the electric field in the wave frame:
    # \begin{align}
    # E'^{x'}(x') &= -B'^z(0,x'), \\ 
    # E'^y(x') &= 0.0, \\ 
    # E'^z(x') &= 1.0  .
    # \end{align}


    EleftpU = ixp.zerorank1()
    EleftpU[0] = -BleftpU[2]
    EleftpU[1] = sp.sympify(0.0)
    EleftpU[2] = sp.sympify(1.0)

    EcenterpU = ixp.zerorank1()
    EcenterpU[0] = -BcenterpU[2]
    EcenterpU[1] = sp.sympify(0.0)
    EcenterpU[2] = sp.sympify(1.0)

    ErightpU = ixp.zerorank1()
    ErightpU[0] = -BrightpU[2]
    ErightpU[1] = sp.sympify(0.0)
    ErightpU[2] = sp.sympify(1.0)


    # Next, we must transform the the fields into the grid frame. We'll do the magnetic fields first.
    # \begin{align}
    #   B^x(0,x) = &\ B'^{x'}(\gamma_\mu x) , \\
    #   B^y(0,x) = &\ \gamma_\mu [ B'^y(\gamma_\mu x) - \mu E'^z(\gamma_\mu x) ] , \\ 
    #   B^z(0,x) = &\ \gamma_\mu [ B'^z(\gamma_\mu x) + \mu E'^y(\gamma_\mu x) ] , 
    # \end{align}
    # 


    BleftU = ixp.zerorank1()
    BleftU[0] = BleftpU[0]
    BleftU[1] = gammamu*(BleftpU[1]-mu_AW*EleftpU[2])
    BleftU[2] = gammamu*(BleftpU[2]+mu_AW*EleftpU[1])

    BcenterU = ixp.zerorank1()
    BcenterU[0] = BcenterpU[0]
    BcenterU[1] = gammamu*(BcenterpU[1]-mu_AW*EcenterpU[2])
    BcenterU[2] = gammamu*(BcenterpU[2]+mu_AW*EcenterpU[1])

    BrightU = ixp.zerorank1()
    BrightU[0] = BrightpU[0]
    BrightU[1] = gammamu*(BrightpU[1]-mu_AW*ErightpU[2])
    BrightU[2] = gammamu*(BrightpU[2]+mu_AW*ErightpU[1])


    # And now the electric fields:
    # \begin{align}
    #   E^x(0,x) = &\ E'^{x'}(\gamma_\mu x) , \\ 
    #   E^y(0,x) = &\ \gamma_\mu [ E'^y(\gamma_\mu x) + \mu B'^z(\gamma_\mu x) ] ,\\ 
    #   E^z(0,x) = &\ \gamma_\mu [ E'^z(\gamma_\mu x) - \mu B'^y(\gamma_\mu x) ],
    # \end{align}
    # 


    EleftU = ixp.zerorank1()
    EleftU[0] = EleftpU[0]
    EleftU[1] = gammamu*(EleftpU[1]+mu_AW*BleftpU[2])
    EleftU[2] = gammamu*(EleftpU[2]-mu_AW*BleftpU[1])

    EcenterU = ixp.zerorank1()
    EcenterU[0] = EcenterpU[0]
    EcenterU[1] = gammamu*(EcenterpU[1]+mu_AW*BcenterpU[2])
    EcenterU[2] = gammamu*(EcenterpU[2]-mu_AW*BcenterpU[1])

    ErightU = ixp.zerorank1()
    ErightU[0] = ErightpU[0]
    ErightU[1] = gammamu*(ErightpU[1]+mu_AW*BrightpU[2])
    ErightU[2] = gammamu*(ErightpU[2]-mu_AW*BrightpU[1])


    # <a id='step3'></a>
    # ### Calculate $v^i$
    # 
    # Now, we calculate $$\mathbf{v} = \frac{\mathbf{E} \times \mathbf{B}}{B^2},$$ which is equivalent to $$v^i = [ijk] \frac{E^j B^k}{B^2},$$ where $[ijk]$ is the Levi-Civita symbol and $B^2 = \gamma_{ij} B^i B^j$ is a trivial dot product in flat space.
    # $$\label{step3}$$


    import WeylScal4NRPy.WeylScalars_Cartesian as weyl
    LeviCivitaSymbolDDD = weyl.define_LeviCivitaSymbol_rank3()

    Bleft2 = BleftU[0]*BleftU[0] + BleftU[1]*BleftU[1] + BleftU[2]*BleftU[2]
    Bcenter2 = BcenterU[0]*BcenterU[0] + BcenterU[1]*BcenterU[1] + BcenterU[2]*BcenterU[2]
    Bright2 = BrightU[0]*BrightU[0] + BrightU[1]*BrightU[1] + BrightU[2]*BrightU[2]

    ValenciavU = ixp.register_gridfunctions_for_single_rank1("AUX","ValenciavU")

    global ValenciavleftU,ValenciavcenterU,ValenciavrightU
    ValenciavleftU = ixp.zerorank1()
    for i in range(DIM):
        for j in range(DIM):
            for k in range(DIM):
                ValenciavleftU[i] = LeviCivitaSymbolDDD[i][j][k] * EleftU[j] * BleftU[k] / Bleft2

    ValenciavcenterU = ixp.zerorank1()
    for i in range(DIM):
        for j in range(DIM):
            for k in range(DIM):
                ValenciavcenterU[i] = LeviCivitaSymbolDDD[i][j][k] * EcenterU[j] * BcenterU[k] / Bcenter2

    ValenciavrightU = ixp.zerorank1()
    for i in range(DIM):
        for j in range(DIM):
            for k in range(DIM):
                ValenciavrightU[i] = LeviCivitaSymbolDDD[i][j][k] * ErightU[j] * BrightU[k] / Bright2

