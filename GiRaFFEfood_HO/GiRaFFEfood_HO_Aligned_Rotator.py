#!/usr/bin/env python
# coding: utf-8

# <a id='top'></a>
# 
# 
# # $\texttt{GiRaFFEfood}$: Initial data for $\texttt{GiRaFFE}$
# 
# ## Aligned Rotator
# 
# $$\label{top}$$
# 
# This module provides another initial data option for $\texttt{GiRaFFE}$. This is a flat-spacetime test with initial data $$A_{\phi} = \frac{\mu \varpi}{r^3},$$ where  $\mu = B_p R_{\rm NS} / 2$, $R_{\rm NS}$ is the neutron star radius, and $\varpi = \sqrt{x^2+y^2}$ is the cylindrical radius. We let $A_r = A_\theta = 0$.
# 
# Additionally, the drift velocity $v^i = \Omega \textbf{e}_z \times \textbf{r} = [ijk] \Omega \textbf{e}^j_z x^k$, where $[ijk]$ is the Levi-Civita permutation symbol and $\textbf{e}^i_z = (0,0,1)$.

# <a id='preliminaries'></a>
# 
# ### Steps 0-1: Preliminaries
# $$\label{preliminaries}$$
# 
# \[Back to [top](#top)\]
# 
# Here, we will import the NRPy+ core modules and set the reference metric to Cartesian, set commonly used NRPy+ parameters, and set C parameters that will be set from outside the code eventually generated from these expressions. We will also set up a parameter to determine what initial data is set up, although it won't do much yet.


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
thismodule = "GiRaFFEfood_HO_Aligned_Rotator"
# Set the spatial dimension parameter to 3.
par.set_parval_from_str("grid::DIM", 3)
DIM = par.parval_from_str("grid::DIM")

B_p_aligned_rotator,R_NS_aligned_rotator = par.Cparameters("REAL",thismodule,["B_p_aligned_rotator","R_NS_aligned_rotator"]) # A constant defining the intensity of the magnetic field and the Neutron star Radius


# <a id='step2'></a>
# 
# ### Step 2: Set the vectors A in Spherical coordinates
# $$\label{step2}$$
# 
# \[Back to [top](#top)\]
# 
# We will first build the fundamental vector $A_i$ in spherical coordinates (see [Table 3](https://arxiv.org/pdf/1704.00599.pdf)). Note that we use reference_metric.py to set $r$ and $\theta$ in terms of Cartesian coordinates; this will save us a step later when we convert to Cartesian coordinates. So, we set 
# \begin{align}
# A_{\phi} &= \frac{\mu \varpi}{r^3}, \\
# \end{align}
# with $\mu = B_p R_{\rm NS} / 2$, $R_{\rm NS}$ is the neutron star radius, and $\varpi = \sqrt{x^2+y^2}$

def GiRaFFEfood_HO_Aligned_Rotator():
    r     = rfm.xxSph[0]
    varpi = sp.sqrt(rfm.xxCart[0]**2 + rfm.xxCart[1]**2)

    mu = B_p_aligned_rotator * R_NS_aligned_rotator**3 / 2

    ASphD = ixp.zerorank1()

    ASphD[2] = mu * varpi**2 / (r**3) # The other components were already declared to be 0.


    # <a id='step3'></a>
    # 
    # ### Step 3: Use the Jacobian matrix to transform the vectors to Cartesian coordinates.
    # $$\label{step3}$$
    # 
    # \[Back to [top](#top)\]
    # 
    # Now, we will use the coordinate transformation definitions provided by reference_metric.py to build the Jacobian 
    # $$ 
    # \frac{\partial x_{\rm Sph}^j}{\partial x_{\rm Cart}^i},
    # $$ 
    # where $x_{\rm Sph}^j \in \{r,\theta,\phi\}$ and $x_{\rm Cart}^i \in \{x,y,z\}$. We would normally compute its inverse, but since none of the quantities we need to transform have upper indices, it is not necessary. Then, since $A_i$ and has one lower index, it will need to be multiplied by the Jacobian:
    # 
    # $$
    # A_i^{\rm Cart} = A_j^{\rm Sph} \frac{\partial x_{\rm Sph}^j}{\partial x_{\rm Cart}^i},
    # $$


    # Step 3: Use the Jacobian matrix to transform the vectors to Cartesian coordinates.
    drrefmetric__dx_0UDmatrix = sp.Matrix([[sp.diff(rfm.xxSph[0],rfm.xx[0]), sp.diff(rfm.xxSph[0],rfm.xx[1]), sp.diff(rfm.xxSph[0],rfm.xx[2])],
                                           [sp.diff(rfm.xxSph[1],rfm.xx[0]), sp.diff(rfm.xxSph[1],rfm.xx[1]), sp.diff(rfm.xxSph[1],rfm.xx[2])],
                                           [sp.diff(rfm.xxSph[2],rfm.xx[0]), sp.diff(rfm.xxSph[2],rfm.xx[1]), sp.diff(rfm.xxSph[2],rfm.xx[2])]])
    #dx__drrefmetric_0UDmatrix = drrefmetric__dx_0UDmatrix.inv() # We don't actually need this in this case.

    global AD
    AD = ixp.register_gridfunctions_for_single_rank1("EVOL","AD")
    for i in range(DIM):
        for j in range(DIM):
            AD[i] = drrefmetric__dx_0UDmatrix[(j,i)]*ASphD[j]


    # <a id='step4'></a>
    # 
    # ### Step 4: Calculate $v^i$
    # $$\label{step4}$$
    # 
    # \[Back to [top](#top)\]
    # 
    # Here, we will calculate the drift velocity $v^i = \Omega \textbf{e}_z \times \textbf{r} = [ijk] \Omega \textbf{e}^j_z x^k$, where $[ijk]$ is the Levi-Civita permutation symbol and $\textbf{e}^i_z = (0,0,1)$. Conveniently, in flat space, the drift velocity reduces to the Valencia velocity because $\alpha = 1$ and $\beta^i = 0$.


    # Step 4: Calculate v^i
    # Here, we build the Levi-Civita tensor from the Levi-Civita symbol.
    import WeylScal4NRPy.WeylScalars_Cartesian as weyl
    LeviCivitaSymbolDDD = weyl.define_LeviCivitaSymbol_rank3()

    Omega_aligned_rotator = par.Cparameters("REAL",thismodule,"Omega_aligned_rotator") # The angular velocity of the NS
    unit_zU = ixp.zerorank1()
    unit_zU[2] = 1

    global ValenciavU
    ValenciavU = ixp.zerorank1()
    for i in range(DIM):
        for j in range(DIM):
            for k in range(DIM):
                ValenciavU[i] += LeviCivitaSymbolDDD[i][j][k] * Omega_aligned_rotator * unit_zU[j] * rfm.xx[k]


    # ### NRPy+ Module Code Validation
    # 
    # \[Back to [top](#top)\]
    # 
    # Here, as a code validation check, we verify agreement in the SymPy expressions for the $\texttt{GiRaFFE}$ Aligned Rotator initial data equations  we intend to use between
    # 1. this tutorial and 
    # 2. the NRPy+ [GiRaFFEfood_HO_Aligned_Rotator.py](../edit/GiRaFFEfood_HO/GiRaFFEfood_HO_Aligned_Rotator.py) module.
    # 




