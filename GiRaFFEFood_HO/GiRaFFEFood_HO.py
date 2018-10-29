#!/usr/bin/env python
# coding: utf-8

# $\newcommand{\giraffe}{\texttt{GiRaFFE}}$
# $\newcommand{\gf}{\texttt{GiRaFFEFood}}$
# # $\gf$: Initial data for $\giraffe$
# 
# With the $\giraffe$ evolution thorn constructed, we now need to "feed" our giraffe with initial data to evolve. While there are several different choices of initial data we can use here, for the moment, we will only be implementing the "Exact Wald" initial data, given by Table 3 in [the original paper](https://arxiv.org/pdf/1704.00599.pdf):
# \begin{align}
# A_{\phi} &= \frac{C_0}{2} r^2 \sin^2 \theta \\
# E_{\phi} &= 2 M C_0 \left( 1+ \frac {2M}{r} \right)^{-1/2} \sin^2 \theta \\
# \end{align}
# (the unspecified components are set to 0). Here, $C_0$ is a constant set to $1$ in our simulations. Now, to use this initial data scheme, we need to transform it into the quantities actually tracked by $\giraffe$ and HydroBase: $A_i$, $B^i$, $S_i$, $v^i$, and $\Phi$. This can be done with eqs. 14, 16, and 18, here given in that same order:
# \begin{align}
# S_\mu &= -n_\nu T^\nu_{{\rm EM} \mu} \\
# v^i &= \alpha \frac{\epsilon^{ijk} E_j B_k}{B^2} -\beta^i \\
# B^i &= \frac{[ijk]}{\sqrt{\gamma}} \partial_j A_k \\
# \end{align}
# -Is this a circular dependency? $\tilde{S}_\mu$ depends on $T^\nu_{{\rm EM} \mu}$ which depends on the velocity, but the velocity depends on $\tilde{S}_\mu$.
# 
# **Zach says:** Check out lines 98-131 in the [original GiRaFFEfood implementation of ExactWald](https://bitbucket.org/zach_etienne/wvuthorns/src/62866e19058cd2477ee88d3a7fb259d5cfc2781a/GiRaFFEfood/src/ExactWald.cc?at=master&fileviewer=file-view-default). You will see that only $A_i$ and $v^i=u^i/u^0$ (not the Valencia 3-velocity that we must specify) are set. Then in the [ID\_Converter\_GiRaFFE](https://bitbucket.org/zach_etienne/wvuthorns/src/62866e19058cd2477ee88d3a7fb259d5cfc2781a/ID_converter_GiRaFFE/src/set_GiRaFFE_metric_GRMHD_variables_based_on_HydroBase_and_ADMBase_variables.C?at=master&fileviewer=file-view-default), $B^i$ and $\tilde{S}_i$ are set based on these quantities. I think that we should have in our version of GiRaFFE a function that computes $\tilde{S}_i$ from $B^i$, $v^i$, and metric quantities.


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

def GiRaFFEFood_HO():
    # Step 1a: Set commonly used parameters.
    thismodule = "GiRaFFEFood_HO"
    # Set the spatial dimension parameter to 3.
    par.set_parval_from_str("grid::DIM", 3)
    DIM = par.parval_from_str("grid::DIM")

    # Create a parameter to control the initial data choice. For now, this will only have Exact Wald as an option.
    par.initialize_param(par.glb_param("char", thismodule, "IDchoice", "Exact_Wald"))

    # Step 1b: Set Cparameters we need to use and the gridfunctions we'll need.
    M,M_PI = par.Cparameters("REAL",thismodule,["M","M_PI"]) # The mass of the black hole, and pi in C
    global StildeD,ValenciavU
    StildeD = ixp.register_gridfunctions_for_single_rank1("AUX","StildeD")
    ValenciavU = ixp.register_gridfunctions_for_single_rank1("AUX","ValenciavU")
    BU = ixp.register_gridfunctions_for_single_rank1("AUX","BU")


    # We will first build the fundamental vectors $A_i$ and $E_i$ in spherical coordinates. Note that we use reference_metric.py to set $r$ and $\theta$ in terms of Cartesian coordinates (see [Table 3](https://arxiv.org/pdf/1704.00599.pdf)).
    # \begin{align}
    # A_{\phi} &= \frac{C_0}{2} r^2 \sin^2 \theta \\
    # E_{\phi} &= 2 M C_0 \left( 1+ \frac {2M}{r} \right)^{-1/2} \sin^2 \theta \\
    # \end{align}
    # While we have $E_i$ set as a variable in NRPy+, note that the final C code won't store these values.


    # Step 2: Set the vectors A and E in Spherical coordinates
    AsphD = ixp.zerorank1()
    EsphD = ixp.zerorank1()

    # The r and theta components (0 and 1) are now 0
    r     = rfm.xxSph[0]
    theta = rfm.xxSph[1]

    IDchoice = par.parval_from_str("IDchoice")

    if IDchoice is "Exact_Wald":
        AsphD[2] = (r * r * sp.sin(theta)**2)/2
        EsphD[2] = 2 * M * sp.sin(theta)**2 / sp.sqrt(1+2*M/r)
    else:
        print("Error: IDchoice == "+par.parval_from_str("IDchoice")+" unsupported!")
        exit(1)


    # Now, we will use the coordinate transformation definitions provided by reference_metric.py to build the Jacobian $$ \frac{\partial x_i}{\partial y_j}, $$ where $x_i \in \{r,\theta,\phi\}$ and $y_i \in \{x,y,z\}$. We will also compute its inverse. Then, since both $A_i$ and $E_i$ have one lower index, both will need to be multiplied by the Jacobian.


    # Step 3: Use the Jacobian matrix to transform the vectors to Cartesian coordinates.
    drrefmetric__dx_0UDmatrix = sp.Matrix([[sp.diff(rfm.xxSph[0],rfm.xx[0]), sp.diff(rfm.xxSph[0],rfm.xx[1]), sp.diff(rfm.xxSph[0],rfm.xx[2])],
                                           [sp.diff(rfm.xxSph[1],rfm.xx[0]), sp.diff(rfm.xxSph[1],rfm.xx[1]), sp.diff(rfm.xxSph[1],rfm.xx[2])],
                                           [sp.diff(rfm.xxSph[2],rfm.xx[0]), sp.diff(rfm.xxSph[2],rfm.xx[1]), sp.diff(rfm.xxSph[2],rfm.xx[2])]])
    #dx__drrefmetric_0UDmatrix = drrefmetric__dx_0UDmatrix.inv() # We don't actually need this in this case.

    global AD
    AD = ixp.register_gridfunctions_for_single_rank1("AUX","AD")
    ED = ixp.register_gridfunctions_for_single_rank1("AUX","ED")

    for i in range(DIM):
        for j in range(DIM):
            AD[i] = drrefmetric__dx_0UDmatrix[(j,i)]*AsphD[j]
            ED[i] = drrefmetric__dx_0UDmatrix[(j,i)]*EsphD[j]

    #Step 4: Register the basic spacetime quantities
    alpha   = gri.register_gridfunctions("AUX","alpha")
    betaU   = ixp.register_gridfunctions_for_single_rank1("AUX","betaU",DIM=3)
    gammaDD = ixp.register_gridfunctions_for_single_rank2("AUX","gammaDD", "sym01",DIM=3)
    gammaUU, gammadet = ixp.symm_matrix_inverter3x3(gammaDD)


    # Now that we have the vector potential and electric fields that we need, we will turn our attention to what other quantities we might need for eqs. 14, 16, and 18 from the [$\giraffe$ paper](https://arxiv.org/pdf/1704.00599.pdf).  We will also need the stress energy tensor $$T^{\mu \nu}_{\rm EM} = b^2 u^{\mu} u^{\nu} + \frac{b^2}{2} g^{\mu \nu} - b^{\mu} b^{\nu}.$$ Note that $T^{\mu \nu}_{\rm EM}$ is in terms of $b^\mu$ and $u^\mu$, provided by $\text{u0_smallb_Poynting__Cartesian}$. 


    # Step 5: Construct the Stress-Energy tensor
    # For TEMUU, we can reuse our code from GiRaFFE_HO.
    import u0_smallb_Poynting__Cartesian.u0_smallb_Poynting__Cartesian as u0b
    u0b.compute_u0_smallb_Poynting__Cartesian(gammaDD,betaU,alpha,ValenciavU,BU)

    # We will now pull in the four metric and its inverse.
    g4DD = ixp.zerorank2(DIM=4)
    g4UU = ixp.zerorank2(DIM=4)
    for mu in range(4):
        for nu in range(4):
            g4DD[mu][nu] = u0b.g4DD[mu][nu]
            g4UU[mu][nu] = u0b.g4UU[mu][nu]

    # We will now pull in the components of the four velocity
    uD = ixp.register_gridfunctions_for_single_rank1("AUX","uD")
    uU = ixp.register_gridfunctions_for_single_rank1("AUX","uU")

    u0 = u0b.u0
    for i in range(DIM):
        uD[i] = u0b.uD[i]
        uU[i] = u0b.uU[i]

    # We will now pull in smallb and related quantities
    smallbU = ixp.zerorank1(DIM=4)
    smallbD = ixp.zerorank1(DIM=4)
    for mu in range(4):
        smallbU[mu] = u0b.smallb4U[mu]
        smallbD[mu] = u0b.smallb4D[mu]

    smallb2 = u0b.smallb2

    TEMUU = ixp.register_gridfunctions_for_single_rank2("AUX","TEMUU","sym01",DIM=4)

    TEMUU[0][0] = smallb2*u0*u0 + smallb2*g4UU[0][0]/2 - smallbU[0]*smallbU[0]
    for mu in range(1,4):
        TEMUU[mu][0] = TEMUU[0][mu] = smallb2*uU[mu-1]*u0 + smallb2*g4UU[mu][0]/2 - smallbU[mu]*smallbU[0]
    for mu in range(1,4):
        for nu in range(1,4):
            TEMUU[mu][nu] = smallb2*uU[mu-1]*uU[nu-1] + smallb2*g4UU[mu][nu]/2 - smallbU[mu]*smallbU[nu]


    # We will now find the densitized Poynting flux given by equation 18 in [the original paper](https://arxiv.org/pdf/1704.00599.pdf), $$S_\mu = -n_\nu T^\nu_{{\rm EM} \mu}$$ and $$\tilde{S}_i = \sqrt{\gamma} S_i, $$ where $n^\mu = (1/\alpha, -\beta^i/\alpha)$.
    # 


    nU = ixp.zerorank1()
    for i in range(DIM):
        nU[i] = -betaU[i] / alpha # we only need the spatial components here.

    # Next, we'll lower the index for both n and T
    nD = ixp.zerorank1()
    TEMUD = ixp.zerorank2()
    for i in range(DIM):
        for j in range(DIM):
            nD[i] = gammaDD[i][j] * nU[j]
            for k in range(DIM):
                TEMUD[i][j] = gammaDD[j][k] * TEMUU[i+1][k+1]

    for i in range(DIM):
        StildeD[i] = sp.sympify(0)
        for j in range(DIM):
            StildeD[i] += -gammadet * nD[j] * TEMUD[j][i]


    # We will now find the magnetic field using equation 18 in [the original paper](https://arxiv.org/pdf/1704.00599.pdf) $$B^i = \frac{[ijk]}{\sqrt{\gamma}} \partial_j A_k. $$ We will need the metric quantites: the lapse $\alpha$, the shift $\beta^i$, and the three-metric $\gamma_{ij}$. We will also need the Levi-Civita symbol, provided by $\text{WeylScal4NRPy}$. 


    # Here, we build the Levi-Civita tensor from the Levi-Civita symbol.
    import WeylScal4NRPy.WeylScalars_Cartesian as weyl
    LeviCivitaDDD = weyl.define_LeviCivitaSymbol_rank3()
    LeviCivitaUUU = ixp.zerorank3()
    for i in range(DIM):
        for j in range(DIM):
            for k in range(DIM):
                LCijk = LeviCivitaDDD[i][j][k]
                LeviCivitaDDD[i][j][k] = LCijk * sp.sqrt(gammadet)
                LeviCivitaUUU[i][j][k] = LCijk / sp.sqrt(gammadet)

    # For the initial data, we can analytically take the derivatives of A_i
    AD_dD = ixp.zerorank2()
    for i in range(DIM):
        for j in range(DIM):
            AD_dD[i][j] = sp.simplify(sp.diff(AD[i],rfm.xxCart[j]))

    BU = ixp.zerorank1()
    for i in range(DIM):
        for j in range(DIM):
            for k in range(DIM):
                BU[i] += LeviCivitaUUU[i][j][k] * AD_dD[k][j]


    # We will now build the drift velocity using equation 152 in [this paper,](https://arxiv.org/pdf/1310.3274v2.pdf) cited in the original $\giraffe$ code: $$ v^i = \alpha \frac{\epsilon^{ijk} E_j B_k}{B^2} -\beta^i. $$ Then, we will need to transform it to the Valencia 3-velocity using the rule $\bar{v}^i = \frac{1}{\alpha} \left(v^i +\beta^i \right)$.
    # 


    # B^2 is an inner product defined in the usual way:
    B2 = sp.sympify(0)
    for i in range(DIM):
        for j in range(DIM):
            B2 += gammaDD[i][j] * BU[i] * BU[j]

    # Lower the index on B^i
    BD = ixp.zerorank1()
    for i in range(DIM):
        for j in range(DIM):
            BD[i] = gammaDD[i][j] * BU[j]

    driftvU = ixp.zerorank1()
    for i in range(DIM):
        driftvU[i] = -betaU[i]
        for j in range(DIM):
            driftvU[i] += alpha*LeviCivitaUUU[i][j][k]*ED[j]*BD[k]/B2

    for i in range(DIM):
        ValenciavU[i] = (driftvU[i] + betaU[i])/alpha

