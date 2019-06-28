# # Computing the 4-Velocity Time-Component $u^0$,
# the Magnetic Field Measured by a Comoving Observer $b^{\mu}$, and the Poynting Vector $S^i$
# 
# ## Part 1 of 2: Computing $u^0$ and $b^{\mu}$
# 
# First some definitions. The spatial components of $b^{\mu}$ are simply the magnetic field as measured by an observer comoving with the plasma $B^{\mu}_{\rm (u)}$, divided by $\sqrt{4\pi}$. In addition, in the ideal MHD limit, $B^{\mu}_{\rm (u)}$ is orthogonal to the plasma 4-velocity $u^\mu$, which sets the $\mu=0$ component. 
# 
# Note also that $B^{\mu}_{\rm (u)}$ is related to the magnetic field as measured by a *normal* observer $B^i$ via a simple projection (Eq 21 in Duez et al), which results in the expressions (Eqs 23 and 24 in Duez et al):
# 
# \begin{align}
# \sqrt{4\pi} b^0 = B^0_{\rm (u)} &= \frac{u_j B^j}{\alpha} \\
# \sqrt{4\pi} b^i = B^i_{\rm (u)} &= \frac{B^i + (u_j B^j) u^i}{\alpha u^0}\\
# \end{align}
# 
# $B^i$ is related to the actual magnetic field evaluated in IllinoisGRMHD, $\tilde{B}^i$ via
# 
# $$B^i = \frac{\tilde{B}^i}{\gamma},$$
# 
# where $\gamma$ is the determinant of the spatial 3-metric.
# 
# The above expressions will require that we compute
# 1. the 4-metric $g_{\mu\nu}$ from the ADM 3+1 variables
# 1. $u^0$ from the Valencia 3-velocity
# 1. $u_j$ from $u^0$, the Valencia 3-velocity, and $g_{\mu\nu}$
# 1. $\gamma$ from the ADM 3+1 variables
# 
# ## Step 1: Compute the 4-metric $g_{\mu\nu}$ from the ADM 3+1 variables
# We are given $\gamma_{ij}$, $\alpha$, and $\beta^i$ from ADMBase, so let's first compute 
# 
# $$
# g_{\mu\nu} = \begin{pmatrix} 
# -\alpha^2 + \beta^k \beta_k & \beta_i \\
# \beta_j & \gamma_{ij}
# \end{pmatrix}.
# $$

import sympy as sp
import NRPy_param_funcs as par
import grid as gri
import indexedexp as ixp
from outputC import *


def compute_u0_smallb_Poynting__Cartesian(gammaDD=None,betaU=None,alpha=None,ValenciavU=None,BU=None):

    if gammaDD==None:
        # Declare these generically if uninitialized.
        gammaDD    = ixp.declarerank2("gammaDD","sym01")
        betaU      = ixp.declarerank1("betaU")
        alpha      = sp.sympify("alpha")
        ValenciavU = ixp.declarerank1("ValenciavU")
        BU         = ixp.declarerank1("BU")
    
    # Set spatial dimension = 3
    DIM=3

    thismodule = "smallbPoynET"

    # To get \gamma_{\mu \nu} = gammabar4DD[mu][nu], we'll need to construct the 4-metric, using Eq. 2.122 in B&S:

    # Eq. 2.121 in B&S
    betaD = ixp.zerorank1()
    for i in range(DIM):
        for j in range(DIM):
            betaD[i] += gammaDD[i][j]*betaU[j]

    # Now compute the beta contraction.
    beta2 = sp.sympify(0)
    for i in range(DIM):
        beta2 += betaU[i]*betaD[i]

    # Eq. 2.122 in B&S
    global g4DD
    g4DD = ixp.zerorank2(DIM=4)
    g4DD[0][0] = -alpha**2 + beta2
    for i in range(DIM):
        g4DD[i+1][0] = g4DD[0][i+1] = betaD[i]
    for i in range(DIM):
        for j in range(DIM):
            g4DD[i+1][j+1] = gammaDD[i][j]

    # ## Step 2: Compute $u^0$ from the Valencia 3-velocity
    #
    # According to Eqs. 9-11 of [the IllinoisGRMHD paper](https://arxiv.org/pdf/1501.07276.pdf), the Valencia 3-velocity $v^i_{(n)}$ is related to the 4-velocity $u^\mu$ via
    #
    # \begin{align}
    # \alpha v^i_{(n)} &= \frac{u^i}{u^0} + \beta^i \\
    # \implies u^i &= u^0 \left(\alpha v^i_{(n)} - \beta^i\right)
    # \end{align}
    #
    # Defining $v^i = \frac{u^i}{u^0}$, we get
    #
    # $$v^i = \alpha v^i_{(n)} - \beta^i,$$
    #
    # and in terms of this variable we get
    #
    # \begin{align}
    # g_{00} \left(u^0\right)^2 + 2 g_{0i} u^0 u^i + g_{ij} u^i u^j &= \left(u^0\right)^2 \left(g_{00} + 2 g_{0i} v^i + g_{ij} v^i v^j\right)\\
    # \implies u^0 &= \pm \sqrt{\frac{-1}{g_{00} + 2 g_{0i} v^i + g_{ij} v^i v^j}} \\
    # &= \pm \sqrt{\frac{-1}{(-\alpha^2 + \beta^2) + 2 \beta_i v^i + \gamma_{ij} v^i v^j}} \\
    # &= \pm \sqrt{\frac{1}{\alpha^2 - \gamma_{ij}\left(\beta^i + v^i\right)\left(\beta^j + v^j\right)}}\\
    # &= \pm \sqrt{\frac{1}{\alpha^2 - \alpha^2 \gamma_{ij}v^i_{(n)}v^j_{(n)}}}\\
    # &= \pm \frac{1}{\alpha}\sqrt{\frac{1}{1 - \gamma_{ij}v^i_{(n)}v^j_{(n)}}}
    # \end{align}
    #
    # Generally speaking, numerical errors will occasionally drive expressions under the radical to either negative values or potentially enormous values (corresponding to enormous Lorentz factors). Thus a reliable approach for computing $u^0$ requires that we first rewrite the above expression in terms of the Lorentz factor squared: $\Gamma^2=\left(\alpha u^0\right)^2$:
    # \begin{align}
    # u^0 &= \pm \frac{1}{\alpha}\sqrt{\frac{1}{1 - \gamma_{ij}v^i_{(n)}v^j_{(n)}}}\\
    # \implies \left(\alpha u^0\right)^2 &= \frac{1}{1 - \gamma_{ij}v^i_{(n)}v^j_{(n)}} \\
    # \implies \gamma_{ij}v^i_{(n)}v^j_{(n)} &= 1 - \frac{1}{\left(\alpha u^0\right)^2}
    # \end{align}
    # In order for the bottom expression to hold true, the left-hand side must be between 0 and 1. Again, this is not guaranteed due to the appearance of numerical errors. In fact, a robust algorithm will not allow $\Gamma^2$ to become too large (which might contribute greatly to the stress-energy of a given gridpoint), so let's define $\Gamma_{\rm max}$, the largest allowed Lorentz factor.
    #
    # Then our algorithm for computing $u^0$ is as follows:
    #
    # If
    # $$R=\gamma_{ij}v^i_{(n)}v^j_{(n)}>1 - \frac{1}{\Gamma_{\rm max}},$$
    # then adjust the 3-velocity $v^i$ as follows:
    #
    # $$v^i_{(n)} = \sqrt{\frac{1 - \frac{1}{\Gamma_{\rm max}}}{R}}v^i_{(n)}.$$
    #
    # After this rescaling, we are then guaranteed that if $R$ is recomputed, it will be set to its ceiling value $R=1 - \frac{1}{\Gamma_{\rm max}}$.
    #
    # Then $u^0$ can be safely computed via
    # u^0 = \frac{1}{\alpha \sqrt{1-R}}.

    # Step 1: Compute R = 1 - 1/max(Gamma)
    R = sp.sympify(0)
    for i in range(DIM):
        for j in range(DIM):
            R += gammaDD[i][j]*ValenciavU[i]*ValenciavU[j]

    GAMMA_SPEED_LIMIT = par.Cparameters("REAL",thismodule,"GAMMA_SPEED_LIMIT")
    Rmax = 1 - 1/GAMMA_SPEED_LIMIT

    rescaledValenciavU = ixp.zerorank1()
    for i in range(DIM):
        rescaledValenciavU[i] = ValenciavU[i]*sp.sqrt(Rmax/R)

    rescaledu0 = 1/(alpha*sp.sqrt(1-Rmax))
    regularu0 =  1/(alpha*sp.sqrt(1-R))

    global computeu0_Cfunction
    computeu0_Cfunction  = "/* Function for computing u^0 from Valencia 3-velocity. */\n"
    computeu0_Cfunction += "/* Inputs: ValenciavU[], alpha, gammaDD[][], GAMMA_SPEED_LIMIT (C parameter) */\n"
    computeu0_Cfunction += "/* Output: u0=u^0 and velocity-limited ValenciavU[] */\n\n"

    computeu0_Cfunction += outputC([R,Rmax],["const double R","const double Rmax"],"returnstring",
                                   params="includebraces=False,CSE_varprefix=tmpR,outCverbose=False")

    computeu0_Cfunction += "if(R <= Rmax) "
    computeu0_Cfunction += outputC(regularu0,"u0","returnstring",
                                   params="includebraces=True,CSE_varprefix=tmpnorescale,outCverbose=False")
    computeu0_Cfunction += " else "
    computeu0_Cfunction += outputC([rescaledValenciavU[0],rescaledValenciavU[1],rescaledValenciavU[2],rescaledu0],
                                   ["ValenciavU0","ValenciavU1","ValenciavU2","u0"],"returnstring",
                                   params="includebraces=True,CSE_varprefix=tmprescale,outCverbose=False")

    # ## Step 3: Compute $u_j$ from $u^0$, the Valencia 3-velocity, and $g_{\mu\nu}$
    # The basic equation is

    # u_j &= g_{\mu j} u^{\mu} \\
    # &= g_{0j} u^0 + g_{ij} u^i \\
    # &= \beta_j u^0 + \gamma_{ij} u^i \\
    # &= \beta_j u^0 + \gamma_{ij} u^0 \left(\alpha v^i_{(n)} - \beta^i\right) \\
    # &= u^0 \left(\beta_j + \gamma_{ij} \left(\alpha v^i_{(n)} - \beta^i\right) \right)\\
    # &= \alpha u^0 \gamma_{ij} v^i_{(n)} \\

    global u0
    u0 = par.Cparameters("REAL",thismodule,"u0")
    global uD
    uD = ixp.zerorank1()
    for i in range(DIM):
        for j in range(DIM):
            uD[j] += alpha*u0*gammaDD[i][j]*ValenciavU[j]

    # ## Step 4: Compute $\gamma=\text{gammaDET}$ from the ADM 3+1 variables
    # This is accomplished simply, using the symmetric matrix inversion function in indexedexp.py:
    gammaUU, gammaDET = ixp.symm_matrix_inverter3x3(gammaDD)


    # ## Step 5: Compute $\beta^\mu$ from above expressions.

    # \sqrt{4\pi} b^0 = B^0_{\rm (u)} &= \frac{u_j B^j}{\alpha} \\
    # \sqrt{4\pi} b^i = B^i_{\rm (u)} &= \frac{B^i + (u_j B^j) u^i}{\alpha u^0}\\

    # $B^i$ is related to the actual magnetic field evaluated in IllinoisGRMHD, $\tilde{B}^i$ via
    #
    # $$B^i = \frac{\tilde{B}^i}{\gamma},$$
    #
    # where $\gamma$ is the determinant of the spatial 3-metric.
    #
    # Pulling this together, we currently have available as input:
    # + $\tilde{B}^i$
    # + $\gamma$
    # + $u_j$
    # + $u^0$,
    #
    # with the goal of outputting now $b^\mu$ and $b^2$:
    M_PI = par.Cparameters("REAL",thismodule,"M_PI")

    # uBcontraction = u_i B^i
    global uBcontraction
    uBcontraction = sp.sympify(0)
    for i in range(DIM):
        uBcontraction += uD[i]*BU[i]

    # uU = 3-vector representing u^i = u^0 \left(\alpha v^i_{(n)} - \beta^i\right)
    global uU
    uU = ixp.zerorank1()
    for i in range(DIM):
        uU[i] = u0*(alpha*ValenciavU[i] - betaU[i])

    global smallb4U
    smallb4U = ixp.zerorank1(DIM=4)
    smallb4U[0] = uBcontraction/(alpha*sp.sqrt(4*M_PI))
    for i in range(DIM):
        smallb4U[1+i] = (BU[i] + uBcontraction*uU[i])/(alpha*u0*sp.sqrt(4*M_PI))

    # ## Part 2 of 2: Computing the Poynting Flux Vector $S^{i}$
    #
    # The Poynting flux is defined in Eq. 11 of [Kelly *et al*](https://arxiv.org/pdf/1710.02132.pdf):
    # S^i = -\alpha T^i_{\rm EM\ 0} = \alpha\left(b^2 u^i u_0 + \frac{1}{2} b^2 g^i{}_0 - b^i b_0\right)

    # ## Part 2, Step 1: Computing $g^{i\nu}$:
    # We have already computed all of these quantities above, except $g^i{}_0$ so let's now construct this object:

    # g^i{}_0 = g^{i\nu} g_{\nu 0},

    # which itself requires $g^{i\nu}$ be defined (Eq. 4.49 in [Gourgoulhon](https://arxiv.org/pdf/gr-qc/0703035.pdf)):

    # g^{\mu\nu} = \begin{pmatrix}
    # -\frac{1}{\alpha^2} & \frac{\beta^i}{\alpha^2} \\
    # \frac{\beta^i}{\alpha^2} & \gamma^{ij} - \frac{\beta^i\beta^j}{\alpha^2}
    # \end{pmatrix},

    # where $\gamma^{ij}$ was defined above where we computed $\text{gammaDET}$.
    global g4UU
    g4UU = ixp.zerorank2(DIM=4)

    g4UU[0][0] = -1 / alpha**2
    for i in range(DIM):
        g4UU[0][i+1] = g4UU[i+1][0] = betaU[i]/alpha**2
    for i in range(DIM):
        for j in range(DIM):
            g4UU[i+1][j+1] = gammaUU[i][j] - betaU[i]*betaU[j]/alpha**2


    # ## Part 2, Step 2: Computing $S^{i}$
    #
    # We start by computing
    # g^\mu{}_\delta = g^{\mu\nu} g_{\nu \delta},

    # and then the rest of the Poynting flux vector can be immediately computed from quantities defined above:
    # S^i = \alpha T^i_{\rm EM\ 0} = -\alpha\left(b^2 u^i u_0 + \frac{1}{2} b^2 g^i{}_0 - b^i b_0\right)

    # Step 2a: compute g^\mu_\delta:
    g4UD = ixp.zerorank2(DIM=4)
    for mu in range(4):
        for delta in range(4):
            for nu in range(4):
                g4UD[mu][delta] += g4UU[mu][nu]*g4DD[nu][delta]

    # Step 2b: compute b_{\mu}
    global smallb4D
    smallb4D = ixp.zerorank1(DIM=4)
    for mu in range(4):
        for nu in range(4):
            smallb4D[mu] += g4DD[mu][nu]*smallb4U[nu]

    # Step 2c: compute u_0 = g_{mu 0} u^{mu} = g4DD[0][0]*u0 + g4DD[i][0]*uU[i]
    u_0 = g4DD[0][0]*u0
    for i in range(DIM):
        u_0 += g4DD[i+1][0]*uU[i]

    # Step 2d: compute b^2
    global smallb2etk
    smallb2etk = sp.sympify(0)
    for mu in range(4):
        smallb2etk += smallb4U[mu]*smallb4D[mu]

    # Step 2d: compute S^i
    global PoynSU
    PoynSU = ixp.zerorank1()
    for i in range(DIM):
        PoynSU[i] = -alpha * (smallb2etk*uU[i]*u_0 + sp.Rational(1,2)*smallb2etk*g4UD[i+1][0] - smallb4U[i+1]*smallb4D[0])
