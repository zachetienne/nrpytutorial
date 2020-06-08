# # Computing the 4-Velocity Time-Component $u^0$,
# the Magnetic Field Measured by a Comoving Observer $b^{\mu}$, and the Poynting Vector $S^i$

# Authors: Zachariah B. Etienne
#          zachetie **at** gmail **dot* com
#          Patrick D. Nelson

# Step 1: Initialize needed Python/NRPy+ modules
import indexedexp as ixp                   # NRPy+: Symbolic indexed expression (e.g., tensors, vectors, etc.) support
from outputC import outputC                # NRPy+: Basic C code output functionality
import NRPy_param_funcs as par             # NRPy+: parameter interface
import sympy as sp                         # SymPy: The Python computer algebra package upon which NRPy+ depends
import BSSN.ADMBSSN_tofrom_4metric as AB4m # NRPy+: ADM/BSSN <-> 4-metric conversions

def compute_u0_smallb_Poynting__Cartesian(gammaDD=None,betaU=None,alpha=None,ValenciavU=None,BU=None):

    if gammaDD is None: # use "is None" instead of "==None", as the former is more correct.
        # Declare these generically if uninitialized.
        gammaDD    = ixp.declarerank2("gammaDD","sym01")
        betaU      = ixp.declarerank1("betaU")
        alpha      = sp.sympify("alpha")
        ValenciavU = ixp.declarerank1("ValenciavU")
        BU         = ixp.declarerank1("BU")

    # Set spatial dimension = 3
    DIM=3

    thismodule = __name__

    # Step 1.a: Compute the 4-metric $g_{\mu\nu}$ and its inverse
    #           $g^{\mu\nu}$ from the ADM 3+1 variables, using the
    #           BSSN.ADMBSSN_tofrom_4metric NRPy+ module
    AB4m.g4DD_ito_BSSN_or_ADM("ADM",gammaDD,betaU,alpha)
    g4DD = AB4m.g4DD
    AB4m.g4UU_ito_BSSN_or_ADM("ADM",gammaDD,betaU,alpha)
    g4UU = AB4m.g4UU

    # Step 1.b: Our algorithm for computing $u^0$ is as follows:
    #
    # Let
    # R = gamma_{ij} v^i_{(n)} v^j_{(n)} > 1 - 1 / Gamma_MAX.
    # Then the velocity exceeds the speed limit (set by the
    # maximum Lorentz Gamma, Gamma_MAX), and adjust the
    # 3-velocity $v^i$ as follows:
    #
    # v^i_{(n)} = \sqrt{(1 - 1/Gamma_MAX)/R} * v^i_{(n)}
    #
    # After this rescaling, we are then guaranteed that if
    # R is recomputed, it will be set to its ceiling value
    #  R = 1 - 1 / Gamma_MAX,
    #
    # Then $u^0$ can be safely computed via
    # u^0 = 1 / (alpha \sqrt{1-R}).

    # Step 1.b.i: Compute R = 1 - 1/max(Gamma)
    R = sp.sympify(0)
    for i in range(DIM):
        for j in range(DIM):
            R += gammaDD[i][j]*ValenciavU[i]*ValenciavU[j]

    # Step 1.b.ii: Output C code for computing u^0
    GAMMA_SPEED_LIMIT = par.Cparameters("REAL",thismodule,"GAMMA_SPEED_LIMIT",10.0) # Default value based on
                                                                                    # IllinoisGRMHD.
                                                                                    # GiRaFFE default = 2000.0
    Rmax = 1 - 1/(GAMMA_SPEED_LIMIT * GAMMA_SPEED_LIMIT)

    rescaledValenciavU = ixp.zerorank1()
    for i in range(DIM):
        rescaledValenciavU[i] = ValenciavU[i]*sp.sqrt(Rmax/R)

    rescaledu0 = 1/(alpha*sp.sqrt(1-Rmax))
    regularu0 =  1/(alpha*sp.sqrt(1-R))

    global computeu0_Cfunction
    computeu0_Cfunction  = """
/* Function for computing u^0 from Valencia 3-velocity. */
/* Inputs: ValenciavU[], alpha, gammaDD[][], GAMMA_SPEED_LIMIT (C parameter) */
/* Output: u0=u^0 and velocity-limited ValenciavU[] */\n\n"""

    computeu0_Cfunction += outputC([R,Rmax],["const double R","const double Rmax"],"returnstring",
                                   params="includebraces=False,CSE_varprefix=tmpR,outCverbose=False")

    computeu0_Cfunction += "if(R <= Rmax) "
    computeu0_Cfunction += outputC(regularu0,"u0","returnstring",
                                   params="includebraces=True,CSE_varprefix=tmpnorescale,outCverbose=False")
    computeu0_Cfunction += " else "
    computeu0_Cfunction += outputC([rescaledValenciavU[0],rescaledValenciavU[1],rescaledValenciavU[2],rescaledu0],
                                   ["ValenciavU0","ValenciavU1","ValenciavU2","u0"],"returnstring",
                                   params="includebraces=True,CSE_varprefix=tmprescale,outCverbose=False")

    # ## Step 1.c: Compute u_j from u^0, the Valencia 3-velocity,
    #    and g_{mu nu}
    # The basic equation is

    # u_j &= g_{\mu j} u^{\mu} \\
    # &= g_{0j} u^0 + g_{ij} u^i \\
    # &= \beta_j u^0 + \gamma_{ij} u^i \\
    # &= \beta_j u^0 + \gamma_{ij} u^0 \left(\alpha v^i_{(n)} - \beta^i\right) \\
    # &= u^0 \left(\beta_j + \gamma_{ij} \left(\alpha v^i_{(n)} - \beta^i\right) \right)\\
    # &= \alpha u^0 \gamma_{ij} v^i_{(n)} \\

    global u0
    u0 = par.Cparameters("REAL",thismodule,"u0",1e300) # Will be overwritten in C code. Set to crazy value to ensure this.
    global uD
    uD = ixp.zerorank1()
    for i in range(DIM):
        for j in range(DIM):
            uD[j] += alpha*u0*gammaDD[i][j]*ValenciavU[i]

    # ## Step 1.d: Compute $b^\mu$ from above expressions.

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

    # with the goal of outputting now $b^\mu$ and $b^2$:
    M_PI = par.Cparameters("#define",thismodule,"M_PI","")

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

    # Step 2: Compute the Poynting flux vector S^i
    #
    # The Poynting flux is defined in Eq. 11 of [Kelly *et al*](https://arxiv.org/pdf/1710.02132.pdf):
    # S^i = -\alpha T^i_{\rm EM\ 0} = \alpha\left(b^2 u^i u_0 + \frac{1}{2} b^2 g^i{}_0 - b^i b_0\right)

    # We start by computing
    # g^\mu{}_\delta = g^{\mu\nu} g_{\nu \delta},

    # and then the rest of the Poynting flux vector can be immediately computed from quantities defined above:
    # S^i = \alpha T^i_{\rm EM\ 0} = -\alpha\left(b^2 u^i u_0 + \frac{1}{2} b^2 g^i{}_0 - b^i b_0\right)

    # Step 2.a.i: compute g^\mu_\delta:
    g4UD = ixp.zerorank2(DIM=4)
    for mu in range(4):
        for delta in range(4):
            for nu in range(4):
                g4UD[mu][delta] += g4UU[mu][nu]*g4DD[nu][delta]

    # Step 2.a.ii: compute b_{\mu}
    global smallb4D
    smallb4D = ixp.zerorank1(DIM=4)
    for mu in range(4):
        for nu in range(4):
            smallb4D[mu] += g4DD[mu][nu]*smallb4U[nu]

    # Step 2.a.iii: compute u_0 = g_{mu 0} u^{mu} = g4DD[0][0]*u0 + g4DD[i][0]*uU[i]
    u_0 = g4DD[0][0]*u0
    for i in range(DIM):
        u_0 += g4DD[i+1][0]*uU[i]

    # Step 2.a.iv: compute b^2, setting b^2 = smallb2etk, as gridfunctions with base names ending in a digit
    #          are forbidden in NRPy+.
    global smallb2etk
    smallb2etk = sp.sympify(0)
    for mu in range(4):
        smallb2etk += smallb4U[mu]*smallb4D[mu]

    # Step 2.a.v: compute S^i
    global PoynSU
    PoynSU = ixp.zerorank1()
    for i in range(DIM):
        PoynSU[i] = -alpha * (smallb2etk*uU[i]*u_0 + sp.Rational(1,2)*smallb2etk*g4UD[i+1][0] - smallb4U[i+1]*smallb4D[0])
