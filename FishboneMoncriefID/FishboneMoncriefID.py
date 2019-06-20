# This module builds the equations that describe the Fishbone-Moncrief initial data.
# More thorough documentation can be found in Tutorial-FishboneMoncriefID.ipynb
# Step 1a: Import needed NRPy+ core modules:
import NRPy_param_funcs as par
import indexedexp as ixp
import grid as gri
import finite_difference as fin
from outputC import *
import loop

import reference_metric as rfm

thismodule = __name__

def FishboneMoncriefID(CoordSystem="Cartesian"):
    par.set_parval_from_str("reference_metric::CoordSystem",CoordSystem)
    rfm.reference_metric()
    #Set the spatial dimension parameter to 3.
    par.set_parval_from_str("grid::DIM", 3)
    DIM = par.parval_from_str("grid::DIM")

    gPhys4UU = ixp.register_gridfunctions_for_single_rank2("AUX","gPhys4UU", "sym01", DIM=4)
    KDD = ixp.register_gridfunctions_for_single_rank2("EVOL","KDD", "sym01")

    # Variables needed for initial data given in spherical basis
    r, th, ph = gri.register_gridfunctions("AUX",["r","th","ph"])

    r_in,r_at_max_density,a,M = par.Cparameters("REAL",thismodule,["r_in","r_at_max_density","a","M"])

    kappa,gamma = par.Cparameters("REAL",thismodule,["kappa","gamma"])

    LorentzFactor = gri.register_gridfunctions("AUX","LorentzFactor")

    def calculate_l_at_r(r):
        l  = sp.sqrt(M/r**3) * (r**4 + r**2*a**2 - 2*M*r*a**2 - a*sp.sqrt(M*r)*(r**2-a**2))
        l /= r**2 - 3*M*r + 2*a*sp.sqrt(M*r)
        return l
    # First compute angular momentum at r_at_max_density, TAKING POSITIVE ROOT. This way disk is co-rotating with black hole
    # Eq 3.8:
    l = calculate_l_at_r(r_at_max_density)

    # Eq 3.6:
    # First compute the radially-independent part of the log of the enthalpy, ln_h_const
    Delta = r**2 - 2*M*r + a**2
    Sigma = r**2 + a**2*sp.cos(th)**2
    A = (r**2 + a**2)**2 - Delta*a**2*sp.sin(th)**2

    # Next compute the radially-dependent part of log(enthalpy), ln_h
    tmp3 = sp.sqrt(1 + 4*l**2*Sigma**2*Delta/(A*sp.sin(th))**2)
    # Term 1 of Eq 3.6
    ln_h  = sp.Rational(1,2)*sp.log( ( 1 + tmp3) / (Sigma*Delta/A))
    # Term 2 of Eq 3.6
    ln_h -= sp.Rational(1,2)*tmp3
    # Term 3 of Eq 3.6
    ln_h -= 2*a*M*r*l/A

    # Next compute the radially-INdependent part of log(enthalpy), ln_h
    # Note that there is some typo in the expression for these terms given in Eq 3.6, so we opt to just evaluate
    #   negative of the first three terms at r=r_in and th=pi/2 (the integration constant), as described in
    #   the text below Eq. 3.6, basically just copying the above lines of code.
    # Delin = Delta_in ; Sigin = Sigma_in ; Ain = A_in .
    Delin = r_in**2 - 2*M*r_in + a**2
    Sigin = r_in**2 + a**2*sp.cos(sp.pi/2)**2
    Ain   = (r_in**2 + a**2)**2 - Delin*a**2*sp.sin(sp.pi/2)**2

    tmp3in = sp.sqrt(1 + 4*l**2*Sigin**2*Delin/(Ain*sp.sin(sp.pi/2))**2)
    # Term 4 of Eq 3.6
    mln_h_in  = -sp.Rational(1,2)*sp.log( ( 1 + tmp3in) / (Sigin*Delin/Ain))
    # Term 5 of Eq 3.6
    mln_h_in += sp.Rational(1,2)*tmp3in
    # Term 6 of Eq 3.6
    mln_h_in += 2*a*M*r_in*l/Ain

    global hm1
    hm1 = sp.exp(ln_h + mln_h_in) - 1

    global rho_initial
    rho_initial,Pressure_initial = gri.register_gridfunctions("AUX",["rho_initial","Pressure_initial"])

    # Python 3.4 + sympy 1.0.0 has a serious problem taking the power here, hangs forever.
    # so instead we use the identity x^{1/y} = exp( [1/y] * log(x) )
    # Original expression (works with Python 2.7 + sympy 0.7.4.1):
    # rho_initial = ( hm1*(gamma-1)/(kappa*gamma) )**(1/(gamma - 1))
    # New expression (workaround):
    rho_initial = sp.exp( (1/(gamma-1)) * sp.log( hm1*(gamma-1)/(kappa*gamma) ))
    Pressure_initial = kappa * rho_initial**gamma
    # Eq 3.3: First compute exp(-2 chi), assuming Boyer-Lindquist coordinates
    #    Eq 2.16: chi = psi - nu, so
    #    Eq 3.5 -> exp(-2 chi) = exp(-2 (psi - nu)) = exp(2 nu)/exp(2 psi)
    exp2nu  = Sigma*Delta / A
    exp2psi = A*sp.sin(th)**2 / Sigma
    expm2chi = exp2nu / exp2psi

    # Eq 3.3: Next compute u_(phi).
    u_pphip = sp.sqrt((-1 + sp.sqrt(1 + 4*l**2*expm2chi))/2)
    # Eq 2.13: Compute u_(t)
    u_ptp   = -sp.sqrt(1 + u_pphip**2)

    # Next compute spatial components of 4-velocity in Boyer-Lindquist coordinates:
    uBL4D = ixp.zerorank1(DIM=4) # Components 1 and 2: u_r = u_theta = 0
    # Eq 2.12 (typo): u_(phi) = e^(-psi) u_phi -> u_phi = e^(psi) u_(phi)
    uBL4D[3] = sp.sqrt(exp2psi)*u_pphip

    # Assumes Boyer-Lindquist coordinates:
    omega = 2*a*M*r/A
    # Eq 2.13: u_(t) = 1/sqrt(exp2nu) * ( u_t + omega*u_phi )
    #     -->  u_t = u_(t) * sqrt(exp2nu) - omega*u_phi
    #     -->  u_t = u_ptp * sqrt(exp2nu) - omega*uBL4D[3]
    uBL4D[0] = u_ptp*sp.sqrt(exp2nu) - omega*uBL4D[3]
    # Eq. 3.5:
    # w = 2*a*M*r/A;
    # Eqs. 3.5 & 2.1:
    # gtt = -Sig*Del/A + w^2*Sin[th]^2*A/Sig;
    # gtp = w*Sin[th]^2*A/Sig;
    # gpp = Sin[th]^2*A/Sig;
    # FullSimplify[Inverse[{{gtt,gtp},{gtp,gpp}}]]
    gPhys4BLUU = ixp.zerorank2(DIM=4)
    gPhys4BLUU[0][0] = -A/(Delta*Sigma)
    # DO NOT NEED TO SET gPhys4BLUU[1][1] or gPhys4BLUU[2][2]!
    gPhys4BLUU[0][3] = gPhys4BLUU[3][0] = -2*a*M*r/(Delta*Sigma)
    gPhys4BLUU[3][3] = -4*a**2*M**2*r**2/(Delta*A*Sigma) + Sigma**2/(A*Sigma*sp.sin(th)**2)

    global uBL4U
    uBL4U = ixp.zerorank1(DIM=4)
    for i in range(4):
        for j in range(4):
            uBL4U[i] += gPhys4BLUU[i][j]*uBL4D[j]

    # https://github.com/atchekho/harmpi/blob/master/init.c
    # Next transform Boyer-Lindquist velocity to Kerr-Schild basis:
    transformBLtoKS = ixp.zerorank2(DIM=4)
    for i in range(4):
        transformBLtoKS[i][i] = 1
    transformBLtoKS[0][1] = 2*r/(r**2 - 2*r + a*a)
    transformBLtoKS[3][1] =   a/(r**2 - 2*r + a*a)
    global uKS4U
    uKS4U = ixp.zerorank1(DIM=4)
    for i in range(4):
        for j in range(4):
            uKS4U[i] += transformBLtoKS[i][j]*uBL4U[j]

    # Adopt the Kerr-Schild metric for Fishbone-Moncrief disks
    # http://gravity.psu.edu/numrel/jclub/jc/Cook___LivRev_2000-5.pdf
    # Alternatively, Appendix of https://arxiv.org/pdf/1704.00599.pdf
    rhoKS2  = r**2 + a**2*sp.cos(th)**2 # Eq 79 of Cook's Living Review article
    DeltaKS = r**2 - 2*M*r + a**2    # Eq 79 of Cook's Living Review article
    alphaKS = 1/sp.sqrt(1 + 2*M*r/rhoKS2)
    betaKSU = ixp.zerorank1()
    betaKSU[0] = alphaKS**2*2*M*r/rhoKS2
    gammaKSDD = ixp.zerorank2()
    gammaKSDD[0][0] = 1 + 2*M*r/rhoKS2
    gammaKSDD[0][2] = gammaKSDD[2][0] = -(1 + 2*M*r/rhoKS2)*a*sp.sin(th)**2
    gammaKSDD[1][1] = rhoKS2
    gammaKSDD[2][2] = (r**2 + a**2 + 2*M*r/rhoKS2 * a**2*sp.sin(th)**2) * sp.sin(th)**2

    AA = a**2 * sp.cos(2*th) + a**2 + 2*r**2
    BB = AA + 4*M*r
    DD = sp.sqrt(2*M*r / (a**2 * sp.cos(th)**2 + r**2) + 1)
    KDD[0][0] = DD*(AA + 2*M*r)/(AA**2*BB) * (4*M*(a**2 * sp.cos(2*th) + a**2 - 2*r**2))
    KDD[0][1] = KDD[1][0] = DD/(AA*BB) * 8*a**2*M*r*sp.sin(th)*sp.cos(th)
    KDD[0][2] = KDD[2][0] = DD/AA**2 * (-2*a*M*sp.sin(th)**2 * (a**2 * sp.cos(2*th) + a**2 - 2*r**2))
    KDD[1][1] = DD/BB * 4*M*r**2
    KDD[1][2] = KDD[2][1] = DD/(AA*BB) * (-8*a**3*M*r*sp.sin(th)**3*sp.cos(th))
    KDD[2][2] = DD/(AA**2*BB) * \
                (2*M*r*sp.sin(th)**2 * (a**4*(r-M)*sp.cos(4*th) + a**4*(M+3*r) +
                 4*a**2*r**2*(2*r-M) + 4*a**2*r*sp.cos(2*th)*(a**2 + r*(M+2*r)) + 8*r**5))

    # For compatibility, we must compute gPhys4UU
    gammaKSUU,gammaKSDET = ixp.symm_matrix_inverter3x3(gammaKSDD)
    # See, e.g., Eq. 4.49 of https://arxiv.org/pdf/gr-qc/0703035.pdf , where N = alpha
    gPhys4UU[0][0] = -1 / alphaKS**2
    for i in range(1,4):
        if i>0:
            # if the quantity does not have a "4", then it is assumed to be a 3D quantity.
            #  E.g., betaKSU[] is a spatial vector, with indices ranging from 0 to 2:
            gPhys4UU[0][i] = gPhys4UU[i][0] = betaKSU[i-1]/alphaKS**2
    for i in range(1,4):
        for j in range(1,4):
            # if the quantity does not have a "4", then it is assumed to be a 3D quantity.
            #  E.g., betaKSU[] is a spatial vector, with indices ranging from 0 to 2,
            #    and gammaKSUU[][] is a spatial tensor, with indices again ranging from 0 to 2.
            gPhys4UU[i][j] = gPhys4UU[j][i] = gammaKSUU[i-1][j-1] - betaKSU[i-1]*betaKSU[j-1]/alphaKS**2

    A_b = par.Cparameters("REAL",thismodule,"A_b")

    A_3vecpotentialD = ixp.zerorank1()
    # Set A_phi = A_b*rho_initial FIXME: why is there a sign error?
    A_3vecpotentialD[2] = -A_b * rho_initial

    BtildeU = ixp.register_gridfunctions_for_single_rank1("EVOL","BtildeU")
    # Eq 15 of https://arxiv.org/pdf/1501.07276.pdf:
    # B = curl A -> B^r = d_th A_ph - d_ph A_th
    BtildeU[0] = sp.diff(A_3vecpotentialD[2],th) - sp.diff(A_3vecpotentialD[1],ph)
    # B = curl A -> B^th = d_ph A_r - d_r A_ph
    BtildeU[1] = sp.diff(A_3vecpotentialD[0],ph) - sp.diff(A_3vecpotentialD[2],r)
    # B = curl A -> B^ph = d_r A_th - d_th A_r
    BtildeU[2] = sp.diff(A_3vecpotentialD[1],r)  - sp.diff(A_3vecpotentialD[0],th)

    # Construct spacetime metric in 3+1 form:
    # See, e.g., Eq. 4.49 of https://arxiv.org/pdf/gr-qc/0703035.pdf , where N = alpha
    alpha = gri.register_gridfunctions("EVOL",["alpha"])
    betaU = ixp.register_gridfunctions_for_single_rank1("EVOL","betaU")

    alpha = sp.sqrt(1/(-gPhys4UU[0][0]))
    betaU = ixp.zerorank1()
    for i in range(3):
        betaU[i] = alpha**2 * gPhys4UU[0][i+1]
    gammaUU = ixp.zerorank2()
    for i in range(3):
        for j in range(3):
            gammaUU[i][j] = gPhys4UU[i+1][j+1] + betaU[i]*betaU[j]/alpha**2

    gammaDD = ixp.register_gridfunctions_for_single_rank2("EVOL","gammaDD","sym01")
    gammaDD,igammaDET = ixp.symm_matrix_inverter3x3(gammaUU)
    gammaDET = 1/igammaDET

    ###############
    # Next compute g_{\alpha \beta} from lower 3-metric, using
    # Eq 4.47 of https://arxiv.org/pdf/gr-qc/0703035.pdf
    betaD = ixp.zerorank1()
    for i in range(3):
        for j in range(3):
            betaD[i] += gammaDD[i][j]*betaU[j]

    beta2 = sp.sympify(0)
    for i in range(3):
        beta2 += betaU[i]*betaD[i]

    gPhys4DD = ixp.zerorank2(DIM=4)
    gPhys4DD[0][0] = -alpha**2 + beta2
    for i in range(3):
        gPhys4DD[0][i+1] = gPhys4DD[i+1][0] = betaD[i]
        for j in range(3):
            gPhys4DD[i+1][j+1] = gammaDD[i][j]

    ###############
    # Next compute b^{\mu} using Eqs 23 and 31 of https://arxiv.org/pdf/astro-ph/0503420.pdf
    uKS4D = ixp.zerorank1(DIM=4)
    for i in range(4):
        for j in range(4):
            uKS4D[i] += gPhys4DD[i][j] * uKS4U[j]

    # Eq 27 of https://arxiv.org/pdf/astro-ph/0503420.pdf
    BU = ixp.zerorank1()
    for i in range(3):
        BU[i] = BtildeU[i]/sp.sqrt(gammaDET)

    # Eq 23 of https://arxiv.org/pdf/astro-ph/0503420.pdf
    BU0_u = sp.sympify(0)
    for i in range(3):
        BU0_u += uKS4D[i+1]*BU[i]/alpha

    smallbU = ixp.zerorank1(DIM=4)
    smallbU[0] = BU0_u   / sp.sqrt(4 * sp.pi)
    # Eqs 24 and 31 of https://arxiv.org/pdf/astro-ph/0503420.pdf
    for i in range(3):
        smallbU[i+1] = (BU[i]/alpha + BU0_u*uKS4U[i+1])/(sp.sqrt(4*sp.pi)*uKS4U[0])

    smallbD = ixp.zerorank1(DIM=4)
    for i in range(4):
        for j in range(4):
            smallbD[i] += gPhys4DD[i][j]*smallbU[j]

    smallb2 = sp.sympify(0)
    for i in range(4):
        smallb2 += smallbU[i]*smallbD[i]

    ###############
    LorentzFactor = alpha * uKS4U[0]
    # Define Valencia 3-velocity v^i_(n), which sets the 3-velocity as measured by normal observers to the spatial slice:
    #  v^i_(n) = u^i/(u^0*alpha) + beta^i/alpha. See eq 11 of https://arxiv.org/pdf/1501.07276.pdf
    Valencia3velocityU = ixp.zerorank1()
    for i in range(3):
        Valencia3velocityU[i] = uKS4U[i + 1] / (alpha * uKS4U[0]) + betaU[i] / alpha

    sqrtgamma4DET = sp.symbols("sqrtgamma4DET")
    sqrtgamma4DET = sp.sqrt(gammaDET)*alpha

    alpha = alpha.subs(r,rfm.xxSph[0]).subs(th,rfm.xxSph[1]).subs(ph,rfm.xxSph[2])
    for i in range(DIM):
        betaU[i] = betaU[i].subs(r,rfm.xxSph[0]).subs(th,rfm.xxSph[1]).subs(ph,rfm.xxSph[2])
        for j in range(DIM):
            gammaDD[i][j] = gammaDD[i][j].subs(r,rfm.xxSph[0]).subs(th,rfm.xxSph[1]).subs(ph,rfm.xxSph[2])
            KDD[i][j]     = KDD[i][j].subs(r,rfm.xxSph[0]).subs(th,rfm.xxSph[1]).subs(ph,rfm.xxSph[2])

    # GRMHD variables:
    # Density and pressure:
    hm1           = hm1.subs(r,rfm.xxSph[0]).subs(th,rfm.xxSph[1]).subs(ph,rfm.xxSph[2])
    rho_initial          = rho_initial.subs(r,rfm.xxSph[0]).subs(th,rfm.xxSph[1]).subs(ph,rfm.xxSph[2])
    Pressure_initial     = Pressure_initial.subs(r,rfm.xxSph[0]).subs(th,rfm.xxSph[1]).subs(ph,rfm.xxSph[2])
    LorentzFactor = LorentzFactor.subs(r,rfm.xxSph[0]).subs(th,rfm.xxSph[1]).subs(ph,rfm.xxSph[2])

    # "Valencia" three-velocity
    for i in range(DIM):
        BtildeU[i] = BtildeU[i].subs(r,rfm.xxSph[0]).subs(th,rfm.xxSph[1]).subs(ph,rfm.xxSph[2])
        uKS4U[i+1] = uKS4U[i+1].subs(r,rfm.xxSph[0]).subs(th,rfm.xxSph[1]).subs(ph,rfm.xxSph[2])
        uBL4U[i+1] = uBL4U[i+1].subs(r,rfm.xxSph[0]).subs(th,rfm.xxSph[1]).subs(ph,rfm.xxSph[2])
        Valencia3velocityU[i] = Valencia3velocityU[i].subs(r,rfm.xxSph[0]).subs(th,rfm.xxSph[1]).subs(ph,rfm.xxSph[2])

    # Transform initial data to our coordinate system:
    # First compute Jacobian and its inverse
    drrefmetric__dx_0UDmatrix = sp.Matrix([[sp.diff( rfm.xxSph[0],rfm.xx[0]), sp.diff( rfm.xxSph[0],rfm.xx[1]), sp.diff( rfm.xxSph[0],rfm.xx[2])],
                                           [sp.diff(rfm.xxSph[1],rfm.xx[0]), sp.diff(rfm.xxSph[1],rfm.xx[1]), sp.diff(rfm.xxSph[1],rfm.xx[2])],
                                           [sp.diff(rfm.xxSph[2],rfm.xx[0]), sp.diff(rfm.xxSph[2],rfm.xx[1]), sp.diff(rfm.xxSph[2],rfm.xx[2])]])
    dx__drrefmetric_0UDmatrix = drrefmetric__dx_0UDmatrix.inv()

    # Declare as gridfunctions the final quantities we will output for the initial data
    global IDalpha,IDgammaDD,IDKDD,IDbetaU,IDValencia3velocityU
    IDalpha = gri.register_gridfunctions("EVOL","IDalpha")
    IDgammaDD = ixp.register_gridfunctions_for_single_rank2("EVOL","IDgammaDD","sym01")
    IDKDD = ixp.register_gridfunctions_for_single_rank2("EVOL","IDKDD","sym01")
    IDbetaU   = ixp.register_gridfunctions_for_single_rank1("EVOL","IDbetaU")
    IDValencia3velocityU = ixp.register_gridfunctions_for_single_rank1("EVOL","IDValencia3velocityU")

    IDalpha = alpha
    for i in range(3):
        IDbetaU[i] = 0
        IDValencia3velocityU[i] = 0
        for j in range(3):
            # Matrices are stored in row, column format, so (i,j) <-> (row,column)
            IDbetaU[i]   += dx__drrefmetric_0UDmatrix[(i,j)]*betaU[j]
            IDValencia3velocityU[i]   += dx__drrefmetric_0UDmatrix[(i,j)]*Valencia3velocityU[j]
            IDgammaDD[i][j] = 0
            IDKDD[i][j] = 0
            for k in range(3):
                for l in range(3):
                    IDgammaDD[i][j] += drrefmetric__dx_0UDmatrix[(k,i)]*drrefmetric__dx_0UDmatrix[(l,j)]*gammaDD[k][l]
                    IDKDD[i][j]     += drrefmetric__dx_0UDmatrix[(k,i)]*drrefmetric__dx_0UDmatrix[(l,j)]*    KDD[k][l]


