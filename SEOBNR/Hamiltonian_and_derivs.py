from outputC import *
import sympy as sp

def simplify_deriv(lhss_deriv,rhss_deriv):
    lhss_deriv_simp = []
    rhss_deriv_simp = []
    for i in range(len(rhss_deriv)):
        lhss_deriv_simp.append(lhss_deriv[i])
        rhss_deriv_simp.append(rhss_deriv[i])
    for i in range(len(rhss_deriv_simp)):
        if rhss_deriv_simp[i] == 0:
            for j in range(i+1,len(rhss_deriv_simp)):
                for var in rhss_deriv_simp[j].free_symbols:
                    if str(var) == str(lhss_deriv_simp[i]):
                        rhss_deriv_simp[j] = rhss_deriv_simp[j].subs(var,0)
    zero_elements_to_remove = []
    for i in range(len(rhss_deriv_simp)):
        if rhss_deriv_simp[i] == sp.sympify(0):
            zero_elements_to_remove.append(i)
    count = 0
    for i in range(len(zero_elements_to_remove)):
        del lhss_deriv_simp[zero_elements_to_remove[i]+count]
        del rhss_deriv_simp[zero_elements_to_remove[i]+count]
        count -= 1
    return lhss_deriv_simp,rhss_deriv_simp

def deriv_onevar(lhss_deriv, rhss_deriv,
                 xprm=0, yprm=0, zprm=0, pxprm=0, pyprm=0, pzprm=0, s1xprm=0, s1yprm=0, s1zprm=0, s2xprm=0, s2yprm=0,
                 s2zprm=0):
    lhss_deriv_new = []
    rhss_deriv_new = []
    for i in range(len(rhss_deriv)):
        lhss_deriv_new.append(lhss_deriv[i])
        rhss_deriv_new.append(rhss_deriv[i])
    for i in range(len(rhss_deriv_new)):
        for var in rhss_deriv_new[i].free_symbols:
            if str(var) == "xprm":
                rhss_deriv_new[i] = rhss_deriv_new[i].subs(var, xprm)
            elif str(var) == "yprm":
                rhss_deriv_new[i] = rhss_deriv_new[i].subs(var, yprm)
            elif str(var) == "zprm":
                rhss_deriv_new[i] = rhss_deriv_new[i].subs(var, zprm)
            elif str(var) == "pxprm":
                rhss_deriv_new[i] = rhss_deriv_new[i].subs(var, pxprm)
            elif str(var) == "pyprm":
                rhss_deriv_new[i] = rhss_deriv_new[i].subs(var, pyprm)
            elif str(var) == "pzprm":
                rhss_deriv_new[i] = rhss_deriv_new[i].subs(var, pzprm)
            elif str(var) == "s1xprm":
                rhss_deriv_new[i] = rhss_deriv_new[i].subs(var, s1xprm)
            elif str(var) == "s1yprm":
                rhss_deriv_new[i] = rhss_deriv_new[i].subs(var, s1yprm)
            elif str(var) == "s1zprm":
                rhss_deriv_new[i] = rhss_deriv_new[i].subs(var, s1zprm)
            elif str(var) == "s2xprm":
                rhss_deriv_new[i] = rhss_deriv_new[i].subs(var, s2xprm)
            elif str(var) == "s2yprm":
                rhss_deriv_new[i] = rhss_deriv_new[i].subs(var, s2yprm)
            elif str(var) == "s2zprm":
                rhss_deriv_new[i] = rhss_deriv_new[i].subs(var, s2zprm)
    lhss_deriv_simp, rhss_deriv_simp = simplify_deriv(lhss_deriv_new, rhss_deriv_new)
    return lhss_deriv_simp, rhss_deriv_simp


def output_H_and_derivs():

    Hamstring = """
sigmaKerr0 = s1x + s2x
sigmaKerr1 = s1y + s2y
sigmaKerr2 = s1z + s2z
invm1m2 = 1/(m1*m2)
m2overm1 = m2*m2*invm1m2
m1overm2 = m1*m1*invm1m2
sigmaStar0 = (m2overm1)*s1x + (m1overm2)*s2x
sigmaStar1 = (m2overm1)*s1y + (m1overm2)*s2y
sigmaStar2 = (m2overm1)*s1z + (m1overm2)*s2z
s1dots1 = s1x*s1x + s1y*s1y + s1z*s1z
s2dots2 = s2x*s2x + s2y*s2y + s2z*s2z
r2 = x*x + y*y + z*z
r = sp.sqrt(r2)
u = 1/r
u2 = u*u
u3 = u2*u
u4 = u2*u2
u5 = u4*u
etau3 = eta*u3
etau4 = eta*u4
nx = x*u
ny = y*u
nz = z*u
sKerrUSCOREx = sigmaKerr0
sKerrUSCOREy = sigmaKerr1
sKerrUSCOREz = sigmaKerr2
sStarUSCOREx = sigmaStar0
sStarUSCOREy = sigmaStar1
sStarUSCOREz = sigmaStar2
a2 = sKerrUSCOREx*sKerrUSCOREx + sKerrUSCOREy*sKerrUSCOREy + sKerrUSCOREz*sKerrUSCOREz
a4 = a2*a2
a = sp.sqrt(a2)
inva = 1/a
m1PlusetaKK = -1 + eta*KK
invm1PlusetaKK = 1/m1PlusetaKK
k0 = KK*(m1PlusetaKK - 1)
k1 = -2*(k0 + KK)*m1PlusetaKK
k2 = c0k2 + c1k2*a2
k3 = c0k3 + c1k3*a2
k4 = c0k4 + c1k4*a2 + c2k4*a4
k5 = c0k5 + c1k5*a2 + c2k5*a4
e3USCOREx = sKerrUSCOREx*inva
e3USCOREy = sKerrUSCOREy*inva
e3USCOREz = sKerrUSCOREz*inva
costheta = e3USCOREx*nx + e3USCOREy*ny + e3USCOREz*nz
xi2 = 1 - costheta*costheta
xiUSCOREx = -e3USCOREz*ny + e3USCOREy*nz
xiUSCOREy =  e3USCOREz*nx - e3USCOREx*nz
xiUSCOREz = -e3USCOREy*nx + e3USCOREx*ny
vx = -nz*xiUSCOREy + ny*xiUSCOREz
vy =  nz*xiUSCOREx - nx*xiUSCOREz
vz = -ny*xiUSCOREx + nx*xiUSCOREy
w2 = r2 + a2
rho2 = r2 + a2*costheta*costheta
bulk = invm1PlusetaKK*(invm1PlusetaKK + 2*u) + a2*u2
logu = sp.log( u )
logarg = k1*u + k2*u2 + k3*u3 + k4*u4 + k5*u5 + k5l*u5*logu
onepluslogarg = (1 + logarg)
invonepluslogarg = 1/onepluslogarg
logTerms = 1 + eta*k0 + eta*sp.log(onepluslogarg)
deltaU = bulk*logTerms
deltaT = r2*deltaU
deltaUUSCOREupt7 = k5 + k5l*logu
deltaUUSCOREupt6 = 4*k4 + 5*deltaUUSCOREupt7*u
deltaUUSCOREupt5 = 3*k3 + u*deltaUUSCOREupt6
deltaUUSCOREupt4 = 2*k2 + u*deltaUUSCOREupt5
deltaUUSCOREupt3 = k1 + u*deltaUUSCOREupt4
deltaUUSCOREupt2 = invm1PlusetaKK + a2*u
deltaUUSCOREupt1 = bulk*eta*deltaUUSCOREupt3
deltaUUSCOREu = 2*deltaUUSCOREupt2*logTerms + deltaUUSCOREupt1*invonepluslogarg
deltaTUSCOREr = 2*r*deltaU - deltaUUSCOREu
Lambda = w2*w2 - a2*deltaT*xi2
rho2xi2Lambda = rho2*xi2*Lambda
invrho2xi2Lambda = 1/rho2xi2Lambda
invrho2 = xi2*Lambda*invrho2xi2Lambda
invxi2 = rho2*Lambda*invrho2xi2Lambda
invLambda = xi2*rho2*invrho2xi2Lambda
invLambdasq = invLambda*invLambda
rho2invLambda = rho2*invLambda
expnu = sp.sqrt(deltaT*rho2invLambda)
expMU = sp.sqrt(rho2)
expMUexpnu = expMU*expnu
expMUsq = expMU*expMU
expnusq = expnu*expnu
expMUsqexpnusq = expMUsq*expnusq
invexpnuexpMU = 1/expMUexpnu
invexpMU = expnu*invexpnuexpMU
invexpMUsq = invexpMU*invexpMU
expnuinvexpMU2 = expnu*invexpMUsq
invexpMUcubinvexpnu = invexpMUsq*invexpnuexpMU
DD = 1 +  sp.log(1 + 6*eta*u2 + 2*(26 - 3*eta)*etau3)
deltaR = deltaT*DD
qq = 2*eta*(4 - 3*eta)
ww = 2*a*r + b3*eta*a2*a*u + bb3*eta*a*u
B = sp.sqrt(deltaT)
sqrtdeltaT = B
sqrtdeltaR = sp.sqrt(deltaR)
deltaTsqrtdeltaR = deltaT*sqrtdeltaR
sqrtdeltaTdeltaTsqrtdeltaR = sqrtdeltaT*deltaTsqrtdeltaR
invdeltaTsqrtdeltaTsqrtdeltaR = 1./sqrtdeltaTdeltaTsqrtdeltaR
invdeltaT = sqrtdeltaT*(sqrtdeltaR*invdeltaTsqrtdeltaTsqrtdeltaR)
invsqrtdeltaT = deltaTsqrtdeltaR*invdeltaTsqrtdeltaTsqrtdeltaR
invsqrtdeltaR = deltaT*sqrtdeltaT*invdeltaTsqrtdeltaTsqrtdeltaR
w = ww*invLambda
LambdaUSCOREr = 4*r*w2 - a2*deltaTUSCOREr*xi2
wwUSCOREr = 2*a - (a2*a*b3*eta)*u2 - bb3*eta*a*u2
BR = (-deltaT*invsqrtdeltaR + deltaTUSCOREr*sp.Rational(1,2))*invsqrtdeltaT
wr = (-LambdaUSCOREr*ww + Lambda*wwUSCOREr)*(invLambdasq)
nurpt2 = w2*(-4.*r*deltaT + w2*deltaTUSCOREr)
nurpt1 = nurpt2*invdeltaT
nur = r*invrho2 + sp.Rational(1,2)*invLambda*nurpt1
mur = r*invrho2 - invsqrtdeltaR
a2costheta = a2*costheta
wcospt2 = deltaT*ww
wcospt1 = invLambdasq*wcospt2
wcos  = -2*(a2costheta)*wcospt1
nucospt3 = invrho2*invLambda
nucospt2 = w2*nucospt3
nucospt1 = a2costheta*nucospt2
nucos = (w2 - deltaT)*nucospt1
mucos = a2costheta*invrho2
etaover12r = eta*sp.Rational(1,12)*u
csi = sp.sqrt(deltaT*deltaR)/w2
csi1  =  1 + (1-sp.Abs(1-tortoise)) * (csi - 1)
#csi2  =  1 + (sp.Rational(1,2)-copysign(0.5, 1.5-tortoise)) * (csi - 1)
csi2  =  1 + (sp.Rational(1,2)-copysignresult) * (csi - 1)
prT = (px*nx + py*ny + pz*nz)*csi2
prTtimesoneminuscsi1inv = prT*(1 - 1/csi1)
tmppx = px - nx*prTtimesoneminuscsi1inv
tmppy = py - ny*prTtimesoneminuscsi1inv
tmppz = pz - nz*prTtimesoneminuscsi1inv
pxir = (tmppx*xiUSCOREx + tmppy*xiUSCOREy + tmppz*xiUSCOREz)*r
pvr  = (tmppx*vx + tmppy*vy + tmppz*vz)*r
pvrsq = pvr*pvr
pn   =  tmppx*nx + tmppy*ny + tmppz*nz
pnsq = pn*pn
pr = pn
prsq = pr*pr
pf = pxir
pxirsq = pxir*pxir
ptheta2 = pvrsq*invxi2
prT4=prT*prT*prT*prT
Hnspt7 = deltaR*invrho2
Hnspt6 = rho2invLambda*invxi2
Hnspt5 = qq*u2
Hnspt4 = (1 + prT4*Hnspt5 + ptheta2*invrho2 + pf*pf*Hnspt6 + prsq*Hnspt7)
Hnspt3 = deltaT*Hnspt4
Hnspt2 = rho2*Hnspt3
Hnspt1 = pf*ww
Hns = sp.sqrt(Hnspt2*invLambda) + invLambda*Hnspt1
Qpt3 = deltaR*invrho2
Qpt2 = rho2invLambda*invxi2
Qpt1 = invrho2*invxi2
Q = 1 + pvrsq*Qpt1 + pxirsq*Qpt2 + pnsq*Qpt3
pn2 = prsq*deltaR*invrho2
pp  = Q - 1
sKerrmultfact = (-8 - 3*r*(12*pn2 - pp))
sStarmultfact = (14 + (- 30*pn2 + 4*pp)*r)
deltaSigmaStarUSCOREx1=etaover12r*(sKerrmultfact*sKerrUSCOREx + sStarmultfact*sStarUSCOREx)
deltaSigmaStarUSCOREy1=etaover12r*(sKerrmultfact*sKerrUSCOREy + sStarmultfact*sStarUSCOREy)
deltaSigmaStarUSCOREz1=etaover12r*(sKerrmultfact*sKerrUSCOREz + sStarmultfact*sStarUSCOREz)
pn2pp = pn2*pp
pp2 = pp*pp
pn2u2 = pn2*u2
ppu2 = pp*u2
pn2ppu2 = pn2pp*u2
sMultiplier1pt6 = -360*pn2*pn2 + 126*pn2pp + 3*pp2
sMultiplier1pt5 = -96*pn2pp + 23*pp2
sMultiplier1pt4 = -120*pp + 324*pn2 + sMultiplier1pt6*r
sMultiplier1pt3 = 206*pp - 282*pn2 + sMultiplier1pt5*r
sMultiplier1pt2 = 54 + sMultiplier1pt4*r
sMultiplier1pt1 = -706 + sMultiplier1pt3*r + sMultiplier1pt2*eta
sMultiplier1 = sMultiplier1pt1*eta*u2*sp.Rational(-1,72)
sMultiplier2pt6 = sp.Rational(45,8)*pn2*pn2u2 - sp.Rational(13,8)*pn2ppu2
sMultiplier2pt5 = pn2ppu2/4 - sp.Rational(5,16)*pp2*u2
sMultiplier2pt4 = sp.Rational(-49,8)*pn2u2 + sp.Rational(17,12)*ppu2 + sMultiplier2pt6*r
sMultiplier2pt3 = sp.Rational(-2,3)*pn2u2 - sp.Rational(109,36)*ppu2 + sMultiplier2pt5*r
sMultiplier2pt2 = sp.Rational(-7,3)*u2 + sMultiplier2pt4*r
sMultiplier2pt1 = sp.Rational(-56,9)*u2 + sMultiplier2pt3*r + sMultiplier2pt2*eta
sMultiplier2 = sMultiplier2pt1*eta
deltaSigmaStarUSCOREx2 = deltaSigmaStarUSCOREx1 + sMultiplier1*sigmaStar0 + sMultiplier2*sigmaKerr0
deltaSigmaStarUSCOREy2 = deltaSigmaStarUSCOREy1 + sMultiplier1*sigmaStar1 + sMultiplier2*sigmaKerr1
deltaSigmaStarUSCOREz2 = deltaSigmaStarUSCOREz1 + sMultiplier1*sigmaStar2 + sMultiplier2*sigmaKerr2
deltaSigmaStarUSCOREx3 = deltaSigmaStarUSCOREx2 + d1*sigmaStar0*etau3
deltaSigmaStarUSCOREy3 = deltaSigmaStarUSCOREy2 + d1*sigmaStar1*etau3
deltaSigmaStarUSCOREz3 = deltaSigmaStarUSCOREz2 + d1*sigmaStar2*etau3
deltaSigmaStarUSCOREx = deltaSigmaStarUSCOREx3 + d1v2*sigmaKerr0*etau3
deltaSigmaStarUSCOREy = deltaSigmaStarUSCOREy3 + d1v2*sigmaKerr1*etau3
deltaSigmaStarUSCOREz = deltaSigmaStarUSCOREz3 + d1v2*sigmaKerr2*etau3
sx = sStarUSCOREx + deltaSigmaStarUSCOREx
sy = sStarUSCOREy + deltaSigmaStarUSCOREy
sz = sStarUSCOREz + deltaSigmaStarUSCOREz
sxi = sx*xiUSCOREx + sy*xiUSCOREy + sz*xiUSCOREz
sv  = sx*vx + sy*vy + sz*vz
sn  = sx*nx + sy*ny + sz*nz
s3 = sx*e3USCOREx + sy*e3USCOREy + sz*e3USCOREz
sqrtQ = sp.sqrt(Q)
oneplus2sqrtQ = 1 + 2*sqrtQ
oneplus1sqrtQ = oneplus2sqrtQ - sqrtQ
twoB1psqrtQsqrtQ = (2*B*oneplus1sqrtQ*sqrtQ)
invtwoB1psqrtQsqrtQ = 1/twoB1psqrtQsqrtQ
expMUsqsqrtQplusQ = (expMUsq)*(sqrtQ + Q)
Hwrpt4a = pxirsq*sv
Hwrpt4 = expMUsqexpnusq*Hwrpt4a
Hwrpt3c = pxir*sxi
Hwrpt3b = pvr*Hwrpt3c
Hwrpt3a = expMUexpnu*Hwrpt3b
Hwrpt3 = B*Hwrpt3a
Hwrpt2g = sv*deltaR
Hwrpt2f = sn*sqrtdeltaR
Hwrpt2e = pvr*Hwrpt2f
Hwrpt2d = pnsq*Hwrpt2g
Hwrpt2c = pn*Hwrpt2e
Hwrpt2b = expMUsqsqrtQplusQ*sv
Hwrpt2a = xi2*(Hwrpt2b + Hwrpt2c - Hwrpt2d)
Hwrpt2 = deltaT*Hwrpt2a
Hwrpt1b = invtwoB1psqrtQsqrtQ*invxi2
Hwrpt1a = sqrtdeltaR*Hwrpt1b
Hwrpt1 = invexpMUcubinvexpnu*Hwrpt1a
Hwr = (Hwrpt4 - Hwrpt3 + Hwrpt2)*Hwrpt1
Hwcospt9 = pxir*sxi
Hwcospt8 = pvr*sv
Hwcospt7 = (B*Hwcospt8 - (expMUexpnu)*Hwcospt9)
Hwcospt6 = sqrtdeltaR*Hwcospt7
Hwcospt5 = (pvrsq - expMUsqsqrtQplusQ*xi2)
Hwcospt4 = pn*Hwcospt6
Hwcospt3 = -(expMUsqexpnusq*pxirsq) + deltaT*Hwcospt5
Hwcospt2 = sn*Hwcospt3 - B*Hwcospt4
Hwcospt1 = invexpMUcubinvexpnu*Hwcospt2
Hwcos = invtwoB1psqrtQsqrtQ*Hwcospt1
deltaTsqrtQ = deltaT*sqrtQ
invdeltatTsqrtQ = 1/deltaTsqrtQ
HSOLpt5 = (-B + (expMUexpnu))*pxir
HSOLpt4 = invexpMU*HSOLpt5
HSOLpt3 = expnusq*HSOLpt4
HSOLpt2 = (HSOLpt3*s3)
HSOLpt1 = HSOLpt2*invxi2
HSOL = HSOLpt1*invdeltatTsqrtQ
deltaTsqrtQplusQ = (deltaT*(sqrtQ + Q))
invdeltaTsqrtQplusQ = 1/deltaTsqrtQplusQ
HSONLmult2 = invxi2*invdeltaTsqrtQplusQ
HSONLmult = expnuinvexpMU2*HSONLmult2
HSONLpt1b = pn*xi2
HSONLpt1a = (mur*pvr - nur*pvr + (-mucos + nucos)*HSONLpt1b)
HSONLpt1 = mur*pvr - (mucos*HSONLpt1b) + sqrtQ*HSONLpt1a
HSONLpt2d = nur*pxir
HSONLpt2c = oneplus2sqrtQ*HSONLpt2d
HSONLpt2b = B*sxi
HSONLpt2a = expMUexpnu*HSONLpt2c
HSONLpt2 = (sv*HSONLpt2a + HSONLpt1*HSONLpt2b)
HSONLpt3c = sv*pxir
HSONLpt3b = oneplus1sqrtQ*HSONLpt3c
HSONLpt3a = expMUexpnu*HSONLpt3b
HSONLpt3 = -BR*HSONLpt3a + B*HSONLpt2
HSONLpt4e = sn*xi2
HSONLpt4d = oneplus2sqrtQ*HSONLpt4e
HSONLpt4c = pxir*HSONLpt4d
HSONLpt4b = nucos*HSONLpt4c
HSONLpt4a = expMUexpnu*HSONLpt4b
HSONLpt4 = (-(B*HSONLpt4a) + HSONLpt3*sqrtdeltaR)
HSONL = HSONLmult*HSONLpt4
Hs = w*s3 + Hwr*wr + Hwcos*wcos + HSOL + HSONL
Hsspt1 = sp.Rational(-1,2)*(sx*sx + sy*sy + sz*sz - 3*sn*sn)
Hss = u3*Hsspt1
sKerrdotsStar = (sKerrUSCOREx*sStarUSCOREx + sKerrUSCOREy*sStarUSCOREy + sKerrUSCOREz*sStarUSCOREz)
Hpt1 = etau4*(s1dots1 + s2dots2)
H = Hns + Hs + Hss + (dheffSS*sKerrdotsStar + dheffSSv2)*Hpt1
Hreal = sp.sqrt(1 + 2*eta*(H - 1))
"""

    f = open("Hamstring.txt", 'r')
    Hamstring1 = str(f.read())
    f.close()
    Hamstring = Hamstring1

    # Split Hamstring by carriage returns:
    Hamterms = Hamstring.splitlines()

    # Create "lr" array, which will store each left-hand side and right-hand side of Hamterms as strings.
    lr = []
    # Loop over each line in Hamstring to separate the left- and right-hand sides.
    for i in range(len(Hamterms)):
        # Ignore lines with 2 or fewer characters and those starting with #
        if len(Hamterms[i]) > 2 and Hamterms[i][0] != "#":
            # Split each line by its equals sign.
            splitHamterms = Hamterms[i].split("=")
            # Append to the "lr" array, removing spaces, "sp." prefixes, and replacing Lambda->Lamb
            #       (Lambda is a protected keyword):
            lr.append(lhrh(lhs=splitHamterms[0].replace(" ", "").replace("Lambda", "Lamb"),
                           rhs=splitHamterms[1].replace(" ", "").replace("sp.", "").replace("Lambda", "Lamb")))

    xx = sp.Symbol('xx')
    func = []
    lhss = []
    rhss = []
    # Affix '(xx)' to each left-hand side as a function designation.
    for i in range(len(lr)):
        func.append(sp.sympify(sp.Function(lr[i].lhs)(xx)))
        lhss.append(sp.sympify(lr[i].lhs))
        rhss.append(sp.sympify(lr[i].rhs))

    # Generate a list of all the "free symbols" in the RHS expressions.
    full_symbol_list_with_dups = []
    for i in range(len(lr)):
        for var in rhss[i].free_symbols:
            full_symbol_list_with_dups.append(var)

    # Remove all duplicated "free symbols" from the RHS expressions.
    full_symbol_list = superfast_uniq(full_symbol_list_with_dups)

    # Declare input constants.
    m1, m2, eta = sp.symbols("m1 m2 eta", real=True)
    c0k2, c1k2, c0k3, c1k3, c0k4, c1k4, c2k4, c0k5, c1k5, c2k5 = sp.symbols(
        "c0k2 c1k2 c0k3 c1k3 c0k4 c1k4 c2k4 c0k5 c1k5 c2k5", real=True)
    KK, k5l, b3, bb3, d1, d1v2, dheffSS, dheffSSv2 = sp.symbols("KK k5l b3 bb3 d1 d1v2 dheffSS dheffSSv2", real=True)
    tortoise, copysignresult = sp.symbols("tortoise copysignresult", real=True)
    input_constants = [m1, m2, eta,
                       c0k2, c1k2, c0k3, c1k3, c0k4, c1k4, c2k4, c0k5, c1k5, c2k5,
                       KK, k5l, b3, bb3, d1, d1v2, dheffSS, dheffSSv2,
                       tortoise, copysignresult]

    # Derivatives of input constants will always be zero, so remove them from the full_symbol_list.
    for inputconst in input_constants:
        for symbol in full_symbol_list:
            if str(symbol) == str(inputconst):
                full_symbol_list.remove(symbol)

    # Add symbols to the function list and replace right-hand side terms with their function equivalent.
    full_function_list = []
    for symb in full_symbol_list:
        func = sp.sympify(sp.Function(str(symb))(xx))
        full_function_list.append(func)
        for i in range(len(rhss)):
            for var in rhss[i].free_symbols:
                if str(var) == str(symb):
                    rhss[i] = rhss[i].subs(var, func)

    # Differentiate with respect to xx, remove '(xx)', and replace xx with 'prm' notation
    lhss_deriv = []
    rhss_deriv = []
    for i in range(len(rhss)):
        lhss_deriv.append(sp.sympify(str(lhss[i]) + "prm"))
        newrhs = sp.sympify(
            str(sp.diff(rhss[i], xx)).replace("(xx)", "").replace(", xx", "prm").replace("Derivative", ""))
        rhss_deriv.append(newrhs)

    lhss_deriv_simp, rhss_deriv_simp = simplify_deriv(lhss_deriv, rhss_deriv)
    lhss_deriv = lhss_deriv_simp
    rhss_deriv = rhss_deriv_simp

    lhss_deriv_x, rhss_deriv_x = deriv_onevar(lhss_deriv, rhss_deriv,
                                              xprm=1, yprm=0, zprm=0, pxprm=0, pyprm=0, pzprm=0, s1xprm=0, s1yprm=0,
                                              s1zprm=0, s2xprm=0, s2yprm=0, s2zprm=0)
    lhss_deriv_y, rhss_deriv_y = deriv_onevar(lhss_deriv, rhss_deriv,
                                              xprm=0, yprm=1, zprm=0, pxprm=0, pyprm=0, pzprm=0, s1xprm=0, s1yprm=0,
                                              s1zprm=0, s2xprm=0, s2yprm=0, s2zprm=0)
    lhss_deriv_z, rhss_deriv_z = deriv_onevar(lhss_deriv, rhss_deriv,
                                              xprm=0, yprm=0, zprm=1, pxprm=0, pyprm=0, pzprm=0, s1xprm=0, s1yprm=0,
                                              s1zprm=0, s2xprm=0, s2yprm=0, s2zprm=0)
    lhss_deriv_px, rhss_deriv_px = deriv_onevar(lhss_deriv, rhss_deriv,
                                                xprm=0, yprm=0, zprm=0, pxprm=1, pyprm=0, pzprm=0, s1xprm=0, s1yprm=0,
                                                s1zprm=0, s2xprm=0, s2yprm=0, s2zprm=0)
    lhss_deriv_py, rhss_deriv_py = deriv_onevar(lhss_deriv, rhss_deriv,
                                                xprm=0, yprm=0, zprm=0, pxprm=0, pyprm=1, pzprm=0, s1xprm=0, s1yprm=0,
                                                s1zprm=0, s2xprm=0, s2yprm=0, s2zprm=0)
    lhss_deriv_pz, rhss_deriv_pz = deriv_onevar(lhss_deriv, rhss_deriv,
                                                xprm=0, yprm=0, zprm=0, pxprm=0, pyprm=1, pzprm=1, s1xprm=0, s1yprm=0,
                                                s1zprm=0, s2xprm=0, s2yprm=0, s2zprm=0)
    lhss_deriv_s1x, rhss_deriv_s1x = deriv_onevar(lhss_deriv, rhss_deriv,
                                                  xprm=0, yprm=0, zprm=0, pxprm=0, pyprm=0, pzprm=0, s1xprm=1, s1yprm=0,
                                                  s1zprm=0, s2xprm=0, s2yprm=0, s2zprm=0)
    lhss_deriv_s1y, rhss_deriv_s1y = deriv_onevar(lhss_deriv, rhss_deriv,
                                                  xprm=0, yprm=0, zprm=0, pxprm=0, pyprm=0, pzprm=0, s1xprm=0, s1yprm=1,
                                                  s1zprm=0, s2xprm=0, s2yprm=0, s2zprm=0)
    lhss_deriv_s1z, rhss_deriv_s1z = deriv_onevar(lhss_deriv, rhss_deriv,
                                                  xprm=0, yprm=0, zprm=0, pxprm=0, pyprm=0, pzprm=0, s1xprm=0, s1yprm=0,
                                                  s1zprm=1, s2xprm=0, s2yprm=0, s2zprm=0)
    lhss_deriv_s2x, rhss_deriv_s2x = deriv_onevar(lhss_deriv, rhss_deriv,
                                                  xprm=0, yprm=0, zprm=0, pxprm=0, pyprm=0, pzprm=0, s1xprm=0, s1yprm=0,
                                                  s1zprm=0, s2xprm=1, s2yprm=0, s2zprm=0)
    lhss_deriv_s2y, rhss_deriv_s2y = deriv_onevar(lhss_deriv, rhss_deriv,
                                                  xprm=0, yprm=0, zprm=0, pxprm=0, pyprm=0, pzprm=0, s1xprm=0, s1yprm=0,
                                                  s1zprm=0, s2xprm=0, s2yprm=1, s2zprm=0)
    lhss_deriv_s2z, rhss_deriv_s2z = deriv_onevar(lhss_deriv, rhss_deriv,
                                                  xprm=0, yprm=0, zprm=0, pxprm=0, pyprm=0, pzprm=0, s1xprm=0, s1yprm=0,
                                                  s1zprm=0, s2xprm=0, s2yprm=0, s2zprm=1)

    outstring = "/* SEOBNR Hamiltonian expression: */\n"
    outstringsp = ""
    outsplhs = []
    outsprhs = []
    for i in range(len(lr)):
        outstring += outputC(sp.sympify(lr[i].rhs), lr[i].lhs, "returnstring",
                             "outCverbose=False,includebraces=False,CSE_enable=False")
        outstringsp += lr[i].lhs + " = " + lr[i].rhs + "\n"
        outsplhs.append(sp.sympify(lr[i].lhs))
        outsprhs.append(sp.sympify(lr[i].rhs))
    outstring += "\n\n\n/* SEOBNR \partial_x H expression: */\n"
    for i in range(len(lhss_deriv_x)):
        outstring += outputC(rhss_deriv_x[i], str(lhss_deriv_x[i]), "returnstring",
                             "outCverbose=False,includebraces=False,CSE_enable=False")
        outstringsp += str(lhss_deriv_x[i]) + " = " + str(rhss_deriv_x[i]) + "\n"
        outsplhs.append(lhss_deriv_x[i])
        outsprhs.append(rhss_deriv_x[i])

    with open("/tmp/sympy_expression.py", "w") as file:
        file.write("""
import sympy as sp
from outputC import *

m1,m2,x,y,z,px,py,pz,s1x,s1y,s1z,s2x,s2y,s2z = sp.symbols("m1 m2 x y z px py pz s1x s1y s1z s2x s2y s2z",real=True)
c0k2,c1k2,c0k3,c1k3,c0k4 = sp.symbols("c0k2 c1k2 c0k3 c1k3 c0k4",real=True)
c1k4,c2k4,c0k5,c1k5,c2k5 = sp.symbols("c1k4 c2k4 c0k5 c1k5 c2k5",real=True)
eta,KK,k5l,b3,bb3,d1,d1v2,dheffSS,dheffSSv2 = sp.symbols("eta KK k5l b3 bb3 d1 d1v2 dheffSS dheffSSv2",real=True)
tortoise = sp.symbols("tortoise",real=True)

""")
        for i in range(len(lr)):
            file.write(lr[i].lhs + " = " + "sp.symbols(\"" + lr[i].lhs + "\")\n")
        file.write("\n")
        for i in range(len(lhss_deriv_x)):
            file.write(str(lhss_deriv_x[i]).replace("prm", "prm_x") + " = " + str(rhss_deriv_x[i]).replace("sqrt(",
                "sp.sqrt(").replace("log(", "sp.log(").replace("Abs(", "sp.Abs(").replace("sign(", "sp.sign(").replace(
                "prm", "prm_x") + "\n")
        for i in range(len(lhss_deriv_y)):
            file.write(str(lhss_deriv_y[i]).replace("prm", "prm_y") + " = " + str(rhss_deriv_y[i]).replace("sqrt(",
                "sp.sqrt(").replace("log(", "sp.log(").replace("Abs(", "sp.Abs(").replace("sign(", "sp.sign(").replace(
                "prm", "prm_y") + "\n")
        for i in range(len(lhss_deriv_z)):
            file.write(str(lhss_deriv_z[i]).replace("prm", "prm_z") + " = " + str(rhss_deriv_z[i]).replace("sqrt(",
                "sp.sqrt(").replace("log(", "sp.log(").replace("Abs(", "sp.Abs(").replace("sign(", "sp.sign(").replace(
                "prm", "prm_z") + "\n")

        for i in range(len(lhss_deriv_px)):
            file.write(str(lhss_deriv_px[i]).replace("prm", "prm_px") + " = " + str(rhss_deriv_px[i]).replace("sqrt(",
                "sp.sqrt(").replace("log(", "sp.log(").replace("Abs(", "sp.Abs(").replace("sign(", "sp.sign(").replace(
                "prm", "prm_px") + "\n")
        for i in range(len(lhss_deriv_py)):
            file.write(str(lhss_deriv_py[i]).replace("prm", "prm_py") + " = " + str(rhss_deriv_py[i]).replace("sqrt(",
                "sp.sqrt(").replace("log(", "sp.log(").replace("Abs(", "sp.Abs(").replace("sign(", "sp.sign(").replace(
                "prm", "prm_py") + "\n")
        for i in range(len(lhss_deriv_pz)):
            file.write(str(lhss_deriv_pz[i]).replace("prm", "prm_pz") + " = " + str(rhss_deriv_pz[i]).replace("sqrt(",
                "sp.sqrt(").replace("log(", "sp.log(").replace("Abs(", "sp.Abs(").replace("sign(", "sp.sign(").replace(
                "prm", "prm_pz") + "\n")

        for i in range(len(lhss_deriv_s1x)):
            file.write(
                str(lhss_deriv_s1x[i]).replace("prm", "prm_s1x") + " = " + str(rhss_deriv_s1x[i]).replace("sqrt(",
                    "sp.sqrt(").replace("log(", "sp.log(").replace("Abs(", "sp.Abs(").replace("sign(",
                    "sp.sign(").replace("prm", "prm_s1x") + "\n")
        for i in range(len(lhss_deriv_s1y)):
            file.write(
                str(lhss_deriv_s1y[i]).replace("prm", "prm_s1y") + " = " + str(rhss_deriv_s1y[i]).replace("sqrt(",
                    "sp.sqrt(").replace("log(", "sp.log(").replace("Abs(", "sp.Abs(").replace("sign(",
                    "sp.sign(").replace("prm", "prm_s1y") + "\n")
        for i in range(len(lhss_deriv_s1z)):
            file.write(
                str(lhss_deriv_s1z[i]).replace("prm", "prm_s1z") + " = " + str(rhss_deriv_s1z[i]).replace("sqrt(",
                    "sp.sqrt(").replace("log(", "sp.log(").replace("Abs(", "sp.Abs(").replace("sign(",
                    "sp.sign(").replace("prm", "prm_s1z") + "\n")

        for i in range(len(lhss_deriv_s2x)):
            file.write(
                str(lhss_deriv_s2x[i]).replace("prm", "prm_s2x") + " = " + str(rhss_deriv_s2x[i]).replace("sqrt(",
                    "sp.sqrt(").replace("log(", "sp.log(").replace("Abs(", "sp.Abs(").replace("sign(",
                    "sp.sign(").replace("prm", "prm_s2x") + "\n")
        for i in range(len(lhss_deriv_s2y)):
            file.write(
                str(lhss_deriv_s2y[i]).replace("prm", "prm_s2y") + " = " + str(rhss_deriv_s2y[i]).replace("sqrt(",
                    "sp.sqrt(").replace("log(", "sp.log(").replace("Abs(", "sp.Abs(").replace("sign(",
                    "sp.sign(").replace("prm", "prm_s2y") + "\n")
        for i in range(len(lhss_deriv_s2z)):
            file.write(
                str(lhss_deriv_s2z[i]).replace("prm", "prm_s2z") + " = " + str(rhss_deriv_s2z[i]).replace("sqrt(",
                    "sp.sqrt(").replace("log(", "sp.log(").replace("Abs(", "sp.Abs(").replace("sign(",
                    "sp.sign(").replace("prm", "prm_s2z") + "\n")
        file.write("""
CSE_results = sp.cse(Hrealprm_x, sp.numbered_symbols("tmp"), order='canonical')
print(CSE_results[0])
print(CSE_results[1])
#for commonsubexpression in CSE_results[0]:
#    print(str(commonsubexpression[0])+" = "+str(commonsubexpression[1]))
#for i,result in enumerate(CSE_results[1]):
#    print(str(CSE_results[i][0])+" = "+str(result))

#outputC([Hreal,Hrealprm_x,Hrealprm_y,Hrealprm_z,Hrealprm_px,Hrealprm_py,Hrealprm_pz,Hrealprm_s1x,Hrealprm_s1y,Hrealprm_s1z,Hrealprm_s2x,Hrealprm_s2y,Hrealprm_s2z],
#["Hreal","Hrealprm_x","Hrealprm_y","Hrealprm_z","Hrealprm_px","Hrealprm_py","Hrealprm_pz",
#"Hrealprm_s1x","Hrealprm_s1y","Hrealprm_s1z","Hrealprm_s2x","Hrealprm_s2y","Hrealprm_s2z"],
#"/tmp/outC.h","outCverbose=False,includebraces=False")
""")

    return 0

output_H_and_derivs()
