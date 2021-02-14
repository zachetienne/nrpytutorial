# As documented in the NRPy+ tutorial module
#   Tutorial-SEOBNR_v3_Hamiltonian.ipynb,
#   this module will compute the numerical
#   value of the real Hamiltonian.  Note that
#   terms are in the reverse order of the
#   notebook.

# Author: Tyler Knowles & Zachariah B. Etienne
#         zachetie **at** gmail **dot* com

# Import all needed modules from NRPy+/Python:
import numpy as np

def compute_Hreal(m1=23., m2=10., EMgamma=0.577215664901532860606512090082402431, tortoise=1,
                  x=2.129681018601393e+01, y=0.000000000000000e+00, z=0.000000000000000e+00,
                  p1=0.000000000000000e+00, p2=2.335391115580442e-01, p3=-4.235164736271502e-22,
                  S1x=4.857667584940312e-03, S1y=9.715161660389764e-03, S1z=-1.457311842632286e-02,
                  S2x=3.673094582185491e-03, S2y=-4.591302628615413e-03, S2z=5.509696538546906e-03):
    # Fundamental Quantities
    M = m1 + m2
    mu=m1*m2/M
    eta=mu/M
    r=np.sqrt(x*x+y*y+z*z)
    u=1/r

    # Spin Combinations
    sigmastar3=m2/m1*S1z+m1/m2*S2z
    sigmastar2 = m2/m1*S1y + m1/m2*S2y
    sigmastar1 = m2/m1*S1x + m1/m2*S2x
    sigma3=S1z+S2z
    sigma2 = S1y + S2y
    sigma1 = S1x + S2x
    Skerr3=sigma3
    Skerr2 = sigma2
    Skerr1 = sigma1
    Skerrmag=np.sqrt(Skerr1*Skerr1+Skerr2*Skerr2+Skerr3*Skerr3)
    Skerrhat3=Skerr3/Skerrmag
    Skerrhat2 = Skerr2/Skerrmag
    Skerrhat1 = Skerr1/Skerrmag
    a=Skerrmag

    # Important Vectors
    n3=z/r
    n2 = y/r
    n1 = x/r
    e33=Skerrhat3
    e32 = Skerrhat2
    e31 = Skerrhat1
    xi3=e31*n2-e32*n1
    xi2 = e31*n3 + e33*n1
    xi1 = e32*n3 - e33*n2
    v3=n1*xi2-n2*xi1
    v2 = n3*xi1 - n1*xi3
    v1 = n2*xi3 - n3*xi2

    # Terms Dependent on Coordinates
    costheta=e31*n1+e32*n2+e33*n3
    sin2theta=1-costheta*costheta
    xisq = sin2theta
    w2=a*a+r*r
    Sigma=r*r+a*a*costheta*costheta

    # Metric Terms
    Dinv=1+np.log(1+6*eta*u*u+2*(26-3*eta)*eta*u*u*u)
    omegatilde=2*a*r
    K=1.712-1.803949138004582*eta-39.77229225266885*eta*eta+103.16588921239249*eta*eta*eta
    etaKminus1 = eta*K - 1
    Delta0=K*(eta*K-2)
    Delta1 = -2*etaKminus1*(K + Delta0)
    Delta2 = np.divide(1,2)*Delta1*(Delta1 - 4*etaKminus1) - a*a*etaKminus1*etaKminus1*Delta0
    Delta3=-np.divide(1,3)*Delta1*Delta1*Delta1+etaKminus1*Delta1*Delta1+Delta2*Delta1-2*etaKminus1*(Delta2-etaKminus1)-a*a*etaKminus1*etaKminus1*Delta1
    Delta4=np.divide(1,12)*(6*a*a*(Delta1*Delta1-2*Delta2)*etaKminus1*etaKminus1+3*Delta1*Delta1*Delta1*Delta1-8*etaKminus1*Delta1*Delta1*Delta1-12*Delta2*Delta1*Delta1+12*(2*etaKminus1*Delta2+Delta3)*Delta1+12*(np.divide(94,3)-np.divide(41,32)*np.pi*np.pi)*etaKminus1*etaKminus1+6*(Delta2*Delta2-4*Delta3*etaKminus1))
    Delta5=etaKminus1*etaKminus1*((np.divide(-4237,60)+np.divide(128,5)*EMgamma+np.divide(2275,512)*np.pi*np.pi-np.divide(1,3)*a*a*(Delta1*Delta1*Delta1-3*Delta1*Delta2+3*Delta3)-(Delta1*Delta1*Delta1*Delta1*Delta1-5*Delta1*Delta1*Delta1*Delta2+5*Delta1*Delta2*Delta2+5*Delta1*Delta1*Delta3-5*Delta2*Delta3-5*Delta1*Delta4)/(5*etaKminus1*etaKminus1)+(Delta1*Delta1*Delta1*Delta1-4*Delta1*Delta1*Delta2+2*Delta2*Delta2+4*Delta1*Delta3-4*Delta4)/(2*etaKminus1)+np.divide(256,5)*np.log(2)))
    Delta5l = etaKminus1*etaKminus1*np.divide(64,5)
    logarg=u*(Delta1+u*(Delta2+u*(Delta3+u*(Delta4+u*(Delta5+Delta5l*np.log(u))))))
    Deltaucalib = 1 + eta*(Delta0 + np.log(1 + logarg))
    Deltaucalibprm=-eta*u*u*(Delta1+u*(2*Delta2+u*(3*Delta3+u*(4*Delta4+u*(5*(Delta5+Delta5l*np.log(u)))))))/(1+logarg)
    Deltaubar=a*a*u*u+(2*u+1/etaKminus1)/etaKminus1
    Deltaubarprm = -2*a*a*u*u*u - 2*u*u/(etaKminus1)
    Deltau=Deltaubar*Deltaucalib
    Deltauprm = Deltaubarprm*Deltaucalib + Deltaubar*Deltaucalibprm
    Deltatprm=2*r*Deltau+r*r*Deltauprm
    Deltat=r*r*Deltau
    Deltar=Deltat*Dinv
    Lambdat=w2*w2-a*a*Deltat*sin2theta

    # Tortoise terms
    csi=np.sqrt(Deltar*Deltat)/w2
    csi1=1+(1-np.abs(1-tortoise))*(csi-1)
    csi2=1+(np.divide(1,2)-np.divide(1,2)*np.sign(np.divide(3,2)-tortoise))*(csi-1)
    prT=csi2*(p1*n1+p2*n2+p3*n3)
    phat3=p3+prT*(1-1/csi1)*n3
    phat2 = p2 + prT*(1 - 1/csi1)*n2
    phat1 = p1 + prT*(1 - 1/csi1)*n1
    pdotxir=(phat1*xi1+phat2*xi2+phat3*xi3)*r
    pdotn=phat1*n1+phat2*n2+phat3*n3
    pdotvr=(phat1*v1+phat2*v2+phat3*v3)*r
    pphi=pdotxir
    Qcoeff2=1/(Sigma*sin2theta)
    Qcoeff1=Sigma/(Lambdat*sin2theta)
    DrSipn2=Deltar*pdotn*pdotn/Sigma

    # The Deformed and Rescaled Metric Potentials
    Q=1+DrSipn2+Qcoeff1*pdotxir*pdotxir+Qcoeff2*pdotvr*pdotvr
    Qminus1 = Q - 1
    Jtilde=np.sqrt(Deltar)
    exp2mu=Sigma
    expmu = np.sqrt(exp2mu)
    Brtilde=(np.sqrt(Deltar)*Deltatprm-2*Deltat)/(2*np.sqrt(Deltar*Deltat))
    Btilde=np.sqrt(Deltat)
    exp2nu=Deltat*Sigma/Lambdat
    expnu = np.sqrt(exp2nu)
    omega=omegatilde/Lambdat

    # Derivatives of the Metric Potential
    omegatildeprm=2*a
    Lambdatprm=4*(a*a+r*r)*r-2*a*a*Deltatprm*sin2theta
    mucostheta=a*a*costheta/Sigma
    nucostheta=a*a*w2*costheta*(w2-Deltat)/(Lambdat*Sigma)
    omegacostheta=-2*a*a*costheta*Deltat*omegatilde/(Lambdat*Lambdat)
    mur=r/Sigma-1/np.sqrt(Deltar)
    nur=r/Sigma+w2*(w2*Deltatprm-4*r*Deltat)/(2*Lambdat*Deltat)
    omegar=(Lambdat*omegatildeprm-Lambdatprm*omegatilde)/(Lambdat*Lambdat)
    dSO=-74.71-156.*eta+627.5*eta*eta

    # Spin Combinations
    sigmacoeffTerm3 = eta*dSO*u*u*u
    sigmacoeffTerm2=eta/(144*r*r)*(-896+r*(-436*Qminus1-96*DrSipn2+r*(-45*Qminus1*Qminus1+36*Qminus1*DrSipn2))+eta*(-336+r*(204*Qminus1-882*DrSipn2+r*(810*DrSipn2*DrSipn2-234*Qminus1*DrSipn2))))
    sigmacoeffTerm1=eta/12*(-8/r+3*Qminus1-36*DrSipn2)
    sigmacoeff=sigmacoeffTerm1+sigmacoeffTerm2+sigmacoeffTerm3
    sigmastarcoeffTerm2=eta/(72*r*r)*(706+r*(-206*Qminus1+282*DrSipn2+r*Qminus1*(96*DrSipn2-23*Qminus1))+eta*(-54+r*(120*Qminus1-324*DrSipn2+r*(360*DrSipn2*DrSipn2+Qminus1*(-126*DrSipn2-3*Qminus1)))))
    sigmastarcoeffTerm1=eta/12*(14/r+4*Qminus1-30*DrSipn2)
    sigmastarcoeff=sigmastarcoeffTerm1+sigmastarcoeffTerm2
    Deltasigmastar3=sigmastar3*sigmastarcoeff+sigma3*sigmacoeff
    Deltasigmastar2 = sigmastar2*sigmastarcoeff + sigma2*sigmacoeff
    Deltasigmastar1 = sigmastar1*sigmastarcoeff + sigma1*sigmacoeff
    Sstar3=sigmastar3+Deltasigmastar3
    Sstar2 = sigmastar2 + Deltasigmastar2
    Sstar1 = sigmastar1 + Deltasigmastar1
    S3 = Sstar3
    S2 = Sstar2
    S1 = Sstar1

    # Common Dot Products
    Sstardotn=Sstar1*n1+Sstar2*n2+Sstar3*n3
    SdotSkerrhat=S1*Skerrhat1+S2*Skerrhat2+S3*Skerrhat3
    Sdotn=S1*n1+S2*n2+S3*n3
    Sdotv=S1*v1+S2*v2+S3*v3
    Sdotxi=S1*xi1+S2*xi2+S3*xi3

    # The $H_{\rm D}$ Terms
    HdsumTerm2=3*Sstardotn*Sstardotn
    HdsumTerm1=Sstar1*Sstar1+Sstar2*Sstar2+Sstar3*Sstar3
    Hdsum=HdsumTerm1-HdsumTerm2
    Hdcoeff=np.divide(1,2)/(r*r*r)
    Q4=2*prT*prT*prT*prT*u*u*(4-3*eta)*eta
    gammappsum=Deltar/Sigma*pdotn*pdotn+1/Sigma*pdotvr*pdotvr/sin2theta+Sigma/Lambdat/sin2theta*pdotxir*pdotxir

    # The $H_{\rm NS}$ Terms
    Hnsradicand=1+gammappsum+Q4
    alpha=np.sqrt(Deltat*Sigma/Lambdat)
    betapsum=omegatilde*pphi/Lambdat

    # The Spin-Spin Term $H_{\rm SS}$
    HssTerm3=expmu*expnu*pdotxir*(Jtilde*pdotn*Sdotxi*Btilde-expmu*expnu*pdotxir*Sdotn)+(pdotvr*(Sdotn*pdotvr-Jtilde*pdotn*Sdotv)-exp2mu*(np.sqrt(Q)+Q)*Sdotn*xisq)*Btilde*Btilde
    HssTerm3coeff=omegacostheta/(2*exp2mu*expmu*expnu*Btilde*(Q+np.sqrt(Q)))
    HssTerm2=expmu*pdotxir*(expmu*exp2nu*pdotxir*Sdotv-expnu*pdotvr*Sdotxi*Btilde)+xisq*Btilde*Btilde*(exp2mu*(np.sqrt(Q)+Q)*Sdotv+Jtilde*pdotn*(pdotvr*Sdotn-Jtilde*pdotn*Sdotv))
    HssTerm2coeff=Jtilde*omegar/(2*exp2mu*expmu*expnu*Btilde*(Q+np.sqrt(Q))*xisq)
    HssTerm1=omega*SdotSkerrhat
    Hss=HssTerm1+HssTerm2coeff*HssTerm2+HssTerm3coeff*HssTerm3

    # The Spin-Orbit Term $H_{\rm SO}$
    HsoTerm2c=Jtilde*Brtilde*expmu*expnu*pdotxir*(np.sqrt(Q)+1)*Sdotv
    HsoTerm2b=expmu*expnu*pdotxir*(2*np.sqrt(Q)+1)*(Jtilde*nur*Sdotv-nucostheta*Sdotn*xisq)*Btilde
    HsoTerm2a=Sdotxi*Jtilde*(mur*pdotvr*(np.sqrt(Q)+1)-mucostheta*pdotn*xisq-np.sqrt(Q)*(nur*pdotvr+(mucostheta-nucostheta)*pdotn*xisq))*Btilde*Btilde
    HsoTerm2=HsoTerm2a+HsoTerm2b-HsoTerm2c
    HsoTerm2coeff=expnu/(exp2mu*Btilde*Btilde*(Q+np.sqrt(Q))*xisq)
    HsoTerm1=exp2nu*(expmu*expnu-Btilde)*pdotxir*SdotSkerrhat/(expmu*Btilde*Btilde*np.sqrt(Q)*xisq)
    Hso=HsoTerm1+HsoTerm2coeff*HsoTerm2

    # Terms of $H_{\rm eff}$
    Hd=Hdcoeff*Hdsum
    Hns=betapsum+alpha*np.sqrt(Hnsradicand)
    Hs=Hso+Hss
    dSS=8.127-154.2*eta+830.8*eta*eta

    # The Effective Hamiltonian $H_{\rm eff}$
    Heff = Hs + Hns - Hd + dSS*eta*u*u*u*u*(S1x*S1x + S1y*S1y + S1z*S1z + S2x*S2x + S2y*S2y + S2z*S2z)

    # The Real Hamiltonian $H_{\rm real}$
    Hreal=np.sqrt(1+2*eta*(Heff-1))

    return Hreal
