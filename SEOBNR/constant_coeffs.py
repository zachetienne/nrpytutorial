# -*- coding: utf-8 -*-
# # The Spinning Effective One-Body Constant Coefficients
# ## Author: Tyler Knowles
#
# ## This module documents the reduced spinning effective one-body constant coefficients (both from numerical relativity and Post-Newtonian expansions) as implemented in LALSuite's SEOBNRv3 gravitational waveform approximant.  <font color='red'><b>FIXME: add a primary source.</b></font>
#
#
# **Module Status:** <font color='red'><b> In progress </b></font>
#
# **Validation Notes:** This module is under active development -- do ***not*** use the resulting code for scientific applications.  In the future, this module will be validated against the LALSuite [SEOBNRv3/SEOBNRv3_opt code]( https://git.ligo.org/lscsoft/lalsuite.) that was reviewed and approved for LIGO parameter estimation by the LIGO Scientific Collaboration.
#
#
# ## Introduction
# ### The Physical System of Interest
#
# Consider two compact objects (e.g. black holes or neutron stars) with masses $m_{1}$, $m_{2}$ (in solar masses) and spin angular momenta ${\bf S}_{1}$, ${\bf S}_{2}$ in a binary system.  The spinning effective one-body ("SEOB") Hamiltonian $H_{\rm real}$ (see [BB2010](https://arxiv.org/abs/0912.3517) Equation (5.69)) describes the dynamics of this system.
#
# The constant coefficients from Post-Newtonian expansions and fits to Numerical Relativity we need to compute the Hamiltonian of this system rely on
# 1. the symmetric mass ratio $\eta$.
#
# Besides $\eta$, we also need the [Euler–Mascheroni constant](https://en.wikipedia.org/wiki/Euler%E2%80%93Mascheroni_constant) $\gamma$ which is hard-coded in LALSuite with the significant digits shown below.  (The following link directly to the appropriate LALSuite documentation: [$\gamma$](https://lscsoft.docs.ligo.org/lalsuite/lal/group___l_a_l_constants__h.html#gac6af32574ff8acaeeafe8bf422281e98).)
#
# Please note that throughout this notebook we adpot the following conventions:
# 1. $c = G = 1$ where $c$ is the speed of light in a vacuum and $G$ is Newton's gravitational constant, and
# 1. $m_{1} \ge m_{2}$.
#
# LALSuite line numbers are taken from Git commit bba40f2 (see [LALSuite's GitLab page](https://git.ligo.org/lscsoft/lalsuite)).

# +
import numpy as np

def compute_const_coeffs(eta, gamma, a):
    
    # We need to find a better primary source for these coefficeints than is listed in the LALSuite code
    # The paper Taracchini, Pan, et. al (2012) (https://arxiv.org/pdf/1202.0790.pdf) does not have as
    # many significant digits as are listed here
    c0  = 1.446
    c1  = -1.7152360250654402
    c2  = -3.246255899738242

    # There is no source for these coefficents listed in the LALSuite code
    c20 = 1.712
    c21 = -1.803949138004582
    c22 = -39.77229225266885
    c23 = 103.16588921239249

    KK = c20 + c21*eta + c22*eta*eta + c23*eta*eta*eta;                                                                                                                                                 
    EtaKKm1 = eta*KK - 1.

    # Tmp variables from Mathematica CSE on the coefficient expressions
    tmp4 = EtaKKm1*EtaKKm1*EtaKKm1
    tmp10 = EtaKKm1*EtaKKm1
    tmp7 = EtaKKm1*EtaKKm1*EtaKKm1*EtaKKm1
    tmp6 = KK*KK
    tmp16 = EtaKKm1**5
    tmp19 = EtaKKm1**6
    tmp18 = KK*KK*KK
    tmp23 = 4.*KK*tmp7
    tmp34 = EtaKKm1**7
    tmp28 = np.pi*np.pi
    tmp32 = 16.*KK*tmp16
    tmp37 = EtaKKm1**8
    tmp36 = KK*KK*KK*KK
    tmp64 = EtaKKm1**9

    kC0 = 0. + 4.*KK*tmp4 + 2.*tmp6*tmp7
    kC1 = 1.*KK*tmp10 - KK*tmp4
    kC2 = 0. + 2.*tmp10 - 1.3333333333333335*tmp18*tmp19 - 8.*tmp16*tmp6-8.*KK*tmp7
    kC3 = tmp23-2.*KK*tmp4 + 2.*tmp16*tmp6 - 2.*tmp6*tmp7
    kC4 = 0. + 31.333333333333332*tmp10 - 1.28125*tmp10*tmp28 + tmp32 + 8.*tmp18*tmp34 + 0.6666666666666666*tmp36*tmp37 - 4.*tmp4 + 24.*tmp19*tmp6 - 4.*KK*tmp7
    kC5 = -12.*KK*tmp16 + 2.*tmp18*tmp19 + tmp23 - 2.*tmp18*tmp34 + 8.*tmp16*tmp6 - 12.*tmp19*tmp6
    kC6 = 1.*KK*tmp16 - tmp16*tmp6 + 0.5*tmp19*tmp6 - KK*tmp7 + 0.5*tmp6*tmp7
    kC7 = -35.12753102199746*tmp10 + 25.6*gamma*tmp10 - 32.*KK*tmp19 + 4.443359375*tmp10*tmp28 + tmp32 - 32.*tmp18*tmp37 - 62.666666666666664*tmp4 + 2.5625*tmp28*tmp4 + 4.*tmp19*tmp6 - 64.*tmp34*tmp6 - 5.333333333333334*tmp36*tmp64 + 8.*tmp7 - 62.666666666666664*KK*tmp7 + 2.5625*KK*tmp28*tmp7 - 0.2666666666666661*KK**5*EtaKKm1**10
    kC8 = -10.*KK*tmp16 + 32.*KK*tmp19 - 12.*tmp18*tmp34 + 16.*tmp18*tmp37 - 1.3333333333333337*tmp36*tmp37 - 24.*tmp19*tmp6 + 48.*tmp34*tmp6 + 1.3333333333333337*tmp36*tmp64 - 2.*tmp7 + 2.*KK*tmp7
    kC9 = 4.*KK*tmp16 - 6.*KK*tmp19 - tmp18*tmp19 + 2.*tmp18*tmp34 - tmp18*tmp37 - 2.*tmp16*tmp6 + 8.*tmp19*tmp6 - 6.*tmp34*tmp6

    asq = a*a
    aft = asq*asq

    k0 = KK*(EtaKKm1 - 1)
    k1 = -2*(k0 + KK)*EtaKKm1
    k2 = kC0 + kC1*asq
    k3 = kC2 + kC3*asq
    k4 = kC4 + kC5*asq + kC6*aft
    k5 = kC7 + kC8*asq + kC9*aft
    k5l = EtaKKm1*EtaKKm1*64./5.

    d1 = -69.5
    d1v2 = -74.71 - 156.*eta + 627.5*eta*eta
    dheffSS = 2.75
    dheffSSv2 = 8.127 - 154.2*eta + 830.8*eta*eta

    return KK, kC0, kC1, kC2, kC3, kC4, kC5, kC6, kC7, kC8, kC9, k0, k1, k2, k3, k4, k5, k5l, d1, d1v2, dheffSS, dheffSSv2
