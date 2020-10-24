# As documented in the NRPyPN notebook
# PN-p_r.ipynb, this Python script
# generates the expression for the radial
# component of momentum p_r up to and
# including terms at 3.5PN order.
# It largely follows the technique of
#  Ramos-Buades, Husa, and Pratten (2018)
#    https://arxiv.org/abs/1810.00036

# Core functions:
# f_Htot_xyplane_binary(m1,m2, n12U,n21U, S1U, S2U, p1U,p2U, q)
#        Given standard input parameters, compute
#        the Hamiltonian for a binary system
#        orbiting instantaneously on the xy plane,
#        and store to the global variable
#        Htot_xyplane_binary

# f_dr_dt(Htot_xyplane_binary, m1,m2, n12U,n21U, chi1U,chi2U, S1U,S2U, p1U,p2U, r)
#         Given Htot_xyplane_binary (computed
#         above) and other standard input
#         parameters, compute
#         dr_dt = dr/dt and store to global
#         variable of the same name.

# f_p_r(m1,m2, chi1U,chi2U, r)
#       Compute p_r and store to
#       global variable of the same name.

# f_p_r_fullHam(m1,m2, n12U,n21U, chi1U,chi2U, S1U,S2U, p1U,p2U, r)
#       Compute p_r using the full
#       Hamiltonian, without truncating
#       higher-order terms self-consistently.

# Author:  Zach Etienne
#          zachetie **at** gmail **dot* com

# Step 0: Add NRPy's directory to the path
# https://stackoverflow.com/questions/16780014/import-file-from-parent-directory
import sympy as sp                        # SymPy: The Python computer algebra package upon which NRPy+ depends
from NRPyPN_shortcuts import Pt,Pr,nU,div # NRPyPN: shortcuts for e.g., vector operations

#################################
#################################
# Basic equation for p_r can be written in the
#   form:
# p_r \approx [dr/dt - (partial_H/partial_{p_r})|_{p_r=0}] * [(partial^2_{H}/partial_{p_r^2})|_{p_r=0}]^{-1},
#  where
# dr/dt = [dE_{\rm GW}/dt + dM/dt] * [dH_{circ} / dr]^{-1},
#  and
# H_{circ} = Htot_xyplane_binary|_{p_r=0}
#  -> [dH_{circ}(r,p_t(r)) / dr] = partial_{H(p_r=0)}/partial_r
#            + partial_{H(p_r=0)}/partial_{p_t} partial_{p_t}/partial_r.
#  Here,
#  * the expression for p_t is given by PN_p_t.py
#  * the expression for [dE_{\rm GW}/dt + dM/dt] is given by
#            PN_dE_GW_dt_and_dM_dt.py
#    + Since [dE_{\rm GW}/dt + dM/dt] is a function of MOmega,
#        we also need input from the PN_MOmega.py Python module.

# Step 1: Construct full Hamiltonian
#         expression for a binary instantaneously
#         orbiting on the xy plane, store
#         result to Htot_xyplane_binary
def f_Htot_xyplane_binary(m1, m2, n12U, n21U, S1U, S2U, p1U, p2U, r):
    def make_replacements(expr):
        zero = sp.sympify(0)
        one = sp.sympify(1)
        return expr.subs(p1U[1], Pt).subs(p2U[1], -Pt).subs(p1U[2], zero).subs(p2U[2], zero).subs(p1U[0], -Pr).subs(
            p2U[0], Pr) \
            .subs(nU[0], one).subs(nU[1], zero).subs(nU[2], zero)

    import PN_Hamiltonian_NS as H_NS
    H_NS.f_H_Newt__H_NS_1PN__H_NS_2PN(m1, m2, p1U, n12U, r)
    H_NS.f_H_NS_3PN(m1, m2, p1U, n12U, r)

    global Htot_xyplane_binary
    Htot_xyplane_binary = make_replacements(+H_NS.H_Newt
                                            + H_NS.H_NS_1PN
                                            + H_NS.H_NS_2PN
                                            + H_NS.H_NS_3PN)

    import PN_Hamiltonian_SO as H_SO
    H_SO.f_H_SO_1p5PN(m1, m2, n12U, n21U, S1U, S2U, p1U, p2U, r)
    H_SO.f_H_SO_2p5PN(m1, m2, n12U, n21U, S1U, S2U, p1U, p2U, r)
    H_SO.f_H_SO_3p5PN(m1, m2, n12U, n21U, S1U, S2U, p1U, p2U, r)
    Htot_xyplane_binary += make_replacements(+H_SO.H_SO_1p5PN
                                             + H_SO.H_SO_2p5PN
                                             + H_SO.H_SO_3p5PN)

    import PN_Hamiltonian_SS as H_SS
    H_SS.f_H_SS_2PN(m1, m2, S1U, S2U, nU, r)
    H_SS.f_H_SS_S1S2_3PN(m1, m2, n12U, S1U, S2U, p1U, p2U, r)
    H_SS.f_H_SS_S1sq_S2sq_3PN(m1, m2, n12U, n21U, S1U, S2U, p1U, p2U, r)
    Htot_xyplane_binary += make_replacements(+H_SS.H_SS_2PN
                                             + H_SS.H_SS_S1S2_3PN
                                             + H_SS.H_SS_S1sq_S2sq_3PN)

    import PN_Hamiltonian_SSS as H_SSS
    H_SSS.f_H_SSS_3PN(m1, m2, n12U, n21U, S1U, S2U, p1U, p2U, r)
    Htot_xyplane_binary += make_replacements(+H_SSS.H_SSS_3PN)

# Function for computing dr/dt
def f_dr_dt(Htot_xyplane_binary, m1,m2, n12U, chi1U,chi2U, S1U,S2U, r):
    # First compute p_t
    import PN_p_t as pt
    pt.f_p_t(m1,m2, chi1U,chi2U, r)

    # Then compute dH_{circ}/dr = partial_H(p_r=0)/partial_r
    #                                  + partial_H(p_r=0)/partial_{p_t} partial_{p_t}/partial_r
    dHcirc_dr = (+sp.diff(Htot_xyplane_binary.subs(Pr,sp.sympify(0)),r)
                 +sp.diff(Htot_xyplane_binary.subs(Pr,sp.sympify(0)),Pt)*sp.diff(pt.p_t,r))

    # Then compute M\Omega
    import PN_MOmega as MOm
    MOm.f_MOmega(m1,m2, chi1U,chi2U, r)

    # Next compute dE_GW_dt_plus_dM_dt
    import PN_dE_GW_dt_and_dM_dt as dEdt
    dEdt.f_dE_GW_dt_and_dM_dt(MOm.MOmega, m1,m2, n12U, S1U,S2U)

    # Finally, compute dr/dt
    global dr_dt
    dr_dt = dEdt.dE_GW_dt_plus_dM_dt / dHcirc_dr

# Next we compute p_r as a function of dr_dt (unknown) and known quantities using
# p_r \approx [dr/dt - (partial_H/partial_{p_r})|_{p_r=0}] * [(partial^2_{H}/partial_{p_r^2})|_{p_r=0}]^{-1}
def f_p_r_fullHam(m1,m2, n12U,n21U, chi1U,chi2U, S1U,S2U, p1U,p2U, r):
    f_Htot_xyplane_binary(m1, m2, n12U, n21U, S1U, S2U, p1U, p2U, r)
    f_dr_dt(Htot_xyplane_binary, m1,m2, n12U, chi1U,chi2U, S1U,S2U, r)

    dHdpr_przero   = sp.diff(Htot_xyplane_binary,Pr).subs(Pr,sp.sympify(0))
    d2Hdpr2_przero = sp.diff(sp.diff(Htot_xyplane_binary,Pr),Pr).subs(Pr,sp.sympify(0))
    global p_r_fullHam
    p_r_fullHam = (dr_dt - dHdpr_przero)/(d2Hdpr2_przero)

# Here's the Ramos-Buades, Husa, and Pratten (2018)
#  approach for computing p_r.
# Transcribed from Eq 2.18 of
# Ramos-Buades, Husa, and Pratten (2018),
#   https://arxiv.org/abs/1810.00036
def f_p_r(m1,m2, n12U,n21U, chi1U,chi2U, S1U,S2U, p1U,p2U, r):
    q = m2/m1 # It is assumed that q >= 1, so m2 >= m1.
    f_Htot_xyplane_binary(m1, m2, n12U, n21U, S1U, S2U, p1U, p2U, r)
    f_dr_dt(Htot_xyplane_binary, m1,m2, n12U, chi1U,chi2U, S1U,S2U, r)
    chi1x = chi1U[0]
    chi1y = chi1U[1]
    chi1z = chi1U[2]
    chi2x = chi2U[0]
    chi2y = chi2U[1]
    chi2z = chi2U[2]
    p_r_num = (-dr_dt
               +(-(6*q+13)*q**2*chi1x*chi2y/(4*(q+1)**4)
                 -(6*q+ 1)*q**2*chi2x*chi2y/(4*(q+1)**4)
                 +chi1y*(-q*(   q+6)*chi1x/(4*(q+1)**4)
                         -q*(13*q+6)*chi2x/(4*(q+1)**4)))/r**div(7,2)
               +(+chi1z*(+3*q   *(5*q+2)*chi1x*chi2y/(2*(q+1)**4)
                         -3*q**2*(2*q+5)*chi2x*chi2y/(2*(q+1)**4))
                 +chi1y*chi2z*(+3*q**2*(2*q+5)*chi2x/(2*(q+1)**4)
                               -3*q   *(5*q+2)*chi1x/(2*(q+1)**4)))/r**4)
    p_r_den = (-(q+1)**2/q - (-7*q**2-15*q-7)/(2*q*r)
               -(47*q**4 + 229*q**3 + 363*q**2 + 229*q + 47)/(8*q*(q+1)**2*r**2)
               -(+( 4*q**2 + 11*q + 12)*chi1z/(4*q*(q+1))
                 +(12*q**2 + 11*q +  4)*chi2z/(4*  (q+1)))/r**div(5,2)
               -(+(- 53*q**5 - 357*q**4 - 1097*q**3 - 1486*q**2 - 842*q - 144)*chi1z/(16*q*(q+1)**4)
                 +(-144*q**5 - 842*q**4 - 1486*q**3 - 1097*q**2 - 357*q -  53)*chi2z/(16  *(q+1)**4))/r**div(7,2)
               -(+(  q**2 + 9*q + 9)*chi1x**2/(2*q*(q+1)**2)
                 +(3*q**2 + 5*q + 3)*chi2x*chi1x/((q+1)**2)
                 +(3*q**2 + 8*q + 3)*chi1y*chi2y/(2*(q+1)**2)
                 -9*q**2*chi2y**2/(4*(q+1))
                 +(3*q**2 + 8*q + 3)*chi1z*chi2z/(2*(q+1)**2)
                 -9*q**2*chi2z**2/(4*(q+1))
                 +(9*q**3 + 9*q**2 + q)*chi2x**2/(2*(q+1)**2)
                 +(-363*q**6 - 2608*q**5 - 7324*q**4 - 10161*q**3 - 7324*q**2 - 2608*q - 363)/(48*q*(q+1)**4)
                 -9*chi1y**2/(4*q*(q+1))
                 -9*chi1z**2/(4*q*(q+1)) - sp.pi**2/16)/r**3)
    global p_r
    p_r = p_r_num/p_r_den
