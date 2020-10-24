# As documented in the NRPyPN notebook
# PN-dE_GW_dt.ipynb, this Python script
# generates dE_GW/dt at highest known
# post-Newtonian order (as of 2015, at
# least).

# Core functions:
# dE_GW_dt_OBKPSS2015_consts(m1,m2, n12U, S1U,S2U):
#       Define constants used in the dE_GW/dt expression.
# f_dE_GW_dt(mOmega, m1,m2, n12U, S1U,S2U):
#       Compute dE_GW_dt and store to global variable of the same name.

# Author:  Zach Etienne
#          zachetie **at** gmail **dot* com

# Step 0: Add NRPy's directory to the path
# https://stackoverflow.com/questions/16780014/import-file-from-parent-directory
import sys                    # Standard Python modules for multiplatform OS-level functions
import sympy as sp               # SymPy: The Python computer algebra package upon which NRPy+ depends
import indexedexpNRPyPN as ixp         # NRPy+: Symbolic indexed expression (e.g., tensors, vectors, etc.) support
from NRPyPN_shortcuts import div,dot,gamma_EulerMascheroni  # NRPyPN: shortcuts for e.g., vector operations

#################################
#################################

# Constants given in Eqs A1-13 of https://arxiv.org/abs/1502.01747
def dE_GW_dt_OBKPSS2015_consts(m1,m2, _n12U, S1U,S2U): # _n12U unused.
    # define scalars:
    m  = (m1+m2)
    nu = m1*m2/m**2
    delta = (m1-m2)/m
    # define vectors:
    Stot = ixp.zerorank1()
    Sigma= ixp.zerorank1()
    l    = ixp.zerorank1()
    l[2] = sp.sympify(1)
    chi1U = ixp.zerorank1()
    chi2U = ixp.zerorank1()
    chi_s = ixp.zerorank1()
    chi_a = ixp.zerorank1()
    for i in range(3):
        Stot[i] = S1U[i] + S2U[i]
        Sigma[i] = (m1+m2)/m2*S2U[i] - (m1+m2)/m1*S1U[i]
        chi1U[i] = S1U[i]/m1**2
        chi2U[i] = S2U[i]/m2**2
        chi_s[i] = div(1,2) * (chi1U[i] + chi2U[i])
        chi_a[i] = div(1,2) * (chi1U[i] - chi2U[i])
    # define scalars that depend on vectors
    s_l = dot(Stot,l)   /m**2
    # s_n = dot(Stot,n12U)/m**2
    sigma_l = dot(Sigma,l)/m**2
    # sigma_n = dot(Sigma,n12U)/m**2
    return nu,delta,  l,chi_a,chi_s,  s_l,sigma_l

#################################
#################################

# Based on Eqs A22-28 of https://arxiv.org/abs/1502.01747, with
#       Eq A.14 of https://arxiv.org/abs/0709.0093 for Mdot
#       and correction on b[7] term by comparison with
#  https://link.springer.com/content/pdf/10.12942/lrr-2014-2.pdf
def f_dE_GW_dt_and_dM_dt(mOmega, m1,m2, n12U, S1U,S2U):
    def f_compute_quantities(mOmega, m1,m2, n12U, S1U,S2U, which_quantity):
        if not which_quantity in ('dM_dt', 'dE_GW_dt', 'dE_GW_dt_plus_dM_dt'):
            print("which_quantity == "+str(which_quantity)+" not supported!")
            sys.exit(1)

        nu,delta,  l,chi_a,chi_s,  s_l,sigma_l = dE_GW_dt_OBKPSS2015_consts(m1,m2, n12U, S1U,S2U)
        x = (mOmega)**div(2,3)

        # Compute b_5_Mdot:
        b_5_Mdot = (-div(1,4)*(+(1-3*nu)*dot(chi_s,l)*(1+3*dot(chi_s,l)**2+9*dot(chi_a,l)**2)
                            +(1-  nu)*delta*dot(chi_a,l)*(1+3*dot(chi_a,l)**2+9*dot(chi_s,l)**2)))
        if which_quantity == "dM_dt":
            return div(32,5)*nu**2*x**5*b_5_Mdot*x**div(5,2)

        b = ixp.zerorank1(DIM=10)
        b[2] = -div(1247,336) - div(35,12)*nu
        b[3] = +4*sp.pi - 4*s_l - div(5,4)*delta*sigma_l
        b[4] =(-div(44711,9072) + div(9271,504)*nu + div(65,18)*nu**2
               +(+div(287,96) + div( 1,24)*nu)*dot(chi_s,l)**2
               -(+div( 89,96) + div( 7,24)*nu)*dot(chi_s,chi_s)
               +(+div(287,96) -         12*nu)*dot(chi_a,l)**2
               +(-div( 89,96) +          4*nu)*dot(chi_a,chi_a)
               +div(287,48)*delta*dot(chi_s,l)*dot(chi_a,l) - div(89,48)*delta*dot(chi_s,chi_a))
        b[5] =(-div(8191,672)*sp.pi - div(9,2)*s_l - div(13,16)*delta*sigma_l
               +nu*(-div(583,24)*sp.pi + div(272,9)*s_l + div(43,4)*delta*sigma_l))
        if which_quantity == "dE_GW_dt_plus_dM_dt":
            b[5]+= b_5_Mdot
        b[6] =(+div(6643739519,69854400) + div(16,3)*sp.pi**2 - div(1712,105)*gamma_EulerMascheroni
               -div(856,105)*sp.log(16*x) + (-div(134543,7776) + div(41,48)*sp.pi**2)*nu
               -div(94403,3024)*nu**2 - div(775,324)*nu**3 - 16*sp.pi*s_l - div(31,6)*sp.pi*delta*sigma_l)
        b[7] =(+(+div(476645,6804) + div(6172,189)*nu - div(2810,27)*nu**2)*s_l
               +(+div(9535,336) + div(1849,126)*nu - div(1501,36)*nu**2)*delta*sigma_l
               +(-div(16285,504) + div(214745,1728)*nu + div(193385,3024)*nu**2)*sp.pi)
        b[8] =(+(-div(3485,96)*sp.pi + div(13879,72)*sp.pi*nu)*s_l
               +(-div(7163,672)*sp.pi + div(130583,2016)*sp.pi*nu)*delta*sigma_l)
        b_sum = sp.sympify(1)
        for k in range(9):
            b_sum += b[k]*x**div(k,2)
        return div(32,5)*nu**2*x**5*b_sum

    global dE_GW_dt_plus_dM_dt, dE_GW_dt, dM_dt
    dE_GW_dt_plus_dM_dt = \
               f_compute_quantities(mOmega, m1,m2, n12U, S1U,S2U, "dE_GW_dt_plus_dM_dt")
    dE_GW_dt = f_compute_quantities(mOmega, m1,m2, n12U, S1U,S2U, "dE_GW_dt")
    dM_dt    = f_compute_quantities(mOmega, m1,m2, n12U, S1U,S2U, "dM_dt")
