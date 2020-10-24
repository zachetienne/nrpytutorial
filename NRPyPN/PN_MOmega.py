# As documented in the NRPyPN notebook
# PN-MOmega.ipynb, this Python script
# generates the expression for the orbital
# angular frequency M*Omega up to and
# including terms at 3.5PN order.
# It largely follows the expression of
#  Ramos-Buades, Husa, and Pratten (2018)
#    https://arxiv.org/abs/1810.00036
# but fixes obvious typos in the given
# equation.

# Core functions:
# f_MOmega(m1,m2, chi1U,chi2U, r)
#       Compute MOmega and store to
#       global variable of the same name.

# Author:  Zach Etienne
#          zachetie **at** gmail **dot* com

# Step 0: Add NRPy's directory to the path
# https://stackoverflow.com/questions/16780014/import-file-from-parent-directory
import sympy as sp               # SymPy: The Python computer algebra package upon which NRPy+ depends
import indexedexpNRPyPN as ixp         # NRPy+: Symbolic indexed expression (e.g., tensors, vectors, etc.) support
from NRPyPN_shortcuts import div # NRPyPN: shortcuts for e.g., vector operations

#################################
#################################
# Basic equation for MOmega can be written in the
#   form:
# MOmega = 1/r^(3/2) * (1 + \sum_{k=2}^7 (a_k/r^{k/2}))
# where we construct the a_k terms in the sum below:

# Step 1: Construct terms a_2, a_3, and a_4, from
#  Eq A2 of Ramos-Buades, Husa, and Pratten (2018)
#    https://arxiv.org/abs/1810.00036
#  These terms have been independently validated
#    against the same terms in Eq 6 of
#  Healy, Lousto, Nakano, and Zlochower (2017)
#    https://arxiv.org/abs/1702.00872
def MOmega__a_2_thru_a_4(m1,m2, chi1x,chi1y,chi1z, chi2x,chi2y,chi2z):
    q = m2/m1 # It is assumed that q >= 1, so m2 >= m1.
    global a_2,a_3,a_4
    a_2 = -((3*q**2+5*q+3)/(2*(q+1)**2))
    a_3 = (-(3*q+4)*chi1z/(4*(q+1)**2) - q*(4*q+3)*chi2z/(4*(q+1)**2))
    a_4 = (-3*q**2*chi2x**2/(2*(q+1)**2)
           +3*q**2*chi2y**2/(4*(q+1)**2)
           +3*q**2*chi2z**2/(4*(q+1)**2)
           +(+24*q**4 + 103*q**3 + 164*q**2 + 103*q + 24)/(16*(q+1)**4)
           -3*chi1x**2/(2*(q+1)**2)
           -3*q*chi1x*chi2x/(q+1)**2
           +3*chi1y**2/(4*(q+1)**2)
           +3*q*chi1y*chi2y/(2*(q+1)**2)
           +3*chi1z**2/(4*(q+1)**2)
           +3*q*chi1z*chi2z/(2*(q+1)**2))

# Construct terms a_5 and a_6, from
#  Eq A1 of Ramos-Buades, Husa, and Pratten (2018)
#    https://arxiv.org/abs/1810.00036
#  These terms have been independently validated
#    against the same terms in Eq 6 of
#  Healy, Lousto, Nakano, and Zlochower (2017)
#    https://arxiv.org/abs/1702.00872
def MOmega__a_5_thru_a_6(m1,m2, chi1x,chi1y,chi1z, chi2x,chi2y,chi2z):
    q = m2/m1 # It is assumed that q >= 1, so m2 >= m1.
    global a_5,a_6
    a_5 = (+3*  (13*q**3 + 34*q**2 + 30*q + 16)*chi1z/(16*(q+1)**4)
           +3*q*(16*q**3 + 30*q**2 + 34*q + 13)*chi2z/(16*(q+1)**4))
    a_6 = (+(+155*q**2 + 180*q + 76)*chi1x**2/(16*(q+1)**4)
           +q*(+120*q**2 + 187*q + 120)*chi1x*chi2x/(8*(q+1)**4)
           -(+55*q**2 + 85*q + 43)*chi1y**2/(8*(q+1)**4)
           -q*(+54*q**2 + 95*q + 54)*chi1y*chi2y/( 4*(q+1)**4)
           -q*(+96*q**2 +127*q + 96)*chi1z*chi2z/(16*(q+1)**4)
           +q**2*(+76*q**2 + 180*q + 155)*chi2x**2/(16*(q+1)**4)
           -q**2*(+43*q**2 +  85*q +  55)*chi2y**2/( 8*(q+1)**4)
           -q**2*(+2*q+5)*(+14*q+27)*chi2z**2/(32*(q+1)**4)
           -     (+5*q+2)*(+27*q+14)*chi1z**2/(32*(q+1)**4)
           +(+501*sp.pi**2*q*(q+1)**4
             -4*(120*q**6 + 2744*q**5 + 10049*q**4 + 14820*q**3 + 10049*q**2 + 2744*q + 120))/(384*(q+1)**6))

# Construct term a_7, from Eq A1 of
#  Ramos-Buades, Husa, and Pratten (2018)
#    https://arxiv.org/abs/1810.00036
def MOmega__a_7(m1,m2, chi1x,chi1y,chi1z, chi2x,chi2y,chi2z):
    q = m2/m1 # It is assumed that q >= 1, so m2 >= m1.
    global a_7
    a_7 = (+3*(4*q+1)*q**3*chi2x**2*chi2z/(2*(q+1)**4)
           -3*(4*q+1)*q**3*chi2y**2*chi2z/(8*(q+1)**4)
           -3*(4*q+1)*q**3*chi2z**3      /(8*(q+1)**4)
           +chi1x*(+9*(2*q+1)*q**2*chi2x*chi2z/(4*(q+1)**4)
                   +9*(1*q+2)*q   *chi2x*chi1z/(4*(q+1)**4))
           +chi1y*(+9*q**2*chi2y*chi1z/(4*(q+1)**4)
                   +9*q**2*chi2y*chi2z/(4*(q+1)**4))
           +chi1z*(+9*q**2*(2*q+3)*chi2x**2/(4*(q+1)**4)
                   -9*q**2*(  q+2)*chi2y**2/(4*(q+1)**4)
                   -9*q**2        *chi2z**2/(4*(q+1)**3)
                   -(135*q**5 + 385*q**4 + 363*q**3 + 377*q**2 + 387*q + 168)/(32*(q+1)**6))
           -(+168*q**5 + 387*q**4 + 377*q**3 + 363*q**2 + 385*q + 135)*q*chi2z/(32*(q+1)**6)
           +chi1x**2*(+3*(q+4)*chi1z/(2*(q+1)**4)
                      +9*q*(3*q+2)*chi2z/(4*(q+1)**4))
           +chi1y**2*(-3*(q+4)*chi1z/(8*(q+1)**4)
                      -9*q*(2*q+1)*chi2z/(4*(q+1)**4))
           -9*q*chi1z**2*chi2z/(4*(q+1)**3)
           -3*(q+4)*chi1z**3/(8*(q+1)**4))

# Finally, sum the expressions for a_k to construct p_t as prescribed:
# MOmega = 1/r^(3/2) * (1 + \sum_{k=2}^7 (a_k/r^{k/2}))
def f_MOmega(m1,m2, chi1U,chi2U, r):
    a = ixp.zerorank1(DIM=10)
    MOmega__a_2_thru_a_4(m1,m2, chi1U[0],chi1U[1],chi1U[2], chi2U[0],chi2U[1],chi2U[2])
    a[2] = a_2
    a[3] = a_3
    a[4] = a_4
    MOmega__a_5_thru_a_6(m1,m2, chi1U[0],chi1U[1],chi1U[2], chi2U[0],chi2U[1],chi2U[2])
    a[5] = a_5
    a[6] = a_6
    MOmega__a_7(         m1,m2, chi1U[0],chi1U[1],chi1U[2], chi2U[0],chi2U[1],chi2U[2])
    a[7] = a_7
    global MOmega
    MOmega = 1 # Term prior to the sum in parentheses
    for k in range(8):
        MOmega += a[k]/r**div(k,2)
    MOmega *= 1/r**div(3,2)
