# As documented in the NRPyPN notebook
# PN-p_t.ipynb, this Python script
# generates the expression for the transverse
# component of momentum p_t up to and
# including terms at 3.5PN order.
# This is an implementation of the equations of
#  Ramos-Buades, Husa, and Pratten (2018)
#    https://arxiv.org/abs/1810.00036
# but validates against the relevant equation
#   in Healy, Lousto, Nakano, and Zlochower (2017)
#    https://arxiv.org/abs/1702.00872

# Core functions:
# f_p_t(m1,m2, chi1U,chi2U, r)
#       Compute p_t and store to
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
# Basic equation for p_t can be written in the
#   form:
# p_t = q/(sqrt(r)*(1+q)^2) (1 + \sum_{k=2}^7 (a_k/r^{k/2})
# where we construct the a_k terms in the sum below:

# Step 1: Construct terms a_2, a_3, and a_4, from
#  Eq A2 of Ramos-Buades, Husa, and Pratten (2018)
#    https://arxiv.org/abs/1810.00036
#  These terms have been independently validated
#    against the same terms in Eq 7 of
#  Healy, Lousto, Nakano, and Zlochower (2017)
#    https://arxiv.org/abs/1702.00872
def p_t__a_2_thru_a_4(m1,m2, chi1x,chi1y,chi1z, chi2x,chi2y,chi2z):
    q = m2/m1 # It is assumed that q >= 1, so m2 >= m1.
    global a_2,a_3,a_4
    a_2 = 2
    a_3 = (-3*(4*q**2+3*q)*chi2z/(4*(q+1)**2) - 3*(3*q+4)*chi1z/(4*(q+1)**2))
    a_4 = (-3*q**2*chi2x**2/(2*(q+1)**2)
           +3*q**2*chi2y**2/(4*(q+1)**2)
           +3*q**2*chi2z**2/(4*(q+1)**2)
           +(+42*q**2 + 41*q + 42)/(16*(q+1)**2)
           -3*chi1x**2/(2*(q+1)**2)
           -3*q*chi1x*chi2x/(q+1)**2
           +3*chi1y**2/(4*(q+1)**2)
           +3*q*chi1y*chi2y/(2*(q+1)**2)
           +3*chi1z**2/(4*(q+1)**2)
           +3*q*chi1z*chi2z/(2*(q+1)**2))

# Construct terms a_5 and a_6, from
#  Eq A2 of Ramos-Buades, Husa, and Pratten (2018)
#    https://arxiv.org/abs/1810.00036
#  These terms have been independently validated
#    against the same terms in Eq 7 of
#  Healy, Lousto, Nakano, and Zlochower (2017)
#    https://arxiv.org/abs/1702.00872
#  and a sign error was corrected in the a_5
#  expression.
def p_t__a_5_thru_a_6(m1,m2, chi1x,chi1y,chi1z, chi2x,chi2y,chi2z, FixSignError=True):
    SignFix = sp.sympify(-1)
    if FixSignError == False:
        SignFix = sp.sympify(+1)
    q = m2/m1 # It is assumed that q >= 1, so m2 >= m1.
    global a_5,a_6
    a_5 = (SignFix*(13*q**3 + 60*q**2 + 116*q + 72)*chi1z/(16*(q+1)**4)
           +(-72*q**4 - 116*q**3 - 60*q**2 - 13*q)*chi2z/(16*(q+1)**4))
    a_6 = (+(+472*q**2 - 640)*chi1x**2/(128*(q+1)**4)
           +(-512*q**2 - 640*q - 64)*chi1y**2/(128*(q+1)**4)
           +(-108*q**2 + 224*q +512)*chi1z**2/(128*(q+1)**4)
           +(+472*q**2 - 640*q**4)*chi2x**2/(128*(q+1)**4)
           +(+192*q**3 + 560*q**2 + 192*q)*chi1x*chi2x/(128*(q+1)**4)
           +(-864*q**3 -1856*q**2 - 864*q)*chi1y*chi2y/(128*(q+1)**4)
           +(+480*q**3 +1064*q**2 + 480*q)*chi1z*chi2z/(128*(q+1)**4)
           +( -64*q**4 - 640*q**3 - 512*q**2)*chi2y**2/(128*(q+1)**4)
           +(+512*q**4 + 224*q**3 - 108*q**2)*chi2z**2/(128*(q+1)**4)
           +(+480*q**4 + 163*sp.pi**2*q**3 - 2636*q**3 + 326*sp.pi**2*q**2 - 6128*q**2 + 163*sp.pi**2*q-2636*q+480)
            /(128*(q+1)**4))

# Construct term a_7, from Eq A2 of
#  Ramos-Buades, Husa, and Pratten (2018)
#    https://arxiv.org/abs/1810.00036
def p_t__a_7(m1,m2, chi1x,chi1y,chi1z, chi2x,chi2y,chi2z):
    q = m2/m1 # It is assumed that q >= 1, so m2 >= m1.
    global a_7
    a_7 = (+5*(4*q+1)*q**3*chi2x**2*chi2z/(2*(q+1)**4)
           -5*(4*q+1)*q**3*chi2y**2*chi2z/(8*(q+1)**4)
           -5*(4*q+1)*q**3*chi2z**3      /(8*(q+1)**4)
           +chi1x*(+15*(2*q+1)*q**2*chi2x*chi2z/(4*(q+1)**4)
                   +15*(1*q+2)*q   *chi2x*chi1z/(4*(q+1)**4))
           +chi1y*(+15*q**2*chi2y*chi1z/(4*(q+1)**4)
                   +15*q**2*chi2y*chi2z/(4*(q+1)**4))
           +chi1z*(+15*q**2*(2*q+3)*chi2x**2/(4*(q+1)**4)
                   -15*q**2*(  q+2)*chi2y**2/(4*(q+1)**4)
                   -15*q**2        *chi2z**2/(4*(q+1)**3)
                   -(103*q**5 + 145*q**4 - 27*q**3 + 252*q**2 + 670*q + 348)/(32*(q+1)**6))
           -(+348*q**5 + 670*q**4 + 252*q**3 - 27*q**2 + 145*q + 103)*q*chi2z/(32*(q+1)**6)
           +chi1x**2*(+5*(q+4)*chi1z/(2*(q+1)**4)
                      +15*q*(3*q+2)*chi2z/(4*(q+1)**4))
           +chi1y**2*(-5*(q+4)*chi1z/(8*(q+1)**4)
                      -15*q*(2*q+1)*chi2z/(4*(q+1)**4))
           -15*q*chi1z**2*chi2z/(4*(q+1)**3)
           -5*(q+4)*chi1z**3/(8*(q+1)**4))

# Finally, sum the expressions for a_k to construct p_t as prescribed:
# p_t = q/(sqrt(r)*(1+q)^2) (1 + \sum_{k=2}^7 (a_k/r^{k/2}))
def f_p_t(m1,m2, chi1U,chi2U, r):
    q = m2/m1 # It is assumed that q >= 1, so m2 >= m1.
    a = ixp.zerorank1(DIM=10)
    p_t__a_2_thru_a_4(m1,m2, chi1U[0],chi1U[1],chi1U[2], chi2U[0],chi2U[1],chi2U[2])
    a[2] = a_2
    a[3] = a_3
    a[4] = a_4
    p_t__a_5_thru_a_6(m1,m2, chi1U[0],chi1U[1],chi1U[2], chi2U[0],chi2U[1],chi2U[2])
    a[5] = a_5
    a[6] = a_6
    p_t__a_7(         m1,m2, chi1U[0],chi1U[1],chi1U[2], chi2U[0],chi2U[1],chi2U[2])
    a[7] = a_7
    global p_t
    p_t = 1 # Term prior to the sum in parentheses
    for k in range(8):
        p_t += a[k]/r**div(k,2)
    p_t *= q / (1+q)**2 * 1/r**div(1,2)
