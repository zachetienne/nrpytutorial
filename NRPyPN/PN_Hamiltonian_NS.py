# As documented in the NRPyPN notebook
# PN-Hamiltonian-Nonspinning.ipynb, this Python script 
# generates nonspinning pieces of the post-Newtonian (PN)
# Hamiltonian, up to and including third PN order.

# Basic functions:
# f_H_Newt__H_NS_1PN__H_NS_2PN(m1,m2, pU, nU, q): Compute H_Newt, 
#                                                 H_NS_1PN, and H_NS_2PN
#                                                 and store to global 
#                                                 variables of the same
#                                                 names.
# f_H_NS_3PN: Compute H_NS_3PN, and store to global variable of same name

# Author:  Zach Etienne
#          zachetie **at** gmail **dot* com

# Step 0: Add NRPy's directory to the path
# https://stackoverflow.com/questions/16780014/import-file-from-parent-directory
import os,sys                    # Standard Python modules for multiplatform OS-level functions
nrpy_dir_path = os.path.join("..")
if nrpy_dir_path not in sys.path:
    sys.path.append(nrpy_dir_path)
import sympy as sp               # SymPy: The Python computer algebra package upon which NRPy+ depends
from outputC import *            # NRPy+: Core C code output module
import indexedexp as ixp         # NRPy+: Symbolic indexed expression (e.g., tensors, vectors, etc.) support
from NRPyPN_shortcuts import *   # NRPyPN: shortcuts for e.g., vector operations

def f_H_Newt__H_NS_1PN__H_NS_2PN(m1,m2, pU, nU, q):
    mu = m1*m2/(m1+m2)
    PU = ixp.zerorank1()
    for i in range(3):
        PU[i] = pU[i]/mu
    eta = m1*m2/(m1+m2)**2
    P_dot_P = dot(PU,PU)
    n_dot_P = dot(nU,PU)

    # H_{\rm Newt} = \frac{p^i p^i}{2} - \frac{1}{q}
    global H_Newt,H_NS_1PN,H_NS_2PN
    H_Newt    = mu*(div(1,2)*P_dot_P - 1/q)
    H_NS_1PN  = mu*(div(1,8)*(3*eta-1)*P_dot_P**2 - \
                 div(1,2)*((3+eta)*P_dot_P + eta*n_dot_P**2)/q + 1/(2*q**2))
    H_NS_2PN  = mu*(div(1,16)*(1 - 5*eta + 5*eta**2)*P_dot_P**3 +
                 div(1,8)*((5 - 20*eta - 3*eta**2)*P_dot_P**2 
                           - 2*eta**2*n_dot_P**2*P_dot_P - 3*eta**2*n_dot_P**4)/q +
                 div(1,2)*((5 + 8*eta)*P_dot_P + 3*eta*n_dot_P**2)/q**2 -
                 div(1,4)*(1 + 3*eta)/q**3)

def f_H_NS_3PN(m1,m2, pU, nU, q):
    mu = m1*m2/(m1+m2)
    PU = ixp.zerorank1()
    for i in range(3):
        PU[i] = pU[i]/mu
    eta = m1*m2/(m1+m2)**2
    P_dot_P = dot(PU,PU)
    n_dot_P = dot(nU,PU)

    global H_NS_3PN
    # The following is simply by-hand search/replaced from the above LaTeX to minimize error
    H_NS_3PN = \
    mu*( div(1,128)*(-5+35*eta-70*eta**2+35*eta**3)*P_dot_P**4 +
         div(1,16)* ( (-7+42*eta-53*eta**2-5*eta**3)*P_dot_P**3
                     +(2-3*eta)*eta**2*n_dot_P**2*P_dot_P**2 +
                     +3*(1-eta)*eta**2*n_dot_P**4*P_dot_P - 5*eta**3*n_dot_P**6 )/(q) +
 (  div(1,16)*(-27+136*eta+109*eta**2)*P_dot_P**2
 + div(1,16)*(17+30*eta)*eta*n_dot_P**2*P_dot_P + div(1,12)*(5+43*eta)*eta*n_dot_P**4)/(q**2) +
 ( ( -div(25,8) + (div(1,64)*sp.pi**2-div(335,48))*eta 
 - div(23,8)*eta**2 )*P_dot_P
 + (-div(85,16)-div(3,64)*sp.pi**2-div(7,4)*eta)*eta*n_dot_P**2)/(q**3) +
 ( div(1,8) + (div(109,12)-div(21,32)*sp.pi**2)*eta)/(q**4) )