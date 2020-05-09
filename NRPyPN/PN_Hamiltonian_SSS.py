# As documented in the NRPyPN notebook
# PN-Hamiltonian-SSS.ipynb, this Python script
# generates spin-spin-spin coupling pieces of the
# post-Newtonian (PN) Hamiltonian, up to and
# including 3PN order.

# Core functions:
# f_H_SSS_3PN(m1,m2, n12U,n21U, S1U,S2U, p1U,p2U, q)
#       Compute the complete H_SSS_3PN term and store to
#                     global variable of the same name.

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
from NRPyPN_shortcuts import *   # NRPyNR: shortcuts for e.g., vector operations

#################################
#################################
# Step 1: 3PN spin-spin-spin term, from Eq. 3.12 of
#        Levi and Steinhoff (2015):
#     https://arxiv.org/abs/1410.2601
def f_H_SSS_3PN(m1,m2, n12U,n21U, S1U,S2U, p1U,p2U, q):
    def SHS2015_HSSS_3PN_pt(m1,m2, nU, S1U, p1U,p2U, q):
        p2_minus_m2_over_4m1_p1 = ixp.zerorank1()
        for i in range(3):
            p2_minus_m2_over_4m1_p1[i] = p2U[i] - m2/(4*m1)*p1U[i]
        H = ( div(3,2)*(+  dot(S1U,S1U)*dot(S2U,cross(nU,p1U))   # line 1
                        +  dot(S1U,nU)*dot(S2U,cross(S1U,p1U))   # line 1
                        -5*dot(S1U,nU)**2*dot(S2U,cross(nU,p1U)) # line 1
                        +dot(nU,cross(S1U,S2U))*(dot(S1U,p1U)-5*dot(S1U,nU)*dot(p1U,nU)) # line 2
                        -div(3,2)*m1/m2*(+  dot(S1U,S1U)*dot(S2U,cross(nU,p2U))     # line 3
                                         +2*dot(S1U,nU)*dot(S2U,cross(S1U,p2U))     # line 3
                                         -5*dot(S1U,nU)**2*dot(S2U,cross(nU,p2U)))) # line 3
              -dot(cross(S1U,nU), p2_minus_m2_over_4m1_p1)*(dot(S1U,S1U) - 5*dot(S1U,nU)**2) )/q**4
        return H

    global H_SSS_3PN
    H_SSS_3PN = (+SHS2015_HSSS_3PN_pt(m1,m2, n12U, S1U, p1U,p2U, q)
                 +SHS2015_HSSS_3PN_pt(m2,m1, n21U, S2U, p2U,p1U, q))