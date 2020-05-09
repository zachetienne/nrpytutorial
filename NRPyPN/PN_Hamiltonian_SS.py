# As documented in the NRPyPN notebook
# PN-Hamiltonian-Spin-Spin.ipynb, this Python script
# generates spin-spin coupling pieces of the
# post-Newtonian (PN) Hamiltonian, up to and
# including 3PN order.

# Core functions:
# f_H_SS_2PN(m1,m2, S1U,S2U, nU, q):
#       Compute the complete H_SS_2PN term and store to
#                     global variable of the same name.
# f_HS1S2_3PN(m1,m2, n12U, S1U,S2U, p1U,p2U, q)):
#       Compute HS1S2_3PN and store to global variable
#                     of the same name.
# f_H_SS_S1sq_S2sq_3PN(m1,m2, n12U,n21U, S1U,S2U, p1U,p2U, q):
#       Compute H_SS_S1sq_S2sq_3PN and store to global
#                     variable of the same name.

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

#################################
#################################

# 2PN spin-spin term, from Eqs. 2.18 and 2.19 of
#      Buonanno, Chen, and Damour (2006):
#     https://arxiv.org/abs/gr-qc/0508067
def f_H_SS_2PN(m1,m2, S1U,S2U, nU, q):
    S_0U = ixp.zerorank1()
    for i in range(3):
        S_0U[i] = (1 + m2/m1)*S1U[i] + (1 + m1/m2)*S2U[i]
    mu = m1*m2 / (m1+m2)
    global H_SS_2PN
    H_SS_2PN = div(1,2)*mu/(m1+m2)*( 3*dot(S_0U,nU)**2 - dot(S_0U,S_0U) )/q**3

#################################
#################################

# 3PN spin-spin S_1,S_2 coupling term, from Eq. 2.11 of
#       Steinhoff, Hergt, and Schäfer (2008a)
#         https://arxiv.org/abs/0712.1716
def f_H_SS_S1S2_3PN(m1,m2, n12U, S1U,S2U, p1U,p2U, q):
    def SHS2008a_HS1S2_3PN_pt1(m1,m2, n12U, S1U,S2U, p1U,p2U, q):
        Hpt1 = ( div(3,2)*(dot(cross(p1U,S1U),n12U)*dot(cross(p2U,S2U),n12U))  # line 1
                 +6 *dot(cross(p2U,S1U),n12U)*dot(cross(p1U,S2U),n12U)         # line 1
                 -15*dot(S1U,n12U)*dot(S2U,n12U)*dot(p1U,n12U)*dot(p2U,n12U)   # line 2
                 -3*dot(S1U,n12U)*dot(S2U,n12U)*dot(p1U,p2U)                   # line 2
                 +3*dot(S1U,p2U)*dot(S2U,n12U)*dot(p1U,n12U) # line 3
                 +3*dot(S2U,p1U)*dot(S1U,n12U)*dot(p2U,n12U) # line 3
                 +3*dot(S1U,p1U)*dot(S2U,n12U)*dot(p1U,n12U) # line 3
                 +3*dot(S2U,p2U)*dot(S1U,n12U)*dot(p1U,n12U)                     # line 4
                 -div(1,2)*dot(S1U,p2U)*dot(S2U,p1U) + dot(S1U,p1U)*dot(S2U,p2U) # line 4
                 -3*dot(S1U,S2U)*dot(p1U,n12U)*dot(p2U,n12U)            # line 5
                 +div(1,2)*dot(S1U,S2U)*dot(p1U,p2U) )/(2*m1*m2*q**3)   # line 5
        return Hpt1
    def SHS2008a_HS1S2_3PN_pt2(m1,m2, n12U, S1U,S2U, p1U,p2U, q):
        Hpt2 = ( -dot(cross(p1U,S1U),n12U)*dot(cross(p1U,S2U),n12U)                # line 6
                 +dot(S1U,S2U)*dot(p1U,n12U)**2                                    # line 6
                 -dot(S1U,n12U)*dot(S2U,p1U)*dot(p1U,n12U) )*div(3,2)/(m1**2*q**3) # line 6
        return Hpt2
    def SHS2008a_HS1S2_3PN_pt3(m1,m2, n12U, S1U,S2U, p1U,p2U, q):
        Hpt3 = ( -dot(cross(p2U,S2U),n12U)*dot(cross(p2U,S2U),n12U)                # line 7
                 +dot(S1U,S2U)*dot(p2U,n12U)**2                                    # line 7
                 -dot(S2U,n12U)*dot(S1U,p1U)*dot(p2U,n12U) )*div(3,2)/(m2**2*q**3) # line 7
        return Hpt3
    def SHS2008a_HS1S2_3PN_pt4(m1,m2, n12U, S1U,S2U, p1U,p2U, q):
        Hpt4 = ( dot(S1U,S2U) - 2*dot(S1U,n12U)*dot(S2U,n12U) ) * 6*(m1+m2)/q**4   # line 8
        return Hpt4
    global H_SS_S1S2_3PN
    H_SS_S1S2_3PN = ( +SHS2008a_HS1S2_3PN_pt1(m1,m2, n12U, S1U,S2U, p1U,p2U, q)
                      +SHS2008a_HS1S2_3PN_pt2(m1,m2, n12U, S1U,S2U, p1U,p2U, q)
                      +SHS2008a_HS1S2_3PN_pt3(m1,m2, n12U, S1U,S2U, p1U,p2U, q)
                      +SHS2008a_HS1S2_3PN_pt4(m1,m2, n12U, S1U,S2U, p1U,p2U, q) )

#################################
#################################
# 3PN spin-orbit coupling term, from Eq. 9 of
#    Steinhoff, Hergt, and Schäfer (2008b)
#       https://arxiv.org/abs/0809.2200
def f_H_SS_S1sq_S2sq_3PN(m1,m2, n12U,n21U, S1U,S2U, p1U,p2U, q):
    def SHS2008b_HSsq_3PN_pt(m1,m2, n12U, S1U, p1U,p2U, q):
        H = ( +div(1,4)*m2/m1**3*dot(p1U,S1U)**2 + div(3,8)*m2/m1**3*dot(p1U,n12U)**2*dot(S1U,S1U) # line 1
              -div(3,8)*m2/m1**3*dot(p1U,p1U)*dot(S1U,n12U)**2                                     # line 1
              -div(3,4)*m2/m1**3*dot(p1U,n12U)*dot(S1U,n12U)*dot(p1U,S1U) # line 2
              -div(3,4)/(m1*m2) *dot(p2U,p2U)*dot(S1U,S1U)                # line 2
              +div(9,4)/(m1*m2) *dot(p2U,p2U)*dot(S1U,n12U)**2 # line 3
              +div(3,4)/m1**2   *dot(p1U,p2U)*dot(S1U,S1U)     # line 3
              -div(3,4)/m1**2   *dot(p1U,p2U)*dot(S1U,n12U)**2 # line 3
              -div(3,2)/m1**2   *dot(p1U,n12U)*dot(p2U,S1U)*dot(S1U,n12U) # line 4
              +       3/m1**2   *dot(p2U,n12U)*dot(p2U,S1U)*dot(S1U,n12U) # line 4
              +div(3,4)/m1**2   *dot(p1U,n12U)*dot(p2U,n12U)*dot(S1U,S1U)            # line 5
              -div(15,4)/m1**2  *dot(p1U,n12U)*dot(p2U,n12U)*dot(S1U,n12U)**2 )/q**3 \
            -( +div(9,2)*dot(S1U,n12U)**2          # line 6
               -div(5,2)*dot(S1U,S1U)              # line 6
               +       7*m2/m1*dot(S1U,n12U)**2    # line 6
               -       3*m2/m1*dot(S1U,S1U) )/q**4 # line 6
        return H

    global H_SS_S1sq_S2sq_3PN
    H_SS_S1sq_S2sq_3PN = ( +SHS2008b_HSsq_3PN_pt(m1,m2, n12U, S1U, p1U,p2U, q)   # S_1^2 term
                           +SHS2008b_HSsq_3PN_pt(m2,m1, n21U, S2U, p2U,p1U, q) ) # S_2^2 term