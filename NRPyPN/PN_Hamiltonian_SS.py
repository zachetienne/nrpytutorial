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
import indexedexpNRPyPN as ixp                   # NRPy+: Symbolic indexed expression (e.g., tensors, vectors, etc.) support
from NRPyPN_shortcuts import div,dot,cross # NRPyPN: shortcuts for e.g., vector operations

#################################
#################################
# 2PN spin-spin term, from Eqs. 2.18 and 2.19 of
#      Buonanno, Chen, and Damour (2006):
#     https://arxiv.org/abs/gr-qc/0508067
def f_H_SS_2PN(m1,m2, S1U,S2U, nU, q):
    S0U = ixp.zerorank1()
    for i in range(3):
        S0U[i] = (1 + m2/m1)*S1U[i] + (1 + m1/m2)*S2U[i]
    global H_SS_2PN
    mu = m1*m2 / (m1 + m2)
    H_SS_2PN = mu/(m1 + m2) * (3*dot(S0U,nU)**2 - dot(S0U,S0U)) / (2*q**3)

#################################
#################################
# 3PN spin-spin S_1,S_2 coupling term, from Eq. 2.11 of
#       Steinhoff, Hergt, and Sch\"afer (2008a)
#         https://arxiv.org/abs/0712.1716
def f_H_SS_S1S2_3PN(m1,m2, n12U, S1U,S2U, p1U,p2U, r12):
    global H_SS_S1S2_3PN
    H_SS_S1S2_3PN = (+div(3,2)*(dot(cross(p1U,S1U),n12U)*dot(cross(p2U,S2U),n12U))
                     +       6*(dot(cross(p2U,S1U),n12U)*dot(cross(p1U,S2U),n12U))
                     -15*dot(S1U,n12U)*dot(S2U,n12U)*dot(p1U,n12U)*dot(p2U,n12U)
                     -3*dot(S1U,n12U)*dot(S2U,n12U)*dot(p1U,p2U)
                     +3*dot(S1U,p2U)*dot(S2U,n12U)*dot(p1U,n12U)
                     +3*dot(S2U,p1U)*dot(S1U,n12U)*dot(p2U,n12U)
                     +3*dot(S1U,p1U)*dot(S2U,n12U)*dot(p2U,n12U)
                     +3*dot(S2U,p2U)*dot(S1U,n12U)*dot(p1U,n12U)
                     -div(1,2)*dot(S1U,p2U)*dot(S2U,p1U)
                     +dot(S1U,p1U)*dot(S2U,p2U)
                     -3*dot(S1U,S2U)*dot(p1U,n12U)*dot(p2U,n12U)
                     +div(1,2)*dot(S1U,S2U)*dot(p1U,p2U))/(2*m1*m2*r12**3)
    H_SS_S1S2_3PN+= (-dot(cross(p1U,S1U),n12U)*dot(cross(p1U,S2U),n12U)
                     +dot(S1U,S2U)*dot(p1U,n12U)**2
                     -dot(S1U,n12U)*dot(S2U,p1U)*dot(p1U,n12U))*3/(2*m1**2*r12**3)
    H_SS_S1S2_3PN+= (-dot(cross(p2U,S2U),n12U)*dot(cross(p2U,S1U),n12U)
                     +dot(S1U,S2U)*dot(p2U,n12U)**2
                     -dot(S2U,n12U)*dot(S1U,p1U)*dot(p2U,n12U))*3/(2*m2**2*r12**3)
    H_SS_S1S2_3PN+= (+dot(S1U,S2U)-2*dot(S1U,n12U)*dot(S2U,n12U))*6*(m1+m2)/r12**4

#################################
#################################
# 3PN spin-orbit coupling term, from Eq. 9 of
#    Steinhoff, Hergt, and Sch\"afer (2008b)
#       https://arxiv.org/abs/0809.2200
def f_H_SS_S1sq_S2sq_3PN(m1,m2, n12U,n21U, S1U,S2U, p1U,p2U, r12):
    def f_H_SS_particle(m1,m2, n12U, S1U,_S2U, p1U,p2U, r12): # _S2U unused.
        H_SS_S1sq_S2sq_3PN_particle = (
            +  m2/(4*m1**3)*dot(p1U,S1U)**2
            +3*m2/(8*m1**3)*dot(p1U,n12U)**2*dot(S1U,S1U)
            -3*m2/(8*m1**3)*dot(p1U,p1U)*dot(S1U,n12U)**2
            -3*m2/(4*m1**3)*dot(p1U,n12U)*dot(S1U,n12U)*dot(p1U,S1U)
            -3/(4*m1*m2)*dot(p2U,p2U)*dot(S1U,S1U)
            +9/(4*m1*m2)*dot(p2U,p2U)*dot(S1U,n12U)**2
            +3/(4*m1**2)*dot(p1U,p2U)*dot(S1U,S1U)
            -9/(4*m1**2)*dot(p1U,p2U)*dot(S1U,n12U)**2
            -3/(2*m1**2)*dot(p1U,n12U)*dot(p2U,S1U)*dot(S1U,n12U)
            +3/(m1**2)  *dot(p2U,n12U)*dot(p1U,S1U)*dot(S1U,n12U)
            +3/(4*m1**2)*dot(p1U,n12U)*dot(p2U,n12U)*dot(S1U,S1U)
            -15/(4*m1**2)*dot(p1U,n12U)*dot(p2U,n12U)*dot(S1U,n12U)**2)/r12**3
        H_SS_S1sq_S2sq_3PN_particle+= -(+div(9,2)*dot(S1U,n12U)**2
                                         -div(5,2)*dot(S1U,S1U)
                                         +7*m2/m1*dot(S1U,n12U)**2
                                         -3*m2/m1*dot(S1U,S1U))*m2/r12**4
        return H_SS_S1sq_S2sq_3PN_particle
    global H_SS_S1sq_S2sq_3PN
    H_SS_S1sq_S2sq_3PN = (+f_H_SS_particle(m1,m2, n12U, S1U,S2U, p1U,p2U, r12)
                          +f_H_SS_particle(m2,m1, n21U, S2U,S1U, p2U,p1U, r12))
