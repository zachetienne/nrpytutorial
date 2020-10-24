# As documented in the NRPyPN notebook
# PN-Hamiltonian-Spin-Orbit.ipynb, this Python script
# generates spin-orbit coupling pieces of the
# post-Newtonian (PN) Hamiltonian, up to and
# including 3.5PN order.

# Core functions:
# f_H_SO_2p5PN(m1,m2, n12U,n21U, S1U, S2U, p1U,p2U, r12):
#       Compute H_SO_1p5PN and store to global variable
#                     of the same name.
# f_H_SO_2p5PN(m1,m2, n12U,n21U, S1U, S2U, p1U,p2U, r12):
#       Compute H_SO_2p5PN and store to global variable
#                     of the same name.
# f_H_SO_3p5PN(m1,m2, n12U,n21U, S1U, S2U, p1U,p2U, r12):
#       Compute H_SO_3p5PN and store to global variable
#                     of the same name.

# Author:  Zach Etienne
#          zachetie **at** gmail **dot* com

# Step 0: Add NRPy's directory to the path
# https://stackoverflow.com/questions/16780014/import-file-from-parent-directory
import indexedexpNRPyPN as ixp                    # NRPy+: Symbolic indexed expression (e.g., tensors, vectors, etc.) support
from NRPyPN_shortcuts import div,dot,cross  # NRPyPN: shortcuts for e.g., vector operations

#################################
#################################

# 1.5PN spin-orbit coupling term, from Eq. 4.11a of
#    Damour, Jaranowski, and Sch\"afer (2008)
#      https://arxiv.org/abs/0711.1048
def f_H_SO_1p5PN(m1,m2, n12U,n21U, S1U, S2U, p1U,p2U, r12):
    def f_Omega1(m1,m2, n12U, p1U,p2U, r12):
        Omega1 = ixp.zerorank1()
        for i in range(3):
            Omega1[i] = (div(3,2)*m2/m1 * cross(n12U,p1U)[i] - 2*cross(n12U,p2U)[i])/r12**2
        return Omega1
    global H_SO_1p5PN
    Omega1 = f_Omega1(m1,m2, n12U, p1U,p2U, r12)
    Omega2 = f_Omega1(m2,m1, n21U, p2U,p1U, r12)
    H_SO_1p5PN = dot(Omega1,S1U) + dot(Omega2,S2U)

#################################
#################################

# 2.5PN spin-orbit coupling term, from Eq. 4.11b of
#    Damour, Jaranowski, and Sch\"afer (2008)
#      https://arxiv.org/abs/0711.1048
def f_H_SO_2p5PN(m1, m2, n12U, n21U, S1U, S2U, p1U, p2U, r12):
    def f_Omega_SO_2p5PN(m1, m2, n12U, p1U, p2U, r12):
        Omega1 = ixp.zerorank1()
        for i in range(3):
            Omega1[i] = (+(+(-div(11, 2) * m2 - 5 * m2 ** 2 / m1) * cross(n12U, p1U)[i]
                           + (6 * m1 + div(15, 2) * m2) * cross(n12U, p2U)[i]) / r12 ** 3
                         + (+(-div(5, 8) * m2 * dot(p1U, p1U) / m1 ** 3
                              - div(3, 4) * dot(p1U, p2U) / m1 ** 2
                              + div(3, 4) * dot(p2U, p2U) / (m1 * m2)
                              - div(3, 4) * dot(n12U, p1U) * dot(n12U, p2U) / m1 ** 2
                              - div(3, 2) * dot(n12U, p2U) ** 2 / (m1 * m2)) * cross(n12U, p1U)[i]
                            + (dot(p1U, p2U) / (m1 * m2) + 3 * dot(n12U, p1U) * dot(n12U, p2U) / (m1 * m2)) *
                            cross(n12U, p2U)[i]
                            + (div(3, 4) * dot(n12U, p1U) / m1 ** 2 - 2 * dot(n12U, p2U) / (m1 * m2)) * cross(p1U, p2U)[
                                i]) / r12 ** 2)
        return Omega1

    Omega1_2p5PNU = f_Omega_SO_2p5PN(m1, m2, n12U, p1U, p2U, r12)
    Omega2_2p5PNU = f_Omega_SO_2p5PN(m2, m1, n21U, p2U, p1U, r12)

    global H_SO_2p5PN
    H_SO_2p5PN = dot(Omega1_2p5PNU, S1U) + dot(Omega2_2p5PNU, S2U)

#################################
#################################
# 3.5PN spin-orbit coupling term, from Eq. 5 of
#    Hartung and Steinhoff (2011)
#   https://arxiv.org/abs/1104.3079

# 3.5PN H_SO:  Omega_1, part 1:
def HS2011_Omega_SO_3p5PN_pt1(m1,m2, n12U, p1U,p2U, r12):
    Omega1 = ixp.zerorank1()
    for i in range(3):
        Omega1[i] = ((+7*m2*dot(p1U,p1U)**2/(16*m1**5)
                      +9*dot(n12U,p1U)*dot(n12U,p2U)*dot(p1U,p1U)/(16*m1**4)
                      +3*dot(p1U,p1U)*dot(n12U,p2U)**2/(4*m1**3*m2)
                      +45*dot(n12U,p1U)*dot(n12U,p2U)**3/(16*m1**2*m2**2)
                      +9*dot(p1U,p1U)*dot(p1U,p2U)/(16*m1**4)
                      -3*dot(n12U,p2U)**2*dot(p1U,p2U)/(16*m1**2*m2**2)
                      -3*dot(p1U,p1U)*dot(p2U,p2U)/(16*m1**3*m2)
                      -15*dot(n12U,p1U)*dot(n12U,p2U)*dot(p2U,p2U)/(16*m1**2*m2**2)
                      +3*dot(n12U,p2U)**2*dot(p2U,p2U)/(4*m1*m2**3)
                      -3*dot(p1U,p2U)*dot(p2U,p2U)/(16*m1**2*m2**2)
                      -3*dot(p2U,p2U)**2/(16*m1*m2**3))*cross(n12U,p1U)[i])/r12**2
    return Omega1

# 3.5PN H_SO:  Omega_1, part 2:
def HS2011_Omega_SO_3p5PN_pt2(m1,m2, n12U, p1U,p2U, r12):
    Omega1 = ixp.zerorank1()
    for i in range(3):
        Omega1[i] = (+(-3*dot(n12U,p1U)*dot(n12U,p2U)*dot(p1U,p1U)/(2*m1**3*m2)
                       -15*dot(n12U,p1U)**2*dot(n12U,p2U)**2/(4*m1**2*m2**2)
                       +3*dot(p1U,p1U)*dot(n12U,p2U)**2/(4*m1**2*m2**2)
                       -dot(p1U,p1U)*dot(p1U,p2U)/(2*m1**3*m2)
                       +dot(p1U,p2U)**2/(2*m1**2*m2**2)
                       +3*dot(n12U,p1U)**2*dot(p2U,p2U)/(4*m1**2*m2**2)
                       -dot(p1U,p1U)*dot(p2U,p2U)/(4*m1**2*m2**2)
                       -3*dot(n12U,p1U)*dot(n12U,p2U)*dot(p2U,p2U)/(2*m1*m2**3)
                       -dot(p1U,p2U)*dot(p2U,p2U)/(2*m1*m2**3))*cross(n12U,p2U)[i])/r12**2
    return Omega1

# 3.5PN H_SO:  Omega_1, part 3:
def HS2011_Omega_SO_3p5PN_pt3(m1,m2, n12U, p1U,p2U, r12):
    Omega1 = ixp.zerorank1()
    for i in range(3):
        Omega1[i] = (+(-9*dot(n12U,p1U)*dot(p1U,p1U)/(16*m1**4)
                       +dot(p1U,p1U)*dot(n12U,p2U)/(m1**3*m2)
                       +27*dot(n12U,p1U)*dot(n12U,p2U)**2/(16*m1**2*m2**2)
                       -dot(n12U,p2U)*dot(p1U,p2U)/(8*m1**2*m2**2)
                       -5*dot(n12U,p1U)*dot(p2U,p2U)/(16*m1**2*m2**2))*cross(p1U,p2U)[i])/r12**2
    return Omega1

# 3.5PN H_SO:  Omega_1, part 4:
def HS2011_Omega_SO_3p5PN_pt4(m1,m2, n12U, p1U,p2U, r12):
    Omega1 = ixp.zerorank1()
    for i in range(3):
        Omega1[i] = (+(-3*m2*dot(n12U,p1U)**2/(2*m1**2)
                       +((-3*m2)/(2*m1**2) + 27*m2**2/(8*m1**3))*dot(p1U,p1U)
                       +(177/(16*m1) + 11/m2)*dot(n12U,p2U)**2
                       +(11/(2*m1) + 9*m2/(2*m1**2))*dot(n12U,p1U)*dot(n12U,p2U)
                       +(23/(4*m1) + 9*m2/(2*m1**2))*dot(p1U,p2U)
                       -(159/(16*m1) + 37/(8*m2))*dot(p2U,p2U))*cross(n12U,p1U)[i])/r12**3
    return Omega1

# 3.5PN H_SO:  Omega_1, part 5:
def HS2011_Omega_SO_3p5PN_pt5(m1,m2, n12U, p1U,p2U, r12):
    Omega1 = ixp.zerorank1()
    for i in range(3):
        Omega1[i] = (+(+4*dot(n12U,p1U)**2/m1
                       +13*dot(p1U,p1U)/(2*m1)
                       +5*dot(n12U,p2U)**2/m2
                       +53*dot(p2U,p2U)/(8*m2)
                       -(211/(8*m1) + 22/m2)*dot(n12U,p1U)*dot(n12U,p2U)
                       -(47/(8*m1) + 5/m2)*dot(p1U,p2U))*cross(n12U,p2U)[i])/r12**3
    return Omega1

# 3.5PN H_SO:  Omega_1, part 6:
def HS2011_Omega_SO_3p5PN_pt6(m1,m2, n12U, p1U,p2U, r12):
    Omega1 = ixp.zerorank1()
    for i in range(3):
        Omega1[i] = (+(-(8/m1 + 9*m2/(2*m1**2))*dot(n12U,p1U)
                       +(59/(4*m1) + 27/(2*m2))*dot(n12U,p2U))*cross(p1U,p2U)[i])/r12**3
    return Omega1

# 3.5PN H_SO:  Omega_1, part 7:
def HS2011_Omega_SO_3p5PN_pt7(m1,m2, n12U, p1U,p2U, r12):
    Omega1 = ixp.zerorank1()
    for i in range(3):
        Omega1[i] = (+(181*m1*m2/16 + 95*m2**2/4 + 75*m2**3/(8*m1))*cross(n12U,p1U)[i]
                     -(21*m1**2/2 + 473*m1*m2/16 + 63*m2**2/4)*cross(n12U,p2U)[i])/r12**4
    return Omega1

# 3.5PN H_SO: Combining all the above Omega terms
#             into the full 3.5PN S-O Hamiltonian
#             expression
def f_H_SO_3p5PN(m1,m2, n12U,n21U, S1U, S2U, p1U,p2U, r12):
    Omega1_3p5PNU = ixp.zerorank1()
    Omega2_3p5PNU = ixp.zerorank1()
    for i in range(3):
        Omega1_3p5PNU[i] = HS2011_Omega_SO_3p5PN_pt1(m1,m2, n12U, p1U,p2U, r12)[i]
        Omega1_3p5PNU[i]+= HS2011_Omega_SO_3p5PN_pt2(m1,m2, n12U, p1U,p2U, r12)[i]
        Omega1_3p5PNU[i]+= HS2011_Omega_SO_3p5PN_pt3(m1,m2, n12U, p1U,p2U, r12)[i]
        Omega1_3p5PNU[i]+= HS2011_Omega_SO_3p5PN_pt4(m1,m2, n12U, p1U,p2U, r12)[i]
        Omega1_3p5PNU[i]+= HS2011_Omega_SO_3p5PN_pt5(m1,m2, n12U, p1U,p2U, r12)[i]
        Omega1_3p5PNU[i]+= HS2011_Omega_SO_3p5PN_pt6(m1,m2, n12U, p1U,p2U, r12)[i]
        Omega1_3p5PNU[i]+= HS2011_Omega_SO_3p5PN_pt7(m1,m2, n12U, p1U,p2U, r12)[i]

        Omega2_3p5PNU[i] = HS2011_Omega_SO_3p5PN_pt1(m2,m1, n21U, p2U,p1U, r12)[i]
        Omega2_3p5PNU[i]+= HS2011_Omega_SO_3p5PN_pt2(m2,m1, n21U, p2U,p1U, r12)[i]
        Omega2_3p5PNU[i]+= HS2011_Omega_SO_3p5PN_pt3(m2,m1, n21U, p2U,p1U, r12)[i]
        Omega2_3p5PNU[i]+= HS2011_Omega_SO_3p5PN_pt4(m2,m1, n21U, p2U,p1U, r12)[i]
        Omega2_3p5PNU[i]+= HS2011_Omega_SO_3p5PN_pt5(m2,m1, n21U, p2U,p1U, r12)[i]
        Omega2_3p5PNU[i]+= HS2011_Omega_SO_3p5PN_pt6(m2,m1, n21U, p2U,p1U, r12)[i]
        Omega2_3p5PNU[i]+= HS2011_Omega_SO_3p5PN_pt7(m2,m1, n21U, p2U,p1U, r12)[i]

    global H_SO_3p5PN
    H_SO_3p5PN = dot(Omega1_3p5PNU,S1U) + dot(Omega2_3p5PNU,S2U)
