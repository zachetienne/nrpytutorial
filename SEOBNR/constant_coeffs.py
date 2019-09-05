#Some of these terms are in the Hamiltonian tutorial.  Move them out so we can use them in the initial conditions
#and not recompute them!

import numpy as np

def compute_const_coeffs(eta, gamma, a):

    asq = a*a
    pisq = np.pi*np.pi

    c20 = 1.712
    c21 = -1.803949138004582
    c22 = -39.77229225266885
    c23 = 103.16588921239249

    K = c20 + c21*eta + c22*eta*eta + c23*eta*eta*eta
    EtaKm1 = eta*K - 1.
    EtaKm1sq = EtaKm1*EtaKm1

    Delta0 = K*(EtaKm1 - 1.)
    Delta1 = -2.*(Delta0 + K)*EtaKm1

    Delta1sq = Delta1*Delta1
    Delta1cu = Delta1*Delta1sq
    Delta1ft = Delta1cu*Delta1

    Delta2 = 0.5*Delta1*(Delta1 - 4.*EtaKm1) - asq*EtaKm1sq*Delta0
    
    Delta2sq = Delta2*Delta2

    Delta3 = -Delta1cu/3. + Delta1*Delta2 + Delta1sq*EtaKm1 - 2.*(Delta2 - EtaKm1)*EtaKm1 - asq*Delta1*EtaKm1sq

    Delta4 = 1./12.*(6*asq*(Delta1sq - 2*Delta2)*EtaKm1sq + 3*Delta1ft - 8*EtaKm1*Delta1cu - 12*Delta2*Delta1sq
                     + 12*(2*EtaKm1*Delta2 + Delta3)*Delta1 + 12*(94./3. - 41./32.*pisq)*EtaKm1sq
                     + 6*(Delta2*Delta2 - 4*Delta3*EtaKm1))

    Delta5 = EtaKm1sq*(-4237./60. + 128./5.*gamma + 2275./512.*pisq - asq*(Delta1cu - 3.*Delta1*Delta2 + 3.*Delta3)/3.
                     - (Delta1ft*Delta1 - 5.*Delta1cu*Delta2 + 5.*Delta1*Delta2sq + 5.*Delta1sq*Delta3
                        - 5.*Delta2*Delta3 - 5.*Delta1*Delta4)/(5.*EtaKm1sq) + (Delta1ft - 4.*(Delta1sq*Delta2)
                        + 2.*Delta2sq + 4.*Delta1*Delta3 - 4.*Delta4)/(2*EtaKm1) + (256./5.)*np.log(2))
    Delta5l = (64./5.)*EtaKm1sq

    d1 = 0.
    d1v2 = -74.71 - 156.*eta + 627.5*eta*eta
    dheffSS = 0.
    dheffSSv2 = 8.127 - 154.2*eta + 830.8*eta*eta

    return K, Delta0, Delta1, Delta2, Delta3, Delta4, Delta5, Delta5l, d1, d1v2, dheffSS, dheffSSv2
