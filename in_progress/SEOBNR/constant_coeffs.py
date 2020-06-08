# Constants of fit to numerical relativity for the spinning effective one-body formulation

# Import necessary NumPy, SymPy, and SEOBNR modules
import numpy as np

# Compute fits to numerical relativity
def compute_const_coeffs(eta, gamma, a):

    # Define frequently-used constants
    asq = a*a
    pisq = np.pi*np.pi

    # Define constants that determine the fitting parameter K
    # See the discussion in https://arxiv.org/pdf/1311.2544.pdf between Equations (3) and (4)
    K0 = 1.712
    K1 = -1.803949138004582
    K2 = -39.77229225266885
    K3 = 103.16588921239249

    # Compute the fitting parameter K
    # See https://arxiv.org/abs/0912.3517 Equation (5.67) and the discussion following Equation 6.9
    # as well as https://arxiv.org/pdf/1311.2544.pdf
    K = K0 + K1*eta + K2*eta*eta + K3*eta*eta*eta

    # Define more frequently-used constants
    EtaKm1 = eta*K - 1.
    EtaKm1sq = EtaKm1*EtaKm1

    # Compute the Post-Newtonian coefficients
    # See https://arxiv.org/abs/0912.3517 Equations (5.77) to (5.81) and
    # https://arxiv.org/pdf/1311.2544.pdf Equation (2)
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
                        - 5.*Delta2*Delta3 - 5.*Delta1*Delta4)/(5.*EtaKm1sq) + (Delta1ft - 4.*Delta1sq*Delta2
                        + 2.*Delta2sq + 4.*Delta1*Delta3 - 4.*Delta4)/(2*EtaKm1) + (256./5.)*np.log(2))
    Delta5l = (64./5.)*EtaKm1sq

    #Add comment here
    dSO = -74.71 - 156.*eta + 627.5*eta*eta
    dSS = 8.127 - 154.2*eta + 830.8*eta*eta

    return K, Delta0, Delta1, Delta2, Delta3, Delta4, Delta5, Delta5l, dSO, dSS
