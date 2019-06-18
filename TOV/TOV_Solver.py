## TOV SOLVER FOR SIMPLE POLYTROPES.
## Authors: Phil Chang, Zachariah B. Etienne

# Full documentation for this module may be found in the NRPy+ tutorial Jupyter notebook:
#  Tutorial-Start_to_Finish-BSSNCurvilinear-Setting_up_TOV_initial_data.ipynb

# Inputs:
# * Output data file name
# * rho_baryon_central, the central density of the TOV star.
# * n, the polytropic equation of state index. n=1 models cold, degenerate neutron star matter.
# * K_Polytrope, the polytropic constant.
# * Verbose output toggle (default = True)

# Output: An initial data file (default file name = "outputTOVpolytrope.txt") that well
#         samples the (spherically symmetric) solution both inside and outside the star.
#         It is up to the initial data module to perform the 1D interpolation to generate
#         the solution at arbitrary radius. The file has the following columns:
# Column 1: Schwarzschild radius
# Column 2: rho(r), *total* mass-energy density (as opposed to baryonic rest-mass density)
# Column 3: P(r), Pressure
# Column 4: m(r), mass enclosed
# Column 5: e^{nu(r)}, g_{tt}(r)
# Column 6: e^{4 phi(r)}, conformal factor g_{rr}(r)
# Column 7: rbar(r), Isotropic radius

# rbar refers to the isotropic radius, and
# R_Schw refers to the Schwarzschild radius

# Step 1: Import needed Python/NRPy+ modules
import numpy as np
import scipy.integrate as si
import math
import sys

# Step 2: The TOV equations
def TOV_Solver(outfile = "outputTOVpolytrope.txt",
               rho_baryon_central=0.129285, n_Polytrope=1.0, K_Polytrope=1.0,
               verbose = True ):
    gamma = 1. + 1. / n_Polytrope
    gam1 = gamma - 1.

    def TOV_pressure(rho_baryon):
        return K_Polytrope * rho_baryon ** gamma

    def TOV_rhs(r_Schw, y):
        P    = y[0]
        m    = y[1]
        nu   = y[2]
        rbar = y[3]

        dPdrSchw = 0.
        drbardrSchw = 0.

        rho_baryon = (P / K_Polytrope) ** (1. / gamma)
        rho = rho_baryon + P / gam1  # rho is the *total* mass-energy density!
        if (r_Schw < 1e-4 or m <= 0.):
            m = 4 * math.pi / 3. * rho * r_Schw ** 3
            dPdrSchw = -(rho + P) * (4. * math.pi / 3. * r_Schw * rho + 4. * math.pi * r_Schw * P) / (
                        1. - 8. * math.pi * rho * r_Schw * r_Schw)
            drbardrSchw = 1. / (1. - 8. * math.pi * rho * r_Schw * r_Schw) ** 0.5
        else:
            dPdrSchw = -(rho + P) * (m + 4. * math.pi * r_Schw ** 3 * P) / (r_Schw * r_Schw * (1. - 2. * m / r_Schw))
            drbardrSchw = 1. / (1. - 2. * m / r_Schw) ** 0.5 * rbar / r_Schw

        dmdrSchw = 4. * math.pi * r_Schw * r_Schw * rho
        dnudrSchw = -2. / (P + rho) * dPdrSchw
        return [dPdrSchw, dmdrSchw, dnudrSchw, drbardrSchw]

    def integrateStar(P, dumpData=False):
        integrator = si.ode(TOV_rhs).set_integrator('dop853')
        y0 = [P, 0., 0., 0.]
        integrator.set_initial_value(y0, 0.)
        dr_Schw = 1e-5
        P = y0[0]

        PArr = []
        r_SchwArr = []
        mArr = []
        nuArr = []
        rbarArr = []

        r_Schw = 0.

        while integrator.successful() and P > 1e-9 * y0[0]:
            P, m, nu, rbar = integrator.integrate(r_Schw + dr_Schw)
            r_Schw = integrator.t

            dPdrSchw, dmdrSchw, dnudrSchw, drbardrSchw = TOV_rhs(r_Schw + dr_Schw, [P, m, nu, rbar])
            dr_Schw = 0.1 * min(abs(P / dPdrSchw), abs(m / dmdrSchw))
            dr_Schw = min(dr_Schw, 1e-2)
            PArr.append(P)
            r_SchwArr.append(r_Schw)
            mArr.append(m)
            nuArr.append(nu)
            rbarArr.append(rbar)

        M = mArr[-1]
        R_Schw = r_SchwArr[-1]

        # Apply integration constant to ensure rbar is continuous across TOV surface
        for ii in range(len(rbarArr)):
            rbarArr[ii] *= 0.5 * (np.sqrt(R_Schw * (R_Schw - 2.0 * M)) + R_Schw - M) / rbarArr[-1]

        nuArr_np = np.array(nuArr)
        # Rescale solution to nu so that it satisfies BC: exp(nu(R))=exp(nutilde-nu(r=R)) * (1 - 2m(R)/R)
        #   Thus, nu(R) = (nutilde - nu(r=R)) + log(1 - 2*m(R)/R)
        nuArr_np = nuArr_np - nuArr_np[-1] + math.log(1. - 2. * mArr[-1] / r_SchwArr[-1])

        r_SchwArrExtend_np = 10. ** (np.arange(0.01, 5.0, 0.01)) * r_SchwArr[-1]

        r_SchwArr.extend(r_SchwArrExtend_np)
        mArr.extend(r_SchwArrExtend_np * 0. + M)
        PArr.extend(r_SchwArrExtend_np * 0.)
        exp2phiArr_np = np.append(np.exp(nuArr_np), 1. - 2. * M / r_SchwArrExtend_np)
        nuArr.extend(np.log(1. - 2. * M / r_SchwArrExtend_np))
        rbarArr.extend(0.5 * (np.sqrt(r_SchwArrExtend_np ** 2 - 2. * M * r_SchwArrExtend_np) + r_SchwArrExtend_np - M))

        # Appending to a Python array does what one would reasonably expect.
        #   Appending to a numpy array allocates space for a new array with size+1,
        #   then copies the data over... over and over... super inefficient.
        r_SchwArr_np = np.array(r_SchwArr)
        PArr_np = np.array(PArr)
        rho_baryonArr_np = (PArr_np / K_Polytrope) ** (1. / gamma)
        mArr_np = np.array(mArr)
        rbarArr_np = np.array(rbarArr)
        confFactor_exp4phi_np = (r_SchwArr_np / rbarArr_np) ** 2

        # Compute the *total* mass-energy density (as opposed to the *baryonic* mass density)
        rhoArr_np = []
        for i in range(len(rho_baryonArr_np)):
            rhoArr_np.append(rho_baryonArr_np[i] + PArr_np[i] / (gamma - 1.))

        #print(len(r_SchwArr_np), len(rhoArr_np), len(PArr_np), len(mArr_np), len(exp2phiArr_np))
        # Special thanks to Leonardo Werneck for pointing out this issue with zip()
        if sys.version_info[0] < 3:
            np.savetxt(outfile,zip(r_SchwArr_np,rhoArr_np,PArr_np,mArr_np,exp2phiArr_np,confFactor_exp4phi_np,rbarArr_np), 
                       fmt="%.15e")
        else:
            np.savetxt(outfile,
                       list(zip(r_SchwArr_np,rhoArr_np,PArr_np,mArr_np,exp2phiArr_np,confFactor_exp4phi_np,rbarArr_np)), 
                       fmt="%.15e")
        return R_Schw, M

    R_Schw_TOV, M_TOV = integrateStar(TOV_pressure(rho_baryon_central), True)
    if verbose:
        print("Just generated a TOV star with R_Schw = " + str(R_Schw_TOV) + " , M = " + str(M_TOV) + " , M/R_Schw = "
              + str(M_TOV / R_Schw_TOV) + " .")