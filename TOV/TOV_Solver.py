## TOV SOLVER FOR SIMPLE POLYTROPES.
## Authors: Phil Chang, Zachariah B. Etienne, Leo Werneck

# Full documentation for this module may be found in the NRPy+ tutorial Jupyter notebook:
#  Tutorial-Start_to_Finish-BSSNCurvilinear-Setting_up_TOV_initial_data.ipynb


##############
# TOV SOLVER #
##############

# Inputs:
# * Output data file name
# * eos, a named tuple containing all EOS parameters.
# * rho_baryon_central, the central density of the TOV star.
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
import numpy as np                  # NumPy: A numerical methods module for Python
import scipy.integrate as si        # SciPy: Python module for mathematics, science, and engineering applications
import math, sys                    # Standard Python modules for math; multiplatform OS-level functions
import TOV.Polytropic_EOSs as ppeos # NRPy+: Piecewise polytrope equation of state support

# Step 2: The TOV equations
def TOV_Solver(eos,
               outfile = "outputTOVpolytrope.txt",
               rho_baryon_central = 0.129285,
               verbose = True,
               return_M_RSchw_and_Riso = False,
               accuracy = "medium",
               integrator_type = "default",
               no_output_File = False):

    def TOV_rhs(r_Schw, y) :
    # In \tilde units
    #
        P    = y[0]
        m    = y[1]
        # nu   = y[2] # nu is not needed as input into TOV_rhs
        rbar = y[3]

        # Compute rho_b and eps_cold, to be used below
        # to compute rho_(total)
        rho_baryon, eps_cold = ppeos.Polytrope_EOS__compute_rhob_and_eps_cold_from_P_cold(eos,P)

#         with open("rhob_P_cold_and_eps_cold.dat","a+") as file:
#             file.write(str(r_Schw).format("%.15e")+"  "+str(rho_baryon).format("%.15e")+"  "+str(P).format("%.15e")+"  "+str(eps_cold).format("%.15e")+"\n")

        # Compute rho, the *total* mass-energy density:
        # .------------------------------.
        # | rho = (1 + eps)*rho_(baryon) |
        # .------------------------------.
        # with eps = eps_cold, for the initial data.
        rho = (1.0 + eps_cold)*rho_baryon

#         m = 4*math.pi/3. * rho*r_Schw**3
        if( r_Schw < 1e-4 or m <= 0.):
            # From https://github.com/natj/tov/blob/master/tov.py#L33:
            # dPdr = -cgs.G*(eden + P/cgs.c**2)*(m + 4.0*pi*r**3*P/cgs.c**2)
            # dPdr = dPdr/(r*(r - 2.0*cgs.G*m/cgs.c**2))
            dPdrSchw = -(rho + P)*(4.*math.pi/3.*r_Schw*rho + 4.*math.pi*r_Schw*P)/(1.-8.*math.pi*rho*r_Schw*r_Schw)
            drbardrSchw = 1./(1. - 8.*math.pi*rho*r_Schw*r_Schw)**0.5
        else:
            dPdrSchw = -(rho + P)*(m + 4.*math.pi*r_Schw**3*P)/(r_Schw*r_Schw*(1.-2.*m/r_Schw))
            drbardrSchw = 1./(1. - 2.*m/r_Schw)**0.5*rbar/r_Schw

        dmdrSchw  =  4.*math.pi*r_Schw*r_Schw*rho
        dnudrSchw = -2./(P + rho)*dPdrSchw
        return [dPdrSchw, dmdrSchw, dnudrSchw, drbardrSchw]

    def integrateStar( eos, P, dumpData = False ):

        if accuracy == "medium":
            min_step_size = 1e-5
            max_step_size = 1e-2
            integrator    = 'dop853'
        elif accuracy == "low":
            min_step_size = 1e-3
            max_step_size = 1e-1
            integrator    = 'dopri5'
        elif accuracy == "verylow":
            min_step_size = 1e-1
            max_step_size = 5e-1
            integrator    = 'dopri5'
        elif accuracy == "high":
            min_step_size = 1e-5
            max_step_size = 1e-5
            integrator    = 'dop853'
        elif accuracy == "veryhigh":
            min_step_size = 1e-7
            max_step_size = 1e-6
            integrator    = 'dop853'
        else:
            print("Unknown accuracy option: "+str(accuracy))

        if integrator_type == "default":
            pass
        else:
            integrator = integrator_type

        integrator = si.ode(TOV_rhs).set_integrator(integrator)#,rtol=1e-4,atol=1e-4)
#         integrator = si.ode(TOV_rhs).set_integrator('dopri5',rtol=1e-4)
        y0 = [P, 0., 0., 0.]
        r_Schw = 0.
        integrator.set_initial_value(y0,r_Schw)
        dr_Schw = min_step_size
        P = y0[0]

        PArr      = []
        r_SchwArr = []
        mArr      = []
        nuArr     = []
        rbarArr   = []

        while integrator.successful() and P > 1e-19*y0[0] :
            P, m, nu, rbar = integrator.integrate(r_Schw + dr_Schw)
            # Update the value of r_Schw to the latest integrated value
            r_Schw += dr_Schw

            dPdrSchw, dmdrSchw, dnudrSchw, drbardrSchw = TOV_rhs( r_Schw, [P,m,nu,rbar])
            dr_Schw = 0.1*min(abs(P/dPdrSchw), abs(m/dmdrSchw))
            dr_Schw = min(dr_Schw, max_step_size)
            PArr.append(P)
            r_SchwArr.append(r_Schw)
            mArr.append(m)
            nuArr.append(nu)
            rbarArr.append(rbar)

        M      = mArr[-1]
        R_Schw = r_SchwArr[-1]
        # Must multiply integration constant to get R_iso to ensure rbar(==R_iso) continuous across TOV surface
        R_iso = rbarArr[-1]*0.5*(np.sqrt(R_Schw*(R_Schw - 2.0*M)) + R_Schw - M) / rbarArr[-1]

        # Apply integration constant to ensure rbar is continuous across TOV surface
        for ii in range(len(rbarArr)):
            rbarArr[ii] *= 0.5*(np.sqrt(R_Schw*(R_Schw - 2.0*M)) + R_Schw - M) / rbarArr[-1]

        if no_output_File == True:
            return M, R_Schw, R_iso

        nuArr_np = np.array(nuArr)
        # Rescale solution to nu so that it satisfies BC: exp(nu(R))=exp(nutilde-nu(r=R)) * (1 - 2m(R)/R)
        #   Thus, nu(R) = (nutilde - nu(r=R)) + log(1 - 2*m(R)/R)
        nuArr_np = nuArr_np - nuArr_np[-1] + math.log(1.-2.*mArr[-1]/r_SchwArr[-1])

        r_SchwArrExtend_np = 10.**(np.arange(0.01,5.0,0.01))*r_SchwArr[-1]

        r_SchwArr.extend(r_SchwArrExtend_np)
        mArr.extend(r_SchwArrExtend_np*0. + M)
        PArr.extend(r_SchwArrExtend_np*0.)
        exp2phiArr_np = np.append( np.exp(nuArr_np), 1. - 2.*M/r_SchwArrExtend_np)
        nuArr.extend(np.log(1. - 2.*M/r_SchwArrExtend_np))
        rbarArr.extend( 0.5*(np.sqrt(r_SchwArrExtend_np**2 - 2.*M*r_SchwArrExtend_np) + r_SchwArrExtend_np - M) )

        # Appending to a Python array does what one would reasonably expect.
        #   Appending to a numpy array allocates space for a new array with size+1,
        #   then copies the data over... over and over... super inefficient.
        r_SchwArr_np     = np.array(r_SchwArr)
        PArr_np          = np.array(PArr)
        rho_baryonArr_np = np.array(PArr) # This is just to initialize the array
        for j in range(len(PArr_np)):
            # Compute rho_b from P
            rho_baryonArr_np[j] = ppeos.Polytrope_EOS__compute_rhob_from_P_cold(eos,PArr_np[j])

        mArr_np               = np.array(mArr)
        rbarArr_np            = np.array(rbarArr)
        confFactor_exp4phi_np = (r_SchwArr_np/rbarArr_np)**2

        # Compute the *total* mass-energy density (as opposed to the *baryonic* mass density)
        rhoArr_np = []
        for i in range(len(PArr)):
            rho_baryon, eps_cold = ppeos.Polytrope_EOS__compute_rhob_and_eps_cold_from_P_cold(eos,PArr[i])
            rho = (1.0 + eps_cold ) * rho_baryon
            rhoArr_np.append(rho)

        if verbose:
            print(len(r_SchwArr_np),len(rhoArr_np),len(rho_baryonArr_np),len(PArr_np),len(mArr_np),len(exp2phiArr_np))

        # Special thanks to Leonardo Werneck for pointing out this issue with zip()
        if sys.version_info[0] < 3:
            np.savetxt(outfile, zip(r_SchwArr_np,rhoArr_np,rho_baryonArr_np,PArr_np,mArr_np,exp2phiArr_np,confFactor_exp4phi_np,rbarArr_np),
                       fmt="%.15e")
        else:
            np.savetxt(outfile, list(zip(r_SchwArr_np,rhoArr_np,rho_baryonArr_np,PArr_np,mArr_np,exp2phiArr_np,confFactor_exp4phi_np,rbarArr_np)),
                       fmt="%.15e")

        return M, R_Schw, R_iso

    # Set initial condition from rho_baryon_central
    P_initial_condition = ppeos.Polytrope_EOS__compute_P_cold_from_rhob(eos, rho_baryon_central)

    # Integrate the initial condition
    M_TOV, R_Schw_TOV, R_iso_TOV = integrateStar(eos, P_initial_condition, True)
    if verbose:
        print("""Just generated a TOV star with
* M        = %.15e ,
* R_Schw   = %.15e ,
* R_iso    = %.15e ,
* M/R_Schw = %.15e \n""" %(M_TOV,R_Schw_TOV,R_iso_TOV,(M_TOV / R_Schw_TOV)))

    if return_M_RSchw_and_Riso:
        return M_TOV, R_Schw_TOV, R_iso_TOV
