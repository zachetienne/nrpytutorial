## TOV SOLVER FOR SIMPLE POLYTROPES.
## Authors: Phil Chang, Zachariah B. Etienne, Leo Werneck

# Full documentation for this module may be found in the NRPy+ tutorial Jupyter notebook:
#  Tutorial-Start_to_Finish-BSSNCurvilinear-Setting_up_TOV_initial_data.ipynb


##############
# TOV SOLVER #
##############

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
               rho_baryon_central=0.129285,
               rho_poly_tab=[],Gamma_poly_tab=[2.0], K_poly_tab0=1.0,
               verbose = True ):

    #####################################
    # Polytropic EOS lowlevel functions #
    #####################################

    # Function     : impose_continuity_on_P_cold()
    # Author(s)    : Leo Werneck
    # Description  : This function populates the array K_poly_tab
    #                by demanding that P_cold be everywhere continuous
    # Dependencies : none
    #
    # Inputs       : eos              - named tuple containing the following:
    #                  neos           - number of EOSs to be used (single polytrope = 1)
    #                  rho_poly_tab   - values of rho distinguish one EOS from the
    #                                   other (not required for a single polytrope)
    #                  Gamma_poly_tab - values of Gamma to be used within each EOS
    #                  K_poly_tab     - uninitialized, see output variable below
    #                  P_poly_tab     - uninitialized, see function
    #                                   compute_P_poly_tab() below
    #                K_poly_tab0      - value of K_poly_tab[0], for the first EOS
    #
    # Outputs      : eos.K_poly_tab   - values of K to be used within each EOS, determined
    #                                   by imposing that P_cold be everywhere continuous

    def impose_continuity_on_P_cold(eos,K_poly_tab0):

        # A piecewise polytropic EOS is given by
        # .--------------------------------------------------------------------------.
        # |      /     K_0 * rho^(Gamma_0)     ,                rho < rho_0 ;        |
        # |      |     K_1 * rho^(Gamma_1)     ,        rho_0 < rho < rho_1 ;        |
        # |      |          ...                                 ...                  |
        # | P = <      K_j * rho^(Gamma_j)     ,    rho_(j-1) < rho < rho_j ;        |
        # |      |          ...                                 ...                  |
        # |      | K_(n-2) * rho^(Gamma_(n-2)) , rho_(neos-3) < rho < rho_(neos-2) ; |
        # |      \ K_(n-1) * rho^(Gamma_(n-1)) ,                rho > rho_(neos-2) . |
        # .--------------------------------------------------------------------------.
        # Notice that the case of a single polytropic EOS corresponds to
        # the first EOS in the boxed equation above, with no condition on
        # rho. Thus we need only return K_poly_tab0.
        eos.K_poly_tab[0] = K_poly_tab0
        if eos.neos==1:
            return

        # For the case of a piecewise polytropic EOS, emanding that P_cold
        # be everywhere continuous results in the relation:
        # .-----------------------------------------------------.
        # | K_j = K_(j-1) * rho_(j-1)^( Gamma_(j-1) - Gamma_j ) |
        # .-----------------------------------------------------.
        for j in range(1,eos.neos):
            eos.K_poly_tab[j] = eos.K_poly_tab[j-1]*eos.rho_poly_tab[j-1]**(eos.Gamma_poly_tab[j-1]-eos.Gamma_poly_tab[j])

        return

    # Function     : compute_P_poly_tab()
    # Author(s)    : Leo Werneck
    # Description  : This function populates the array eos.P_poly_tab,
    #                used to distinguish which EOS we are using in the
    #                case of a piecewise polytropic EOS
    # Dependencies : none
    #
    # Inputs       : eos              - named tuple containing the following:
    #                  neos           - number of EOSs to be used (single polytrope = 1)
    #                  rho_poly_tab   - values of rho used to distinguish one EOS from
    #                                   the other (not required for a single polytrope)
    #                  Gamma_poly_tab - values of Gamma to be used within each EOS
    #                  K_poly_tab     - value of K to be used within each EOS
    #                  P_poly_tab     - uninitialized, see output variable below
    #
    # Outputs      : eos.P_poly_tab   - values of P used to distinguish one EOS from
    #                                   the other (not required for a single polytrope)

    def compute_P_poly_tab(eos):

        # We now compute the values of P_poly_tab that are used
        # to find the appropriate polytropic index and, thus,
        # EOS we must use.
        # First, if we have a single polytrope EOS, we need to
        # do nothing.
        if eos.neos==1:
            return

        # For the case of a piecewise polytropic EOS, we have
        # .---------------------------.
        # | P_j = K_j*rho_j^(Gamma_j) |
        # .---------------------------.
        for j in range(eos.neos-1):
            eos.P_poly_tab[j] = eos.K_poly_tab[j]*rho_poly_tab[j]**(Gamma_poly_tab[j])

        return

    # Function     : set_single_or_piecewise_polytrope_EOS_parameters()
    # Author(s)    : Leo Werneck
    # Description  : This function determine all polytropic related
    #                parameters from user input
    # Dependencies : impose_continuity_on_P_cold()
    #                compute_P_poly_tab()
    #
    # Inputs       : neos             - number of EOSs to be used (single polytrope = 1)
    #                rho_poly_tab     - values of rho distinguish one EOS from the
    #                                   other (not required for a single polytrope)
    #                Gamma_poly_tab   - values of Gamma to be used within each EOS
    #                K_poly_tab0      - value of K_poly_tab[0], for the first EOS
    #
    # Outputs      : eos              - named tuple containing the following:
    #                  neos           - number of EOSs to be used (single polytrope = 1)
    #                  rho_poly_tab   - values of rho used to distinguish one EOS from
    #                                   the other (not required for a single polytrope)
    #                  Gamma_poly_tab - values of Gamma to be used within each EOS
    #                  K_poly_tab     - value of K to be used within each EOS
    #                  P_poly_tab     - values of P used to distinguish one EOS from
    #                                   the other (not required for a single polytrope)

    def set_single_or_piecewise_polytrope_EOS_parameters(neos,rho_poly_tab,Gamma_poly_tab,K_poly_tab0):

        # Error check #1: Verify if the correct number of rho_poly_tab has been given by the user
        if (neos == 1):
            pass
        elif len(rho_poly_tab) != neos-1:
            print("Error: neos="+str(neos)+". Expected "+str(neos-1)+" values of rho_poly_tab, but "+str(len(rho_poly_tab))+" values were given.")
            sys.exit(1)

        # Error check #2: Verify if the correct number of Gamma_poly_tab has been given by the user
        if len(Gamma_poly_tab) != neos:
            print("Error: neos="+str(neos)+". Expected "+str(neos)+" values of Gamma_poly_tab, but "+str(len(Gamma_poly_tab))+" values were given.")
            sys.exit(2)

        # Create the arrays to store the values of K_poly_tab and eps_integ_const_tab
        K_poly_tab = [0 for i in range(neos)]
        P_poly_tab = [0 for i in range(neos-1)]

        # Create the EOS "struct" (named tuple)
        from collections import namedtuple
        eos_struct = namedtuple("eos_struct","neos rho_poly_tab Gamma_poly_tab K_poly_tab P_poly_tab")
        eos = eos_struct(neos,rho_poly_tab,Gamma_poly_tab,K_poly_tab,P_poly_tab)

        # Step 1: Determine K_poly_tab. For the details, please see the implementation
        #         of the function impose_continuity_on_P_cold() below.
        impose_continuity_on_P_cold(eos,K_poly_tab0)

        # Step 2: Determine eps_integ_const_tab. For the details, please see the
        #         implementation of the function impose_continuity_on_eps_cold() below.
        compute_P_poly_tab(eos)

        return eos

    # Function     : Polytrope_EOS__compute_P_cold_from_rhob()
    # Author(s)    : Leo Werneck
    # Description  : This function computes P_cold for a polytropic EOS
    # Dependencies : polytropic_index_from_rhob()
    #
    # Inputs       : eos              - named tuple containing the following:
    #                  neos           - number of EOSs to be used (single polytrope = 1)
    #                  rho_poly_tab   - values of rho distinguish one EOS from the
    #                                   other (not required for a single polytrope)
    #                  Gamma_poly_tab - values of Gamma to be used within each EOS
    #                  K_poly_tab     - value of K to be used within each EOS
    #                  P_poly_tab     - values of P used to distinguish one EOS from
    #                                   the other (not required for a single polytrope)
    #                rho_baryon       - the value of rho for which we want to
    #                                   compute P_cold
    #
    # Outputs      : P_cold           - for a single or piecewise polytropic EOS

    def Polytrope_EOS__compute_P_cold_from_rhob(eos, rho_baryon):

        # Compute the polytropic index from rho_baryon
        j = polytropic_index_from_rhob(eos, rho_baryon)

        # Return the value of P_cold for a polytropic EOS
        # .--------------------------------.
        # | P_cold = K_j * rho_b^(Gamma_j) |
        # .--------------------------------.
        return eos.K_poly_tab[j]*rho_baryon**eos.Gamma_poly_tab[j]

    # Function     : Polytrope_EOS__compute_rhob_from_P_cold()
    # Author(s)    : Leo Werneck
    # Description  : This function computes rho_b for a polytropic EOS
    # Dependencies : polytropic_index_from_P()
    #
    # Inputs       : eos              - named tuple containing the following:
    #                  neos           - number of EOSs to be used (single polytrope = 1)
    #                  rho_poly_tab   - values of rho distinguish one EOS from the
    #                                   other (not required for a single polytrope)
    #                  Gamma_poly_tab - values of Gamma to be used within each EOS
    #                  K_poly_tab     - value of K to be used within each EOS
    #                  P_poly_tab     - values of P used to distinguish one EOS from
    #                                   the other (not required for a single polytrope)
    #                P                - the value of P for which we want to
    #                                   compute rho_b
    #
    # Outputs      : rho_baryon       - for a single or piecewise polytropic EOS

    def Polytrope_EOS__compute_rhob_from_P_cold(eos,P):

        # Compute the polytropic index from P
        j = polytropic_index_from_P(eos,P)

        # Return the value of rho_b for a polytropic EOS
        # .----------------------------------.
        # | rho_b = (P_cold/K_j)^(1/Gamma_j) |
        # .----------------------------------.
        return (P/eos.K_poly_tab[j])**(1.0/eos.Gamma_poly_tab[j])

    # Function     : polytropic_index_from_rhob()
    # Author(s)    : Leo Werneck and Zach Etienne
    # Description  : This function computes P_cold for a polytropic EOS
    # Dependencies : none
    #
    # Input(s)     : eos              - named tuple containing the following:
    #                  neos           - number of EOSs to be used (single polytrope = 1)
    #                  rho_poly_tab   - values of rho distinguish one EOS from the
    #                                   other (not required for a single polytrope)
    #                  Gamma_poly_tab - values of Gamma to be used within each EOS
    #                  K_poly_tab     - value of K to be used within each EOS
    #                  P_poly_tab     - values of P used to distinguish one EOS from
    #                                   the other (not required for a single polytrope)
    #                rho_in           - value of rho for which we compute the
    #                                   polytropic index
    #
    # Output(s)    : polytropic index computed from rho_in

    def polytropic_index_from_rhob(eos, rho_in):

        # Returns the value of the polytropic index based on rho_in
        polytropic_index = 0
        if not (eos.neos==1):
            for j in range(eos.neos-1):
                polytropic_index += (rho_in > eos.rho_poly_tab[j])

        return polytropic_index

    # Function     : polytropic_index_from_P()
    # Author(s)    : Leo Werneck and Zach Etienne
    # Description  : This function computes P_cold for a polytropic EOS
    # Dependencies : none
    #
    # Input(s)     : eos              - named tuple containing the following:
    #                  neos           - number of EOSs to be used (single polytrope = 1)
    #                  rho_poly_tab   - values of rho distinguish one EOS from the
    #                                   other (not required for a single polytrope)
    #                  Gamma_poly_tab - values of Gamma to be used within each EOS
    #                  K_poly_tab     - value of K to be used within each EOS
    #                  P_poly_tab     - values of P used to distinguish one EOS from
    #                                   the other (not required for a single polytrope)
    #                P_in             - value of P for which we compute the
    #                                   polytropic index
    #
    # Output(s)    : polytropic index computed from P_in

    def polytropic_index_from_P(eos, P_in):

        # Returns the value of the polytropic index based on P_in
        polytropic_index = 0
        if not (eos.neos==1):
            for j in range(eos.neos-1):
                polytropic_index += (P_in > eos.P_poly_tab[j])

        return polytropic_index
    def TOV_rhs(r_Schw, y) : 
    # In \tilde units
    #
        P    = y[0]
        m    = y[1]
        nu   = y[2]
        rbar = y[3]

        j = polytropic_index_from_P(eos,P)
        Gamma = Gamma_poly_tab[j]
        Gam1  = Gamma-1.0

        rho_baryon = Polytrope_EOS__compute_rhob_from_P_cold(eos,P)
        rho = rho_baryon + P/Gam1 # rho is the *total* mass-energy density!

        if( r_Schw < 1e-4 or m <= 0.): 
            m = 4*math.pi/3. * rho*r_Schw**3
            dPdrSchw = -(rho + P)*(4.*math.pi/3.*r_Schw*rho + 4.*math.pi*r_Schw*P)/(1.-8.*math.pi*rho*r_Schw*r_Schw)
            drbardrSchw = 1./(1. - 8.*math.pi*rho*r_Schw*r_Schw)**0.5
        else:
            dPdrSchw = -(rho + P)*(m + 4.*math.pi*r_Schw**3*P)/(r_Schw*r_Schw*(1.-2.*m/r_Schw))
            drbardrSchw = 1./(1. - 2.*m/r_Schw)**0.5*rbar/r_Schw

        dmdrSchw  =  4.*math.pi*r_Schw*r_Schw*rho
        dnudrSchw = -2./(P + rho)*dPdrSchw
        return [dPdrSchw, dmdrSchw, dnudrSchw, drbardrSchw]

    def integrateStar( eos, P, dumpData = False ):
        integrator = si.ode(TOV_rhs).set_integrator('dop853')
        y0 = [P, 0., 0., 0.]
        integrator.set_initial_value(y0,0.)
        dr_Schw = 1e-5
        P = y0[0]

        PArr      = []
        r_SchwArr = []
        mArr      = []
        nuArr     = []
        rbarArr   = []

        r_Schw = 0.

        while integrator.successful() and P > 1e-9*y0[0] : 
            P, m, nu, rbar = integrator.integrate(r_Schw + dr_Schw)
            r_Schw = integrator.t

            dPdrSchw, dmdrSchw, dnudrSchw, drbardrSchw = TOV_rhs( r_Schw+dr_Schw, [P,m,nu,rbar])
            dr_Schw = 0.1*min(abs(P/dPdrSchw), abs(m/dmdrSchw))
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
            rbarArr[ii] *= 0.5*(np.sqrt(R_Schw*(R_Schw - 2.0*M)) + R_Schw - M) / rbarArr[-1]

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
        rho_baryonArr_np = np.array(PArr)
        for j in range(len(PArr_np)):
            # Compute rho_b from P
            rho_baryonArr_np[j] = Polytrope_EOS__compute_rhob_from_P_cold(eos,PArr_np[j])

        mArr_np               = np.array(mArr)
        rbarArr_np            = np.array(rbarArr)
        confFactor_exp4phi_np = (r_SchwArr_np/rbarArr_np)**2

        # Compute the *total* mass-energy density (as opposed to the *baryonic* mass density)
        rhoArr_np = []
        for i in range(len(rho_baryonArr_np)):
            polytropic_index = 0
            if not (eos.neos==1):
                for i in range(eos.neos-1):
                    polytropic_index += (PArr_np[j] > P_poly_tab[i])

            rhoArr_np.append(rho_baryonArr_np[i] + PArr_np[i]/(eos.Gamma_poly_tab[polytropic_index] - 1.))

        print(len(r_SchwArr_np),len(rhoArr_np),len(PArr_np),len(mArr_np),len(exp2phiArr_np))
        # Special thanks to Leonardo Werneck for pointing out this issue with zip()
        if sys.version_info[0] < 3:
            np.savetxt("outputTOVpolytrope.txt", zip(r_SchwArr_np,rhoArr_np,PArr_np,mArr_np,exp2phiArr_np,confFactor_exp4phi_np,rbarArr_np), 
                       fmt="%.15e")
        else:
            np.savetxt("outputTOVpolytrope.txt", list(zip(r_SchwArr_np,rhoArr_np,PArr_np,mArr_np,exp2phiArr_np,confFactor_exp4phi_np,rbarArr_np)), 
                       fmt="%.15e")

        return R_Schw, M

    # Set neos from input variables
    neos = len(Gamma_poly_tab)
    
    # Set polytropic quantities
    eos = set_single_or_piecewise_polytrope_EOS_parameters(neos,rho_poly_tab,Gamma_poly_tab,K_poly_tab0)
    
    # Set initial condition from rho_baryon_central
    P_initial_condition = Polytrope_EOS__compute_P_cold_from_rhob(eos, rho_baryon_central)
    
    # Integrate the initial condition
    R_Schw_TOV, M_TOV = integrateStar(eos, P_initial_condition, True)
    if verbose:
        print("Just generated a TOV star with R_Schw = " + str(R_Schw_TOV) + " , M = " + str(M_TOV) + " , M/R_Schw = "
              + str(M_TOV / R_Schw_TOV) + " .")