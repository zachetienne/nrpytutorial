## Polytropic EOSs Python Module
## Author(s): Leo Werneck and Zach Etienne
## In this NRPy+ module we set up useful "lowlevel" functions that compute
## useful polytropic quantities.

# Full documentation for this module may be found in the NRPy+ tutorial Jupyter notebook:
#  Tutorial-TOV-Piecewise_Polytrope_EOSs.ipynb

# This ensures the availability of the argument end="" in the print() function.
from __future__ import print_function

# Step 0: Import needed Python/NRPy+ modules
import numpy as np                  # NumPy: A numerical methods module for Python
import sys                          # This module is used for system related function calls
from collections import namedtuple  # This module is used to create named tuples


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
        eos.P_poly_tab[j] = eos.K_poly_tab[j]*eos.rho_poly_tab[j]**(eos.Gamma_poly_tab[j])

    return


# Function     : impose_continuity_on_eps_cold()
# Author(s)    : Leo Werneck
# Description  : This function populates the array eps_integ_const_tab
#                by demanding that eps_cold be everywhere continuous
# Dependencies : none
#
# Inputs       : eos                     - named tuple containing the following:
#                  neos                  - number of EOSs to be used (single polytrope = 1)
#                  rho_poly_tab          - values of rho distinguish one EOS from the
#                                          other (not required for a single polytrope)
#                  Gamma_poly_tab        - values of Gamma to be used within each EOS
#                  K_poly_tab            - value of K to be used within each EOS
#                  eps_integ_const_tab   - uninitialized, see output variable below
#
# Outputs      : eos.eps_integ_const_tab - value of C used to compute eps_cold within each EOS,
#                                          determined by imposing that eps_cold be everywhere
#                                          continuous

def impose_continuity_on_eps_cold(eos):

    # Computing eps_cold for the case of a polytropic EOS, we have
    # .------------------------------------------------------------------------------------------------------.
    # |        / C_0     + K_0*rho^(Gamma_0 - 1)/(Gamma_0 - 1)         ,                rho < rho_0 ;        |
    # |        | C_1     + K_1*rho^(Gamma_1 - 1)/(Gamma_1 - 1)         ,        rho_0 < rho < rho_1 ;        |
    # |        |                       ...                                              ...                  |
    # | eps = <  C_j     + K_j*rho^(Gamma_j - 1)/(Gamma_j - 1)         ,    rho_(j-1) < rho < rho_j ;        |
    # |        |                       ...                                              ...                  |
    # |        | C_(n-2) + K_(n-2)*rho^(Gamma_(n-2)-1)/(Gamma_(n-2)-1) , rho_(neos-3) < rho < rho_(neos-2) ; |
    # |        \ C_(n-1) + K_(n-1)*rho^(Gamma_(n-1)-1)/(Gamma_(n-1)-1) ,                rho > rho_(neos-2) . |
    # .------------------------------------------------------------------------------------------------------.
    # By demanding that eps_cold(rho -> 0) = 0, we fix C_0 = 0. Thus, for
    # a single polytrope we need only return this
    if eos.neos==1:
        return

    # For the case of a piecewise polytropic EOS, emanding that eps_cold
    # be everywhere continuous results in the relation:
    # .-----------------------------------------------------------------.
    # | C_j = C_(j-1)                                                   |
    # |     + K_(j-1)*rho_(j-1)^( Gamma_(j-1) - 1 )/( Gamma_(j-1) - 1 ) |
    # |     - K_(j+0)*rho_(j-1)^( Gamma_(j+0) - 1 )/( Gamma_(j+0) - 1 ) |
    # .-----------------------------------------------------------------.
    eos.eps_integ_const_tab[0] = 0.0
    for j in range(1,eos.neos):
        # Second line of the boxed equation above
        aux_jm1 = eos.K_poly_tab[j-1]*eos.rho_poly_tab[j-1]**(eos.Gamma_poly_tab[j-1]-1.0)/(eos.Gamma_poly_tab[j-1]-1)

        # Third line of the boxed equation above
        aux_jp0 = eos.K_poly_tab[j+0]*eos.rho_poly_tab[j-1]**(eos.Gamma_poly_tab[j+0]-1.0)/(eos.Gamma_poly_tab[j+0]-1)

        # Boxed equation above
        eos.eps_integ_const_tab[j] = eos.eps_integ_const_tab[j-1] + aux_jm1 - aux_jp0

    return


# Function     : set_up_EOS_parameters__complete_set_of_input_variables()
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

def set_up_EOS_parameters__complete_set_of_input_variables(neos,rho_poly_tab,Gamma_poly_tab,K_poly_tab0):

    # Error check #1: Verify if the correct number of rho_poly_tab has been given by the user
    if (neos == 1):
        pass
    elif len(rho_poly_tab) != neos-1:
        print("Error: neos="+str(neos)+". Expected "+str(neos-1)+" values of rho_poly_tab, but "+str(len(rho_poly_tab))+" values were given.")
        sys.exit(1)

    # Error check #2: Verify if the correct number of Gamma_poly_tab has been given by the user
    if len(Gamma_poly_tab) != neos:
        print("Error: neos="+str(neos)+". Expected "+str(neos)+" values of Gamma_poly_tab, but "+str(len(Gamma_poly_tab))+" values were given.")
        sys.exit(1)

    # Create the arrays to store the values of K_poly_tab and eps_integ_const_tab
    K_poly_tab          = [0 for i in range(neos)]
    P_poly_tab          = [0 for i in range(neos-1)]
    eps_integ_const_tab = [0 for i in range(neos)]

    # Create the EOS "struct" (named tuple)
    eos_struct = namedtuple("eos_struct","neos rho_poly_tab Gamma_poly_tab K_poly_tab P_poly_tab eps_integ_const_tab")
    eos = eos_struct(neos,rho_poly_tab,Gamma_poly_tab,K_poly_tab,P_poly_tab,eps_integ_const_tab)

    # Step 1: Determine K_poly_tab. For the details, please see the implementation
    #         of the function impose_continuity_on_P_cold() below.
    impose_continuity_on_P_cold(eos,K_poly_tab0)

    # Step 2: Determine eps_integ_const_tab. For the details, please see the
    #         implementation of the function impose_continuity_on_eps_cold() below.
    impose_continuity_on_eps_cold(eos)

    # Step 3: Determine P_poly_tab. For the details, please see the implementation
    #         of the function compute_P_poly_tab() below.
    compute_P_poly_tab(eos)

    return eos


# Function     : set_up_EOS_parameters__Read_et_al_input_variables()
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

def set_up_EOS_parameters__Read_et_al_input_variables(EOSname,units="rescaledensity"):

    # Check if the input units are implemented below
    available_units = ["rescaledensity","geometrized","cgs"]
    if units not in available_units:
        print("ERROR: unknown units ",units)
        print("Available units are: ",end="")
        for unit in available_units:
            print(unit,end=" ")
        print("")
        sys.exit(1)

    # Set up the number of polytropic EOSs, which is
    # fixed at seven for this type of input
    neos = 7

    # Set up input from table II of Read et al. (2008),
    # and from the legend in FIG. 3.
    # Source: https://arxiv.org/pdf/0812.2163.pdf
    # .--------------.---------.-------------.-----.
    # |    rho_j     | Gamma_j |     K_j     | P_j |
    # .--------------.---------.-------------.-----.
    # | 2.440340e+07 | 1.58425 | 6.80110e-09 | P_0 |
    # | 3.78358e+11  | 1.28733 |     K_1     | P_1 |
    # | 2.62780e+12  | 0.62223 |     K_2     | P_2 |
    # |    rho_3     | 1.35692 |     K_3     | P_3 |
    # |  10^(14.7)   | Gamma_4 |     K_4     | P_4 |
    # |  10^(15.0)   | Gamma_5 |     K_5     | P_5 |
    # |       -      | Gamma_6 |     K_6     |  -  |
    # .--------------.---------.-------------.-----.
    #
    # Load up the NRPy+ dictionary containing information about the Read et al.
    # EOS tables (table III in Read et al. - https://arxiv.org/pdf/0812.2163.pdf)
    import TOV.Piecewise_Polytrope__dict as PPdict
    log_of_p4 = PPdict.EOS_Read_et_al_dict[EOSname].log_of_p4
    Gamma4    = PPdict.EOS_Read_et_al_dict[EOSname].Gamma4
    Gamma5    = PPdict.EOS_Read_et_al_dict[EOSname].Gamma5
    Gamma6    = PPdict.EOS_Read_et_al_dict[EOSname].Gamma6

    # Set up the speed of light and change the units of the input pressure
    c = 2.997924580000000e+10 # Speed of light
    G = 6.674299999999999e-08 # Gravitational constant
    M = 1.988409870698051e+33 # Mass of the sun
    log_of_p4 -= 2.0*np.log10(c)

    # Set up tabulated polytropic values following the table above
    # and the user input. All quantities which are still unknown are
    # set to absurd values to make sure they are overwritten
    rho_poly_tab   = [2.440340e+07, 3.78358e+11, 2.62780e+12, -1e30  , 10**(14.7)     , 10**(15.0)]
    P_poly_tab     = [-1e30       , -1e30      , -1e30      , -1e30  , 10**(log_of_p4), -1e30     ]
    Gamma_poly_tab = [1.58425     , 1.28733    , 0.62223    , 1.35692, Gamma4         , Gamma5, Gamma6]
    K_poly_tab     = [6.80110e-09 , -1e30      , -1e30      , -1e30  , -1e30          , -1e30  , -1e30]

    # Compute {K_1,K_2,K_3}, using
    # .-----------------------------------------------------.
    # | K_j = K_(j-1) * rho_(j-1)^( Gamma_(j-1) - Gamma_j ) |
    # .-----------------------------------------------------.
    for j in range(1,4):
        K_poly_tab[j] = K_poly_tab[j-1] * rho_poly_tab[j-1]**(Gamma_poly_tab[j-1] - Gamma_poly_tab[j])

    # Compute {P_0,P_1,P_2}, using
    # .-------------------------------.
    # | P_j = K_j * rho_j^( Gamma_j ) |
    # .-------------------------------.
    for j in range(3):
        P_poly_tab[j] = K_poly_tab[j] * rho_poly_tab[j]**(Gamma_poly_tab[j])


    # Set up auxiliary variables for the evaluation of rho_3
    P4          = P_poly_tab[4]
    K3          = K_poly_tab[3]
    rho4_p_Gam4 = rho_poly_tab[4]**(Gamma_poly_tab[4])
    G3m4        = Gamma_poly_tab[3] - Gamma_poly_tab[4]

    # Compute rho_3 using
    # .----------------------------------------------------------------------.
    # | rho_3 = ( P_4 /( K_3 * rho_4^(Gamma_4) ) )^(1.0/(Gamma_3 - Gamma_4)) |
    # .----------------------------------------------------------------------.
    rho_poly_tab[3] = ( P4/(K3 * rho4_p_Gam4) )**(1.0/G3m4)

    # Compute {P_3,P_4,P_5} and {K_4,K_5,K_6}
    for j in range(3,neos-1):
        P_poly_tab[j]   = K_poly_tab[j] * rho_poly_tab[j]**(Gamma_poly_tab[j])
        K_poly_tab[j+1] = K_poly_tab[j] * rho_poly_tab[j]**(Gamma_poly_tab[j] - Gamma_poly_tab[j+1])

    if units == "rescaledensity":

        # We impose a "ratio preserving rescaling" of rhob:
        #
        # rhob_rescaled[j] / rhob[j] = rhob_rescaled[j-1] / rhob[j-1]
        #
        # which implies the relation
        # .-------------------------------------------------------------.
        # | rhob_rescaled[j-1] = (rhob[j-1]/rhob[j]) * rhob_rescaled[j] |
        # .-------------------------------------------------------------.
        # after setting rhob_nuclear_rescaled = 1
        rhob_rescaled  = [1.0 for i in range(neos-1)]
        for j in range(neos-2,0,-1):
            rhob_rescaled[j-1] = (rho_poly_tab[j-1]/rho_poly_tab[j]) * rhob_rescaled[j]

        # Now because the values of P and rho given by Read et al. are already
        # in the same units, namely (g/cm^3), the ratio P/rho should be invariant
        # under this rescalling procedure. Therefore
        # .---------------------------------------------------------------------.
        # | P_rescaled[j] = (rhob_rescaled[j]/rhob_readetal[j]) * P_readetal[j] |
        # .---------------------------------------------------------------------.
        P_rescaled = [0.0 for i in range(neos-1)]
        for j in range(neos-1):
            P_rescaled[j] = (rhob_rescaled[j]/rho_poly_tab[j]) * P_poly_tab[j]

        rho_poly_tab = rhob_rescaled
        P_poly_tab   = P_rescaled

    elif units == "geometrized" :
        # Now convert to units in which Msun = 1, G = 1, c = 1
        csq               = c**2
        units_of_length   = G * M / csq
        units_of_time     = units_of_length/c
        units_of_mass     = M
        units_of_density  = units_of_mass / units_of_length**3
        units_of_pressure = units_of_mass / units_of_length / units_of_time**2

        for i in range(neos-1):
            rho_poly_tab[i] /= units_of_density
            P_poly_tab[i]   *= csq
            P_poly_tab[i]   /= units_of_pressure

    elif units == "cgs" :

        # Restore P to cgs units
        csq = c*c
        for i in range(neos-1):
            P_poly_tab[i] *= csq

    # Demanding that the pressure be everywhere continuous then imposes
    # .-------------------------------------------------------------------------------------------.
    # | K_dimensionless[j-1] = K_dimensionless[j]/rhob_dimensionless[j-1]^(Gamma[j-1] - Gamma[j]) |
    # .-------------------------------------------------------------------------------------------.
    K_poly_tab[0] = P_poly_tab[0]/rho_poly_tab[0]**(Gamma_poly_tab[0])
    for j in range(1,neos):
        K_poly_tab[j] = K_poly_tab[j-1]*rho_poly_tab[j-1]**(Gamma_poly_tab[j-1]-Gamma_poly_tab[j])

    # Allocate memory for the integration constants of eps_cold
    eps_integ_const_tab = [0 for i in range(neos)]

    # Create the EOS "struct" (named tuple)
    eos_struct = namedtuple("eos_struct","neos rho_poly_tab Gamma_poly_tab K_poly_tab P_poly_tab eps_integ_const_tab")
    eos = eos_struct(neos,rho_poly_tab,Gamma_poly_tab,K_poly_tab,P_poly_tab,eps_integ_const_tab)

    # Populate the integration constants of eps_cold
    impose_continuity_on_eps_cold(eos)

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


# Function     : Polytrope_EOS__compute_eps_cold_from_rhob()
# Author(s)    : Leo Werneck
# Description  : This function computes eps_cold for a polytropic EOS
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
# Outputs      : eps_cold         - for a single or piecewise polytropic EOS

def Polytrope_EOS__compute_eps_cold_from_rhob(eos, rho_baryon):

    if rho_baryon == 0.0:
        return 0.0

    # Compute the polytropic index from rho_baryon
    j = polytropic_index_from_rhob(eos, rho_baryon)

    # Compute P_cold
    P_cold = Polytrope_EOS__compute_P_cold_from_rhob(eos, rho_baryon)

    # Return the value of P_cold for a polytropic EOS
    # .----------------------------------------------.
    # | eps_cold = C_j + P_cold/( rhob*(Gamma_j-1) ) |
    # .----------------------------------------------.
    return ( eos.eps_integ_const_tab[j] + P_cold/(rho_baryon*(eos.Gamma_poly_tab[j] - 1.0)) )

# Function     : Polytrope_EOS__compute_rhob_and_eps_cold_from_P_cold()
# Author(s)    : Leo Werneck
# Description  : This function computes rho_b and eps_cold for a polytropic EOS
# Dependencies : polytropic_index_from_P()
#                Polytrope_EOS__compute_rhob_from_P_cold()
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

def Polytrope_EOS__compute_rhob_and_eps_cold_from_P_cold(eos,P):

    # Compute the polytropic index from P and set Gamma
    j     = polytropic_index_from_P(eos,P)
    Gamma = eos.Gamma_poly_tab[j]
    # Compute the value of rho_b for a polytropic EOS
    # .----------------------------------.
    # | rho_b = (P_cold/K_j)^(1/Gamma_j) |
    # .----------------------------------.
    rho_b = (P/eos.K_poly_tab[j])**(1.0/Gamma)

    return rho_b, Polytrope_EOS__compute_eps_cold_from_rhob(eos, rho_b)


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


# Function     : generate_IllinoisGRMHD_EOS_parameter_file()
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
# Output(s)    : parameter file to be used by IllinoisGRMHD

def generate_IllinoisGRMHD_EOS_parameter_file(EOSname,outfilename, \
                                              Gamma_thermal=None,  \
                                              EOS_struct=None, \
                                              tau_atmosphere=4.876083025795607e-12, \
                                              rho_atmosphere=1.292852735094440e-10, \
                                              K_single_polytrope=1.0, \
                                              Gamma_single_polytrope=2.0):

    with open(outfilename,"w") as file:
        file.write("""
#vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
#
#.-------------------------------------------------------------------------.
#| IllinoisGRMHD Equation of State (EOS) parameter file Generated by NRPy+ |
#|-------------------------------------------------------------------------|
#|         This section of the parameter file has been generated by        |
#|                 the TOV/Polytropic_EOSs.py NRPy+ module                 |
#|-------------------------------------------------------------------------|
#|    Recommended reading: Tutorial-TOV-Piecewise_Polytrope_EOSs.ipynb     |
#|-------------------------------------------------------------------------|
#| NRPy+ repositoryon github: https://github.com/zachetienne/nrpytutorial/ |
#|-------------------------------------------------------------------------|
#| *Warning*: it is highly recommended not to change this section manually |
#.-------------------------------------------------------------------------.
""")
    if EOSname == "single":
        with open(outfilename,"a") as file:
            file.write("""#
#.-------------------------------.
#|  EOS Type: Single Polytrope   |
#.-------------------------------.
#| Required inputs:              |
#|   - K_single_polytrope        |
#|   - Gamma_single_polytrope    |
#|   - tau_atmosphere            |
#|   - rho_atmosphere            |
#|                               |
#| IllinoisGRMHD parameters set: |
#|   - neos                      |
#|   - K_ppoly_tab0              |
#|   - rho_ppoly_tab_in[0]       |
#|   - Gamma_ppoly_tab_in[0]     |
#|   - Gamma_th                  |
#|   - tau_atm                   |
#|   - rho_b_atm                 |
#|                               |
#| NRPyPlusTOVID parameters set: |
#|   - rho_atmosphere            |
#|   - Gamma_atmosphere          |
#|   - K_atmosphere              |
#|-------------------------------|
#| For single polytropes, we     |
#| always assume:                |
#| Gamma_th = Gamma_poly_tab     |
#.-------------------------------.
#
# Set up initial data file name
NRPyPlusTOVID::TOV_filename = "outputTOVpolytrope.txt"

# Set the number of EOSs to 1 (single polytrope)
IllinoisGRMHD::neos = 1

# Set atmospheric value of tau
IllinoisGRMHD::tau_atm = %.15e

# Set atmospheric value of rho
IllinoisGRMHD::rho_b_atm      = %.15e
NRPyPlusTOVID::rho_atmosphere = %.15e

# Set K_ppoly_tab0 and K_atmosphere
IllinoisGRMHD::K_ppoly_tab0 = %.15e
NRPyPlusTOVID::K_atmosphere = %.15e

# Set Gamma_ppoly_tab_in[0]  and Gamma_atmosphere
IllinoisGRMHD::Gamma_ppoly_tab_in[0] = %.15e
NRPyPlusTOVID::Gamma_atmosphere      = %.15e

# Set Gamma_thermal
# (must be the same as Gamma_ppoly_tab for a single polytrope)
IllinoisGRMHD::Gamma_th = %.15e

# Set rho_ppoly_tab_in[0] to zero
# (for a single polytrope this value is not used)
IllinoisGRMHD::rho_ppoly_tab_in[0] = 0.0

#.----------------------.
#| EOS_Omni parameters: |
#|  - n_pieces          |
#|  - hybrid_k0         |
#|  - hybrid_gamma[0]   |
#|  - hybrid_gamma_th   |
#.----------------------.
# Set up the number of polytropic EOSs.
EOS_Omni::n_pieces = 1

# Set hybrid_k0 to K_ppoly_tab0
EOS_Omni::hybrid_k0 = %.15e

# Set hybrid_gamma to Gamma_ppoly_tab_in
EOS_Omni::hybrid_gamma[0] = %.15e

# Set hybrid_gamma_th to Gamma_th
EOS_Omni::hybrid_gamma_th = %.15e

#.--------------------------------.
#| End of NRPy+ generated section |
#.--------------------------------.
#
#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
"""%(tau_atmosphere,          # sets IllinoisGRMHD::tau_atm
     rho_atmosphere,          # sets IllinoisGRMHD::rho_b_atm
     rho_atmosphere,          # sets NRPyPlusTOVID::rho_atmosphere
     K_single_polytrope,      # sets IllinoisGRMHD::K_ppoly_tab0
     K_single_polytrope,      # sets NRPyPlusTOVID::K_atmosphere
     Gamma_single_polytrope,  # sets IllinoisGRMHD::Gamma_ppoly_tab_in[0]
     Gamma_single_polytrope,  # sets NRPyPlusTOVID::Gamma_atmosphere
     Gamma_single_polytrope,  # sets IllinoisGRMHD::Gamma_th
     K_single_polytrope,      # sets EOS_Omni::hybrid_k0
     Gamma_single_polytrope,  # sets EOS_Omni::hybrid_gamma[0]
     Gamma_single_polytrope)) # sets EOS_Omni::hybrid_gamma_th

    elif EOSname == "piecewise":
        if EOS_struct is None: # Use "is None" instead of "==None", as the former is more correct.
            print("Error: Please set the EOS named tuple. Usage:")
            print("generate_IllinoisGRMHD_EOS_parameter_file(\"piecewise\",outfilename,Gamma_thermal=Gamma_th,EOS_struct=eos_named_tuple)")
            sys.exit(1)

        if Gamma_thermal is None: # Use "is None" instead of "==None", as the former is more correct.
            print("Error: Please set Gamma_thermal. Usage:")
            print("generate_IllinoisGRMHD_EOS_parameter_file(\"piecewise\",outfilename,Gamma_thermal=Gamma_th,EOS_struct=eos_named_tuple)")
            sys.exit(1)

        atm_index  = polytropic_index_from_rhob(EOS_struct,rho_atmosphere)
        Gamma_atm  = EOS_struct.Gamma_poly_tab[atm_index]
        Kpoly_atm  = EOS_struct.K_poly_tab[atm_index]
        IDfilename = "outputTOVpolytrope-"+EOSname+".txt"

        with open(outfilename,"a") as file:
            file.write("""#
#.---------------------------------------.
#| EOS Type: Generic Piecewise Polytrope |
#.---------------------------------------.
#| Required parameters:                  |
#|  - EOS_struct                         |
#|  - Gamma_thermal                      |
#|  - tau_atmosphere                     |
#|  - rho_atmosphere                     |
#|                                       |
#| IllinoisGRMHD parameters set:         |
#|  - neos                               |
#|  - K_ppoly_tab0                       |
#|  - rho_ppoly_tab_in[j]   0<=j<=neos-2 |
#|  - Gamma_ppoly_tab_in[j] 0<=j<=neos-1 |
#|  - Gamma_th                           |
#|  - tau_atm                            |
#|  - rho_b_atm                          |
#.---------------------------------------.
#| NRPyPlusTOVID parameters set:         |
#|  - rho_atmosphere                     |
#|  - Gamma_atmosphere                   |
#|  - K_atmosphere                       |
#.---------------------------------------.
#| EOS_Omni parameters set:              |
#|  - n_pieces                           |
#|  - hybrid_k0                          |
#|  - hybrid_rho[j]   0<=j<=neos-2       |
#|  - hybrid_gamma[j] 0<=j<=neos-1       |
#|  - hybrid_gamma_th                    |
#.---------------------------------------.
#
# Set up initial data file name
NRPyPlusTOVID::TOV_filename = \"%s\"

# Set up the number of polytropic EOSs.
IllinoisGRMHD::neos = %d

# Set atmospheric value of tau
IllinoisGRMHD::tau_atm = %.15e

# Set K_ppoly_tab0 and K_atmosphere
IllinoisGRMHD::K_ppoly_tab0 = %.15e
NRPyPlusTOVID::K_atmosphere = %.15e

# Set atmospheric value of rho
IllinoisGRMHD::rho_b_atm      = %.15e
NRPyPlusTOVID::rho_atmosphere = %.15e

# Set rho_ppoly_tab_in""" %(IDfilename,EOS_struct.neos,tau_atmosphere,EOS_struct.K_poly_tab[0],Kpoly_atm,rho_atmosphere,rho_atmosphere))
            for j in range(EOS_struct.neos-1):
                file.write("""
IllinoisGRMHD::rho_ppoly_tab_in[%d] = %.15e""" %(j,EOS_struct.rho_poly_tab[j]))
            file.write("""

# Set Gamma_atmosphere and Gamma_ppoly_tab_in
NRPyPlusTOVID::Gamma_atmosphere      = %.15e""" %(Gamma_atm))
            for j in range(EOS_struct.neos):
                file.write("""
IllinoisGRMHD::Gamma_ppoly_tab_in[%d] = %.15e""" %(j,EOS_struct.Gamma_poly_tab[j]))
            file.write("""

# Set Gamma_th
IllinoisGRMHD::Gamma_th = %.15e

#.---------------------------------.
#| EOS_Omni parameters:            |
#|  - n_pieces                     |
#|  - hybrid_k0                    |
#|  - hybrid_rho[j]   0<=j<=neos-2 |
#|  - hybrid_gamma[j] 0<=j<=neos-1 |
#|  - hybrid_gamma_th              |
#.---------------------------------.
# Set up the number of polytropic EOSs.
EOS_Omni::n_pieces = %d

# Set hybrid_k0 to K_ppoly_tab0
EOS_Omni::hybrid_k0 = %.15e

# Set hybrid_rho to rho_ppoly_tab_in""" %(Gamma_thermal,EOS_struct.neos,EOS_struct.K_poly_tab[0]))
            for j in range(EOS_struct.neos-1):
                file.write("""
EOS_Omni::hybrid_rho[%d] = %.15e""" %(j,EOS_struct.rho_poly_tab[j]))
            file.write("""

# Set hybrid_gamma to Gamma_ppoly_tab_in""")
            for j in range(EOS_struct.neos):
                file.write("""
EOS_Omni::hybrid_gamma[%d] = %.15e""" %(j,EOS_struct.Gamma_poly_tab[j]))
            file.write("""

# Set hybrid_gamma_th to Gamma_th
EOS_Omni::hybrid_gamma_th = %.15e

#.--------------------------------.
#| End of NRPy+ generated section |
#.--------------------------------.
#
#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^""" %(Gamma_thermal))

    else:
        import TOV.Piecewise_Polytrope__dict
        if EOSname not in TOV.Piecewise_Polytrope__dict.EOS_Read_et_al_dict:
            print("ERROR: Unknown EOS name "+EOSname)
            sys.exit(1)

        if Gamma_thermal is None: # Use "is None" instead of "==None", as the former is more correct.
            print("Error: Please set Gamma_thermal. Usage:")
            print("generate_IllinoisGRMHD_EOS_parameter_file(EOSname,outfilename,Gamma_thermal=None)")
            sys.exit(1)

        eos        = set_up_EOS_parameters__Read_et_al_input_variables(EOSname)
        atm_index  = polytropic_index_from_rhob(eos,rho_atmosphere)
        Gamma_atm  = EOS_struct.Gamma_poly_tab[atm_index]
        Kpoly_atm  = EOS_struct.K_poly_tab[atm_index]
        IDfilename = "outputTOVpolytrope-"+EOSname+".txt"

        # This is done for cosmetic purposes, so that parameter files
        # of different EOS names all look the same.
        largest_name_in_EOS_table = 6
        if len(EOSname)==largest_name_in_EOS_table:
            pass
        else:
            for _k in range(largest_name_in_EOS_table - len(EOSname)): # _k is unused.
                EOSname += " "

        with open(outfilename,"a") as file:
            file.write("""#
#.---------------------------------------.
#|     EOS Type: Piecewise Polytrope     |
#.---------------------------------------.
#|           EOS name: """+EOSname+"""            |
#.---------------------------------------.
#|  Reference: Table II and III in       |
#|    Read et al. PRD 79,124032 (2009)   |
#|  https://arxiv.org/pdf/0812.2163.pdf  |
#.---------------------------------------.
#| Note that while we use the values in  |
#| Read et al. (2009), we write them in  |
#| geometrized units where G = 1 = c. We |
#| also normalize the nuclear density to |
#| unity.                                |
#| You can read more about this in the   |
#| following NRPy+ tutorial module:      |
#| Tutorial-TOV-Piecewise_Polytrope_EOSs |
#.---------------------------------------.
#| Required inputs:                      |
#|  - EOS name                           |
#|  - Gamma_thermal                      |
#.---------------------------------------.
#| IllinoisGRMHD parameters:             |
#|  - neos                               |
#|  - K_ppoly_tab0                       |
#|  - rho_ppoly_tab_in[j]   0<=j<=neos-2 |
#|  - Gamma_ppoly_tab_in[j] 0<=j<=neos-1 |
#|  - Gamma_th                           |
#|  - tau_atm                            |
#|  - rho_b_atm                          |
#.---------------------------------------.
# Set up the number of polytropic EOSs.
IllinoisGRMHD::neos = %d

# Set atmospheric value of tau
IllinoisGRMHD::tau_atm = %.15e

# Set K_ppoly_tab0
IllinoisGRMHD::K_ppoly_tab0 = %.15e

# Set atmospheric value of rho
IllinoisGRMHD::rho_b_atm = %.15e

# Set rho_ppoly_tab_in""" %(EOS_struct.neos,tau_atmosphere,EOS_struct.K_poly_tab[0],rho_atmosphere))
            for j in range(EOS_struct.neos-1):
                file.write("""
IllinoisGRMHD::rho_ppoly_tab_in[%d] = %.15e""" %(j,EOS_struct.rho_poly_tab[j]))
            file.write("""

# Set Gamma_ppoly_tab_in""")
            for j in range(EOS_struct.neos):
                file.write("""
IllinoisGRMHD::Gamma_ppoly_tab_in[%d] = %.15e""" %(j,EOS_struct.Gamma_poly_tab[j]))
            file.write("""

# Set Gamma_th
IllinoisGRMHD::Gamma_th = %.15e

#.---------------------------.
#| NRPyPlusTOVID parameters: |
#|  - TOV_filename           |
#|  - rho_atmosphere         |
#|  - Gamma_atmosphere       |
#|  - K_atmosphere           |
#.---------------------------.
# Set up initial data file name
NRPyPlusTOVID::TOV_filename = \"%s\"

# Set atmospheric value of rho
NRPyPlusTOVID::rho_atmosphere = %.15e

# Set Gamma_atmosphere
NRPyPlusTOVID::Gamma_atmosphere = %.15e

# Set K_atmosphere
NRPyPlusTOVID::K_atmosphere = %.15e

#.---------------------------------.
#| EOS_Omni parameters:            |
#|  - n_pieces                     |
#|  - hybrid_k0                    |
#|  - hybrid_rho[j]   0<=j<=neos-2 |
#|  - hybrid_gamma[j] 0<=j<=neos-1 |
#|  - hybrid_gamma_th              |
#.---------------------------------.
# Set up the number of polytropic EOSs.
EOS_Omni::n_pieces = %d

# Set hybrid_k0 to K_ppoly_tab0
EOS_Omni::hybrid_k0 = %.15e

# Set hybrid_rho to rho_ppoly_tab_in""" %(Gamma_thermal,IDfilename,rho_atmosphere,Gamma_atm,Kpoly_atm,EOS_struct.neos,EOS_struct.K_poly_tab[0]))
            for j in range(EOS_struct.neos-1):
                file.write("""
EOS_Omni::hybrid_rho[%d] = %.15e""" %(j,EOS_struct.rho_poly_tab[j]))
            file.write("""

# Set hybrid_gamma to Gamma_ppoly_tab_in""")
            for j in range(EOS_struct.neos):
                file.write("""
EOS_Omni::hybrid_gamma[%d] = %.15e""" %(j,EOS_struct.Gamma_poly_tab[j]))
            file.write("""

# Set hybrid_gamma_th to Gamma_th
EOS_Omni::hybrid_gamma_th = %.15e

#.--------------------------------.
#| End of NRPy+ generated section |
#.--------------------------------.
#
#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^""" %(Gamma_thermal))
