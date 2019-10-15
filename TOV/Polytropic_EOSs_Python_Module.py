## Polytropic EOSs Python Module
## Author(s): Leo Werneck and Zach Etienne
## In this NRPy+ module we set up useful "lowlevel" functions that compute
## useful polytropic quantities.

# Full documentation for this module may be found in the NRPy+ tutorial Jupyter notebook:
#  Tutorial-TOV_Polytropic_EOSs_Python_Module.ipynb

# Step 0: Import needed Python/NRPy+ modules
import numpy as np                  # This module is used for math functions
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

def set_up_EOS_parameters__Read_et_al_input_variables(log_of_P4,Gamma_4,Gamma_5,Gamma_6):

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
    # Set up the speed of light and change the units of the input pressure
    log_of_P4     -= 2.0*np.log10(2.9979e10)
    
    # Set up tabulated polytropic values following the table above
    # and the user input. All quantities which are still unknown are
    # set to absurd values to make sure they are overwritten
    rho_poly_tab   = [2.440340e+07, 3.78358e+11, 2.62780e+12, -1e30  , 10**(14.7)     , 10**(15.0)]
    P_poly_tab     = [-1e30       , -1e30      , -1e30      , -1e30  , 10**(log_of_P4), -1e30     ]
    Gamma_poly_tab = [1.58425     , 1.28733    , 0.62223    , 1.35692, Gamma_4        , Gamma_5, Gamma_6]
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
    
    # Create the EOS "struct" (named tuple)
    eos_struct = namedtuple("eos_struct","neos rho_poly_tab Gamma_poly_tab K_poly_tab P_poly_tab")
    
    return eos_struct(neos,rho_poly_tab,Gamma_poly_tab,K_poly_tab,P_poly_tab)

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

