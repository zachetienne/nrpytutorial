# Load necessary Python modules
import os,sys                                   # Multiplatform OS specific functions
from bisect import bisect_left                  # Bisection algorithm
import numpy as np                              # Support for large, multi-dimensional arrays and matrices and mathematical functions
import h5py as h5                               # HDF5 file system support
from collections import namedtuple              # C-like struct functionality
from scipy.interpolate import interp1d,interp2d # 1d and 2d interpolating functions
import bisect                                   # Bisection algorithms
import astropy.constants as constants

# Bisection algorithm
# Courtesy https://stackoverflow.com/questions/6065697/python-numpy-quickly-find-the-index-in-an-array-closest-to-some-value
def bisection_index_finder(t, x):
    i = bisect_left(t, x)
    if(i == len(t)):
        i-=1

    if t[i] - x < x - t[i]:
        i-=1
    return i

# Lagrange interpolation algorithm
def eos_table_1D_interpolator(interp_stencil_size,y_values,x_values, x_star):

    if x_star < x_values[0]:
        return 0.0

    elif x_star > x_values[-1]:
        return y_values[-1]

    # Find the index closest to x_star
    idx = bisection_index_finder(x_values, x_star)

    # If the index is too close to the beginning
    # of the table, then shift it accordingly
    while(idx - int(interp_stencil_size/2) < 0):
        idx += 1

    # If the index is too close to the end
    # of the table, then shift it accordingly
    while(idx + int(interp_stencil_size/2) >= len(y_values)):
        idx -= 1

    # Allocate memory for the points used by the interpolator
    interp_x = np.zeros(interp_stencil_size)
    interp_y = np.zeros(interp_stencil_size)

    # Find the values of rho and P that
    # will be used in the interpolation
    for i in range(int(interp_stencil_size/2)):

        idx1 = int(i)
        idx2 = int(i + interp_stencil_size/2)

        interp_x[idx1] = x_values[idx+i-int(interp_stencil_size/2)]
        interp_x[idx2] = x_values[idx+i  ]

        interp_y[idx1] = y_values[idx+i-int(interp_stencil_size/2)]
        interp_y[idx2] = y_values[idx+i  ]

    l_i_of_x = np.zeros(interp_stencil_size)
    # Construct the Lagrange polynomial
    for i in range(interp_stencil_size):
        numer = 1.0
        denom = 1.0
        for j in range(i):
            numer *= x_star - interp_x[j]
            denom *= interp_x[i] - interp_x[j]

        for j in range(i+1,interp_stencil_size):
            numer *= x_star - interp_x[j]
            denom *= interp_x[i] - interp_x[j]

        l_i_of_x[i] = numer/denom

    # Interpolate to y_star
    y_star = 0.0
    for i in range(interp_stencil_size):
        y_star += l_i_of_x[i] * interp_y[i]

    return y_star

def initialize_eos_struct(eos_table_path=os.path.join('IllinoisGRMHD','doc','eostable.h5'),idx_Ye=0,idx_T=0):
    # Unit conversion: cgs to dimensionless
    # Set physical constants
    G    = constants.G.cgs.value     # cm^3/g/s^2
    c    = constants.c.cgs.value     # cm/s
    Msun = constants.M_sun.cgs.value # g

    # Set the length factor:
    # .----------------------.
    # | [Msun * G / c^2] = L |
    # .----------------------.
    units_of_length = Msun * G / c**2

    # Set the time factor:
    # .----------------------.
    # | [Msun * G / c^3] = T |
    # .----------------------.
    units_of_time = units_of_length/c

    # Pressure has units of M * L^(-1) * T^(-2), so
    units_of_pressure = Msun * units_of_length**(-1) * units_of_time**(-2)

    # Density has units of M * L^(-3), so
    units_of_density = Msun * units_of_length**(-3)

    # Specific internal energy has units of L^(2) * T^(-2),so
    units_of_spec_int_energy = c**2

    # Load the relevant quantities from the EOS table
    eos_table           = h5.File(eos_table_path,'r')
    interp_stencil_size = 8
    tab_P               = 10**np.array(eos_table["logpress" ][:])
    tab_rho             = 10**np.array(eos_table["logrho"   ][:])
    tab_eps             = 10**np.array(eos_table["logenergy"][:])
    tab_P_of_rho        = tab_P[idx_Ye,idx_T,:]
    tab_eps_of_rho      = tab_eps[idx_Ye,idx_T,:]

    # Define the eos_struct named tuple
    eos_struct = namedtuple("eos_struct","tab_P_of_rho tab_eps_of_rho tab_rho")

    # Set the eos "struct" with dimensionless quantities
    eos = eos_struct(tab_P_of_rho/units_of_pressure,
                     tab_eps_of_rho/units_of_spec_int_energy,
                     tab_rho/units_of_density)

    return eos

def Read_EOS_table_and_sliced_table(eos_table_path,sliced_table_path):
    # Read in the EOS table
    eos_table = h5.File(eos_table_path,'r')
    logrho    = np.array(eos_table["logrho"])

    # Read in the sliced table values
    tempindex, rhoindex, yeindex, P, S, dummy, dummy, dummy, dummy, eps = np.loadtxt(os.path.join(sliced_table_path,"eos.thermo")).T
    yeb = np.loadtxt(os.path.join(sliced_table_path,"eos.yqb")).T

    # Convert P to cgs units
    # pressconv = 17.980714384531204 # math.log10(m_n * 1.0E-3 * 1.0E-4 / mev_si)  #useful constant for pressure conversion
    # P = 10**(np.log10(P) + logrho + pressconv)

    # Set rho
    rho = 10**logrho

    # Check that the length of all quantities match
    length_check = len(logrho)
    entries = [rhoindex,yeindex,P,rho,yeb]
    for entry in entries:
        if( len(entry) != length_check ):
            print("Error! Incorrect length found!")
            sys.exit(1)
    print("Loaded EOS quantities successfully")

    # Convert from cgs to dimensionless units
    P,rho,eps = Convert_table_to_code_units(P,rho,eps)

    print("Converted P,rho,eps from cgs to dimensionless units")

    # Set the eos namedtuple
    eos_tuple = namedtuple("eos_tuple","rho P yeindex yeb eps")
    eos = eos_tuple(rho,P,yeindex,yeb,eps)

    print("Finished without errors")

    return eos

def read_and_slice_EOS_table_at_given_temperature(eos_file_path,T_in_MeV):
    eos_file = h5.File(eos_file_path,"r")

    # Read in the ye, P, eps, T, and rho
    tab_Ye   =     np.array(eos_file['ye'       ])[:]
    tab_P    = 10**np.array(eos_file['logpress' ])[:]
    tab_eps  = 10**np.array(eos_file['logenergy'])[:] - eos_file['energy_shift'][0]
    tab_temp = 10**np.array(eos_file['logtemp'  ])[:]
    tab_rho  = 10**np.array(eos_file['logrho'   ])[:]
    tab_ent  =     np.array(eos_file['entropy'  ])[:]

    # Read in the chemical potentials
#     munu     =     np.array(eos_file['munu'   ])[:]
    mu_e     =     np.array(eos_file['mu_e'  ])[:]
    mu_hat   =     np.array(eos_file['muhat' ])[:]
    mu_l     = mu_e - mu_hat

    # Find temperature index
    T_idx = bisect.bisect_left(tab_temp,T_in_MeV)
    print("Slicing table at index %d. T[%d] = %.3lf MeV | Input temperature: %.3lf MeV"%(T_idx,T_idx,tab_temp[T_idx],T_in_MeV))

    # Set 2d interpolators for P and eps
    # as functions of rho and Ye
    P_interp   = interp2d(tab_rho,tab_Ye,tab_P[  :,T_idx,:])
    eps_interp = interp2d(tab_rho,tab_Ye,tab_eps[:,T_idx,:])
    ent_interp = interp2d(tab_rho,tab_Ye,tab_ent[:,T_idx,:])

    # Now, loop over the density and set the array
    Ye_be   = np.empty(len(tab_rho))
    P_be    = np.empty(len(tab_rho))
    eps_be  = np.empty(len(tab_rho))
    ent_be  = np.empty(len(tab_rho))
    Y_e_min = 1.1*tab_Ye.min()
    for rho_idx in range(len(tab_rho)):
        # Set the interpolator
#         Ye_ito_munu_of_T_and_rho = interp1d(munu[:,T_idx,rho_idx],tab_Ye)
        Ye_ito_munu_of_T_and_rho = interp1d(mu_l[:,T_idx,rho_idx],tab_Ye)
        # Then find the value of Ye which sets munu to zero
        try:
            YeL = Ye_ito_munu_of_T_and_rho(0.0)
        except:
            YeL = Y_e_min
        # Set Ye in beta equilibrium
        Ye_be[ rho_idx] = YeL
        # Finally, determine P and eps in beta equilibrium
        P_be[  rho_idx] = P_interp(  tab_rho[rho_idx],Ye_be[rho_idx])
        eps_be[rho_idx] = eps_interp(tab_rho[rho_idx],Ye_be[rho_idx])
        ent_be[rho_idx] = ent_interp(tab_rho[rho_idx],Ye_be[rho_idx])

    # Convert table to code units
    # Set physical constants in cgs units
    G    = 6.6738480e-8   # cm^3/g/s^2
    c    = 2.99792458e10  # cm/s
    Msun = 1.9884e33      # g

    # Set the length factor:
    # .----------------------.
    # | [Msun * G / c^2] = L |
    # .----------------------.
    units_of_length = Msun * G / c**2

    # Set the time factor:
    # .----------------------.
    # | [Msun * G / c^3] = T |
    # .----------------------.
    units_of_time = units_of_length/c

    # Pressure has units of M * L^(-1) * T^(-2), so
    units_of_pressure = Msun * units_of_length**(-1) * units_of_time**(-2)

    # Density has units of M * L^(-3), so
    units_of_density = Msun * units_of_length**(-3)

    # Specific internal energy has units of L^(2) * T^(-2),so
    units_of_spec_int_energy = c**2

    # Convert the table to dimensionless units
    P_be    /= units_of_pressure
    tab_rho /= units_of_density
    eps_be  /= units_of_spec_int_energy

    # Functions to interpolate P and eps
    rho_of_P   = interp1d(P_be   ,tab_rho)
    P_of_rho   = interp1d(tab_rho,P_be)
    eps_of_rho = interp1d(tab_rho,eps_be)
    Ye_of_rho  = interp1d(tab_rho,Ye_be)
    ent_of_rho = interp1d(tab_rho,ent_be)

    # Set the EOS tuple
    eos_tuple = namedtuple("eos_tuple","rho Ye T P eps Temperature rho_of_P_interpolator P_of_rho_interpolator eps_of_rho_interpolator Ye_of_rho_interpolator units_of_density units_of_pressure units_of_spec_int_energy ent_of_rho_interpolator")
    beta_equilibrated_eos_slice = eos_tuple(tab_rho,tab_Ye,tab_temp,
                                            tab_P,tab_eps,
                                            tab_temp[T_idx],
                                            rho_of_P,
                                            P_of_rho,
                                            eps_of_rho,
                                            Ye_of_rho,
                                            units_of_density,
                                            units_of_pressure,
                                            units_of_spec_int_energy,
                                            ent_of_rho)

    return beta_equilibrated_eos_slice
