# Step 0: Add NRPy's directory to the path
# https://stackoverflow.com/questions/16780014/import-file-from-parent-directory
import os,sys
nrpy_dir_path = os.path.join("..")
if nrpy_dir_path not in sys.path:
    sys.path.append(nrpy_dir_path)

# Step 0.a: Import the NRPy+ core modules and set the reference metric to Cartesian
import sympy as sp               # SymPy: The Python computer algebra package upon which NRPy+ depends
from sympy import Rational as rl
import NRPy_param_funcs as par   # NRPy+: Parameter interface
import indexedexp as ixp         # NRPy+: Symbolic indexed expression (e.g., tensors, vectors, etc.) support
import reference_metric as rfm   # NRPy+: Reference metric support
import Min_Max_and_Piecewise_Expressions as noif

thismodule = __name__

def ConstantDensitySphere(r):

    ConstantDensitySphere_sphere_radius =  par.Cparameters("REAL", thismodule, 
                                                ["ConstantDensitySphere_sphere_radius"],1.0)

    ConstantDensitySphere_rho_interior  =   par.Cparameters("REAL", thismodule, 
                                                ["ConstantDensitySphere_rho_interior"],1.586831634477515e-04) # 9.8e13 g/cm3

    ConstantDensitySphere_Y_e_interior  =   par.Cparameters("REAL", thismodule, 
                                                ["ConstantDensitySphere_Y_e_interior"],0.1)

    ConstantDensitySphere_T_interior    =   par.Cparameters("REAL", thismodule, 
                                                ["ConstantDensitySphere_T_interior"],8.0)

    ConstantDensitySphere_rho_exterior  =   par.Cparameters("REAL", thismodule, 
                                                ["ConstantDensitySphere_rho_exterior"],9.715295721290909e-11) # 6.0e07 g/cm3

    ConstantDensitySphere_Y_e_exterior  =   par.Cparameters("REAL", thismodule, 
                                                ["ConstantDensitySphere_Y_e_exterior"],0.5)

    ConstantDensitySphere_T_exterior    =   par.Cparameters("REAL", thismodule, 
                                                ["ConstantDensitySphere_T_exterior"],0.01 )


    r_cutoff = ConstantDensitySphere_sphere_radius

    rho_in  = ConstantDensitySphere_rho_interior
    rho_out = ConstantDensitySphere_rho_exterior

    rho = noif.coord_leq_bound(    r, r_cutoff)*rho_in\
         +noif.coord_greater_bound(r, r_cutoff)*rho_out


    T_in  = ConstantDensitySphere_T_interior
    T_out = ConstantDensitySphere_T_exterior

    T = noif.coord_leq_bound(    r, r_cutoff)*T_in\
       +noif.coord_greater_bound(r, r_cutoff)*T_out


    Ye_in  = ConstantDensitySphere_Y_e_interior
    Ye_out = ConstantDensitySphere_Y_e_exterior

    Ye = noif.coord_leq_bound(   r, r_cutoff)*Ye_in\
       +noif.coord_greater_bound(r, r_cutoff)*Ye_out

    return rho, T, Ye


def IsotropicGas():

    IsotropicGasID_rho         = par.Cparameters("REAL", thismodule, ["IsotropicGasID_rho"],         1e-12)
    IsotropicGasID_Y_e         = par.Cparameters("REAL", thismodule, ["IsotropicGasID_Y_e"],         0.5)
    IsotropicGasID_temperature = par.Cparameters("REAL", thismodule, ["IsotropicGasID_temperature"], 1.0)

    rho = IsotropicGasID_rho
    T   = IsotropicGasID_temperature
    Ye  = IsotropicGasID_Y_e
    return rho, T, Ye


