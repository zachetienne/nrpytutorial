# Step 0: Add NRPy's directory to the path
# https://stackoverflow.com/questions/16780014/import-file-from-parent-directory
import os,sys
nrpy_dir_path = os.path.join("..")
if nrpy_dir_path not in sys.path:
    sys.path.append(nrpy_dir_path)
giraffefood_dir_path = os.path.join("GiRaFFEfood_NRPy")
if giraffefood_dir_path not in sys.path:
    sys.path.append(giraffefood_dir_path)

# Step 0: Import the NRPy+ core modules and set the reference metric to Cartesian
import NRPy_param_funcs as par   # NRPy+: Parameter interface
import reference_metric as rfm   # NRPy+: Reference metric support
import GiRaFFEfood_NRPy_Common_Functions as gfcf # Some useful functions for GiRaFFE initial data.

par.set_parval_from_str("reference_metric::CoordSystem","Cartesian")
rfm.reference_metric()

# We import all ID modules ahead of time so that options can be changed *before* generating the functions.
import GiRaFFEfood_NRPy_Exact_Wald as gfew
import GiRaFFEfood_NRPy_Split_Monopole as gfsm

# Step 1a: Set commonly used parameters.
thismodule = __name__

def GiRaFFEfood_NRPy_generate_initial_data(ID_type = "DegenAlfvenWave", stagger_enable = False,**params):
    global AD, ValenciavU
#     if ID_type == "DegenAlfvenWave":
#         mu_DAW = par.Cparameters("REAL",thismodule,["mu_DAW"], -0.5) # The wave speed
#         AD = Axyz_func(Ax_DAW, Ay_DAW, Az_DAW, stagger_enable,
#                        gammamu=sp.sympify(1)/sp.sqrt(sp.sympify(1)-mu_AW**2))
#         ValenciaVU = ValenciavU_DAW(mu_DAW=mu_DAW, gammamu=gammamu)
#     elif ID_type == "FastWave":
#         AD = Axyz_func(Ax_FW, Ay_FW, Az_FW, stagger_enable)
#         ValenciaVU = ValenciavU_FW()
    if ID_type == "ExactWald":
        AD = gfcf.Axyz_func_spherical(gfew.Ar_EW,gfew.Ath_EW,gfew.Aph_EW,stagger_enable,**params)
        ValenciavU = gfew.ValenciavU_func_EW(**params)
    elif ID_type == "SplitMonopole":
        AD = gfcf.Axyz_func_spherical(gfsm.Ar_SM,gfsm.Ath_SM,gfsm.Aph_SM,stagger_enable,**params)
        ValenciavU = gfsm.ValenciavU_func_SM(**params)