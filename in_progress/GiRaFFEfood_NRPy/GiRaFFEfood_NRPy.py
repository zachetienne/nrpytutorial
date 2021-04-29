# Step 0: Add NRPy's directory to the path
# https://stackoverflow.com/questions/16780014/import-file-from-parent-directory
import os,sys
nrpy_dir_path = os.path.join("..")
if nrpy_dir_path not in sys.path:
    sys.path.append(nrpy_dir_path)

# Step 0: Import the NRPy+ core modules and set the reference metric to Cartesian
import NRPy_param_funcs as par   # NRPy+: Parameter interface
import reference_metric as rfm   # NRPy+: Reference metric support
import GiRaFFEfood_NRPy.GiRaFFEfood_NRPy_Common_Functions as gfcf # Some useful functions for GiRaFFE initial data.

par.set_parval_from_str("reference_metric::CoordSystem","Cartesian")
rfm.reference_metric()

# Step 1a: Set commonly used parameters.
thismodule = __name__

def GiRaFFEfood_NRPy_generate_initial_data(ID_type = "DegenAlfvenWave", stagger_enable = False,**params):
    global AD, Valencia_vU
#     if ID_type == "DegenAlfvenWave":
#         mu_DAW = par.Cparameters("REAL",thismodule,["mu_DAW"], -0.5) # The wave speed
#         AD = Axyz_func(Ax_DAW, Ay_DAW, Az_DAW, stagger_enable,
#                        gammamu=sp.sympify(1)/sp.sqrt(sp.sympify(1)-mu_AW**2))
#         ValenciaVU = ValenciavU_DAW(mu_DAW=mu_DAW, gammamu=gammamu)
#     elif ID_type == "FastWave":
#         AD = Axyz_func(Ax_FW, Ay_FW, Az_FW, stagger_enable)
#         ValenciaVU = ValenciavU_FW()
    if ID_type == "SplitMonopole":
        import GiRaFFEfood_NRPy.GiRaFFEfood_NRPy_Split_Monopole as gfsm
        AD = gfcf.Axyz_func_spherical(gfsm.Ar_SM,gfsm.Ath_SM,gfsm.Aph_SM,stagger_enable,**params)
        Valencia_vU = gfsm.ValenciavU_func_SM(**params)