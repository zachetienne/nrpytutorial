# Step 0: Import the NRPy+ core modules and set the reference metric to Cartesian
import NRPy_param_funcs as par
import grid as gri               # NRPy+: Functions having to do with numerical grids
import indexedexp as ixp
import sympy as sp               # SymPy: The Python computer algebra package upon which NRPy+ depends
import reference_metric as rfm
par.set_parval_from_str("reference_metric::CoordSystem","Cartesian")
rfm.reference_metric()

# Step 1a: Set commonly used parameters.
thismodule = __name__

def GiRaFFEfood_NRPy_generate_initial_data(ID_type = "DegenAlfvenWave", stagger_enable = False):
    global AD, ValenciaVU
#     if ID_type == "DegenAlfvenWave":
#         mu_DAW = par.Cparameters("REAL",thismodule,["mu_DAW"], -0.5) # The wave speed
#         AD = Axyz_func(Ax_DAW, Ay_DAW, Az_DAW, stagger_enable,
#                        gammamu=sp.sympify(1)/sp.sqrt(sp.sympify(1)-mu_AW**2))
#         ValenciaVU = ValenciavU_DAW(mu_DAW=mu_DAW, gammamu=gammamu)
#     elif ID_type == "FastWave":
#         AD = Axyz_func(Ax_FW, Ay_FW, Az_FW, stagger_enable)
#         ValenciaVU = ValenciavU_FW()
    if ID_type == "SplitMonopole":
        import GiRaFFEfood_NRPy_Split_Monopole
        AD = Axyz_func(Ax_SM,Ay_SM,Az_SM,stagger_enable)