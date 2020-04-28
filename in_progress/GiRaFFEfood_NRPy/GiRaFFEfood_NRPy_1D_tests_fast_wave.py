# Step 0: Add NRPy's directory to the path
# https://stackoverflow.com/questions/16780014/import-file-from-parent-directory
import os,sys
nrpy_dir_path = os.path.join("..")
if nrpy_dir_path not in sys.path:
    sys.path.append(nrpy_dir_path)

# Step 0.a: Import the NRPy+ core modules and set the reference metric to Cartesian
from outputC import *            # NRPy+: Core C code output module
import finite_difference as fin  # NRPy+: Finite difference C code generation module
import NRPy_param_funcs as par   # NRPy+: Parameter interface
import grid as gri               # NRPy+: Functions having to do with numerical grids
import indexedexp as ixp         # NRPy+: Symbolic indexed expression (e.g., tensors, vectors, etc.) support

import reference_metric as rfm   # NRPy+: Reference metric support
par.set_parval_from_str("reference_metric::CoordSystem","Cartesian")
rfm.reference_metric()

# Step 1a: Set commonly used parameters.
thismodule = "GiRaFFEfood_NRPy_1D"

def GiRaFFEfood_NRPy_1D_tests_fast_wave():
    # We'll use reference_metric.py to define x and y
    x = rfm.xxCart[0]
    y = rfm.xxCart[1]

    global AD
    AD = ixp.zerorank1(DIM=3)

    import Min_Max_and_Piecewise_Expressions as noif
    bound = sp.Rational(1,10)

    # A_x = 0, A_y = 0
    # A_z = y+ (-x-0.0075) if x <= -0.1
    #          (0.75x^2 - 0.85x) if -0.1 < x <= 0.1 
    #          (-0.7x-0.0075) if x > 0.1 

    Azleft = y - x - sp.Rational(75,10000)
    Azcenter = y + sp.Rational(75,100)*x*x - sp.Rational(85,100)*x
    Azright = y - sp.Rational(7,10)*x - sp.Rational(75,10000)

    AD[0] = sp.sympify(0)
    AD[1] = sp.sympify(0)
    AD[2] = noif.coord_leq_bound(x,-bound)*Azleft\
           +noif.coord_greater_bound(x,-bound)*noif.coord_leq_bound(x,bound)*Azcenter\
           +noif.coord_greater_bound(x,bound)*Azright
    
    # B^x(0,x) = 1.0
    # B^y(0,x) = 1.0 if x <= -0.1
    #            1.0-1.5(x+0.1) if -0.1 < x <= 0.1
    #            0.7 if x > 0.1
    # B^z(0,x) = 0

    Byleft = sp.sympify(1)
    Bycenter = sp.sympify(1) - sp.Rational(15,10)*(x+sp.Rational(1,10))
    Byright = sp.Rational(7,10)

    BU = ixp.zerorank1()
    BU[0] = sp.sympify(1)
    BU[1] = noif.coord_leq_bound(x,-bound)*Byleft\
            +noif.coord_greater_bound(x,-bound)*noif.coord_leq_bound(x,bound)*Bycenter\
            +noif.coord_greater_bound(x,bound)*Byright
    BU[2] = 0
    
    # E^x(0,x) = 0.0 , E^y(x) = 0.0 , E^z(x) = -B^y(0,x)
    EU = ixp.zerorank1()
    EU[0] = sp.sympify(0)
    EU[1] = sp.sympify(0)
    EU[2] = -BU[1]
    
    # Define Levi-Civita symbol
    def define_LeviCivitaSymbol_rank3(DIM=-1):
        if DIM == -1:
            DIM = par.parval_from_str("DIM")

        LeviCivitaSymbol = ixp.zerorank3()

        for i in range(DIM):
            for j in range(DIM):
                for k in range(DIM):
                    # From https://codegolf.stackexchange.com/questions/160359/levi-civita-symbol :
                    LeviCivitaSymbol[i][j][k] = (i - j) * (j - k) * (k - i) * sp.Rational(1,2)
        return LeviCivitaSymbol
    # Here, we import the Levi-Civita tensor and compute the tensor with lower indices
    LeviCivitaSymbolDDD = define_LeviCivitaSymbol_rank3()

    B2 = sp.sympify(0)
    for i in range(3):
        # In flat spacetime, gamma_{ij} is just a Kronecker delta
        B2 += BU[i]**2 # This is trivial to extend to curved spacetime

    # v^i = [ijk] (E^j B^k) / (B^2)
    global ValenciavU
    ValenciavU = ixp.zerorank1()
    for i in range(3):
        for j in range(3):
            for k in range(3):
                ValenciavU[i] += LeviCivitaSymbolDDD[i][j][k] * EU[j] * BU[k] / B2   

    