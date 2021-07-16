# Step 0: Import the NRPy+ core modules and set the reference metric to Cartesian
# Step 0: Add NRPy's directory to the path
# https://stackoverflow.com/questions/16780014/import-file-from-parent-directory
import os,sys
nrpy_dir_path = os.path.join("..")
if nrpy_dir_path not in sys.path:
    sys.path.append(nrpy_dir_path)
nrpy_dir_path = os.path.join("../..")
if nrpy_dir_path not in sys.path:
    sys.path.append(nrpy_dir_path)
giraffefood_dir_path = os.path.join("GiRaFFEfood_NRPy")
if giraffefood_dir_path not in sys.path:
    sys.path.append(giraffefood_dir_path)

# Step 0.a: Import the NRPy+ core modules and set the reference metric to Cartesian
import sympy as sp               # SymPy: The Python computer algebra package upon which NRPy+ depends
import NRPy_param_funcs as par   # NRPy+: Parameter interface
import indexedexp as ixp         # NRPy+: Symbolic indexed expression (e.g., tensors, vectors, etc.) support
import GiRaFFEfood_NRPy_Common_Functions as gfcf # Some useful functions for GiRaFFE initial data.

import reference_metric as rfm
par.set_parval_from_str("reference_metric::CoordSystem","Cartesian")
rfm.reference_metric()

# Step 1a: Set commonly used parameters.
thismodule = __name__

# Step 2: Set the vectors A and E in Spherical coordinates
def Ar_EW(r,theta,phi, **params):
    return sp.sympify(0)

def Ath_EW(r,theta,phi, **params):
    return sp.sympify(0)

def Aph_EW(r,theta,phi, **params):
    # 1/2 r^2 \sin^2 \theta
    return (r * r * sp.sin(theta)**2)/2

#Step 3: Compute v^i from A_i and E_i
def ValenciavU_func_EW(**params):
    M = params["M"]
    gammaDD = params["gammaDD"] # Note that this must use a Cartesian basis!
    sqrtgammaDET = params["sqrtgammaDET"]
    KerrSchild_radial_shift = params["KerrSchild_radial_shift"]
    r     = rfm.xxSph[0] + KerrSchild_radial_shift # We are setting the data up in Shifted Kerr-Schild coordinates
    theta = rfm.xxSph[1]

    LeviCivitaTensorUUU = ixp.LeviCivitaTensorUUU_dim3_rank3(sqrtgammaDET)

    AD = gfcf.Axyz_func_spherical(Ar_EW,Ath_EW,Aph_EW,False,**params)
    # For the initial data, we can analytically take the derivatives of A_i
    ADdD = ixp.zerorank2()
    for i in range(3):
        for j in range(3):
            ADdD[i][j] = sp.simplify(sp.diff(AD[i],rfm.xx_to_Cart[j]))

    BU = ixp.zerorank1()
    for i in range(3):
        for j in range(3):
            for k in range(3):
                BU[i] += LeviCivitaTensorUUU[i][j][k] * ADdD[k][j]

    EsphD = ixp.zerorank1()
    # 2 M ( 1+ 2M/r )^{-1/2} \sin^2 \theta
    EsphD[2] = 2 * M * sp.sin(theta)**2 / sp.sqrt(1+2*M/r)

    ED = gfcf.change_basis_spherical_to_Cartesian_D(EsphD)

    return gfcf.compute_ValenciavU_from_ED_and_BU(ED, BU, gammaDD)
