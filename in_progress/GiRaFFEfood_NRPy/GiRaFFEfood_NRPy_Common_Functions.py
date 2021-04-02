# Step 0: Add NRPy's directory to the path
# https://stackoverflow.com/questions/16780014/import-file-from-parent-directory
import os,sys
nrpy_dir_path = os.path.join("..")
if nrpy_dir_path not in sys.path:
    sys.path.append(nrpy_dir_path)

# Step 0: Import the NRPy+ core modules and set the reference metric to Cartesian
import NRPy_param_funcs as par   # NRPy+: Parameter interface
import grid as gri               # NRPy+: Functions having to do with numerical grids
import indexedexp as ixp         # NRPy+: Symbolic indexed expression (e.g., tensors, vectors, etc.) support
import sympy as sp               # SymPy: The Python computer algebra package upon which NRPy+ depends
import reference_metric as rfm   # NRPy+: Reference metric support
par.set_parval_from_str("reference_metric::CoordSystem","Cartesian")
rfm.reference_metric()

# Step 1a: Set commonly used parameters.
thismodule = __name__

# Generic function for all 1D tests: Compute Ax,Ay,Az
def Axyz_func_Cartesian(Ax_func,Ay_func,Az_func, stagger_enable, **params):
    x = rfm.xx_to_Cart[0]
    y = rfm.xx_to_Cart[1]
    z = rfm.xx_to_Cart[2]
    AD = ixp.zerorank1()
    # First Ax
    if stagger_enable:
        y += sp.Rational(1,2)*gri.dxx[1]
        z += sp.Rational(1,2)*gri.dxx[2]
    AD[0] = Ax_func(x,y,z, **params)
    # Then Ay
    if stagger_enable:
        x += sp.Rational(1,2)*gri.dxx[0]
        y -= sp.Rational(1,2)*gri.dxx[1]
        z += sp.Rational(1,2)*gri.dxx[2]
    AD[1] = Ay_func(x,y,z, **params)
    # Finally Az
    if stagger_enable:
        x += sp.Rational(1,2)*gri.dxx[0]
        y += sp.Rational(1,2)*gri.dxx[1]
        z -= sp.Rational(1,2)*gri.dxx[2]
    AD[2] = Az_func(x,y,z, **params)

    return AD

# Generic function to convert contravariant vectors from a spherical to Cartesian basis.
def change_basis_spherical_to_Cartesian(somevector_sphD):
    # Use the Jacobian matrix to transform the vectors to Cartesian coordinates.
    drrefmetric__dx_0UDmatrix = sp.Matrix([[sp.diff(rfm.xxSph[0],rfm.xx[0]), sp.diff(rfm.xxSph[0],rfm.xx[1]), sp.diff(rfm.xxSph[0],rfm.xx[2])],
                                           [sp.diff(rfm.xxSph[1],rfm.xx[0]), sp.diff(rfm.xxSph[1],rfm.xx[1]), sp.diff(rfm.xxSph[1],rfm.xx[2])],
                                           [sp.diff(rfm.xxSph[2],rfm.xx[0]), sp.diff(rfm.xxSph[2],rfm.xx[1]), sp.diff(rfm.xxSph[2],rfm.xx[2])]])

    somevectorD = ixp.zerorank1()

    for i in range(3):
        for j in range(3):
            somevectorD[i] = drrefmetric__dx_0UDmatrix[(j,i)]*somevector_sphD[j]

    return somevectorD

# Generic function for all 1D tests: Compute Ax,Ay,Az
def Axyz_func_spherical(Ar_func,At_func,Ap_func, stagger_enable, **params):
    KerrSchild_radial_shift = params["KerrSchild_radial_shift"]
    r     = rfm.xxSph[0] + KerrSchild_radial_shift # We are setting the data up in Shifted Kerr-Schild coordinates
    theta = rfm.xxSph[1]
    phi   = rfm.xxSph[2]
    AsphD = ixp.zerorank1()
    if stagger_enable:
        # First Ax
        AsphD[0] = Ar_func(r,theta,phi, **params)
        # Then Ay
        AsphD[1] = At_func(r,theta,phi, **params)
        # Finally Az
        AsphD[2] = Ap_func(r,theta,phi, **params)
    else:
        # First Ax
        AsphD[0] = Ar_func(r.subs(rfm.xx[1],rfm.xx[1]+sp.Rational(1,2)*gri.dxx[1]).subs(rfm.xx[2],rfm.xx[2]+sp.Rational(1,2)*gri.dxx[2]),
                           theta.subs(rfm.xx[1],rfm.xx[1]+sp.Rational(1,2)*gri.dxx[1]).subs(rfm.xx[2],rfm.xx[2]+sp.Rational(1,2)*gri.dxx[2]),
                           phi.subs(rfm.xx[1],rfm.xx[1]+sp.Rational(1,2)*gri.dxx[1]).subs(rfm.xx[2],rfm.xx[2]+sp.Rational(1,2)*gri.dxx[2]), **params)
        # Then Ay
        AsphD[1] = At_func(r.subs(rfm.xx[0],rfm.xx[0]+sp.Rational(1,2)*gri.dxx[0]).subs(rfm.xx[2],rfm.xx[2]+sp.Rational(1,2)*gri.dxx[2]),
                           theta.subs(rfm.xx[0],rfm.xx[0]+sp.Rational(1,2)*gri.dxx[0]).subs(rfm.xx[2],rfm.xx[2]+sp.Rational(1,2)*gri.dxx[2]),
                           phi.subs(rfm.xx[0],rfm.xx[0]+sp.Rational(1,2)*gri.dxx[0]).subs(rfm.xx[2],rfm.xx[2]+sp.Rational(1,2)*gri.dxx[2]), **params)
        # Finally Az
        AsphD[2] = Ap_func(r.subs(rfm.xx[0],rfm.xx[0]+sp.Rational(1,2)*gri.dxx[0]).subs(rfm.xx[1],rfm.xx[1]+sp.Rational(1,2)*gri.dxx[1]),
                           theta.subs(rfm.xx[0],rfm.xx[0]+sp.Rational(1,2)*gri.dxx[0]).subs(rfm.xx[1],rfm.xx[1]+sp.Rational(1,2)*gri.dxx[1]),
                           phi.subs(rfm.xx[0],rfm.xx[0]+sp.Rational(1,2)*gri.dxx[0]).subs(rfm.xx[1],rfm.xx[1]+sp.Rational(1,2)*gri.dxx[1]), **params)

    # Use the Jacobian matrix to transform the vectors to Cartesian coordinates.
    AD = change_basis_spherical_to_Cartesian(AsphD)
    return AD

# Generic function for all 1D tests: Valencia 3-velocity from EU and BU
def compute_ValenciavU_from_EU_and_BU(EU, BU, gammaDD=None):
    # Now, we calculate v^i = ([ijk] E^j B^k) / B^2,
    # where [ijk] is the Levi-Civita symbol and B^2 = \gamma_{ij} B^i B^j$ is a trivial dot product in flat space.
    LeviCivitaSymbolDDD = ixp.LeviCivitaSymbol_dim3_rank3()

    B2 = sp.sympify(0)
    # In flat spacetime, use the Minkowski metric; otherwise, use the input metric.
    if gammaDD==None:
        gammaDD = ixp.zerorank2()
        for i in range(3):
            gammaDD[i][i] = sp.sympify(1)
    for i in range(3):
        for j in range(3):
            B2 += gammaDD[i][j] * BU[i] * BU[j]

    ValenciavU = ixp.zerorank1()
    for i in range(3):
        for j in range(3):
            for k in range(3):
                ValenciavU[i] += LeviCivitaSymbolDDD[i][j][k] * EU[j] * BU[k] / B2

    return ValenciavU