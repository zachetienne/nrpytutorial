# Author: Terrence Pierre Jacques
#         terrencepierrej **at** gmail **dot* com

# License: BSD 2-Clause

# COMPLETE DOCUMENTATION (JUPYTER NOTEBOOKS):
# START PAGE (start here!):  ../NRPy+_Tutorial.ipynb
# THIS MODULE: ../Tutorial-VacuumMaxwell_Flat_ID.ipynb

# Step P1: Import needed NRPy+ core modules:
import NRPy_param_funcs as par                # NRPy+: Parameter interface
import sympy as sp                            # SymPy: The Python computer algebra package upon which NRPy+ depends
import indexedexp as ixp                      # NRPy+: Symbolic indexed expression (e.g., tensors, vectors, etc.) support
import reference_metric as rfm                # NRPy+: Reference metric support
import sys                                    # Standard Python module for multiplatform OS-level functions
from Maxwell.CommonParams import amp, lam, time # NRPy+: Common parameters for all VacuumMaxwell modules (defines amp, lam, and time)


# The name of this module ("InitialData") is given by __name__:
thismodule = __name__

# define parameter for which system to use
par.initialize_param(par.glb_param("char",thismodule,"System_to_use","System_I"))

# Initial data for toroidal dipole field
def Toroidal():
    system = par.parval_from_str(thismodule+"::System_to_use")
    DIM = par.parval_from_str("grid::DIM")

    dst_basis = par.parval_from_str("reference_metric::CoordSystem")

    # Set coordinate system to Cartesian
    par.set_parval_from_str("reference_metric::CoordSystem","Cartesian")
    rfm.reference_metric()

    global AidU, EidU, psi_ID

    x = rfm.xx_to_Cart[0]
    y = rfm.xx_to_Cart[1]
    z = rfm.xx_to_Cart[2]

    AidD_Sph = ixp.zerorank1()
    # Set coordinate transformations:
    r = sp.sqrt(x*x + y*y + z*z)
    sin_theta = z / r

    u = time + r
    v = time - r
    e_lam_u = sp.exp(-lam*u**2)
    e_lam_v = sp.exp(-lam*v**2)

    # Equation 16 from https://arxiv.org/abs/gr-qc/0201051
    AD_phi_hat = (amp*sin_theta)*( ((e_lam_v - e_lam_u)/r**2) - \
                            2*lam*(v*e_lam_v + u*e_lam_u)/r )

    AidD_Sph[2] = AD_phi_hat/(r*sin_theta)

    # Coordinate transformation from spherical to Cartesian
    AidU_Cart = ixp.zerorank1()
    Jac_dxSphU_dxCartD = ixp.zerorank2()
    for i in range(DIM):
        for j in range(DIM):
            Jac_dxSphU_dxCartD[i][j] = sp.diff(rfm.xxSph[i],rfm.xx_to_Cart[j])

    #         Jac_dxCartU_dxSphD[i][j] = sp.diff(rfm.xx_to_Cart[i],rfm.xx[j])
    Jac_dxCartU_dxSphD,dummy = ixp.generic_matrix_inverter3x3(Jac_dxSphU_dxCartD)

    for i in range(DIM):
        for j in range(DIM):
            AidU_Cart[i] += Jac_dxCartU_dxSphD[i][j]*AidD_Sph[j]
    for i in range(DIM):
        AidU_Cart[i] = sp.simplify(AidU_Cart[i])

    # rfm is still defined in Cartesian coordinates
    cart_xx = ixp.declarerank1("cart_xx")
    for i in range(3):
        for k in range(3):
            AidU_Cart[i] = AidU_Cart[i].subs(rfm.xx[k], cart_xx[k])

    # Set coordinate system to dst_basis
    par.set_parval_from_str("reference_metric::CoordSystem", dst_basis)
    rfm.reference_metric()

    for i in range(3):
        for k in range(3):
            AidU_Cart[i] = AidU_Cart[i].subs(cart_xx[k], rfm.xx_to_Cart[k])

#     if radial_like_dst_xx0:
#         for j in range(3):
#             AidU_Cart[j] =  sp.refine(sp.simplify(AidU_Cart[j]), sp.Q.positive(rfm.xx[0]))

        # Step 3: Transform BSSN tensors in Cartesian basis to destination grid basis, using center of dest. grid as origin

    # Step 3.a: Next construct Jacobian and inverse Jacobian matrices:
    # _Jac_dUCart_dDrfmUD is unused.
    _Jac_dUCart_dDrfmUD, Jac_dUrfm_dDCartUD = rfm.compute_Jacobian_and_inverseJacobian_tofrom_Cartesian()

    # Step 3.b: Convert basis of all BSSN *vectors* from Cartesian to destination basis
    AidU = rfm.basis_transform_vectorU_from_Cartesian_to_rfmbasis(Jac_dUrfm_dDCartUD, AidU_Cart)

    # Define electric field --> E^i = -\partial_t A^i
    EidU = ixp.zerorank1()
    for j in range(DIM):
        EidU[j] = -sp.diff(AidU[j], time)

    psi_ID = sp.sympify(0)

    if system == "System_II":
        global Gamma_ID
        Gamma_ID = sp.sympify(0)
        print('Currently using ' + system + ' initial data')
    elif system == "System_I":
        print('Currently using ' + system + ' initial data')
    else:
        print("Invalid choice of system: System_to_use must be either System_I or System_II")
        sys.exit(1)

# Initial data driver routine
def InitialData(Type="ToroidalDipole"):
    if Type=="ToroidalDipole":
        Toroidal()
    else:
        print("Error: Vacuum Maxwell initial data Type="+str(Type)+" not supported.")
        sys.exit(1)
