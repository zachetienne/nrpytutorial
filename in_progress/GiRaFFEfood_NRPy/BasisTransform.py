import os, sys           # Standard Python modules for multiplatform OS-level functions
# First, we'll add the parent directory to the list of directories Python will check for modules.
nrpy_dir_path = os.path.join("..")
if nrpy_dir_path not in sys.path:
    sys.path.append(nrpy_dir_path)

# Import needed Python modules
import NRPy_param_funcs as par   # NRPy+: Parameter interface
import indexedexp as ixp         # NRPy+: Symbolic indexed expression (e.g., tensors, vectors, etc.) support
import reference_metric as rfm   # NRPy+: Reference metric support
import sympy as sp                # SymPy: The Python computer algebra package upon which NRPy+ depends

#Step 0: Set the spatial dimension parameter to 3.
par.set_parval_from_str("grid::DIM", 3)
DIM = par.parval_from_str("grid::DIM")

def basis_transform(CoordSystem, AD, ValenciavU, BU=None):
    global AD_dst, ValenciavU_dst, BU_dst

    # Set coordinate system to dst_basis
    par.set_parval_from_str("reference_metric::CoordSystem",CoordSystem)
    rfm.reference_metric()

#         We define Jacobians relative to the center of the destination grid, at a point $x^j_{\rm dst}=$(`xx0,xx1,xx2`)${}_{\rm dst}$ on the destination grid:
#     $$
#     {\rm Jac\_dUCart\_dDdstUD[i][j]} = \frac{\partial x^i_{\rm Cart}}{\partial x^j_{\rm dst}},
#     $$

#     via exact differentiation (courtesy SymPy), and the inverse Jacobian
#     $$
#     {\rm Jac\_dUdst\_dDCartUD[i][j]} = \frac{\partial x^i_{\rm dst}}{\partial x^j_{\rm Cart}},
#     $$

#     using NRPy+'s `generic_matrix_inverter3x3()` function. In terms of these, the transformation of BSSN tensors from Cartesian to the destination grid's `"reference_metric::CoordSystem"` coordinates may be written:

#     $$
#     B^i_{\rm dst} = \frac{\partial x^i_{\rm dst}}{\partial x^\ell_{\rm Cart}} B^\ell_{\rm Cart},
#     $$

#     while for lowered indices we have

#     $$
#     A^{\rm dst}_{i} =
#     \frac{\partial x^\ell_{\rm Cart}}{\partial x^i_{\rm dst}} A^{\rm Cart}_{\ell}\\
#     $$

    # Step 3.a: Next construct Jacobian and inverse Jacobian matrices:
    Jac_dUCart_dDrfmUD,Jac_dUrfm_dDCartUD = rfm.compute_Jacobian_and_inverseJacobian_tofrom_Cartesian()

    # Step 3.b: Convert basis of all BSSN *vectors* from Cartesian to destination basis
    ValenciavU_dst = rfm.basis_transform_vectorU_from_Cartesian_to_rfmbasis(Jac_dUrfm_dDCartUD, ValenciavU)

    # Note that the below the function should really be "...basis_transform_vectorUDfrom_Cartesian_to_rfmbasis.."
    AD_dst = rfm.basis_transform_vectorU_from_Cartesian_to_rfmbasis(Jac_dUCart_dDrfmUD, AD)

    BU_dst = ixp.zerorank1()
    if BU != None:
        BU_dst = rfm.basis_transform_vectorU_from_Cartesian_to_rfmbasis(Jac_dUrfm_dDCartUD, BU)
