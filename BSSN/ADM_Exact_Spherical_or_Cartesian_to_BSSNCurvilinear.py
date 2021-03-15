# This module performs the conversion between ADM
# spacetime variables in Spherical or Cartesian coordinates
# given as *exact* SymPy expressions (i.e., direct functions
# of r,th,ph or x,y,z), to rescaled BSSN-in-curvilinear
# coordinate quantities, as defined in BSSN_RHSs.py

# Author: Zachariah B. Etienne
#         zachetie **at** gmail **dot* com

# Step P0: Import needed Python/NRPy+ modules
import sympy as sp                # SymPy: The Python computer algebra package upon which NRPy+ depends
import indexedexp as ixp          # NRPy+: Symbolic indexed expression (e.g., tensors, vectors, etc.) support
import reference_metric as rfm    # NRPy+: Reference metric support
import sys                        # Standard Python module for multiplatform OS-level functions

def Convert_Spherical_or_Cartesian_ADM_to_BSSN_curvilinear(CoordType_in, Sph_r_th_ph_or_Cart_xyz,
        gammaDD_inSphorCart, KDD_inSphorCart, alpha_inSphorCart, betaU_inSphorCart, BU_inSphorCart):
    # This routine converts the ADM variables
    #    $$\left\{\gamma_{ij}, K_{ij}, \alpha, \beta^i\right\}$$
    #    in Spherical or Cartesian basis+coordinates, first to the BSSN variables
    #    in the chosen reference_metric::CoordSystem coordinate system+basis:
    # $$\left\{\bar{\gamma}_{i j},\bar{A}_{i j},\phi, K, \bar{\Lambda}^{i}, \alpha, \beta^i, B^i\right\},$$
    #    and then to the rescaled variables:
    # $$\left\{h_{i j},a_{i j},\phi, K, \lambda^{i}, \alpha, \mathcal{V}^i, \mathcal{B}^i\right\}.$$

    # The ADM & BSSN formalisms only work in 3D; they are 3+1 decompositions of Einstein's equations.
    #    To implement axisymmetry or spherical symmetry, simply set all spatial derivatives in
    #    the relevant angular directions to zero; DO NOT SET DIM TO ANYTHING BUT 3.
    # Step P1: Set spatial dimension (must be 3 for BSSN)
    DIM = 3

    # Step P2: Copy gammaSphDD_in to gammaSphDD, KSphDD_in to KSphDD, etc.
    #    This ensures that the input arrays are not modified below;
    #    modifying them would result in unexpected effects outside
    #    this function.
    alphaSphorCart = alpha_inSphorCart
    betaSphorCartU   = ixp.zerorank1()
    BSphorCartU      = ixp.zerorank1()
    gammaSphorCartDD = ixp.zerorank2()
    KSphorCartDD     = ixp.zerorank2()
    for i in range(DIM):
        betaSphorCartU[i] = betaU_inSphorCart[i]
        BSphorCartU[i]    = BU_inSphorCart[i]
        for j in range(DIM):
            gammaSphorCartDD[i][j] = gammaDD_inSphorCart[i][j]
            KSphorCartDD[i][j]     = KDD_inSphorCart[i][j]

    # Make sure that rfm.reference_metric() has been called.
    #    We'll need the variables it defines throughout this module.
    if rfm.have_already_called_reference_metric_function == False:
        print("Error. Called Convert_Spherical_ADM_to_BSSN_curvilinear() without")
        print("       first setting up reference metric, by calling rfm.reference_metric().")
        sys.exit(1)

    # Step 1: All input quantitiefs are in terms of r,th,ph or x,y,z. We want them in terms
    #         of xx0,xx1,xx2, so here we call sympify_integers__replace_rthph() to replace
    #         r,th,ph or x,y,z, respectively, with the appropriate functions of xx0,xx1,xx2
    #         as defined for this particular reference metric in reference_metric.py's
    #         xxSph[] or xx_to_Cart[], respectively:

    # Note that substitution only works when the variable is not an integer. Hence the
    #         if isinstance(...,...) stuff:
    def sympify_integers__replace_rthph_or_Cartxyz(obj, rthph_or_xyz, rthph_or_xyz_of_xx):
        if isinstance(obj, int):
            return sp.sympify(obj)
        return obj.subs(rthph_or_xyz[0], rthph_or_xyz_of_xx[0]).\
            subs(rthph_or_xyz[1], rthph_or_xyz_of_xx[1]).\
            subs(rthph_or_xyz[2], rthph_or_xyz_of_xx[2])

    r_th_ph_or_Cart_xyz_of_xx = []
    if CoordType_in == "Spherical":
        r_th_ph_or_Cart_xyz_of_xx = rfm.xxSph
    elif CoordType_in == "Cartesian":
        r_th_ph_or_Cart_xyz_of_xx = rfm.xx_to_Cart
    else:
        print("Error: Can only convert ADM Cartesian or Spherical initial data to BSSN Curvilinear coords.")
        sys.exit(1)

    alphaSphorCart = sympify_integers__replace_rthph_or_Cartxyz(
        alphaSphorCart, Sph_r_th_ph_or_Cart_xyz, r_th_ph_or_Cart_xyz_of_xx)
    for i in range(DIM):
        betaSphorCartU[i] = sympify_integers__replace_rthph_or_Cartxyz(
            betaSphorCartU[i], Sph_r_th_ph_or_Cart_xyz, r_th_ph_or_Cart_xyz_of_xx)
        BSphorCartU[i]    = sympify_integers__replace_rthph_or_Cartxyz(
            BSphorCartU[i],    Sph_r_th_ph_or_Cart_xyz, r_th_ph_or_Cart_xyz_of_xx)
        for j in range(DIM):
            gammaSphorCartDD[i][j] = sympify_integers__replace_rthph_or_Cartxyz(
                gammaSphorCartDD[i][j], Sph_r_th_ph_or_Cart_xyz, r_th_ph_or_Cart_xyz_of_xx)
            KSphorCartDD[i][j]     = sympify_integers__replace_rthph_or_Cartxyz(
                KSphorCartDD[i][j],     Sph_r_th_ph_or_Cart_xyz, r_th_ph_or_Cart_xyz_of_xx)

    # Step 2: All ADM initial data quantities are now functions of xx0,xx1,xx2, but
    #         they are still in the Spherical or Cartesian basis. We can now directly apply
    #         Jacobian transformations to get them in the correct xx0,xx1,xx2 basis:

    # alpha is a scalar, so no Jacobian transformation is necessary.
    alpha = alphaSphorCart

    Jac_dUSphorCart_dDrfmUD = ixp.zerorank2()
    for i in range(DIM):
        for j in range(DIM):
            Jac_dUSphorCart_dDrfmUD[i][j] = sp.diff(r_th_ph_or_Cart_xyz_of_xx[i],rfm.xx[j])

    Jac_dUrfm_dDSphorCartUD, dummyDET = ixp.generic_matrix_inverter3x3(Jac_dUSphorCart_dDrfmUD)

    betaU   = ixp.zerorank1()
    BU      = ixp.zerorank1()
    gammaDD = ixp.zerorank2()
    KDD     = ixp.zerorank2()
    for i in range(DIM):
        for j in range(DIM):
            betaU[i] += Jac_dUrfm_dDSphorCartUD[i][j] * betaSphorCartU[j]
            BU[i]    += Jac_dUrfm_dDSphorCartUD[i][j] * BSphorCartU[j]
            for k in range(DIM):
                for l in range(DIM):
                    gammaDD[i][j] += Jac_dUSphorCart_dDrfmUD[k][i]*Jac_dUSphorCart_dDrfmUD[l][j] * gammaSphorCartDD[k][l]
                    KDD[i][j]     += Jac_dUSphorCart_dDrfmUD[k][i]*Jac_dUSphorCart_dDrfmUD[l][j] *     KSphorCartDD[k][l]


    # Step 3: All ADM quantities were input into this function in the Spherical or Cartesian
    #         basis, as functions of r,th,ph or x,y,z, respectively. In Steps 1 and 2 above,
    #         we converted them to the xx0,xx1,xx2 basis, and as functions of xx0,xx1,xx2.
    #         Here we convert ADM quantities in the "rfm" basis to their BSSN Curvilinear
    #         counterparts:

    import BSSN.BSSN_in_terms_of_ADM as BitoA
    BitoA.gammabarDD_hDD(gammaDD)
    BitoA.trK_AbarDD_aDD(gammaDD, KDD)
    BitoA.LambdabarU_lambdaU__exact_gammaDD(gammaDD)
    BitoA.cf_from_gammaDD(gammaDD)
    BitoA.betU_vetU(betaU, BU)

    # Step 4: Return the BSSN Curvilinear variables in the desired xx0,xx1,xx2
    #         basis, and as functions of the consistent xx0,xx1,xx2 coordinates.
    return BitoA.cf, BitoA.hDD, BitoA.lambdaU, BitoA.aDD, BitoA.trK, alpha, BitoA.vetU, BitoA.betU
