# This module provides functions that convert BSSN tensors and vectors
#   from one coordinate basis to another

# Author: Zachariah B. Etienne
#         zachetie **at** gmail **dot* com

# Step P1: Import needed NRPy+ core modules:
import indexedexp as ixp          # NRPy+: Symbolic indexed expression (e.g., tensors, vectors, etc.) support
import reference_metric as rfm    # NRPy+: Reference metric support
import NRPy_param_funcs as par    # NRPy+: Parameter interface

def BSSN_basis_transform(src_basis,src_xx, dst_basis,dst_xx,
                         src_hDD,src_aDD,src_lambdaU,src_vetU,src_betU):

    # Step 0: Store the current reference_metric::CoordSystem,
    #         so we can restore it at the end of this function
    origCoordSystem = par.parval_from_str("reference_metric::CoordSystem")

    # Step 1: Unrescale all BSSN variables

    par.set_parval_from_str("reference_metric::CoordSystem",src_basis)
    rfm.reference_metric()

    # STOLEN FROM BSSN/BSSN_quantities.py:
    # Step 1.a: gammabarDD and AbarDD:
    src_gammabarDD = ixp.zerorank2()
    src_AbarDD = ixp.zerorank2()
    for i in range(3):
        for j in range(3):
            # gammabar_{ij}  = h_{ij}*ReDD[i][j] + gammahat_{ij}
            src_gammabarDD[i][j] = src_hDD[i][j] * rfm.ReDD[i][j] + rfm.ghatDD[i][j]
            # Abar_{ij}      = a_{ij}*ReDD[i][j]
            src_AbarDD[i][j]     = src_aDD[i][j] * rfm.ReDD[i][j]

    # Step 1.b: LambdabarU, betaU, and BU:
    src_LambdabarU = ixp.zerorank1()
    src_betaU = ixp.zerorank1()
    src_BU = ixp.zerorank1()
    for i in range(3):
        src_LambdabarU[i] = src_lambdaU[i] * rfm.ReU[i]
        src_betaU[i]      =    src_vetU[i] * rfm.ReU[i]
        src_BU[i]         =    src_betU[i] * rfm.ReU[i]


    # Step 2: Transform source grid basis to Cartesian, using center of source grid as origin

    # Step 2.a: Construct Jacobian & Inverse Jacobians:
    Jac_dUCart_dDrfmUD,Jac_dUrfm_dDCartUD = rfm.compute_Jacobian_and_inverseJacobian_tofrom_Cartesian()

    # Step 2.b: Convert basis of all BSSN *vectors* to Cartesian
    CartLambdabarU = rfm.basis_transform_vectorU_from_rfmbasis_to_Cartesian(Jac_dUCart_dDrfmUD, src_LambdabarU)
    CartbetaU      = rfm.basis_transform_vectorU_from_rfmbasis_to_Cartesian(Jac_dUCart_dDrfmUD, src_betaU)
    CartBU         = rfm.basis_transform_vectorU_from_rfmbasis_to_Cartesian(Jac_dUCart_dDrfmUD, src_BU)

    # Step 2.c: Convert basis of all BSSN *tensors* to Cartesian
    CartgammabarDD = rfm.basis_transform_tensorDD_from_rfmbasis_to_Cartesian(Jac_dUrfm_dDCartUD, src_gammabarDD)
    CartAbarDD     = rfm.basis_transform_tensorDD_from_rfmbasis_to_Cartesian(Jac_dUrfm_dDCartUD, src_AbarDD)

    # Step 2.d: All BSSN tensor/vector quantities are written in terms of
    #           rescaled quantities and (xx0,xx1,xx2) on the SOURCE grid.
    #           To avoid confusion with (xx0,xx1,xx2) on the DESTINATION grid,
    #           we replace (xx0,xx1,xx2) with (src_xx0,src_xx1,src_xx2) here:
    for i in range(3):
        for k in range(3):
            CartLambdabarU[i] = CartLambdabarU[i].subs(rfm.xx[k], src_xx[k])
            CartbetaU[i]      =      CartbetaU[i].subs(rfm.xx[k], src_xx[k])
            CartBU[i]         =         CartBU[i].subs(rfm.xx[k], src_xx[k])
    for i in range(3):
        for j in range(3):
            for k in range(3):
                CartgammabarDD[i][j] = CartgammabarDD[i][j].subs(rfm.xx[k], src_xx[k])
                CartAbarDD[i][j]     =     CartAbarDD[i][j].subs(rfm.xx[k], src_xx[k])

    # Step 3: Transform BSSN tensors in Cartesian basis to destination grid basis, using center of dest. grid as origin

    # Step 3.a: Set up destination grid coordinate system
    par.set_parval_from_str("reference_metric::CoordSystem",dst_basis)
    rfm.reference_metric()

    # Step 3.b: Next construct Jacobian and inverse Jacobian matrices:
    Jac_dUCart_dDrfmUD,Jac_dUrfm_dDCartUD = rfm.compute_Jacobian_and_inverseJacobian_tofrom_Cartesian()

    # Step 3.c: Convert basis of all BSSN *vectors* from Cartesian to destination basis
    dst_LambdabarU = rfm.basis_transform_vectorU_from_Cartesian_to_rfmbasis(Jac_dUrfm_dDCartUD, CartLambdabarU)
    dst_betaU      = rfm.basis_transform_vectorU_from_Cartesian_to_rfmbasis(Jac_dUrfm_dDCartUD, CartbetaU)
    dst_BU         = rfm.basis_transform_vectorU_from_Cartesian_to_rfmbasis(Jac_dUrfm_dDCartUD, CartBU)

    # Step 3.d: Convert basis of all BSSN *tensors* from Cartesian to destination basis
    dst_gammabarDD = rfm.basis_transform_tensorDD_from_Cartesian_to_rfmbasis(Jac_dUCart_dDrfmUD, CartgammabarDD)
    dst_AbarDD     = rfm.basis_transform_tensorDD_from_Cartesian_to_rfmbasis(Jac_dUCart_dDrfmUD, CartAbarDD)

    # Step 4: Rescale all BSSN quantities

    # BASED ON BSSN/BSSN_quantities.py:
    # Step 4.a: hDD and aDD:
    global dst_hDD, dst_aDD
    dst_hDD = ixp.zerorank2()
    dst_aDD = ixp.zerorank2()
    for i in range(3):
        for j in range(3):
            # gammabar_{ij}  = h_{ij}*ReDD[i][j] + gammahat_{ij}
            # ==>     h_{ij} = (gammabar_{ij} - gammahat_{ij}) / ReDD[i][j]
            dst_hDD[i][j] = (dst_gammabarDD[i][j] - rfm.ghatDD[i][j]) / rfm.ReDD[i][j]
            # Abar_{ij}      = a_{ij}*ReDD[i][j]
            # ==>     a_{ij} = Abar_{ij}/ReDD[i][j]
            dst_aDD[i][j] = dst_AbarDD[i][j] / rfm.ReDD[i][j]

    # Step 4.b: lambdaU, vetU, and betU:
    global dst_lambdaU, dst_vetU, dst_betU
    dst_lambdaU = ixp.zerorank1()
    dst_vetU    = ixp.zerorank1()
    dst_betU    = ixp.zerorank1()
    for i in range(3):
        # Lambdabar^i = \lambda^i * ReU[i]
        # ==>  \lambda^i = Lambdabar^i / ReU[i]
        dst_lambdaU[i] = dst_LambdabarU[i] / rfm.ReU[i]
        dst_vetU[i]    =      dst_betaU[i] / rfm.ReU[i]
        dst_betU[i]    =         dst_BU[i] / rfm.ReU[i]

    # Step 4.c: All BSSN tensor/vector quantities are written in terms of
    #           rescaled quantities and (xx0,xx1,xx2) on the DESTINATION grid.
    #           To avoid confusion with (xx0,xx1,xx2) on the SOURCE grid,
    #           we replace (xx0,xx1,xx2) with (dst_xx0,dst_xx1,dst_xx2) here:
    for i in range(3):
        for k in range(3):
            dst_lambdaU[i] = dst_lambdaU[i].subs(rfm.xx[k],dst_xx[k])
            dst_vetU[i]    =    dst_vetU[i].subs(rfm.xx[k],dst_xx[k])
            dst_betU[i]    =    dst_betU[i].subs(rfm.xx[k],dst_xx[k])
    for i in range(3):
        for j in range(3):
            for k in range(3):
                dst_hDD[i][j] = dst_hDD[i][j].subs(rfm.xx[k],dst_xx[k])
                dst_aDD[i][j] = dst_aDD[i][j].subs(rfm.xx[k],dst_xx[k])

    # Step 5: Restore reference_metric::CoordSystem to its value when
    #         this function was called.
    par.set_parval_from_str("reference_metric::CoordSystem",origCoordSystem)
    rfm.reference_metric()
