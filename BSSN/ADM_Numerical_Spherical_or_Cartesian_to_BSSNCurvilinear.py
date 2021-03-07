# This module performs the conversion between ADM
# spacetime variables in Spherical or Cartesian coordinates
# given as *numerical* expressions (i.e., ADM quantities
# are only given to floating-point precision; e.g., in the
# case of an initial data solver), to rescaled BSSN-in-curvilinear
# coordinate quantities, as defined in BSSN_RHSs.py

# Author: Zachariah B. Etienne
#         zachetie **at** gmail **dot* com
# Step P1: Initialize core Python/NRPy+ modules
from outputC import lhrh,outCfunction,outputC  # NRPy+: Core C code output module
import sympy as sp                # SymPy: The Python computer algebra package upon which NRPy+ depends
import finite_difference as fin   # NRPy+: Finite difference C code generation module
import grid as gri                # NRPy+: Functions having to do with numerical grids
import indexedexp as ixp          # NRPy+: Symbolic indexed expression (e.g., tensors, vectors, etc.) support
import reference_metric as rfm    # NRPy+: Reference metric support
import BSSN.BSSN_quantities as Bq # NRPy+: Computes useful BSSN quantities; e.g., gammabarUU & GammabarUDD needed below
import os, sys                    # Standard Python modules for multiplatform OS-level functions

def Convert_Spherical_or_Cartesian_ADM_to_BSSN_curvilinear(CoordType_in, ADM_input_function_name,
                                                          Ccodesdir = "BSSN", pointer_to_ID_inputs=False,loopopts=",oldloops"):
    # The ADM & BSSN formalisms only work in 3D; they are 3+1 decompositions of Einstein's equations.
    #    To implement axisymmetry or spherical symmetry, simply set all spatial derivatives in
    #    the relevant angular directions to zero; DO NOT SET DIM TO ANYTHING BUT 3.

    # Step 0: Set spatial dimension (must be 3 for BSSN)
    DIM = 3

    # Step 1: All ADM initial data quantities are now functions of xx0,xx1,xx2, but
    #         they are still in the Spherical or Cartesian basis. We can now directly apply
    #         Jacobian transformations to get them in the correct xx0,xx1,xx2 basis:

    #         All input quantities are in terms of r,th,ph or x,y,z. We want them in terms
    #         of xx0,xx1,xx2, so here we call sympify_integers__replace_rthph() to replace
    #         r,th,ph or x,y,z, respectively, with the appropriate functions of xx0,xx1,xx2
    #         as defined for this particular reference metric in reference_metric.py's
    #         xxSph[] or xx_to_Cart[], respectively:

    # Define the input variables:
    gammaSphorCartDD = ixp.declarerank2("gammaSphorCartDD", "sym01")
    KSphorCartDD = ixp.declarerank2("KSphorCartDD", "sym01")
    alphaSphorCart = sp.symbols("alphaSphorCart")
    betaSphorCartU = ixp.declarerank1("betaSphorCartU")
    BSphorCartU = ixp.declarerank1("BSphorCartU")

    # Make sure that rfm.reference_metric() has been called.
    #    We'll need the variables it defines throughout this module.
    if rfm.have_already_called_reference_metric_function == False:
        print("Error. Called Convert_Spherical_ADM_to_BSSN_curvilinear() without")
        print("       first setting up reference metric, by calling rfm.reference_metric().")
        sys.exit(1)

    r_th_ph_or_Cart_xyz_oID_xx = []
    if CoordType_in == "Spherical":
        r_th_ph_or_Cart_xyz_oID_xx = rfm.xxSph
    elif CoordType_in == "Cartesian":
        r_th_ph_or_Cart_xyz_oID_xx = rfm.xx_to_Cart
    else:
        print("Error: Can only convert ADM Cartesian or Spherical initial data to BSSN Curvilinear coords.")
        sys.exit(1)

    # Step 2: All ADM initial data quantities are now functions of xx0,xx1,xx2, but
    #         they are still in the Spherical or Cartesian basis. We can now directly apply
    #         Jacobian transformations to get them in the correct xx0,xx1,xx2 basis:

    # alpha is a scalar, so no Jacobian transformation is necessary.
    alpha = alphaSphorCart

    Jac_dUSphorCart_dDrfmUD = ixp.zerorank2()
    for i in range(DIM):
        for j in range(DIM):
            Jac_dUSphorCart_dDrfmUD[i][j] = sp.diff(r_th_ph_or_Cart_xyz_oID_xx[i], rfm.xx[j])

    Jac_dUrfm_dDSphorCartUD, dummyDET = ixp.generic_matrix_inverter3x3(Jac_dUSphorCart_dDrfmUD)

    betaU = ixp.zerorank1()
    BU = ixp.zerorank1()
    gammaDD = ixp.zerorank2()
    KDD = ixp.zerorank2()
    for i in range(DIM):
        for j in range(DIM):
            betaU[i] += Jac_dUrfm_dDSphorCartUD[i][j] * betaSphorCartU[j]
            BU[i] += Jac_dUrfm_dDSphorCartUD[i][j] * BSphorCartU[j]
            for k in range(DIM):
                for l in range(DIM):
                    gammaDD[i][j] += Jac_dUSphorCart_dDrfmUD[k][i] * Jac_dUSphorCart_dDrfmUD[l][j] * \
                                     gammaSphorCartDD[k][l]
                    KDD[i][j] += Jac_dUSphorCart_dDrfmUD[k][i] * Jac_dUSphorCart_dDrfmUD[l][j] * KSphorCartDD[k][l]

    # Step 3: All ADM quantities were input into this function in the Spherical or Cartesian
    #         basis, as functions of r,th,ph or x,y,z, respectively. In Steps 1 and 2 above,
    #         we converted them to the xx0,xx1,xx2 basis, and as functions of xx0,xx1,xx2.
    #         Here we convert ADM quantities in the "rfm" basis to their BSSN Curvilinear
    #         counterparts, for all BSSN quantities *except* lambda^i:
    import BSSN.BSSN_in_terms_of_ADM as BitoA
    BitoA.gammabarDD_hDD(gammaDD)
    BitoA.trK_AbarDD_aDD(gammaDD, KDD)
    BitoA.cf_from_gammaDD(gammaDD)
    BitoA.betU_vetU(betaU, BU)
    hDD = BitoA.hDD
    trK = BitoA.trK
    aDD = BitoA.aDD
    cf = BitoA.cf
    vetU = BitoA.vetU
    betU = BitoA.betU

    # Step 4: Compute $\bar{\Lambda}^i$ (Eqs. 4 and 5 of
    #         [Baumgarte *et al.*](https://arxiv.org/pdf/1211.6632.pdf)),
    #         from finite-difference derivatives of rescaled metric
    #         quantities $h_{ij}$:

    # \bar{\Lambda}^i = \bar{\gamma}^{jk}\left(\bar{\Gamma}^i_{jk} - \hat{\Gamma}^i_{jk}\right).

    # The reference_metric.py module provides us with analytic expressions for
    #         $\hat{\Gamma}^i_{jk}$, so here we need only compute
    #         finite-difference expressions for $\bar{\Gamma}^i_{jk}$, based on
    #         the values for $h_{ij}$ provided in the initial data. Once
    #         $\bar{\Lambda}^i$ has been computed, we apply the usual rescaling
    #         procedure:

    # \lambda^i = \bar{\Lambda}^i/\text{ReU[i]},

    # and then output the result to a C file using the NRPy+
    #         finite-difference C output routine.

    # We will need all BSSN gridfunctions to be defined, as well as
    #     expressions for gammabarDD_dD in terms of exact derivatives of
    #     the rescaling matrix and finite-difference derivatives of
    #     hDD's. This functionality is provided by BSSN.BSSN_unrescaled_and_barred_vars,
    #     which we call here to overwrite above definitions of gammabarDD,gammabarUU, etc.
    Bq.gammabar__inverse_and_derivs() # Provides gammabarUU and GammabarUDD
    gammabarUU    = Bq.gammabarUU
    GammabarUDD   = Bq.GammabarUDD

    # Next evaluate \bar{\Lambda}^i, based on GammabarUDD above and GammahatUDD
    #       (from the reference metric):
    LambdabarU = ixp.zerorank1()
    for i in range(DIM):
        for j in range(DIM):
            for k in range(DIM):
                LambdabarU[i] += gammabarUU[j][k] * (GammabarUDD[i][j][k] - rfm.GammahatUDD[i][j][k])

    # Finally apply rescaling:
    # lambda^i = Lambdabar^i/\text{ReU[i]}
    lambdaU = ixp.zerorank1()
    for i in range(DIM):
        lambdaU[i] = LambdabarU[i] / rfm.ReU[i]

    if ADM_input_function_name == "DoNotOutputADMInputFunction":
        return hDD,aDD,trK,vetU,betU,alpha,cf,lambdaU

    # Step 5.A: Output files containing finite-differenced lambdas.
    outCparams = "preindent=1,outCfileaccess=a,outCverbose=False,includebraces=False"
    lambdaU_expressions = [lhrh(lhs=gri.gfaccess("in_gfs", "lambdaU0"), rhs=lambdaU[0]),
                           lhrh(lhs=gri.gfaccess("in_gfs", "lambdaU1"), rhs=lambdaU[1]),
                           lhrh(lhs=gri.gfaccess("in_gfs", "lambdaU2"), rhs=lambdaU[2])]

    desc = "Output lambdaU[i] for BSSN, built using finite-difference derivatives."
    name = "ID_BSSN_lambdas"
    params = "const paramstruct *restrict params,REAL *restrict xx[3],REAL *restrict in_gfs"
    preloop = ""
    opts = ""
    idx4replace = "IDX4S"
    if "oldloops" in loopopts:
        params = "const int Nxx[3],const int Nxx_plus_2NGHOSTS[3],REAL *xx[3],const REAL dxx[3],REAL *in_gfs"
        opts = "DisableCparameters"
        idx4replace = "IDX4"
        preloop = """
const REAL invdx0 = 1.0/dxx[0];
const REAL invdx1 = 1.0/dxx[1];
const REAL invdx2 = 1.0/dxx[2];
"""
    outCfunction(
        outfile=os.path.join(Ccodesdir, name + ".h"), desc=desc, name=name, params=params,
        preloop=preloop,
        body=fin.FD_outputC("returnstring", lambdaU_expressions, outCparams).replace("IDX4",idx4replace),
        loopopts="InteriorPoints,Read_xxs"+loopopts, opts=opts)

    # Step 5: Output all ADM-to-BSSN expressions to a C function. This function
    #         must first call the ID_ADM_SphorCart() defined above. Using these
    #         Spherical or Cartesian data, it sets up all quantities needed for
    #         BSSNCurvilinear initial data, *except* $\lambda^i$, which must be
    #         computed from numerical data using finite-difference derivatives.
    ID_inputs_param = "ID_inputs other_inputs,"
    if pointer_to_ID_inputs == True:
        ID_inputs_param = "ID_inputs *other_inputs,"

    desc = "Write BSSN variables in terms of ADM variables at a given point xx0,xx1,xx2"
    name = "ID_ADM_xx0xx1xx2_to_BSSN_xx0xx1xx2__ALL_BUT_LAMBDAs"
    opts = ""
    params = "const paramstruct *restrict params, "
    if "oldloops" in loopopts:
        opts = "DisableCparameters"
        params = ""
    params += "const REAL xx0xx1xx2[3]," + ID_inputs_param + """
                    REAL *hDD00,REAL *hDD01,REAL *hDD02,REAL *hDD11,REAL *hDD12,REAL *hDD22,
                    REAL *aDD00,REAL *aDD01,REAL *aDD02,REAL *aDD11,REAL *aDD12,REAL *aDD22,
                    REAL *trK,
                    REAL *vetU0,REAL *vetU1,REAL *vetU2,
                    REAL *betU0,REAL *betU1,REAL *betU2,
                    REAL *alpha,  REAL *cf"""
    outCparams = "preindent=1,outCverbose=False,includebraces=False"
    outCfunction(
        outfile=os.path.join(Ccodesdir, name + ".h"), desc=desc, name=name, params=params,
        body="""
      REAL gammaSphorCartDD00,gammaSphorCartDD01,gammaSphorCartDD02,
           gammaSphorCartDD11,gammaSphorCartDD12,gammaSphorCartDD22;
      REAL KSphorCartDD00,KSphorCartDD01,KSphorCartDD02,
           KSphorCartDD11,KSphorCartDD12,KSphorCartDD22;
      REAL alphaSphorCart,betaSphorCartU0,betaSphorCartU1,betaSphorCartU2;
      REAL BSphorCartU0,BSphorCartU1,BSphorCartU2;
      const REAL xx0 = xx0xx1xx2[0];
      const REAL xx1 = xx0xx1xx2[1];
      const REAL xx2 = xx0xx1xx2[2];
      REAL xyz_or_rthph[3];\n""" +
             outputC(r_th_ph_or_Cart_xyz_oID_xx[0:3], ["xyz_or_rthph[0]", "xyz_or_rthph[1]", "xyz_or_rthph[2]"],
                     "returnstring",
                     outCparams + ",CSE_enable=False") + "      " + ADM_input_function_name + """(xyz_or_rthph, other_inputs,
                       &gammaSphorCartDD00,&gammaSphorCartDD01,&gammaSphorCartDD02,
                       &gammaSphorCartDD11,&gammaSphorCartDD12,&gammaSphorCartDD22,
                       &KSphorCartDD00,&KSphorCartDD01,&KSphorCartDD02,
                       &KSphorCartDD11,&KSphorCartDD12,&KSphorCartDD22,
                       &alphaSphorCart,&betaSphorCartU0,&betaSphorCartU1,&betaSphorCartU2,
                       &BSphorCartU0,&BSphorCartU1,&BSphorCartU2);
      // Next compute all rescaled BSSN curvilinear quantities:\n""" +
             outputC([hDD[0][0], hDD[0][1], hDD[0][2], hDD[1][1], hDD[1][2], hDD[2][2],
                      aDD[0][0], aDD[0][1], aDD[0][2], aDD[1][1], aDD[1][2], aDD[2][2],
                      trK, vetU[0], vetU[1], vetU[2], betU[0], betU[1], betU[2],
                      alpha, cf],
                     ["*hDD00", "*hDD01", "*hDD02", "*hDD11", "*hDD12", "*hDD22",
                      "*aDD00", "*aDD01", "*aDD02", "*aDD11", "*aDD12", "*aDD22",
                      "*trK", "*vetU0", "*vetU1", "*vetU2", "*betU0", "*betU1", "*betU2",
                      "*alpha", "*cf"], "returnstring", params=outCparams),
        opts = opts)

    # Step 5.a: Output the driver function for the above
    #           function ID_ADM_xx0xx1xx2_to_BSSN_xx0xx1xx2__ALL_BUT_LAMBDAs()
    # Next write the driver function for ID_ADM_xx0xx1xx2_to_BSSN_xx0xx1xx2__ALL_BUT_LAMBDAs():
    desc = """Driver function for ID_ADM_xx0xx1xx2_to_BSSN_xx0xx1xx2__ALL_BUT_LAMBDAs(),
which writes BSSN variables in terms of ADM variables at a given point xx0,xx1,xx2"""
    name = "ID_BSSN__ALL_BUT_LAMBDAs"
    params = "const paramstruct *restrict params,REAL *restrict xx[3]," + ID_inputs_param + "REAL *in_gfs"
    opts = ""
    funccallparams = "params, "
    idx3replace   = "IDX3S"
    idx4ptreplace = "IDX4ptS"
    if "oldloops" in loopopts:
        params = "const int Nxx_plus_2NGHOSTS[3],REAL *xx[3]," + ID_inputs_param + "REAL *in_gfs"
        opts = "DisableCparameters"
        funccallparams = ""
        idx3replace   = "IDX3"
        idx4ptreplace = "IDX4pt"
    outCfunction(
        outfile=os.path.join(Ccodesdir, name + ".h"), desc=desc, name=name, params=params,
        body="""
const int idx = IDX3(i0,i1,i2);
const REAL xx0xx1xx2[3] = {xx0,xx1,xx2};
ID_ADM_xx0xx1xx2_to_BSSN_xx0xx1xx2__ALL_BUT_LAMBDAs(""".replace("IDX3",idx3replace)+funccallparams+"""xx0xx1xx2,other_inputs,
                    &in_gfs[IDX4pt(HDD00GF,idx)],&in_gfs[IDX4pt(HDD01GF,idx)],&in_gfs[IDX4pt(HDD02GF,idx)],
                    &in_gfs[IDX4pt(HDD11GF,idx)],&in_gfs[IDX4pt(HDD12GF,idx)],&in_gfs[IDX4pt(HDD22GF,idx)],
                    &in_gfs[IDX4pt(ADD00GF,idx)],&in_gfs[IDX4pt(ADD01GF,idx)],&in_gfs[IDX4pt(ADD02GF,idx)],
                    &in_gfs[IDX4pt(ADD11GF,idx)],&in_gfs[IDX4pt(ADD12GF,idx)],&in_gfs[IDX4pt(ADD22GF,idx)],
                    &in_gfs[IDX4pt(TRKGF,idx)],
                    &in_gfs[IDX4pt(VETU0GF,idx)],&in_gfs[IDX4pt(VETU1GF,idx)],&in_gfs[IDX4pt(VETU2GF,idx)],
                    &in_gfs[IDX4pt(BETU0GF,idx)],&in_gfs[IDX4pt(BETU1GF,idx)],&in_gfs[IDX4pt(BETU2GF,idx)],
                    &in_gfs[IDX4pt(ALPHAGF,idx)],&in_gfs[IDX4pt(CFGF,idx)]);
""".replace("IDX4pt",idx4ptreplace),
        loopopts="AllPoints,Read_xxs"+loopopts, opts=opts)
