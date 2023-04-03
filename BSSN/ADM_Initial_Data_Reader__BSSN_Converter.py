# This module performs the conversion from ADM
# spacetime variables in Spherical or Cartesian coordinates,
# to BSSN quantities in any basis supported by NRPy+.

# Author: Zachariah B. Etienne
#         zachetie **at** gmail **dot* com

# Initialize core Python/NRPy+ modules
# Step 1: Initialize core Python/NRPy+ modules
from outputC import outputC,lhrh,add_to_Cfunction_dict, Cfunction # NRPy+: Core C code output module
from outputC import outC_NRPy_basic_defines_h_dict
import NRPy_param_funcs as par    # NRPy+: Parameter interface
import finite_difference as fin   # NRPy+: Finite difference C code generation module
import grid as gri                # NRPy+: Functions having to do with numerical grids
import indexedexp as ixp          # NRPy+: Symbolic indexed expression (e.g., tensors, vectors, etc.) support
import reference_metric as rfm    # NRPy+: Reference metric support
import BSSN.BSSN_quantities as Bq # NRPy+: Computes useful BSSN quantities; e.g., gammabarUU & GammabarUDD needed below
from pickling import pickle_NRPy_env # NRPy+: Pickle/unpickle NRPy+ environment, for parallel codegen
import sympy as sp                # SymPy: The Python computer algebra package upon which NRPy+ depends
import os, sys                    # Standard Python modules for multiplatform OS-level functions


def add_to_Cfunction_dict_exact_ADM_ID_function(IDtype, IDCoordSystem, alpha, betaU, BU, gammaDD, KDD):
    includes = ["NRPy_basic_defines.h", "NRPy_function_prototypes.h"]
    desc = IDtype + " initial data"
    c_type = "void"
    name = IDtype
    params = "const paramstruct *params, const REAL xCart[3], const ID_persist_struct *restrict ID_persist, initial_data_struct *restrict initial_data"
    desired_rfm_coord = par.parval_from_str("reference_metric::CoordSystem")
    par.set_parval_from_str("reference_metric::CoordSystem", IDCoordSystem)
    rfm.reference_metric()
    body = ""
    if IDCoordSystem == "Spherical":
        body += r"""
  REAL xx0,xx1,xx2 __attribute__((unused));  // xx2 might be unused in the case of axisymmetric initial data.
  {
    int unused_Cart_to_i0i1i2[3];
    REAL xx[3];
    Cart_to_xx_and_nearest_i0i1i2(params, xCart, xx, unused_Cart_to_i0i1i2);
    xx0=xx[0];  xx1=xx[1];  xx2=xx[2];
  }
  const REAL r  = xx0; // Some ID only specify r,th,ph.
  const REAL th = xx1;
  const REAL ph = xx2;
"""
    elif IDCoordSystem == "Cartesian":
        body += r"""  const REAL Cartxyz0=xCart[0], Cartxyz1=xCart[1], Cartxyz2=xCart[2];
"""
    else:
        print("add_to_Cfunction_dict_exact_ADM_ID_function() Error: IDCoordSystem == " + IDCoordSystem + " unsupported")
        sys.exit(1)
    list_of_output_exprs = [alpha]
    list_of_output_varnames = ["initial_data->alpha"]
    for i in range(3):
        list_of_output_exprs += [betaU[i]]
        list_of_output_varnames += ["initial_data->betaSphorCartU" + str(i)]
        list_of_output_exprs += [BU[i]]
        list_of_output_varnames += ["initial_data->BSphorCartU" + str(i)]
        for j in range(i, 3):
            list_of_output_exprs += [gammaDD[i][j]]
            list_of_output_varnames += ["initial_data->gammaSphorCartDD" + str(i) + str(j)]
            list_of_output_exprs += [KDD[i][j]]
            list_of_output_varnames += ["initial_data->KSphorCartDD" + str(i) + str(j)]
    # Sort the outputs before calling outputC()
    # https://stackoverflow.com/questions/9764298/is-it-possible-to-sort-two-listswhich-reference-each-other-in-the-exact-same-w
    list_of_output_varnames, list_of_output_exprs = (list(t) for t in zip(*sorted(zip(list_of_output_varnames, list_of_output_exprs))))

    body += outputC(list_of_output_exprs, list_of_output_varnames,
                    filename="returnstring", params="outCverbose=False,includebraces=False,preindent=1")

    # Restore CoordSystem:
    par.set_parval_from_str("reference_metric::CoordSystem", desired_rfm_coord)
    rfm.reference_metric()
    add_to_Cfunction_dict(
        includes=includes,
        desc=desc, c_type=c_type, name=name, params=params,
        body=body,
        enableCparameters=True)
    return pickle_NRPy_env()


def Cfunction_ADM_SphorCart_to_Cart(input_Coord="Spherical", include_T4UU=False, extra_gfs=None):

    def get_gf_string(name, rank, dim):
        if rank==0:
            return "const REAL "+name+" = initial_data->"+name+";\n"
        else:
            gf_string = ""
            for i in range(dim):
                gf_string += get_gf_string(name+str(i), rank-1, dim)
            return gf_string
    
    includes = []

    desired_rfm_basis = par.parval_from_str("reference_metric::CoordSystem")

    desc = "Convert ADM variables from the spherical or Cartesian basis to the Cartesian basis"
    c_type = "static void"
    name = "ADM_SphorCart_to_Cart"
    params = """paramstruct *restrict params, const REAL xCart[3], const initial_data_struct *restrict initial_data,
                                  ADM_Cart_basis_struct *restrict ADM_Cart_basis"""

    body = r"""
  // Unpack initial_data for ADM vectors/tensors
"""
    gfs = ["betaU", "BU", "gammaDD", "KDD"]

    if extra_gfs is not None:
        if not isinstance(extra_gfs, list):
            extra_gfs = list(extra_gfs)
        gfs += extra_gfs

    if include_T4UU:
        gfs += ["T4UU"]

    string = ""
    for gf in gfs:
        rank = 0
        last = ""
        while gf[-1] == "U" or gf[-1] == "D":
            last  = gf[-1] + last
            gf    = gf[:-1]
            rank += 1

        if rank == 0:
            continue

        dim = 3 if "4" not in gf else 4
        string += get_gf_string(gf+"SphorCart"+last, rank, dim)

    body += string
    
    # for i in ["betaSphorCartU", "BSphorCartU"]:
    #     for j in range(3):
    #         varname = i + str(j)
    #         body += "  const REAL " + varname + " = initial_data->" + varname + ";\n"
    #     body += "\n"
    # for i in ["gammaSphorCartDD", "KSphorCartDD"]:
    #     for j in range(3):
    #         for k in range(j, 3):
    #             varname = i + str(j) + str(k)
    #             body += "  const REAL " + varname + " = initial_data->" + varname + ";\n"
    #     body += "\n"
    # # Read stress-energy tensor in spherical or Cartesian basis if desired.
    # if include_T4UU:
    #     for mu in range(4):
    #         for nu in range(mu, 4):
    #             varname = "T4SphorCartUU" + str(mu) + str(nu)
    #             body += "  const REAL " + varname + " = initial_data->" + varname + ";\n"
    #     body += "\n"

    if input_Coord == "Spherical":
        body += r"""
  // Perform the basis transform on ADM vectors/tensors from """ + input_Coord + """ to Cartesian:

  // Set destination xx[3] based on desired xCart[3]
  REAL xx0,xx1,xx2;
  {
    int unused_Cart_to_i0i1i2[3];
    REAL xx[3];
    Cart_to_xx_and_nearest_i0i1i2(params, xCart, xx, unused_Cart_to_i0i1i2);
    xx0=xx[0];  xx1=xx[1];  xx2=xx[2];
  }
"""
    # Set reference_metric to the input_Coord
    par.set_parval_from_str("reference_metric::CoordSystem", input_Coord)
    rfm.reference_metric()

    # Define the input variables:
    gammaSphorCartDD = ixp.declarerank2("gammaSphorCartDD", "sym01")
    KSphorCartDD     = ixp.declarerank2("KSphorCartDD", "sym01")
    betaSphorCartU = ixp.declarerank1("betaSphorCartU")
    BSphorCartU    = ixp.declarerank1("BSphorCartU")
    T4SphorCartUU = ixp.declarerank2("T4SphorCartUU", "sym01", DIM=4)

    # Compute Jacobian to convert to Cartesian coordinates
    Jac_dUCart_dDrfmUD, Jac_dUrfm_dDCartUD = rfm.compute_Jacobian_and_inverseJacobian_tofrom_Cartesian()

    gammaCartDD = rfm.basis_transform_tensorDD_from_rfmbasis_to_Cartesian(Jac_dUrfm_dDCartUD, gammaSphorCartDD)
    KCartDD = rfm.basis_transform_tensorDD_from_rfmbasis_to_Cartesian(Jac_dUrfm_dDCartUD, KSphorCartDD)
    betaCartU = rfm.basis_transform_vectorU_from_rfmbasis_to_Cartesian(Jac_dUCart_dDrfmUD, betaSphorCartU)
    BCartU = rfm.basis_transform_vectorU_from_rfmbasis_to_Cartesian(Jac_dUCart_dDrfmUD, BSphorCartU)
    T4CartUU = ixp.zerorank2(DIM=4)
    if include_T4UU:
        T4CartUU = rfm.basis_transform_4tensorUU_from_time_indep_rfmbasis_to_Cartesian(Jac_dUCart_dDrfmUD,
                                                                                       T4SphorCartUU)

    alpha = sp.symbols("initial_data->alpha", real=True)
    list_of_output_exprs    = [alpha]
    list_of_output_varnames = ["ADM_Cart_basis->alpha"]
    for i in range(3):
        list_of_output_exprs += [betaCartU[i]]
        list_of_output_varnames += ["ADM_Cart_basis->betaU" + str(i)]
        list_of_output_exprs += [BCartU[i]]
        list_of_output_varnames += ["ADM_Cart_basis->BU" + str(i)]
        for j in range(i, 3):
            list_of_output_exprs += [gammaCartDD[i][j]]
            list_of_output_varnames += ["ADM_Cart_basis->gammaDD" + str(i) + str(j)]
            list_of_output_exprs += [KCartDD[i][j]]
            list_of_output_varnames += ["ADM_Cart_basis->KDD" + str(i) + str(j)]
    if include_T4UU:
        for mu in range(4):
            for nu in range(mu, 4):
                list_of_output_exprs += [T4CartUU[mu][nu]]
                list_of_output_varnames += ["ADM_Cart_basis->T4UU" + str(mu) + str(nu)]

    # Sort the outputs before calling outputC()
    # https://stackoverflow.com/questions/9764298/is-it-possible-to-sort-two-listswhich-reference-each-other-in-the-exact-same-w
    list_of_output_varnames, list_of_output_exprs = (list(t) for t in zip(*sorted(zip(list_of_output_varnames, list_of_output_exprs))))

    body += outputC(list_of_output_exprs, list_of_output_varnames,
                    filename="returnstring", params="outCverbose=False,includebraces=False,preindent=1")

    # Restore reference metric globals to coordsystem on grid.
    par.set_parval_from_str("reference_metric::CoordSystem", desired_rfm_basis)
    rfm.reference_metric()

    _func_prototype, func = Cfunction(
        includes=includes,
        desc=desc,
        c_type=c_type, name=name, params=params,
        body=body,
        enableCparameters=False)
    return func


def Cfunction_ADM_Cart_to_BSSN_Cart(include_T4UU=False):
    includes = []

    desired_rfm_basis = par.parval_from_str("reference_metric::CoordSystem")

    desc = "Convert ADM variables in the Cartesian basis to BSSN variables in the Cartesian basis"
    c_type = "static void"
    name = "ADM_Cart_to_BSSN_Cart"
    params = """paramstruct *restrict params, const REAL xCart[3], const ADM_Cart_basis_struct *restrict ADM_Cart_basis,
                                  BSSN_Cart_basis_struct *restrict BSSN_Cart_basis"""

    # Extract desired rfm basis from reference_metric::CoordSystem
    desired_rfm_basis = par.parval_from_str("reference_metric::CoordSystem")

    # Set CoordSystem to Cartesian
    par.set_parval_from_str("reference_metric::CoordSystem", "Cartesian")
    rfm.reference_metric()

    gammaCartDD = ixp.declarerank2("ADM_Cart_basis->gammaDD", "sym01")
    KCartDD     = ixp.declarerank2("ADM_Cart_basis->KDD", "sym01")

    import BSSN.BSSN_in_terms_of_ADM as BitoA
    BitoA.trK_AbarDD_aDD(gammaCartDD, KCartDD)
    BitoA.gammabarDD_hDD(gammaCartDD)
    BitoA.cf_from_gammaDD(gammaCartDD)

    body = r"""
  // *In the Cartesian basis*, convert ADM quantities gammaDD & KDD
  //   into BSSN gammabarDD, AbarDD, cf, and trK.
  BSSN_Cart_basis->alpha = ADM_Cart_basis->alpha;
  BSSN_Cart_basis->betaU0 = ADM_Cart_basis->betaU0;
  BSSN_Cart_basis->betaU1 = ADM_Cart_basis->betaU1;
  BSSN_Cart_basis->betaU2 = ADM_Cart_basis->betaU2;
  BSSN_Cart_basis->BU0 = ADM_Cart_basis->BU0;
  BSSN_Cart_basis->BU1 = ADM_Cart_basis->BU1;
  BSSN_Cart_basis->BU2 = ADM_Cart_basis->BU2;
"""
    list_of_output_exprs    = [BitoA.cf, BitoA.trK]
    list_of_output_varnames = ["BSSN_Cart_basis->cf", "BSSN_Cart_basis->trK"]
    for i in range(3):
        for j in range(i, 3):
            list_of_output_exprs += [BitoA.gammabarDD[i][j]]
            list_of_output_varnames += ["BSSN_Cart_basis->gammabarDD" + str(i) + str(j)]
            list_of_output_exprs += [BitoA.AbarDD[i][j]]
            list_of_output_varnames += ["BSSN_Cart_basis->AbarDD" + str(i) + str(j)]
    if include_T4UU:
        T4CartUU = ixp.declarerank2("ADM_Cart_basis->T4UU", "sym01", DIM=4)
        for mu in range(4):
            for nu in range(mu, 4):
                list_of_output_exprs += [T4CartUU[mu][nu]]
                list_of_output_varnames += ["BSSN_Cart_basis->T4UU" + str(mu) + str(nu)]

    # Sort the outputs before calling outputC()
    # https://stackoverflow.com/questions/9764298/is-it-possible-to-sort-two-listswhich-reference-each-other-in-the-exact-same-w
    list_of_output_varnames, list_of_output_exprs = (list(t) for t in zip(*sorted(zip(list_of_output_varnames, list_of_output_exprs))))
    body += outputC(list_of_output_exprs, list_of_output_varnames,
                    filename="returnstring", params="outCverbose=False,includebraces=False,preindent=1")

    # Restore reference metric globals to desired reference metric.
    par.set_parval_from_str("reference_metric::CoordSystem", desired_rfm_basis)
    rfm.reference_metric()

    _func_prototype, func = Cfunction(
        includes=includes,
        desc=desc,
        c_type=c_type, name=name, params=params,
        body=body,
        enableCparameters=False)
    return func


def Cfunction_BSSN_Cart_to_rescaled_BSSN_rfm(rel_path_to_Cparams=os.path.join("."), include_T4UU=False):
    includes = []

    desc = r"""Convert Cartesian-basis BSSN vectors/tensors *except* lambda^i,
to the basis specified by `reference_metric::CoordSystem`, then rescale these BSSN quantities"""
    c_type = "static void"
    name = "BSSN_Cart_to_rescaled_BSSN_rfm"
    params = """paramstruct *restrict params, const REAL xCart[3],
                                           const BSSN_Cart_basis_struct *restrict BSSN_Cart_basis,
                                           rescaled_BSSN_rfm_basis_struct *restrict rescaled_BSSN_rfm_basis"""

    body = r"""
  REAL xx0,xx1,xx2 __attribute__((unused));  // xx2 might be unused in the case of axisymmetric initial data.
  {
    int unused_Cart_to_i0i1i2[3];
    REAL xx[3];
    Cart_to_xx_and_nearest_i0i1i2(params, xCart, xx, unused_Cart_to_i0i1i2);
    xx0=xx[0];  xx1=xx[1];  xx2=xx[2];
  }
"""

    # Define the input variables:
    gammabarCartDD = ixp.declarerank2("BSSN_Cart_basis->gammabarDD", "sym01")
    AbarCartDD     = ixp.declarerank2("BSSN_Cart_basis->AbarDD", "sym01")
    betaCartU = ixp.declarerank1("BSSN_Cart_basis->betaU")
    BCartU    = ixp.declarerank1("BSSN_Cart_basis->BU")

    # Compute Jacobian to convert to Cartesian coordinates
    Jac_dUCart_dDrfmUD, Jac_dUrfm_dDCartUD = rfm.compute_Jacobian_and_inverseJacobian_tofrom_Cartesian()

    gammabarDD = rfm.basis_transform_tensorDD_from_Cartesian_to_rfmbasis(Jac_dUCart_dDrfmUD, gammabarCartDD)
    AbarDD = rfm.basis_transform_tensorDD_from_Cartesian_to_rfmbasis(Jac_dUCart_dDrfmUD, AbarCartDD)
    betaU = rfm.basis_transform_vectorU_from_Cartesian_to_rfmbasis(Jac_dUrfm_dDCartUD, betaCartU)
    BU = rfm.basis_transform_vectorU_from_Cartesian_to_rfmbasis(Jac_dUrfm_dDCartUD, BCartU)

    # Next rescale:
    vetU = ixp.zerorank1()
    betU = ixp.zerorank1()
    hDD  = ixp.zerorank2()
    aDD  = ixp.zerorank2()
    for i in range(3):
        vetU[i] = betaU[i] / rfm.ReU[i]
        betU[i] =    BU[i] / rfm.ReU[i]
        for j in range(3):
            hDD[i][j] = (gammabarDD[i][j] - rfm.ghatDD[i][j]) / rfm.ReDD[i][j]
            aDD[i][j] = AbarDD[i][j] / rfm.ReDD[i][j]

    if include_T4UU:
        T4CartUU = ixp.declarerank2("BSSN_Cart_basis->T4UU", "sym01", DIM=4)
        T4UU = rfm.basis_transform_4tensorUU_from_Cartesian_to_time_indep_rfmbasis(Jac_dUrfm_dDCartUD, T4CartUU)

    alpha, cf, trK = sp.symbols('BSSN_Cart_basis->alpha BSSN_Cart_basis->cf BSSN_Cart_basis->trK', real=True)

    list_of_output_exprs    = [alpha, cf, trK]
    list_of_output_varnames = ["rescaled_BSSN_rfm_basis->alpha",
                               "rescaled_BSSN_rfm_basis->cf",
                               "rescaled_BSSN_rfm_basis->trK"]
    for i in range(3):
        list_of_output_exprs += [vetU[i]]
        list_of_output_varnames += ["rescaled_BSSN_rfm_basis->vetU" + str(i)]
        list_of_output_exprs += [betU[i]]
        list_of_output_varnames += ["rescaled_BSSN_rfm_basis->betU" + str(i)]
        for j in range(i, 3):
            list_of_output_exprs += [hDD[i][j]]
            list_of_output_varnames += ["rescaled_BSSN_rfm_basis->hDD" + str(i) + str(j)]
            list_of_output_exprs += [aDD[i][j]]
            list_of_output_varnames += ["rescaled_BSSN_rfm_basis->aDD" + str(i) + str(j)]
    if include_T4UU:
        for mu in range(4):
            for nu in range(mu, 4):
                # T4UU IS ASSUMED NOT RESCALED; RESCALINGS ARE HANDLED WITHIN BSSN RHSs, etc.
                list_of_output_exprs += [T4UU[mu][nu]]
                list_of_output_varnames += ["rescaled_BSSN_rfm_basis->T4UU" + str(mu) + str(nu)]

    # Sort the outputs before calling outputC()
    # https://stackoverflow.com/questions/9764298/is-it-possible-to-sort-two-listswhich-reference-each-other-in-the-exact-same-w
    list_of_output_varnames, list_of_output_exprs = (list(t) for t in zip(*sorted(zip(list_of_output_varnames, list_of_output_exprs))))

    body += outputC(list_of_output_exprs, list_of_output_varnames,
                    filename="returnstring", params="outCverbose=False,includebraces=False,preindent=1")

    _func_prototype, func = Cfunction(
        includes=includes,
        desc=desc,
        c_type=c_type, name=name, params=params,
        body=body,
        enableCparameters=True, rel_path_to_Cparams=rel_path_to_Cparams)
    return func


# initial_data_lambdaU_grid_interior() computes lambdaU from
#  finite-difference derivatives of rescaled metric quantities
def Cfunction_initial_data_lambdaU_grid_interior():
    includes = []
    c_type = "static void"

    output_Coord = par.parval_from_str("reference_metric::CoordSystem")
    desc = "Compute lambdaU in " + output_Coord + " coordinates"
    name = "initial_data_lambdaU_grid_interior"
    params = """const paramstruct *restrict params, REAL *restrict xx[3], REAL *restrict in_gfs"""
    # Step 7: Compute $\bar{\Lambda}^i$ from finite-difference derivatives of rescaled metric quantities

    # We will need all BSSN gridfunctions to be defined, as well as
    #     expressions for gammabarDD_dD in terms of exact derivatives of
    #     the rescaling matrix and finite-difference derivatives of
    #     hDD's. This functionality is provided by BSSN.BSSN_unrescaled_and_barred_vars,
    #     which we call here to overwrite above definitions of gammabarDD,gammabarUU, etc.
    Bq.gammabar__inverse_and_derivs()  # Provides gammabarUU and GammabarUDD
    gammabarUU    = Bq.gammabarUU
    GammabarUDD   = Bq.GammabarUDD

    # Next evaluate \bar{\Lambda}^i, based on GammabarUDD above and GammahatUDD
    #       (from the reference metric):
    LambdabarU = ixp.zerorank1()
    for i in range(3):
        for j in range(3):
            for k in range(3):
                LambdabarU[i] += gammabarUU[j][k] * (GammabarUDD[i][j][k] - rfm.GammahatUDD[i][j][k])

    # Finally apply rescaling:
    # lambda^i = Lambdabar^i/\text{ReU[i]}
    lambdaU = ixp.zerorank1()
    for i in range(3):
        lambdaU[i] = LambdabarU[i] / rfm.ReU[i]

    lambdaU_expressions = [lhrh(lhs=gri.gfaccess("in_gfs", "lambdaU0"), rhs=lambdaU[0]),
                           lhrh(lhs=gri.gfaccess("in_gfs", "lambdaU1"), rhs=lambdaU[1]),
                           lhrh(lhs=gri.gfaccess("in_gfs", "lambdaU2"), rhs=lambdaU[2])]
    body = fin.FD_outputC("returnstring", lambdaU_expressions,
                           params="outCverbose=False,includebraces=False,preindent=0")
    _func_prototype, func = Cfunction(
        includes=includes,
        desc=desc,
        c_type=c_type, name=name, params=params,
        body=body,
        loopopts="InteriorPoints,Read_xxs",
        enableCparameters=True)
    return func


def add_to_Cfunction_dict_initial_data_reader__convert_ADM_Sph_or_Cart_to_BSSN(addl_includes=None,
                                                                               rel_path_to_Cparams=os.path.join("."),
                                                                               input_Coord="Spherical",
                                                                               include_T4UU=False,
                                                                               extra_adm_gfs=None,
                                                                               extra_bssn_gfs=None):
    includes = ["NRPy_basic_defines.h", "NRPy_function_prototypes.h"]
    if par.parval_from_str("finite_difference::enable_FD_functions"):
        includes += ["finite_difference_functions.h"]
    if addl_includes is not None:
        if not isinstance(addl_includes, list):
            print("Error: addl_includes must be a list.")
            sys.exit(1)
        includes += addl_includes

    extra_adm_gfs_str = ""
    if extra_adm_gfs is not None:
        if not isinstance(extra_adm_gfs, list):
            extra_adm_gfs = list(extra_adm_gfs)
        extra_adm_gfs_str = "REAL "+", ".join(extra_adm_gfs)+";\n"

    extra_bssn_gfs_str = ""
    if extra_adm_gfs is not None:
        if not isinstance(extra_adm_gfs, list):
            extra_bssn_gfs = list(extra_bssn_gfs)
        extra_bssn_gfs_str = "REAL "+", ".join(extra_bssn_gfs)+";\n"

    def T4UU_prettyprint():
        return r"""
  REAL T4UU00,T4UU01,T4UU02,T4UU03;
  REAL        T4UU11,T4UU12,T4UU13;
  REAL               T4UU22,T4UU23;
  REAL                      T4UU33;
"""
    prefunc = """
// ADM variables in the Cartesian basis:
typedef struct __ADM_Cart_basis_struct__ {
  REAL alpha, betaU0,betaU1,betaU2, BU0,BU1,BU2;
  REAL gammaDD00,gammaDD01,gammaDD02,gammaDD11,gammaDD12,gammaDD22;
  REAL KDD00,KDD01,KDD02,KDD11,KDD12,KDD22;
"""
    if extra_adm_gfs is not None:
        prefunc += extra_adm_gfs_str

    if include_T4UU:
        prefunc += T4UU_prettyprint()
    prefunc += "} ADM_Cart_basis_struct;\n"
    ##############
    prefunc += """
// BSSN variables in the Cartesian basis:
typedef struct __BSSN_Cart_basis_struct__ {
  REAL alpha, betaU0,betaU1,betaU2, BU0,BU1,BU2;
  REAL cf, trK;
  REAL gammabarDD00,gammabarDD01,gammabarDD02,gammabarDD11,gammabarDD12,gammabarDD22;
  REAL AbarDD00,AbarDD01,AbarDD02,AbarDD11,AbarDD12,AbarDD22;
"""
    if extra_bssn_gfs is not None:
        prefunc += extra_bssn_gfs_str

    if include_T4UU:
        prefunc += T4UU_prettyprint()
    prefunc += "} BSSN_Cart_basis_struct;\n"
    ##############
    prefunc += """
// Rescaled BSSN variables in the rfm basis:
typedef struct __rescaled_BSSN_rfm_basis_struct__ {
  REAL alpha, vetU0,vetU1,vetU2, betU0,betU1,betU2;
  REAL cf, trK;
  REAL hDD00,hDD01,hDD02,hDD11,hDD12,hDD22;
  REAL aDD00,aDD01,aDD02,aDD11,aDD12,aDD22;
"""
    if extra_bssn_gfs is not None:
        prefunc += extra_bssn_gfs_str

    if include_T4UU:
        prefunc += T4UU_prettyprint()
    prefunc += "} rescaled_BSSN_rfm_basis_struct;\n"
    ##############
    ##############
    prefunc += Cfunction_ADM_SphorCart_to_Cart(input_Coord=input_Coord, include_T4UU=include_T4UU,)
    prefunc += Cfunction_ADM_Cart_to_BSSN_Cart(                         include_T4UU=include_T4UU)
    prefunc += Cfunction_BSSN_Cart_to_rescaled_BSSN_rfm(rel_path_to_Cparams=rel_path_to_Cparams  ,
                                                        include_T4UU=include_T4UU)
    prefunc += Cfunction_initial_data_lambdaU_grid_interior()

    output_Coord = par.parval_from_str("reference_metric::CoordSystem")
    desc = "Read in ADM initial data in the " + input_Coord + " basis, and convert to BSSN data in " + output_Coord + " coordinates"
    c_type = "void"
    name = "initial_data_reader__convert_ADM_" + input_Coord + "_to_BSSN"
    params = """griddata_struct *restrict griddata, ID_persist_struct *restrict ID_persist,
                                                             void ID_function(const paramstruct *params, const REAL xCart[3],
                                                                              const ID_persist_struct *restrict ID_persist,
                                                                              initial_data_struct *restrict initial_data)"""

    body = r"""
  const int Nxx_plus_2NGHOSTS0 = griddata->params.Nxx_plus_2NGHOSTS0;
  const int Nxx_plus_2NGHOSTS1 = griddata->params.Nxx_plus_2NGHOSTS1;
  const int Nxx_plus_2NGHOSTS2 = griddata->params.Nxx_plus_2NGHOSTS2;

  LOOP_OMP("omp parallel for", i0,0,Nxx_plus_2NGHOSTS0, i1,0,Nxx_plus_2NGHOSTS1, i2,0,Nxx_plus_2NGHOSTS2) {
    // xCart is the global Cartesian coordinate, which accounts for any grid offsets from the origin.
    REAL xCart[3];  xx_to_Cart(&griddata->params, griddata->xx, i0,i1,i2, xCart);

    // Read or compute initial data at destination point xCart
    initial_data_struct initial_data;
    ID_function(&griddata->params, xCart, ID_persist, &initial_data);

    ADM_Cart_basis_struct ADM_Cart_basis;
    ADM_SphorCart_to_Cart(&griddata->params, xCart, &initial_data, &ADM_Cart_basis);

    BSSN_Cart_basis_struct BSSN_Cart_basis;
    ADM_Cart_to_BSSN_Cart(&griddata->params, xCart, &ADM_Cart_basis, &BSSN_Cart_basis);

    rescaled_BSSN_rfm_basis_struct rescaled_BSSN_rfm_basis;
    BSSN_Cart_to_rescaled_BSSN_rfm(&griddata->params, xCart, &BSSN_Cart_basis, &rescaled_BSSN_rfm_basis);

    const int idx3 = IDX3S(i0,i1,i2);

    // Output data to gridfunctions
"""
    gf_list = ["alpha", "trK", "cf"]

    if extra_bssn_gfs is not None:
        gf_list += extra_bssn_gfs
    
    for i in range(3):
        gf_list += ["vetU"+str(i), "betU"+str(i)]
        for j in range(i, 3):
            gf_list += ["hDD"+str(i)+str(j), "aDD"+str(i)+str(j)]
    for gf in sorted(gf_list):
        body += "    griddata->gridfuncs.y_n_gfs[IDX4ptS("+gf.upper()+"GF, idx3)] = rescaled_BSSN_rfm_basis."+gf+";\n"
    if include_T4UU:
        for mu in range(4):
            for nu in range(mu, 4):
                gf = "T4UU" + str(mu) + str(nu)
                body += "    griddata->gridfuncs.auxevol_gfs[IDX4ptS("+gf.upper()+"GF, idx3)] = rescaled_BSSN_rfm_basis."+gf+";\n"
    body += """
  } // END LOOP over all gridpoints on given grid

  initial_data_lambdaU_grid_interior(&griddata->params, griddata->xx, griddata->gridfuncs.y_n_gfs);
"""

    add_to_Cfunction_dict(
        includes=includes,
        prefunc=prefunc,
        desc=desc,
        c_type=c_type, name=name, params=params,
        body=body,
        enableCparameters=False)
    return pickle_NRPy_env()


# Other than its core use as a means to store ADM input quantities,
# `initial_data_struct` is designed to be extensible. For example, it may be
# used to store e.g., pseudospectral coefficients for TwoPunctures,
# initial data gridfunctions from NRPyElliptic, pointers to TOV 1D data
# from the TOV solver, etc.
def register_NRPy_basic_defines(ID_persist_struct_contents_str="", include_T4UU=False, extra_ID_gfs=None):
    Nbd = r"""typedef struct __initial_data_struct__ {
  REAL alpha;

  REAL betaSphorCartU0, betaSphorCartU1, betaSphorCartU2;
  REAL BSphorCartU0, BSphorCartU1, BSphorCartU2;

  REAL gammaSphorCartDD00, gammaSphorCartDD01, gammaSphorCartDD02;
  REAL gammaSphorCartDD11, gammaSphorCartDD12, gammaSphorCartDD22;

  REAL KSphorCartDD00, KSphorCartDD01, KSphorCartDD02;
  REAL KSphorCartDD11, KSphorCartDD12, KSphorCartDD22;
"""
    if extra_ID_gfs is not None:
        if not isinstance(extra_ID_gfs, list):
            extra_gfs = list(extra_ID_gfs)
        Nbd += "REAL " + ", ".join(extra_ID_gfs) + ";\n"
        
    
    if include_T4UU:
        Nbd += """
  REAL T4SphorCartUU00,T4SphorCartUU01,T4SphorCartUU02,T4SphorCartUU03;
  REAL                 T4SphorCartUU11,T4SphorCartUU12,T4SphorCartUU13;
  REAL                                 T4SphorCartUU22,T4SphorCartUU23;
  REAL                                                 T4SphorCartUU33;
"""
    Nbd += """
} initial_data_struct;
"""

    Nbd += "typedef struct __ID_persist_struct__ {\n"
    Nbd += ID_persist_struct_contents_str + "\n"
    Nbd += "} ID_persist_struct;\n"
    outC_NRPy_basic_defines_h_dict["BSSN_initial_data"] = Nbd
