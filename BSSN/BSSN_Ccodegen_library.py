# Generating full C library for solving
#  the BSSN equations in
#  ***curvilinear*** coordinates, using a
#  reference-metric formalism
#
# Author: Zachariah B. Etienne
#         zachetie **at** gmail **dot* com
#
# License: BSD 2-Clause

# COMPLETE DOCUMENTATION (JUPYTER NOTEBOOKS):
# START PAGE (start here!):  ../NRPy+_Tutorial.ipynb
# THIS MODULE: ../Tutorial-BSSN_time_evolution-C_codegen_library.ipynb

# RULES FOR ADDING FUNCTIONS TO THIS ROUTINE:
# 1. The function must be runnable from a multiprocessing environment,
#    which means that the function
# 1.a: cannot depend on previous function calls.
# 1.b: cannot create directories (this is not multiproc friendly)


# Step P1: Import needed NRPy+ core modules:
from outputC import lhrh, add_to_Cfunction_dict, outputC  # NRPy+: Core C code output module
import finite_difference as fin  # NRPy+: Finite difference C code generation module
import NRPy_param_funcs as par   # NRPy+: Parameter interface
import grid as gri               # NRPy+: Functions having to do with numerical grids
import indexedexp as ixp         # NRPy+: Symbolic indexed expression (e.g., tensors, vectors, etc.) support
import reference_metric as rfm   # NRPy+: Reference metric support
from pickling import pickle_NRPy_env   # NRPy+: Pickle/unpickle NRPy+ environment, for parallel codegen
import os, time, sys             # Standard Python modules for multiplatform OS-level functions, benchmarking
import sympy as sp               # SymPy: The Python computer algebra package upon which NRPy+ depends
import BSSN.BSSN_RHSs as rhs
import BSSN.BSSN_gauge_RHSs as gaugerhs
import BSSN.Enforce_Detgammahat_Constraint as EGC
import BSSN.Psi4 as psi4
import BSSN.Psi4_tetrads as psi4tet
import SpinWeight_minus2_SphHarmonics.SpinWeight_minus2_SphHarmonics as SWm2SH
import loop as lp

###############################################
# Helper Python functions for C code generation
# print_msg_with_timing() gives the user an idea of what's going on/taking so long. Also outputs timing info.
def print_msg_with_timing(desc, msg="Symbolic", startstop="start", starttime=0.0):
    CoordSystem = par.parval_from_str("reference_metric::CoordSystem")
    elapsed = time.time()-starttime
    if msg == "Symbolic":
        if startstop == "start":
            print("Generating symbolic expressions for " + desc + " (%s coords)..." % CoordSystem)
            return time.time()
        else:
            print("Finished generating symbolic expressions for "+desc+
                  " (%s coords) in %.1f seconds. Next up: C codegen..." % (CoordSystem, elapsed))
    elif msg == "Ccodegen":
        if startstop == "start":
            print("Generating C code for "+desc+" (%s coords)..." % CoordSystem)
            return time.time()
        else:
            print("Finished generating C code for "+desc+" (%s coords) in %.1f seconds." % (CoordSystem, elapsed))


# get_loopopts() sets up options for NRPy+'s loop module
def get_loopopts(points_to_update, enable_SIMD, enable_rfm_precompute, OMP_pragma_on, enable_xxs=True):
    loopopts = points_to_update + ",includebraces=False"
    if enable_SIMD:
        loopopts += ",enable_SIMD"
    if enable_rfm_precompute:
        loopopts += ",enable_rfm_precompute"
    elif not enable_xxs:
        pass
    else:
        loopopts += ",Read_xxs"
    if OMP_pragma_on != "i2":
        loopopts += ",pragma_on_"+OMP_pragma_on
    return loopopts


# register_stress_energy_source_terms_return_T4UU() registers gridfunctions
#        for T4UU if needed and not yet registered.
def register_stress_energy_source_terms_return_T4UU(enable_stress_energy_source_terms):
    if enable_stress_energy_source_terms:
        registered_already = False
        for i in range(len(gri.glb_gridfcs_list)):
            if gri.glb_gridfcs_list[i].name == "T4UU00":
                registered_already = True
        if not registered_already:
            return ixp.register_gridfunctions_for_single_rank2("AUXEVOL", "T4UU", "sym01", DIM=4)
        else:
            return ixp.declarerank2("T4UU", "sym01", DIM=4)
    return None

def EinsteinToolkit_keep_param__return_type(paramtuple):
    keep_param = True # We'll not set some parameters in param.ccl;
                      #   e.g., those that should be #define'd like M_PI.
    typestring = ""
    # Separate thorns within the ETK take care of grid/coordinate parameters;
    #   thus we ignore NRPy+ grid/coordinate parameters:
    if paramtuple.module == "grid" or paramtuple.module == "reference_metric":
        keep_param = False

    partype = paramtuple.type
    if partype == "bool":
        typestring += "BOOLEAN "
    elif partype == "REAL":
        if paramtuple.defaultval != 1e300: # 1e300 is a magic value indicating that the C parameter should be mutable
            typestring += "CCTK_REAL "
        else:
            keep_param = False
    elif partype == "int":
        typestring += "CCTK_INT "
    elif partype == "#define":
        keep_param = False
    elif partype == "char":
        # FIXME: char/string parameter types should in principle be supported
        print("Error: parameter "+paramtuple.module+"::"+paramtuple.parname+
              " has unsupported type: \""+ paramtuple.type + "\"")
        sys.exit(1)
    else:
        print("Error: parameter "+paramtuple.module+"::"+paramtuple.parname+
              " has unsupported type: \""+ paramtuple.type + "\"")
        sys.exit(1)
    return keep_param, typestring


def EinsteinToolkit_declare_loop_params():
    return """
  if(cctk_nghostzones[0] != cctk_nghostzones[1] ||
     cctk_nghostzones[0] != cctk_nghostzones[2]) {
    CCTK_ERROR("cctk_nghostzones[i] must be the same in all directions i");
  }

  const CCTK_INT NGHOSTS CCTK_ATTRIBUTE_UNUSED = cctk_nghostzones[0];

  const CCTK_INT Nxx0 CCTK_ATTRIBUTE_UNUSED = cctk_lsh[0] - 2*NGHOSTS;
  const CCTK_INT Nxx1 CCTK_ATTRIBUTE_UNUSED = cctk_lsh[1] - 2*NGHOSTS;
  const CCTK_INT Nxx2 CCTK_ATTRIBUTE_UNUSED = cctk_lsh[2] - 2*NGHOSTS;

  const CCTK_INT Nxx_plus_2NGHOSTS0 CCTK_ATTRIBUTE_UNUSED = cctk_lsh[0];
  const CCTK_INT Nxx_plus_2NGHOSTS1 CCTK_ATTRIBUTE_UNUSED = cctk_lsh[1];
  const CCTK_INT Nxx_plus_2NGHOSTS2 CCTK_ATTRIBUTE_UNUSED = cctk_lsh[2];
"""

def EinsteinToolkit_SIMD_declare_C_params(ThornName):
    SIMD_declare_C_params_str = """

  const CCTK_REAL NOSIMDinvdx0 CCTK_ATTRIBUTE_UNUSED = 1.0/CCTK_DELTA_SPACE(0);
  const REAL_SIMD_ARRAY invdx0 CCTK_ATTRIBUTE_UNUSED = ConstSIMD(NOSIMDinvdx0);
  const CCTK_REAL NOSIMDinvdx1 CCTK_ATTRIBUTE_UNUSED = 1.0/CCTK_DELTA_SPACE(1);
  const REAL_SIMD_ARRAY invdx1 CCTK_ATTRIBUTE_UNUSED = ConstSIMD(NOSIMDinvdx1);
  const CCTK_REAL NOSIMDinvdx2 CCTK_ATTRIBUTE_UNUSED = 1.0/CCTK_DELTA_SPACE(2);
  const REAL_SIMD_ARRAY invdx2 CCTK_ATTRIBUTE_UNUSED = ConstSIMD(NOSIMDinvdx2);
"""

    for i in range(len(par.glb_Cparams_list)):
        # keep_param is a boolean indicating whether we should accept or reject
        #    the parameter. singleparstring will contain the string indicating
        #    the variable type.
        keep_param, singleparstring = EinsteinToolkit_keep_param__return_type(par.glb_Cparams_list[i])

        if (keep_param) and ("CCTK_REAL" in singleparstring):
            parname = par.glb_Cparams_list[i].parname
            SIMD_declare_C_params_str += "  const "+singleparstring + "*NOSIMD"+parname+\
                " CCTK_ATTRIBUTE_UNUSED = CCTK_ParameterGet(\""+parname+"\", \""+ThornName+"\", NULL);\n"
            SIMD_declare_C_params_str += "  const REAL_SIMD_ARRAY "+parname+" CCTK_ATTRIBUTE_UNUSED = ConstSIMD(*NOSIMD"+parname+");\n"
    return SIMD_declare_C_params_str


def set_ETK_func_params_preloop(func_name_suffix, enable_SIMD=True):
    params = "CCTK_ARGUMENTS"
    preloop = """  DECLARE_CCTK_ARGUMENTS;"""
    ThornName = "Baikal"
    if "BaikalVacuum" in func_name_suffix:
        ThornName = "BaikalVacuum"
    if enable_SIMD:
        preloop += EinsteinToolkit_SIMD_declare_C_params(ThornName)
    preloop += EinsteinToolkit_declare_loop_params()
    return params, preloop


#################################
#################################
# Generate symbolic BSSN RHS expressions
def BSSN_RHSs__generate_symbolic_expressions(LapseCondition="OnePlusLog",
                                             ShiftCondition="GammaDriving2ndOrder_Covariant",
                                             enable_KreissOliger_dissipation=True,
                                             enable_stress_energy_source_terms=False,
                                             leave_Ricci_symbolic=True):
    ######################################
    # START: GENERATE SYMBOLIC EXPRESSIONS
    starttime = print_msg_with_timing("BSSN_RHSs", msg="Symbolic", startstop="start")

    # Returns None if enable_stress_energy_source_terms==False; otherwise returns symb expressions for T4UU
    T4UU = register_stress_energy_source_terms_return_T4UU(enable_stress_energy_source_terms)

    # Evaluate BSSN RHSs:
    import BSSN.BSSN_quantities as Bq
    par.set_parval_from_str("BSSN.BSSN_quantities::LeaveRicciSymbolic", str(leave_Ricci_symbolic))
    rhs.BSSN_RHSs()

    if enable_stress_energy_source_terms:
        import BSSN.BSSN_stress_energy_source_terms as Bsest
        Bsest.BSSN_source_terms_for_BSSN_RHSs(T4UU)
        rhs.trK_rhs += Bsest.sourceterm_trK_rhs
        for i in range(3):
            # Needed for Gamma-driving shift RHSs:
            rhs.Lambdabar_rhsU[i] += Bsest.sourceterm_Lambdabar_rhsU[i]
            # Needed for BSSN RHSs:
            rhs.lambda_rhsU[i] += Bsest.sourceterm_lambda_rhsU[i]
            for j in range(3):
                rhs.a_rhsDD[i][j] += Bsest.sourceterm_a_rhsDD[i][j]

    par.set_parval_from_str("BSSN.BSSN_gauge_RHSs::LapseEvolutionOption", LapseCondition)
    par.set_parval_from_str("BSSN.BSSN_gauge_RHSs::ShiftEvolutionOption", ShiftCondition)
    gaugerhs.BSSN_gauge_RHSs()  # Can depend on above RHSs
    # Restore BSSN.BSSN_quantities::LeaveRicciSymbolic to False
    par.set_parval_from_str("BSSN.BSSN_quantities::LeaveRicciSymbolic", "False")

    # Add Kreiss-Oliger dissipation to the BSSN RHSs:
    if enable_KreissOliger_dissipation:
        thismodule = "KO_Dissipation"
        diss_strength = par.Cparameters("REAL", thismodule, "diss_strength", 0.1) # *Bq.cf # *Bq.cf*Bq.cf*Bq.cf # cf**1 is found better than cf**4 over the long term.

        alpha_dKOD = ixp.declarerank1("alpha_dKOD")
        cf_dKOD = ixp.declarerank1("cf_dKOD")
        trK_dKOD = ixp.declarerank1("trK_dKOD")
        betU_dKOD = ixp.declarerank2("betU_dKOD", "nosym")
        vetU_dKOD = ixp.declarerank2("vetU_dKOD", "nosym")
        lambdaU_dKOD = ixp.declarerank2("lambdaU_dKOD", "nosym")
        aDD_dKOD = ixp.declarerank3("aDD_dKOD", "sym01")
        hDD_dKOD = ixp.declarerank3("hDD_dKOD", "sym01")
        for k in range(3):
            gaugerhs.alpha_rhs += diss_strength * alpha_dKOD[k] * rfm.ReU[k]  # ReU[k] = 1/scalefactor_orthog_funcform[k]
            rhs.cf_rhs += diss_strength * cf_dKOD[k] * rfm.ReU[k]  # ReU[k] = 1/scalefactor_orthog_funcform[k]
            rhs.trK_rhs += diss_strength * trK_dKOD[k] * rfm.ReU[k]  # ReU[k] = 1/scalefactor_orthog_funcform[k]
            for i in range(3):
                if "2ndOrder" in ShiftCondition:
                    gaugerhs.bet_rhsU[i] += diss_strength * betU_dKOD[i][k] * rfm.ReU[k]  # ReU[k] = 1/scalefactor_orthog_funcform[k]
                gaugerhs.vet_rhsU[i] += diss_strength * vetU_dKOD[i][k] * rfm.ReU[k]  # ReU[k] = 1/scalefactor_orthog_funcform[k]
                rhs.lambda_rhsU[i] += diss_strength * lambdaU_dKOD[i][k] * rfm.ReU[k]  # ReU[k] = 1/scalefactor_orthog_funcform[k]
                for j in range(3):
                    rhs.a_rhsDD[i][j] += diss_strength * aDD_dKOD[i][j][k] * rfm.ReU[k]  # ReU[k] = 1/scalefactor_orthog_funcform[k]
                    rhs.h_rhsDD[i][j] += diss_strength * hDD_dKOD[i][j][k] * rfm.ReU[k]  # ReU[k] = 1/scalefactor_orthog_funcform[k]

    # We use betaU as our upwinding control vector:
    Bq.BSSN_basic_tensors()
    betaU = Bq.betaU

    # END: GENERATE SYMBOLIC EXPRESSIONS
    ######################################

    lhs_names = ["alpha", "cf", "trK"]
    rhs_exprs = [gaugerhs.alpha_rhs, rhs.cf_rhs, rhs.trK_rhs]
    for i in range(3):
        lhs_names.append("betU" + str(i))
        rhs_exprs.append(gaugerhs.bet_rhsU[i])
        lhs_names.append("lambdaU" + str(i))
        rhs_exprs.append(rhs.lambda_rhsU[i])
        lhs_names.append("vetU" + str(i))
        rhs_exprs.append(gaugerhs.vet_rhsU[i])
        for j in range(i, 3):
            lhs_names.append("aDD" + str(i) + str(j))
            rhs_exprs.append(rhs.a_rhsDD[i][j])
            lhs_names.append("hDD" + str(i) + str(j))
            rhs_exprs.append(rhs.h_rhsDD[i][j])

    # Sort the lhss list alphabetically, and rhss to match.
    #   This ensures the RHSs are evaluated in the same order
    #   they're allocated in memory:
    lhs_names, rhs_exprs = [list(x) for x in zip(*sorted(zip(lhs_names, rhs_exprs), key=lambda pair: pair[0]))]

    # Declare the list of lhrh's
    BSSN_RHSs_SymbExpressions = []
    for var in range(len(lhs_names)):
        BSSN_RHSs_SymbExpressions.append(lhrh(lhs=gri.gfaccess("rhs_gfs", lhs_names[var]), rhs=rhs_exprs[var]))

    print_msg_with_timing("BSSN_RHSs", msg="Symbolic", startstop="stop", starttime=starttime)
    return [betaU, BSSN_RHSs_SymbExpressions]


# Register C code rhs_eval() for BSSN RHS expressions
def add_rhs_eval_to_Cfunction_dict(includes=None, rel_path_to_Cparams=os.path.join("."),
                                   enable_rfm_precompute=True, enable_golden_kernels=False,
                                   enable_SIMD=True, enable_split_for_optimizations_doesnt_help=False,
                                   LapseCondition="OnePlusLog", ShiftCondition="GammaDriving2ndOrder_Covariant",
                                   enable_KreissOliger_dissipation=False, enable_stress_energy_source_terms=False,
                                   leave_Ricci_symbolic=True, OMP_pragma_on="i2",
                                   func_name_suffix=""):
    if includes is None:
        includes = []
    if enable_SIMD:
        includes += [os.path.join("SIMD", "SIMD_intrinsics.h")]
    enable_FD_functions = bool(par.parval_from_str("finite_difference::enable_FD_functions"))
    if enable_FD_functions:
        includes += ["finite_difference_functions.h"]

    # Set up the C function for the BSSN RHSs
    desc = "Evaluate the BSSN RHSs"
    name = "rhs_eval" + func_name_suffix
    params = "const paramstruct *restrict params, "
    if enable_rfm_precompute:
        params += "const rfm_struct *restrict rfmstruct, "
    else:
        params += "REAL *xx[3], "
    params += """
              const REAL *restrict auxevol_gfs,const REAL *restrict in_gfs,REAL *restrict rhs_gfs"""

    betaU, BSSN_RHSs_SymbExpressions = \
        BSSN_RHSs__generate_symbolic_expressions(LapseCondition=LapseCondition, ShiftCondition=ShiftCondition,
                                                 enable_KreissOliger_dissipation=enable_KreissOliger_dissipation,
                                                 enable_stress_energy_source_terms=enable_stress_energy_source_terms,
                                                 leave_Ricci_symbolic=leave_Ricci_symbolic)

    # Construct body:
    preloop=""
    enableCparameters=True
    # Set up preloop in case we're outputting code for the Einstein Toolkit (ETK)
    if par.parval_from_str("grid::GridFuncMemAccess") == "ETK":
        params, preloop = set_ETK_func_params_preloop(func_name_suffix)
        enableCparameters=False

    FD_outCparams = "outCverbose=False,enable_SIMD=" + str(enable_SIMD)
    FD_outCparams += ",GoldenKernelsEnable=" + str(enable_golden_kernels)

    loopopts = get_loopopts("InteriorPoints", enable_SIMD, enable_rfm_precompute, OMP_pragma_on)
    FDorder = par.parval_from_str("finite_difference::FD_CENTDERIVS_ORDER")
    starttime = print_msg_with_timing("BSSN_RHSs (FD order="+str(FDorder)+")", msg="Ccodegen", startstop="start")
    if enable_split_for_optimizations_doesnt_help and FDorder == 6:
        loopopts += ",DisableOpenMP"
        BSSN_RHSs_SymbExpressions_pt1 = []
        BSSN_RHSs_SymbExpressions_pt2 = []
        for lhsrhs in BSSN_RHSs_SymbExpressions:
            if "BETU" in lhsrhs.lhs or "LAMBDAU" in lhsrhs.lhs:
                BSSN_RHSs_SymbExpressions_pt1.append(lhrh(lhs=lhsrhs.lhs, rhs=lhsrhs.rhs))
            else:
                BSSN_RHSs_SymbExpressions_pt2.append(lhrh(lhs=lhsrhs.lhs, rhs=lhsrhs.rhs))
        preloop += """#pragma omp parallel
    {
"""
        preloopbody = fin.FD_outputC("returnstring", BSSN_RHSs_SymbExpressions_pt1,
                                     params=FD_outCparams,
                                     upwindcontrolvec=betaU)
        preloop += "\n#pragma omp for\n" + lp.simple_loop(loopopts, preloopbody)
        preloop += "\n#pragma omp for\n"
        body = fin.FD_outputC("returnstring", BSSN_RHSs_SymbExpressions_pt2,
                              params=FD_outCparams,
                              upwindcontrolvec=betaU)
        postloop = "\n    } // END #pragma omp parallel\n"
    else:
        preloop += ""
        body = fin.FD_outputC("returnstring", BSSN_RHSs_SymbExpressions,
                              params=FD_outCparams,
                              upwindcontrolvec=betaU)
        postloop = ""
    print_msg_with_timing("BSSN_RHSs (FD order="+str(FDorder)+")", msg="Ccodegen", startstop="stop", starttime=starttime)

    add_to_Cfunction_dict(
        includes=includes,
        desc=desc,
        name=name, params=params,
        preloop=preloop, body=body, loopopts=loopopts, postloop=postloop,
        rel_path_to_Cparams=rel_path_to_Cparams, enableCparameters=enableCparameters)
    return pickle_NRPy_env()


# Generate symbolic expressions for 3-Ricci tensor Rbar[i][j]
def Ricci__generate_symbolic_expressions():
    ######################################
    # START: GENERATE SYMBOLIC EXPRESSIONS
    starttime = print_msg_with_timing("3-Ricci tensor", msg="Symbolic", startstop="start")

    # Evaluate 3-Ricci tensor:
    import BSSN.BSSN_quantities as Bq
    par.set_parval_from_str("BSSN.BSSN_quantities::LeaveRicciSymbolic", "False")

    # Register all BSSN gridfunctions if not registered already
    Bq.BSSN_basic_tensors()
    # Next compute Ricci tensor
    Bq.RicciBar__gammabarDD_dHatD__DGammaUDD__DGammaU()
    # END: GENERATE SYMBOLIC EXPRESSIONS
    ######################################
    # Must register RbarDD as gridfunctions, as we're outputting them to gridfunctions here:
    foundit = False
    for i in range(len(gri.glb_gridfcs_list)):
        if "RbarDD00" in gri.glb_gridfcs_list[i].name:
            foundit = True
    if not foundit:
        ixp.register_gridfunctions_for_single_rank2("AUXEVOL", "RbarDD", "sym01")

    Ricci_SymbExpressions = [lhrh(lhs=gri.gfaccess("auxevol_gfs", "RbarDD00"), rhs=Bq.RbarDD[0][0]),
                             lhrh(lhs=gri.gfaccess("auxevol_gfs", "RbarDD01"), rhs=Bq.RbarDD[0][1]),
                             lhrh(lhs=gri.gfaccess("auxevol_gfs", "RbarDD02"), rhs=Bq.RbarDD[0][2]),
                             lhrh(lhs=gri.gfaccess("auxevol_gfs", "RbarDD11"), rhs=Bq.RbarDD[1][1]),
                             lhrh(lhs=gri.gfaccess("auxevol_gfs", "RbarDD12"), rhs=Bq.RbarDD[1][2]),
                             lhrh(lhs=gri.gfaccess("auxevol_gfs", "RbarDD22"), rhs=Bq.RbarDD[2][2])]
    print_msg_with_timing("3-Ricci tensor", msg="Symbolic", startstop="stop", starttime=starttime)
    return Ricci_SymbExpressions


# Register C function Ricci_eval() for evaluating 3-Ricci tensor
def add_Ricci_eval_to_Cfunction_dict(includes=None, rel_path_to_Cparams=os.path.join("."),
                                     enable_rfm_precompute=True, enable_golden_kernels=False, enable_SIMD=True,
                                     enable_split_for_optimizations_doesnt_help=False, OMP_pragma_on="i2",
                                     func_name_suffix=""):
    if includes is None:
        includes = []
    if enable_SIMD:
        includes += [os.path.join("SIMD", "SIMD_intrinsics.h")]
    enable_FD_functions = bool(par.parval_from_str("finite_difference::enable_FD_functions"))
    if enable_FD_functions:
        includes += ["finite_difference_functions.h"]

    # Set up the C function for the 3-Ricci tensor
    desc = "Evaluate the 3-Ricci tensor"
    name = "Ricci_eval" + func_name_suffix
    params = "const paramstruct *restrict params, "
    if enable_rfm_precompute:
        params += "const rfm_struct *restrict rfmstruct, "
    else:
        params += "REAL *xx[3], "
    params += "const REAL *restrict in_gfs, REAL *restrict auxevol_gfs"

    # Construct body:
    Ricci_SymbExpressions = Ricci__generate_symbolic_expressions()
    FD_outCparams = "outCverbose=False,enable_SIMD=" + str(enable_SIMD)
    FD_outCparams += ",GoldenKernelsEnable=" + str(enable_golden_kernels)
    loopopts = get_loopopts("InteriorPoints", enable_SIMD, enable_rfm_precompute, OMP_pragma_on)

    FDorder = par.parval_from_str("finite_difference::FD_CENTDERIVS_ORDER")
    starttime = print_msg_with_timing("3-Ricci tensor (FD order="+str(FDorder)+")", msg="Ccodegen", startstop="start")

    # Construct body:
    preloop=""
    enableCparameters=True
    # Set up preloop in case we're outputting code for the Einstein Toolkit (ETK)
    if par.parval_from_str("grid::GridFuncMemAccess") == "ETK":
        params, preloop = set_ETK_func_params_preloop(func_name_suffix)
        enableCparameters=False

    if enable_split_for_optimizations_doesnt_help and FDorder >= 8:
        loopopts += ",DisableOpenMP"
        Ricci_SymbExpressions_pt1 = []
        Ricci_SymbExpressions_pt2 = []
        for lhsrhs in Ricci_SymbExpressions:
            if "RBARDD00" in lhsrhs.lhs or "RBARDD11" in lhsrhs.lhs or "RBARDD22" in lhsrhs.lhs:
                Ricci_SymbExpressions_pt1.append(lhrh(lhs=lhsrhs.lhs, rhs=lhsrhs.rhs))
            else:
                Ricci_SymbExpressions_pt2.append(lhrh(lhs=lhsrhs.lhs, rhs=lhsrhs.rhs))
        preloop = """#pragma omp parallel
    {
#pragma omp for
"""
        preloopbody = fin.FD_outputC("returnstring", Ricci_SymbExpressions_pt1,
                                     params=FD_outCparams)
        preloop += lp.simple_loop(loopopts, preloopbody)
        preloop += "#pragma omp for\n"
        body = fin.FD_outputC("returnstring", Ricci_SymbExpressions_pt2,
                              params=FD_outCparams)
        postloop = "\n    } // END #pragma omp parallel\n"
    else:
        body = fin.FD_outputC("returnstring", Ricci_SymbExpressions,
                              params=FD_outCparams)
        postloop = ""
    print_msg_with_timing("3-Ricci tensor (FD order="+str(FDorder)+")", msg="Ccodegen", startstop="stop", starttime=starttime)

    add_to_Cfunction_dict(
        includes=includes,
        desc=desc,
        name=name, params=params,
        preloop=preloop, body=body, loopopts=loopopts, postloop=postloop,
        rel_path_to_Cparams=rel_path_to_Cparams, enableCparameters=enableCparameters)
    return pickle_NRPy_env()


# Generate symbolic expressions for BSSN Hamiltonian & momentum constraints
def BSSN_constraints__generate_symbolic_expressions(enable_stress_energy_source_terms=False, leave_Ricci_symbolic=True,
                                                    output_H_only=False):
    ######################################
    # START: GENERATE SYMBOLIC EXPRESSIONS
    starttime = print_msg_with_timing("BSSN constraints", msg="Symbolic", startstop="start")

    # Define the Hamiltonian constraint and output the optimized C code.
    par.set_parval_from_str("BSSN.BSSN_quantities::LeaveRicciSymbolic", str(leave_Ricci_symbolic))
    import BSSN.BSSN_constraints as bssncon

    # Returns None if enable_stress_energy_source_terms==False; otherwise returns symb expressions for T4UU
    T4UU = register_stress_energy_source_terms_return_T4UU(enable_stress_energy_source_terms)

    bssncon.BSSN_constraints(add_T4UUmunu_source_terms=False, output_H_only=output_H_only)  # We'll add them below if desired.
    if enable_stress_energy_source_terms:
        import BSSN.BSSN_stress_energy_source_terms as Bsest
        Bsest.BSSN_source_terms_for_BSSN_constraints(T4UU)
        bssncon.H += Bsest.sourceterm_H
        for i in range(3):
            bssncon.MU[i] += Bsest.sourceterm_MU[i]

    BSSN_constraints_SymbExpressions = [lhrh(lhs=gri.gfaccess("aux_gfs", "H"), rhs=bssncon.H)]
    if not output_H_only:
        BSSN_constraints_SymbExpressions += [lhrh(lhs=gri.gfaccess("aux_gfs", "MU0"), rhs=bssncon.MU[0]),
                                             lhrh(lhs=gri.gfaccess("aux_gfs", "MU1"), rhs=bssncon.MU[1]),
                                             lhrh(lhs=gri.gfaccess("aux_gfs", "MU2"), rhs=bssncon.MU[2])]
    par.set_parval_from_str("BSSN.BSSN_quantities::LeaveRicciSymbolic", "False")
    print_msg_with_timing("BSSN constraints", msg="Symbolic", startstop="stop", starttime=starttime)
    # END: GENERATE SYMBOLIC EXPRESSIONS
    ######################################
    return BSSN_constraints_SymbExpressions


# Register C function BSSN_constraints() for evaluating BSSN Hamiltonian & momentum constraints
def add_BSSN_constraints_to_Cfunction_dict(includes=None, rel_path_to_Cparams=os.path.join("."),
                                           enable_rfm_precompute=True, enable_golden_kernels=False, enable_SIMD=True,
                                           enable_stress_energy_source_terms=False, leave_Ricci_symbolic=True,
                                           output_H_only=False, OMP_pragma_on="i2", func_name_suffix=""):
    if includes is None:
        includes = []
    if enable_SIMD:
        includes += [os.path.join("SIMD", "SIMD_intrinsics.h")]
    enable_FD_functions = bool(par.parval_from_str("finite_difference::enable_FD_functions"))
    if enable_FD_functions:
        includes += ["finite_difference_functions.h"]

    # Set up the C function for the BSSN constraints
    desc = "Evaluate the BSSN constraints"
    name = "BSSN_constraints" + func_name_suffix
    params = "const paramstruct *restrict params, "
    if enable_rfm_precompute:
        params += "const rfm_struct *restrict rfmstruct, "
    else:
        params += "REAL *xx[3], "
    params += """
                 const REAL *restrict in_gfs, const REAL *restrict auxevol_gfs, REAL *restrict aux_gfs"""

    # Construct body:

    BSSN_constraints_SymbExpressions = BSSN_constraints__generate_symbolic_expressions(enable_stress_energy_source_terms,
                                                                                       leave_Ricci_symbolic=leave_Ricci_symbolic,
                                                                                       output_H_only=output_H_only)

    preloop=""
    enableCparameters=True
    # Set up preloop in case we're outputting code for the Einstein Toolkit (ETK)
    if par.parval_from_str("grid::GridFuncMemAccess") == "ETK":
        params, preloop = set_ETK_func_params_preloop(func_name_suffix)
        enableCparameters=False

    FD_outCparams = "outCverbose=False,enable_SIMD=" + str(enable_SIMD)
    FD_outCparams += ",GoldenKernelsEnable=" + str(enable_golden_kernels)
    FDorder = par.parval_from_str("finite_difference::FD_CENTDERIVS_ORDER")
    starttime = print_msg_with_timing("BSSN constraints (FD order="+str(FDorder)+")", msg="Ccodegen", startstop="start")
    body = fin.FD_outputC("returnstring", BSSN_constraints_SymbExpressions,
                          params=FD_outCparams)
    print_msg_with_timing("BSSN constraints (FD order="+str(FDorder)+")", msg="Ccodegen", startstop="stop", starttime=starttime)

    add_to_Cfunction_dict(
        includes=includes,
        desc=desc,
        name=name, params=params,
        preloop=preloop,
        body=body,
        loopopts=get_loopopts("InteriorPoints", enable_SIMD, enable_rfm_precompute, OMP_pragma_on),
        rel_path_to_Cparams=rel_path_to_Cparams, enableCparameters=enableCparameters)
    return pickle_NRPy_env()


# Register C function enforce_detgammahat_constraint for enforcing the conformal 3-metric: det(gammahat)=det(gammabar)
def add_enforce_detgammahat_constraint_to_Cfunction_dict(includes=None, rel_path_to_Cparams=os.path.join("."),
                                                         enable_rfm_precompute=True, enable_golden_kernels=False,
                                                         OMP_pragma_on="i2", func_name_suffix=""):
    # This function disables SIMD, as it includes cbrt() and abs() functions.
    if includes is None:
        includes = []
    # This function does not use finite differences!
    # enable_FD_functions = bool(par.parval_from_str("finite_difference::enable_FD_functions"))
    # if enable_FD_functions:
    #     includes += ["finite_difference_functions.h"]

    # Set up the C function for enforcing the det(gammabar) = det(gammahat) BSSN algebraic constraint
    desc = "Enforce the det(gammabar) = det(gammahat) (algebraic) constraint"
    name = "enforce_detgammahat_constraint" + func_name_suffix
    params = "const paramstruct *restrict params, "
    if enable_rfm_precompute:
        params += "const rfm_struct *restrict rfmstruct, "
    else:
        params += "REAL *xx[3], "
    params += "REAL *restrict in_gfs"

    # Construct body:
    enforce_detg_constraint_symb_expressions = EGC.Enforce_Detgammahat_Constraint_symb_expressions()

    preloop=""
    enableCparameters=True
    # Set up preloop in case we're outputting code for the Einstein Toolkit (ETK)
    if par.parval_from_str("grid::GridFuncMemAccess") == "ETK":
        params, preloop = set_ETK_func_params_preloop(func_name_suffix, enable_SIMD=False)
        enableCparameters=False

    FD_outCparams = "outCverbose=False,enable_SIMD=False"
    FD_outCparams += ",GoldenKernelsEnable=" + str(enable_golden_kernels)
    starttime = print_msg_with_timing("Enforcing det(gammabar)=det(gammahat) constraint", msg="Ccodegen", startstop="start")
    body = fin.FD_outputC("returnstring", enforce_detg_constraint_symb_expressions,
                          params=FD_outCparams)
    print_msg_with_timing("Enforcing det(gammabar)=det(gammahat) constraint", msg="Ccodegen", startstop="stop", starttime=starttime)

    enable_SIMD = False
    add_to_Cfunction_dict(
        includes=includes,
        desc=desc,
        name=name, params=params,
        preloop=preloop,
        body=body,
        loopopts=get_loopopts("AllPoints", enable_SIMD, enable_rfm_precompute, OMP_pragma_on),
        rel_path_to_Cparams=rel_path_to_Cparams, enableCparameters=enableCparameters)
    return pickle_NRPy_env()


# Add psi4 to Cfunction dictionary. psi4 is a really huge expression, so we output
#   it in 3 parts: psi4_part0, psi4_part1, and psi4_part2
def add_psi4_part_to_Cfunction_dict(includes=None, rel_path_to_Cparams=os.path.join("."), whichpart=0,
                                    setPsi4tozero=False, OMP_pragma_on="i2"):
    starttime = print_msg_with_timing("psi4, part " + str(whichpart), msg="Ccodegen", startstop="start")

    # Set up the C function for psi4
    if includes is None:
        includes = []
    includes += ["NRPy_function_prototypes.h"]

    desc = "Compute psi4 at all interior gridpoints, part " + str(whichpart)
    name = "psi4_part" + str(whichpart)
    params = """const paramstruct *restrict params, const REAL *restrict in_gfs, REAL *restrict xx[3], REAL *restrict aux_gfs"""

    body = ""
    gri.register_gridfunctions("AUX", ["psi4_part" + str(whichpart) + "re", "psi4_part" + str(whichpart) + "im"])
    FD_outCparams = "outCverbose=False,enable_SIMD=False,CSE_sorting=none"
    if not setPsi4tozero:
        # Set the body of the function
        # First compute the symbolic expressions
        psi4.Psi4(specify_tetrad=False)

        # We really don't want to store these "Cparameters" permanently; they'll be set via function call...
        #   so we make a copy of the original par.glb_Cparams_list (sans tetrad vectors) and restore it below
        Cparams_list_orig = par.glb_Cparams_list.copy()
        par.Cparameters("REAL", __name__, ["mre4U0", "mre4U1", "mre4U2", "mre4U3"], [0, 0, 0, 0])
        par.Cparameters("REAL", __name__, ["mim4U0", "mim4U1", "mim4U2", "mim4U3"], [0, 0, 0, 0])
        par.Cparameters("REAL", __name__, ["n4U0", "n4U1", "n4U2", "n4U3"], [0, 0, 0, 0])

        body += """
REAL mre4U0,mre4U1,mre4U2,mre4U3,mim4U0,mim4U1,mim4U2,mim4U3,n4U0,n4U1,n4U2,n4U3;
psi4_tetrad(params,
              in_gfs[IDX4S(CFGF, i0,i1,i2)],
              in_gfs[IDX4S(HDD00GF, i0,i1,i2)],
              in_gfs[IDX4S(HDD01GF, i0,i1,i2)],
              in_gfs[IDX4S(HDD02GF, i0,i1,i2)],
              in_gfs[IDX4S(HDD11GF, i0,i1,i2)],
              in_gfs[IDX4S(HDD12GF, i0,i1,i2)],
              in_gfs[IDX4S(HDD22GF, i0,i1,i2)],
              &mre4U0,&mre4U1,&mre4U2,&mre4U3,&mim4U0,&mim4U1,&mim4U2,&mim4U3,&n4U0,&n4U1,&n4U2,&n4U3,
              xx, i0,i1,i2);
"""
        body += "REAL xCart_rel_to_globalgrid_center[3];\n"
        body += "xx_to_Cart(params, xx, i0, i1, i2,  xCart_rel_to_globalgrid_center);\n"
        body += "int ignore_Cart_to_i0i1i2[3];  REAL xx_rel_to_globalgridorigin[3];\n"
        body += "Cart_to_xx_and_nearest_i0i1i2_global_grid_center(params, xCart_rel_to_globalgrid_center,xx_rel_to_globalgridorigin,ignore_Cart_to_i0i1i2);\n"
        for i in range(3):
            body += "const REAL xx" + str(i) + "=xx_rel_to_globalgridorigin[" + str(i) + "];\n"
        body += fin.FD_outputC("returnstring",
                               [lhrh(lhs=gri.gfaccess("in_gfs", "psi4_part" + str(whichpart) + "re"),
                                     rhs=psi4.psi4_re_pt[whichpart]),
                                lhrh(lhs=gri.gfaccess("in_gfs", "psi4_part" + str(whichpart) + "im"),
                                     rhs=psi4.psi4_im_pt[whichpart])],
                               params=FD_outCparams)
        par.glb_Cparams_list = Cparams_list_orig.copy()

    elif setPsi4tozero:
        body += fin.FD_outputC("returnstring",
                               [lhrh(lhs=gri.gfaccess("in_gfs", "psi4_part" + str(whichpart) + "re"),
                                     rhs=sp.sympify(0)),
                                lhrh(lhs=gri.gfaccess("in_gfs", "psi4_part" + str(whichpart) + "im"),
                                     rhs=sp.sympify(0))],
                               params=FD_outCparams)
    enable_SIMD = False
    enable_rfm_precompute = False

    print_msg_with_timing("psi4, part " + str(whichpart), msg="Ccodegen", startstop="stop", starttime=starttime)
    add_to_Cfunction_dict(
        includes=includes,
        desc=desc,
        name=name, params=params,
        body=body,
        loopopts=get_loopopts("InteriorPoints", enable_SIMD, enable_rfm_precompute, OMP_pragma_on,
                              enable_xxs=False),
        rel_path_to_Cparams=rel_path_to_Cparams)
    return pickle_NRPy_env()


def add_psi4_tetrad_to_Cfunction_dict(includes=None, rel_path_to_Cparams=os.path.join("."), setPsi4tozero=False):
    starttime = print_msg_with_timing("psi4 tetrads", msg="Ccodegen", startstop="start")

    # Set up the C function for BSSN basis transformations
    desc = "Compute tetrad for psi4"
    name = "psi4_tetrad"

    # First set up the symbolic expressions (RHSs) and their names (LHSs)
    psi4tet.Psi4_tetrads()
    list_of_varnames = []
    list_of_symbvars = []
    for i in range(4):
        list_of_varnames.append("*mre4U" + str(i))
        list_of_symbvars.append(psi4tet.mre4U[i])
    for i in range(4):
        list_of_varnames.append("*mim4U" + str(i))
        list_of_symbvars.append(psi4tet.mim4U[i])
    for i in range(4):
        list_of_varnames.append("*n4U" + str(i))
        list_of_symbvars.append(psi4tet.n4U[i])

    paramsindent = "                  "
    params = """const paramstruct *restrict params,\n""" + paramsindent
    list_of_metricvarnames = ["cf"]
    for i in range(3):
        for j in range(i, 3):
            list_of_metricvarnames.append("hDD" + str(i) + str(j))
    for var in list_of_metricvarnames:
        params += "const REAL " + var + ","
    params += "\n" + paramsindent
    for var in list_of_varnames:
        params += "REAL " + var + ","
    params += "\n" + paramsindent + "REAL *restrict xx[3], const int i0,const int i1,const int i2"

    # Set the body of the function
    body = ""
    outCparams = "includebraces=False,outCverbose=False,enable_SIMD=False,preindent=1"
    if not setPsi4tozero:
        for i in range(3):
            body += "  const REAL xx" + str(i) + " = xx[" + str(i) + "][i" + str(i) + "];\n"
        body += "  // Compute tetrads:\n"
        body += "  {\n"
        # Sort the lhss list alphabetically, and rhss to match:
        lhss, rhss = [list(x) for x in zip(*sorted(zip(list_of_varnames, list_of_symbvars), key=lambda pair: pair[0]))]
        body += outputC(rhss, lhss, filename="returnstring", params=outCparams)
        body += "  }\n"

    elif setPsi4tozero:
        body += "return;\n"
    loopopts = ""

    print_msg_with_timing("psi4 tetrads", msg="Ccodegen", startstop="stop", starttime=starttime)
    add_to_Cfunction_dict(
        includes=includes,
        desc=desc,
        name=name, params=params,
        body=body,
        loopopts=loopopts,
        rel_path_to_Cparams=rel_path_to_Cparams)
    return pickle_NRPy_env()


def add_SpinWeight_minus2_SphHarmonics_to_Cfunction_dict(includes=None, rel_path_to_Cparams=os.path.join("."),
                                                         maximum_l=8):
    starttime = print_msg_with_timing("Spin-weight s=-2 Spherical Harmonics", msg="Ccodegen", startstop="start")

    # Set up the C function for computing the spin-weight -2 spherical harmonic at theta,phi: Y_{s=-2, l,m}(theta,phi)
    prefunc = r"""// Compute at a single point (th,ph) the spin-weight -2 spherical harmonic Y_{s=-2, l,m}(th,ph)
// Manual "inline void" of this function results in compilation error with clang.
void SpinWeight_minus2_SphHarmonics(const int l, const int m, const REAL th, const REAL ph,
                                           REAL *reYlmswm2_l_m, REAL *imYlmswm2_l_m) {
"""
    # Construct prefunc:
    outCparams = "preindent=1,outCfileaccess=a,outCverbose=False,includebraces=False"
    prefunc += """
  switch(l) {
"""
    for l in range(maximum_l + 1):  # Output values up to and including l=8.
        prefunc += "  case " + str(l) + ":\n"
        prefunc += "    switch(m) {\n"
        for m in range(-l, l + 1):
            prefunc += "    case " + str(m) + ":\n"
            prefunc += "      {\n"
            Y_m2_lm = SWm2SH.Y(-2, l, m, SWm2SH.th, SWm2SH.ph)
            prefunc += outputC([sp.re(Y_m2_lm), sp.im(Y_m2_lm)], ["*reYlmswm2_l_m", "*imYlmswm2_l_m"],
                            "returnstring", outCparams)
            prefunc += "      }\n"
            prefunc += "      return;\n"
        prefunc += "    }  // END switch(m)\n"
    prefunc += "  } // END switch(l)\n"

    prefunc += r"""
  fprintf(stderr, "ERROR: SpinWeight_minus2_SphHarmonics handles only l=[0,"""+str(maximum_l)+r"""] and only m=[-l,+l] is defined.\n");
  fprintf(stderr, "       You chose l=%d and m=%d, which is out of these bounds.\n",l,m);
  exit(1);
}

void lowlevel_decompose_psi4_into_swm2_modes(const int Nxx_plus_2NGHOSTS1,const int Nxx_plus_2NGHOSTS2,
                                             const REAL dxx1, const REAL dxx2,
                                             const REAL curr_time, const REAL R_ext,
                                             const REAL *restrict th_array, const REAL *restrict sinth_array, const REAL *restrict ph_array,
                                             const REAL *restrict psi4r_at_R_ext, const REAL *restrict psi4i_at_R_ext) {
  for(int l=2;l<="""+str(maximum_l)+r""";l++) {  // The maximum l here is set in Python.
    for(int m=-l;m<=l;m++) {
      // Parallelize the integration loop:
      REAL psi4r_l_m = 0.0;
      REAL psi4i_l_m = 0.0;
#pragma omp parallel for reduction(+:psi4r_l_m,psi4i_l_m)
      for(int i1=0;i1<Nxx_plus_2NGHOSTS1-2*NGHOSTS;i1++) {
        const REAL th    = th_array[i1];
        const REAL sinth = sinth_array[i1];
        for(int i2=0;i2<Nxx_plus_2NGHOSTS2-2*NGHOSTS;i2++) {
          const REAL ph = ph_array[i2];
          // Construct integrand for psi4 spin-weight s=-2,l=2,m=0 spherical harmonic
          REAL ReY_sm2_l_m,ImY_sm2_l_m;
          SpinWeight_minus2_SphHarmonics(l,m, th,ph,  &ReY_sm2_l_m,&ImY_sm2_l_m);

          const int idx2d = i1*(Nxx_plus_2NGHOSTS2-2*NGHOSTS)+i2;
          const REAL a = psi4r_at_R_ext[idx2d];
          const REAL b = psi4i_at_R_ext[idx2d];
          const REAL c = ReY_sm2_l_m;
          const REAL d = ImY_sm2_l_m;
          psi4r_l_m += (a*c + b*d) * dxx2  * sinth*dxx1;
          psi4i_l_m += (b*c - a*d) * dxx2  * sinth*dxx1;
        }
      }
      // Step 4: Output the result of the integration to file.
      char filename[100];
      sprintf(filename,"outpsi4_l%d_m%d-r%.2f.txt",l,m, (double)R_ext);
      // If you love "+"'s in filenames by all means enable this (ugh):
      //if(m>=0) sprintf(filename,"outpsi4_l%d_m+%d-r%.2f.txt",l,m, (double)R_ext);
      FILE *outpsi4_l_m;
      // 0 = n*dt when n=0 is exactly represented in double/long double precision,
      //          so no worries about the result being ~1e-16 in double/ld precision
      if(curr_time==0) outpsi4_l_m = fopen(filename, "w");
      else             outpsi4_l_m = fopen(filename, "a");
      fprintf(outpsi4_l_m,"%e %.15e %.15e\n", (double)(curr_time),
              (double)psi4r_l_m,(double)psi4i_l_m);
      fclose(outpsi4_l_m);
    }
  }
}
"""

    desc = ""
    name = "driver__spherlikegrids__psi4_spinweightm2_decomposition"
    params = r"""const paramstruct *restrict params,  REAL *restrict diagnostic_output_gfs,
        const int *restrict list_of_R_ext_idxs, const int num_of_R_ext_idxs, const REAL time,
        REAL *restrict xx[3],void xx_to_Cart(const paramstruct *restrict params, REAL *restrict xx[3],const int i0,const int i1,const int i2, REAL xCart[3])"""

    body = r"""  // Step 1: Allocate memory for 2D arrays used to store psi4, theta, sin(theta), and phi.
  const int sizeof_2Darray = sizeof(REAL)*(Nxx_plus_2NGHOSTS1-2*NGHOSTS)*(Nxx_plus_2NGHOSTS2-2*NGHOSTS);
  REAL *restrict psi4r_at_R_ext = (REAL *restrict)malloc(sizeof_2Darray);
  REAL *restrict psi4i_at_R_ext = (REAL *restrict)malloc(sizeof_2Darray);
  //         ... also store theta, sin(theta), and phi to corresponding 1D arrays.
  REAL *restrict sinth_array = (REAL *restrict)malloc(sizeof(REAL)*(Nxx_plus_2NGHOSTS1-2*NGHOSTS));
  REAL *restrict th_array    = (REAL *restrict)malloc(sizeof(REAL)*(Nxx_plus_2NGHOSTS1-2*NGHOSTS));
  REAL *restrict ph_array    = (REAL *restrict)malloc(sizeof(REAL)*(Nxx_plus_2NGHOSTS2-2*NGHOSTS));

  // Step 2: Loop over all extraction indices:
  for(int ii0=0;ii0<num_of_R_ext_idxs;ii0++) {
    // Step 2.a: Set the extraction radius R_ext based on the radial index R_ext_idx
    REAL R_ext;
    {
      REAL xCart[3];  xx_to_Cart(params,xx,list_of_R_ext_idxs[ii0],1,1,xCart); // values for itheta and iphi don't matter.
      R_ext = sqrt(xCart[0]*xCart[0] + xCart[1]*xCart[1] + xCart[2]*xCart[2]);
    }

    // Step 2.b: Compute psi_4 at this extraction radius and store to a local 2D array.
    const int i0=list_of_R_ext_idxs[ii0];
#pragma omp parallel for
    for(int i1=NGHOSTS;i1<Nxx_plus_2NGHOSTS1-NGHOSTS;i1++) {
      th_array[i1-NGHOSTS]    =     xx[1][i1];
      sinth_array[i1-NGHOSTS] = sin(xx[1][i1]);
      for(int i2=NGHOSTS;i2<Nxx_plus_2NGHOSTS2-NGHOSTS;i2++) {
        ph_array[i2-NGHOSTS] = xx[2][i2];

        // Compute real & imaginary parts of psi_4, output to diagnostic_output_gfs
        const REAL psi4r = (diagnostic_output_gfs[IDX4S(PSI4_PART0REGF, i0,i1,i2)] +
                            diagnostic_output_gfs[IDX4S(PSI4_PART1REGF, i0,i1,i2)] +
                            diagnostic_output_gfs[IDX4S(PSI4_PART2REGF, i0,i1,i2)]);
        const REAL psi4i = (diagnostic_output_gfs[IDX4S(PSI4_PART0IMGF, i0,i1,i2)] +
                            diagnostic_output_gfs[IDX4S(PSI4_PART1IMGF, i0,i1,i2)] +
                            diagnostic_output_gfs[IDX4S(PSI4_PART2IMGF, i0,i1,i2)]);

        // Store result to "2D" array (actually 1D array with 2D storage):
        const int idx2d = (i1-NGHOSTS)*(Nxx_plus_2NGHOSTS2-2*NGHOSTS)+(i2-NGHOSTS);
        psi4r_at_R_ext[idx2d] = psi4r;
        psi4i_at_R_ext[idx2d] = psi4i;
      }
    }
    // Step 3: Perform integrations across all l,m modes from l=2 up to and including L_MAX (global variable):
    lowlevel_decompose_psi4_into_swm2_modes(Nxx_plus_2NGHOSTS1,Nxx_plus_2NGHOSTS2,
                                            dxx1,dxx2,
                                            time, R_ext, th_array, sinth_array, ph_array,
                                            psi4r_at_R_ext,psi4i_at_R_ext);
  }

  // Step 4: Free all allocated memory:
  free(psi4r_at_R_ext); free(psi4i_at_R_ext);
  free(sinth_array); free(th_array); free(ph_array);
"""

    print_msg_with_timing("Spin-weight s=-2 Spherical Harmonics", msg="Ccodegen", startstop="stop", starttime=starttime)

    add_to_Cfunction_dict(
        includes=includes,
        prefunc=prefunc,
        desc=desc,
        name=name, params=params,
        body=body,
        rel_path_to_Cparams=rel_path_to_Cparams)
    return pickle_NRPy_env()
