# BaikalETK: NRPy+-Based BSSN Solver for the Einstein Toolkit

# As documented in the NRPy+ tutorial module
#   Tutorial-BaikalETK.ipynb
#   this module will generate an Einstein
#   Toolkit thorn for solving the BSSN
#   equations with or without T^{mu nu} source
#   terms.

# Author: Zachariah B. Etienne
#         zachetie **at** gmail **dot* com


# Step 1: Import needed core NRPy+ modules
from outputC import *            # NRPy+: Core C code output module
import finite_difference as fin  # NRPy+: Finite difference C code generation module
import NRPy_param_funcs as par   # NRPy+: Parameter interface
import grid as gri               # NRPy+: Functions having to do with numerical grids
import loop as lp                # NRPy+: Generate C code loops
import indexedexp as ixp         # NRPy+: Symbolic indexed expression (e.g., tensors, vectors, etc.) support
import reference_metric as rfm   # NRPy+: Reference metric support
import cmdline_helper as cmd     # NRPy+: Multi-platform Python command-line interface
import shutil, os, sys           # Standard Python modules for multiplatform OS-level functions

def BaikalETK_codegen(outrootdir = "BaikalETK/",
              FD_order=4, # Finite difference order: even numbers only, starting with 2. 12 is generally unstable
              LapseCondition  = "OnePlusLog", # Set the standard 1+log lapse condition
              ShiftCondition  = "GammaDriving2ndOrder_NoCovariant", # Set the standard, second-order advecting-shift,
                                                                    # Gamma-driving shift condition
              add_stress_energy_source_terms = False,  # Enable stress-energy terms?
              default_KO_strength = 0.1 # default Kreiss-Oliger dissipation strength; adjustable within thorn's param.ccl.
              ):
    # Create directory for BaikalETK thorn & subdirectories in case they don't exist.
    cmd.mkdir(os.path.join(outrootdir))
    outdir = os.path.join(outrootdir,"src") # Main C code output directory
    cmd.mkdir(outdir)

    # Set spatial dimension (must be 3 for BSSN)
    DIM = 3
    par.set_parval_from_str("grid::DIM",DIM)

    # Step 2: Set some core parameters, including CoordSystem MoL timestepping algorithm,
    #                                 FD order, floating point precision, and CFL factor:
    # Choices are: Spherical, SinhSpherical, SinhSphericalv2, Cylindrical, SinhCylindrical,
    #              SymTP, SinhSymTP
    # NOTE: Only CoordSystem == Cartesian makes sense here; new
    #       boundary conditions are needed within the ETK for
    #       Spherical, etc. coordinates.
    CoordSystem     = "Cartesian"

    par.set_parval_from_str("reference_metric::CoordSystem",CoordSystem)
    rfm.reference_metric() # Create ReU, ReDD needed for rescaling B-L initial data, generating BSSN RHSs, etc.

    REAL      = "CCTK_REAL" # Set REAL to CCTK_REAL, the ETK data type for
                            # floating point precision (typically `double`)
    # Set finite differencing order:
    par.set_parval_from_str("finite_difference::FD_CENTDERIVS_ORDER", FD_order)

    # Copy SIMD/SIMD_intrinsics.h to $outdir/SIMD/SIMD_intrinsics.h
    cmd.mkdir(os.path.join(outdir,"SIMD"))
    shutil.copy(os.path.join("SIMD/")+"SIMD_intrinsics.h",os.path.join(outdir,"SIMD/"))

    # Set the gridfunction memory access type to ETK-like, so that finite_difference
    #    knows how to read and write gridfunctions from/to memory.
    par.set_parval_from_str("grid::GridFuncMemAccess","ETK")

    # Step 2: Output C code for BSSN spacetime solve
    # Step 2.a: BSSN RHS expressions

    import time  # Standard Python module; useful for benchmarking below expression & code generation.

    import BSSN.BSSN_RHSs as rhs
    import BSSN.BSSN_gauge_RHSs as gaugerhs
    par.set_parval_from_str("BSSN.BSSN_gauge_RHSs::ShiftEvolutionOption", ShiftCondition)
    par.set_parval_from_str("BSSN.BSSN_gauge_RHSs::LapseEvolutionOption", LapseCondition)

    print("Generating symbolic expressions for BSSN RHSs...")
    start = time.time()
    # Enable rfm_precompute infrastructure, which results in
    #   BSSN RHSs that are free of transcendental functions,
    #   even in curvilinear coordinates, so long as
    #   ConformalFactor is set to "W" (default).
    cmd.mkdir(os.path.join(outdir, "rfm_files/"))
    par.set_parval_from_str("reference_metric::enable_rfm_precompute", "True")
    par.set_parval_from_str("reference_metric::rfm_precompute_Ccode_outdir", os.path.join(outdir, "rfm_files/"))

    # Evaluate BSSN + BSSN gauge RHSs with rfm_precompute enabled:
    import BSSN.BSSN_quantities as Bq
    par.set_parval_from_str("BSSN.BSSN_quantities::LeaveRicciSymbolic", "True")

    rhs.BSSN_RHSs()

    if add_stress_energy_source_terms == True:
        T4UU = ixp.register_gridfunctions_for_single_rank2("AUXEVOL", "T4UU", "sym01", DIM=4)
        import BSSN.BSSN_stress_energy_source_terms as Bsest
        Bsest.BSSN_source_terms_for_BSSN_RHSs(T4UU)
        rhs.trK_rhs += Bsest.sourceterm_trK_rhs
        for i in range(DIM):
            # Needed for Gamma-driving shift RHSs:
            rhs.Lambdabar_rhsU[i] += Bsest.sourceterm_Lambdabar_rhsU[i]
            # Needed for BSSN RHSs:
            rhs.lambda_rhsU[i] += Bsest.sourceterm_lambda_rhsU[i]
            for j in range(DIM):
                rhs.a_rhsDD[i][j] += Bsest.sourceterm_a_rhsDD[i][j]

    gaugerhs.BSSN_gauge_RHSs()

    # Add Kreiss-Oliger dissipation to the BSSN RHSs:
    thismodule = "KO_Dissipation"
    diss_strength = par.Cparameters("REAL", thismodule, "diss_strength", default_KO_strength)

    alpha_dKOD = ixp.declarerank1("alpha_dKOD")
    cf_dKOD = ixp.declarerank1("cf_dKOD")
    trK_dKOD = ixp.declarerank1("trK_dKOD")
    betU_dKOD = ixp.declarerank2("betU_dKOD", "nosym")
    vetU_dKOD = ixp.declarerank2("vetU_dKOD", "nosym")
    lambdaU_dKOD = ixp.declarerank2("lambdaU_dKOD", "nosym")
    aDD_dKOD = ixp.declarerank3("aDD_dKOD", "sym01")
    hDD_dKOD = ixp.declarerank3("hDD_dKOD", "sym01")
    for k in range(DIM):
        gaugerhs.alpha_rhs += diss_strength * alpha_dKOD[k] * rfm.ReU[k]  # ReU[k] = 1/scalefactor_orthog_funcform[k]
        rhs.cf_rhs += diss_strength * cf_dKOD[k] * rfm.ReU[k]  # ReU[k] = 1/scalefactor_orthog_funcform[k]
        rhs.trK_rhs += diss_strength * trK_dKOD[k] * rfm.ReU[k]  # ReU[k] = 1/scalefactor_orthog_funcform[k]
        for i in range(DIM):
            if "2ndOrder" in ShiftCondition:
                gaugerhs.bet_rhsU[i] += diss_strength * betU_dKOD[i][k] * rfm.ReU[
                    k]  # ReU[k] = 1/scalefactor_orthog_funcform[k]
            gaugerhs.vet_rhsU[i] += diss_strength * vetU_dKOD[i][k] * rfm.ReU[
                k]  # ReU[k] = 1/scalefactor_orthog_funcform[k]
            rhs.lambda_rhsU[i] += diss_strength * lambdaU_dKOD[i][k] * rfm.ReU[
                k]  # ReU[k] = 1/scalefactor_orthog_funcform[k]
            for j in range(DIM):
                rhs.a_rhsDD[i][j] += diss_strength * aDD_dKOD[i][j][k] * rfm.ReU[
                    k]  # ReU[k] = 1/scalefactor_orthog_funcform[k]
                rhs.h_rhsDD[i][j] += diss_strength * hDD_dKOD[i][j][k] * rfm.ReU[
                    k]  # ReU[k] = 1/scalefactor_orthog_funcform[k]

    # We use betaU as our upwinding control vector:
    Bq.BSSN_basic_tensors()
    betaU = Bq.betaU

    import BSSN.Enforce_Detgammabar_Constraint as EGC
    enforce_detg_constraint_symb_expressions = EGC.Enforce_Detgammabar_Constraint_symb_expressions()

    # Next compute Ricci tensor
    par.set_parval_from_str("BSSN.BSSN_quantities::LeaveRicciSymbolic", "False")
    Bq.RicciBar__gammabarDD_dHatD__DGammaUDD__DGammaU()

    # Now that we are finished with all the rfm hatted
    #           quantities in generic precomputed functional
    #           form, let's restore them to their closed-
    #           form expressions.
    par.set_parval_from_str("reference_metric::enable_rfm_precompute",
                            "False")  # Reset to False to disable rfm_precompute.
    rfm.ref_metric__hatted_quantities()
    end = time.time()
    print("Finished BSSN symbolic expressions in " + str(end - start) + " seconds.")

    def BSSN_RHSs():
        print("Generating C code for BSSN RHSs in " + par.parval_from_str(
            "reference_metric::CoordSystem") + " coordinates.")
        start = time.time()

        BSSN_evol_rhss = [ \
            lhrh(lhs=gri.gfaccess("rhs_gfs", "aDD00"), rhs=rhs.a_rhsDD[0][0]),
            lhrh(lhs=gri.gfaccess("rhs_gfs", "aDD01"), rhs=rhs.a_rhsDD[0][1]),
            lhrh(lhs=gri.gfaccess("rhs_gfs", "aDD02"), rhs=rhs.a_rhsDD[0][2]),
            lhrh(lhs=gri.gfaccess("rhs_gfs", "aDD11"), rhs=rhs.a_rhsDD[1][1]),
            lhrh(lhs=gri.gfaccess("rhs_gfs", "aDD12"), rhs=rhs.a_rhsDD[1][2]),
            lhrh(lhs=gri.gfaccess("rhs_gfs", "aDD22"), rhs=rhs.a_rhsDD[2][2]),
            lhrh(lhs=gri.gfaccess("rhs_gfs", "alpha"), rhs=gaugerhs.alpha_rhs),
            lhrh(lhs=gri.gfaccess("rhs_gfs", "betU0"), rhs=gaugerhs.bet_rhsU[0]),
            lhrh(lhs=gri.gfaccess("rhs_gfs", "betU1"), rhs=gaugerhs.bet_rhsU[1]),
            lhrh(lhs=gri.gfaccess("rhs_gfs", "betU2"), rhs=gaugerhs.bet_rhsU[2]),
            lhrh(lhs=gri.gfaccess("rhs_gfs", "cf"), rhs=rhs.cf_rhs),
            lhrh(lhs=gri.gfaccess("rhs_gfs", "hDD00"), rhs=rhs.h_rhsDD[0][0]),
            lhrh(lhs=gri.gfaccess("rhs_gfs", "hDD01"), rhs=rhs.h_rhsDD[0][1]),
            lhrh(lhs=gri.gfaccess("rhs_gfs", "hDD02"), rhs=rhs.h_rhsDD[0][2]),
            lhrh(lhs=gri.gfaccess("rhs_gfs", "hDD11"), rhs=rhs.h_rhsDD[1][1]),
            lhrh(lhs=gri.gfaccess("rhs_gfs", "hDD12"), rhs=rhs.h_rhsDD[1][2]),
            lhrh(lhs=gri.gfaccess("rhs_gfs", "hDD22"), rhs=rhs.h_rhsDD[2][2]),
            lhrh(lhs=gri.gfaccess("rhs_gfs", "lambdaU0"), rhs=rhs.lambda_rhsU[0]),
            lhrh(lhs=gri.gfaccess("rhs_gfs", "lambdaU1"), rhs=rhs.lambda_rhsU[1]),
            lhrh(lhs=gri.gfaccess("rhs_gfs", "lambdaU2"), rhs=rhs.lambda_rhsU[2]),
            lhrh(lhs=gri.gfaccess("rhs_gfs", "trK"), rhs=rhs.trK_rhs),
            lhrh(lhs=gri.gfaccess("rhs_gfs", "vetU0"), rhs=gaugerhs.vet_rhsU[0]),
            lhrh(lhs=gri.gfaccess("rhs_gfs", "vetU1"), rhs=gaugerhs.vet_rhsU[1]),
            lhrh(lhs=gri.gfaccess("rhs_gfs", "vetU2"), rhs=gaugerhs.vet_rhsU[2])]

        BSSN_RHSs_string = fin.FD_outputC("returnstring", BSSN_evol_rhss, params="outCverbose=False,SIMD_enable=True",
                                          upwindcontrolvec=betaU)

        with open(os.path.join(outdir, "BSSN_RHSs.h"), "w") as file:
            file.write(
                lp.loop(["i2", "i1", "i0"], ["cctk_nghostzones[2]", "cctk_nghostzones[1]", "cctk_nghostzones[0]"],
                        ["cctk_lsh[2]-cctk_nghostzones[2]", "cctk_lsh[1]-cctk_nghostzones[1]",
                         "cctk_lsh[0]-cctk_nghostzones[0]"],
                        ["1", "1", "SIMD_width"],
                        ["#pragma omp parallel for",
                         "#include \"rfm_files/rfm_struct__SIMD_outer_read2.h\"",
                         "#include \"rfm_files/rfm_struct__SIMD_outer_read1.h\""], "",
                        "#include \"rfm_files/rfm_struct__SIMD_inner_read0.h\"\n" + BSSN_RHSs_string))
        end = time.time()
        print("Finished BSSN_RHS C codegen in " + str(end - start) + " seconds.")

    def Ricci():
        print("Generating C code for Ricci tensor in " + par.parval_from_str(
            "reference_metric::CoordSystem") + " coordinates.")
        start = time.time()
        Ricci_string = fin.FD_outputC("returnstring",
                                      [lhrh(lhs=gri.gfaccess("auxevol_gfs", "RbarDD00"), rhs=Bq.RbarDD[0][0]),
                                       lhrh(lhs=gri.gfaccess("auxevol_gfs", "RbarDD01"), rhs=Bq.RbarDD[0][1]),
                                       lhrh(lhs=gri.gfaccess("auxevol_gfs", "RbarDD02"), rhs=Bq.RbarDD[0][2]),
                                       lhrh(lhs=gri.gfaccess("auxevol_gfs", "RbarDD11"), rhs=Bq.RbarDD[1][1]),
                                       lhrh(lhs=gri.gfaccess("auxevol_gfs", "RbarDD12"), rhs=Bq.RbarDD[1][2]),
                                       lhrh(lhs=gri.gfaccess("auxevol_gfs", "RbarDD22"), rhs=Bq.RbarDD[2][2])],
                                      params="outCverbose=False,SIMD_enable=True")
        with open(os.path.join(outdir, "BSSN_Ricci.h"), "w") as file:
            file.write(
                lp.loop(["i2", "i1", "i0"], ["cctk_nghostzones[2]", "cctk_nghostzones[1]", "cctk_nghostzones[0]"],
                        ["cctk_lsh[2]-cctk_nghostzones[2]", "cctk_lsh[1]-cctk_nghostzones[1]",
                         "cctk_lsh[0]-cctk_nghostzones[0]"],
                        ["1", "1", "SIMD_width"],
                        ["#pragma omp parallel for",
                         "#include \"rfm_files/rfm_struct__SIMD_outer_read2.h\"",
                         "#include \"rfm_files/rfm_struct__SIMD_outer_read1.h\""], "",
                        "#include \"rfm_files/rfm_struct__SIMD_inner_read0.h\"\n" + Ricci_string))
        end = time.time()
        print("Finished Ricci C codegen in " + str(end - start) + " seconds.")

    # Step 2.b: Hamiltonian & momentum constraints
    # First register the Hamiltonian as a gridfunction.
    H = gri.register_gridfunctions("AUX", "H")
    MU = ixp.register_gridfunctions_for_single_rank1("AUX", "MU")

    # Then define the Hamiltonian constraint and output the optimized C code.
    import BSSN.BSSN_constraints as bssncon

    def BSSNconstraints():
        bssncon.BSSN_constraints(add_T4UUmunu_source_terms=False)
        if add_stress_energy_source_terms == True:
            #         T4UU = gri.register_gridfunctions_for_single_rank2("AUXEVOL","T4UU","sym01",DIM=4)
            import BSSN.BSSN_stress_energy_source_terms as Bsest
            Bsest.BSSN_source_terms_for_BSSN_constraints(T4UU)
            bssncon.H += Bsest.sourceterm_H
            for i in range(DIM):
                bssncon.MU[i] += Bsest.sourceterm_MU[i]

        start = time.time()
        print("Generating optimized C code for Ham. & mom. constraints. May take a while, depending on CoordSystem.")
        Ham_mom_string = fin.FD_outputC("returnstring",
                                        [lhrh(lhs=gri.gfaccess("aux_gfs", "H"), rhs=bssncon.H),
                                         lhrh(lhs=gri.gfaccess("aux_gfs", "MU0"), rhs=bssncon.MU[0]),
                                         lhrh(lhs=gri.gfaccess("aux_gfs", "MU1"), rhs=bssncon.MU[1]),
                                         lhrh(lhs=gri.gfaccess("aux_gfs", "MU2"), rhs=bssncon.MU[2])],
                                        params="outCverbose=False")

        with open(os.path.join(outdir, "BSSN_constraints.h"), "w") as file:
            file.write(
                lp.loop(["i2", "i1", "i0"], ["cctk_nghostzones[2]", "cctk_nghostzones[1]", "cctk_nghostzones[0]"],
                        ["cctk_lsh[2]-cctk_nghostzones[2]", "cctk_lsh[1]-cctk_nghostzones[1]",
                         "cctk_lsh[0]-cctk_nghostzones[0]"],
                        ["1", "1", "1"], ["#pragma omp parallel for", "", ""], "", Ham_mom_string))
        end = time.time()
        print("Finished Hamiltonian & momentum constraint C codegen in " + str(end - start) + " seconds.")

    # Step 2.c: Enforce det(gammabar) = det(gammahat) constraint (det(gammahat)=1 in Cartesian coordinates)
    def gammadet():
        start = time.time()
        print("Generating optimized C code for gamma constraint. May take a while, depending on CoordSystem.")
        enforce_gammadet_string = fin.FD_outputC("returnstring", enforce_detg_constraint_symb_expressions,
                                                 params="outCverbose=False,preindent=0,includebraces=False")

        with open(os.path.join(outdir, "enforce_detgammabar_constraint.h"), "w") as file:
            file.write(lp.loop(["i2", "i1", "i0"], ["0", "0", "0"],
                               ["cctk_lsh[2]", "cctk_lsh[1]", "cctk_lsh[0]"],
                               ["1", "1", "1"],
                               ["#pragma omp parallel for",
                                "#include \"rfm_files/rfm_struct__read2.h\"",
                                "#include \"rfm_files/rfm_struct__read1.h\""], "",
                               "#include \"rfm_files/rfm_struct__read0.h\"\n" + enforce_gammadet_string))
        end = time.time()
        print("Finished gamma constraint C codegen in " + str(end - start) + " seconds.")

    # Step 2.d: Generate all C codes in parallel
    # Step 2.d.0: Import the multiprocessing module.
    import multiprocessing
    
    # Step 2.d.1: Create a list of functions we wish to evaluate in parallel
    funcs = [BSSN_RHSs,Ricci,Hamiltonian,gammadet]
    # Step 2.d.1.a: Define master function for parallelization.
    #           Note that lambdifying this doesn't work in Python 3
    def master_func(arg):
        funcs[arg]()

    # Step 2.d.2: Evaluate list of functions in parallel if allowed; 
    #         otherwise fallback to serial evaluation:
    try:
#            if __name__ == '__main__':
        pool = multiprocessing.Pool()
        pool.map(master_func,range(len(funcs)))
    except:
        # If multiprocessing didn't work for whatever reason,
        #        evaluate functions in serial.
        for func in funcs:
            func()import multiprocessing

    # Step 3: ETK .ccl file generation
    # Step 3.a: param.ccl: specify free parameters within BaikalETK

    def keep_param__return_type(paramtuple):
        keep_param = True  # We'll not set some parameters in param.ccl;
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
            if paramtuple.defaultval != 1e300:  # 1e300 is a magic value indicating that the C parameter should be mutable
                typestring += "CCTK_REAL "
            else:
                keep_param = False
        elif partype == "int":
            typestring += "CCTK_INT "
        elif partype == "#define":
            keep_param = False
        elif partype == "char":
            # FIXME: char/string parameter types should in principle be supported
            print("Error: parameter " + paramtuple.module + "::" + paramtuple.parname +
                  " has unsupported type: \"" + paramtuple.type + "\"")
            sys.exit(1)
        else:
            print("Error: parameter " + paramtuple.module + "::" + paramtuple.parname +
                  " has unsupported type: \"" + paramtuple.type + "\"")
            sys.exit(1)
        return keep_param, typestring

    with open(os.path.join(outrootdir, "param.ccl"), "w") as file:
        file.write("""
# This param.ccl file was automatically generated by NRPy+. 
#   You are advised against modifying it directly.

shares: ADMBase
USES CCTK_INT lapse_timelevels   # Needed to ensure ADMBase gridfunctions are allocated (see top of schedule.ccl)
USES CCTK_INT shift_timelevels   # Needed to ensure ADMBase gridfunctions are allocated (see top of schedule.ccl)
USES CCTK_INT metric_timelevels  # Needed to ensure ADMBase gridfunctions are allocated (see top of schedule.ccl)

EXTENDS CCTK_KEYWORD evolution_method "evolution_method"
{
  "BaikalETK" :: ""
} 

EXTENDS CCTK_KEYWORD lapse_evolution_method "lapse_evolution_method"
{
  "BaikalETK" :: ""
} 

EXTENDS CCTK_KEYWORD shift_evolution_method "shift_evolution_method"
{
  "BaikalETK" :: ""
} 

EXTENDS CCTK_KEYWORD dtshift_evolution_method "dtshift_evolution_method"
{
  "BaikalETK" :: ""
} 

EXTENDS CCTK_KEYWORD dtlapse_evolution_method "dtlapse_evolution_method"
{
  "BaikalETK" :: ""
} 

restricted:
""")
        paramccl_str = ""
        for i in range(len(par.glb_Cparams_list)):
            # keep_param is a boolean indicating whether we should accept or reject
            #    the parameter. singleparstring will contain the string indicating
            #    the variable type.
            keep_param, singleparstring = keep_param__return_type(par.glb_Cparams_list[i])

            if keep_param:
                parname = par.glb_Cparams_list[i].parname
                partype = par.glb_Cparams_list[i].type
                singleparstring += parname + " \"" + parname + " (see NRPy+ for parameter definition)\"\n"
                singleparstring += "{\n"
                if partype != "bool":
                    singleparstring += " *:* :: \"All values accepted. NRPy+ does not restrict the allowed ranges of parameters yet.\"\n"
                singleparstring += "} " + str(par.glb_Cparams_list[i].defaultval) + "\n\n"

                paramccl_str += singleparstring
        file.write(paramccl_str)

    # Step 3.b: interface.ccl: define needed gridfunctions; provide keywords denoting what this
    #                          thorn provides and what it should inherit from other thorns

    # First construct lists of the basic gridfunctions used in NRPy+.
    #    Each type will be its own group in BaikalETK.
    evol_gfs_list = []
    auxevol_gfs_list = []
    aux_gfs_list = []
    for i in range(len(gri.glb_gridfcs_list)):
        if gri.glb_gridfcs_list[i].gftype == "EVOL":
            evol_gfs_list.append(gri.glb_gridfcs_list[i].name + "GF")

        if gri.glb_gridfcs_list[i].gftype == "AUX":
            aux_gfs_list.append(gri.glb_gridfcs_list[i].name + "GF")

        if gri.glb_gridfcs_list[i].gftype == "AUXEVOL":
            auxevol_gfs_list.append(gri.glb_gridfcs_list[i].name + "GF")

    # NRPy+'s finite-difference code generator assumes gridfunctions
    #    are alphabetized; not sorting may result in unnecessary
    #    cache misses.
    evol_gfs_list.sort()
    aux_gfs_list.sort()
    auxevol_gfs_list.sort()

    with open(os.path.join(outrootdir, "interface.ccl"), "w") as file:
        file.write("""
# With "implements", we give our thorn its unique name.
implements: BaikalETK

# By "inheriting" other thorns, we tell the Toolkit that we 
#   will rely on variables/function that exist within those
#   functions. 
inherits:   ADMBase Boundary Grid MethodofLines\n""")
        if add_stress_energy_source_terms == True:
            file.write("inherits:   TmunuBase")
        file.write("""
# Needed functions and #include's:
USES INCLUDE: Symmetry.h
USES INCLUDE: Boundary.h

# Needed Method of Lines function
CCTK_INT FUNCTION MoLRegisterEvolvedGroup(CCTK_INT IN EvolvedIndex, \
                                          CCTK_INT IN RHSIndex)
REQUIRES FUNCTION MoLRegisterEvolvedGroup

# Needed Boundary Conditions function
CCTK_INT FUNCTION GetBoundarySpecification(CCTK_INT IN size, CCTK_INT OUT ARRAY nboundaryzones, CCTK_INT OUT ARRAY is_internal, CCTK_INT OUT ARRAY is_staggered, CCTK_INT OUT ARRAY shiftout)
USES FUNCTION GetBoundarySpecification

CCTK_INT FUNCTION SymmetryTableHandleForGrid(CCTK_POINTER_TO_CONST IN cctkGH)
USES FUNCTION SymmetryTableHandleForGrid

CCTK_INT FUNCTION Boundary_SelectVarForBC(CCTK_POINTER_TO_CONST IN GH, CCTK_INT IN faces, CCTK_INT IN boundary_width, CCTK_INT IN table_handle, CCTK_STRING IN var_name, CCTK_STRING IN bc_name)
USES FUNCTION Boundary_SelectVarForBC

# Needed for EinsteinEvolve/NewRad outer boundary condition driver:
CCTK_INT FUNCTION                         \\
    NewRad_Apply                          \\
        (CCTK_POINTER_TO_CONST IN cctkGH, \\
         CCTK_REAL ARRAY IN var,          \\
         CCTK_REAL ARRAY INOUT rhs,       \\
         CCTK_REAL IN var0,               \\
         CCTK_REAL IN v0,                 \\
         CCTK_INT IN radpower)
REQUIRES FUNCTION NewRad_Apply

# Needed to convert ADM initial data into BSSN initial data (gamma extrapolation)
CCTK_INT FUNCTION                         \\
    ExtrapolateGammas                     \\
        (CCTK_POINTER_TO_CONST IN cctkGH, \\
         CCTK_REAL ARRAY INOUT var)
REQUIRES FUNCTION ExtrapolateGammas

# Tell the Toolkit that we want all gridfunctions
#    to be visible to other thorns by using 
#    the keyword "public". Note that declaring these 
#    gridfunctions *does not* allocate memory for them;
#    that is done by the schedule.ccl file.

# FIXME: add info for symmetry conditions: 
#    https://einsteintoolkit.org/thornguide/CactusBase/SymBase/documentation.html
public:
""")

        # Next we declare gridfunctions based on their corresponding gridfunction groups as registered within NRPy+

        def output_list_of_gfs(gfs_list, description="User did not provide description"):
            gfsstr = "    "
            for i in range(len(gfs_list)):
                gfsstr += gfs_list[i]
                if i != len(gfs_list) - 1:
                    gfsstr += ","  # This is a comma-separated list of gridfunctions
                else:
                    gfsstr += "\n} \"" + description + "\"\n\n"
            return gfsstr

        # First EVOL type:
        file.write("CCTK_REAL evol_variables type = GF Timelevels=3\n{\n")
        file.write(output_list_of_gfs(evol_gfs_list, "BSSN evolved gridfunctions"))
        # Second EVOL right-hand-sides
        file.write(
            "CCTK_REAL evol_variables_rhs type = GF Timelevels=1 TAGS=\'InterpNumTimelevels=1 prolongation=\"none\"\'\n{\n")
        rhs_list = []
        for gf in evol_gfs_list:
            rhs_list.append(gf.replace("GF", "") + "_rhsGF")
        file.write(output_list_of_gfs(rhs_list, "right-hand-side storage for BSSN evolved gridfunctions"))
        # Then AUX type:
        file.write("CCTK_REAL aux_variables type = GF Timelevels=3\n{\n")
        file.write(output_list_of_gfs(aux_gfs_list, "Auxiliary gridfunctions for BSSN diagnostics"))
        # Finally, AUXEVOL type:
        file.write(
            "CCTK_REAL auxevol_variables type = GF Timelevels=1 TAGS=\'InterpNumTimelevels=1 prolongation=\"none\"\'\n{\n")
        file.write(output_list_of_gfs(auxevol_gfs_list, "Auxiliary gridfunctions needed for evaluating the BSSN RHSs"))

        # Step 3.c: schedule.ccl: schedule all functions used within BaikalETK,
        #                         specify data dependencies within said functions,
        #                         and allocate memory for gridfunctions
        with open(os.path.join(outrootdir, "schedule.ccl"), "w") as file:
            file.write("""
# First allocate storage for all ADMBase gridfunctions, which are needed by NRPy+
STORAGE: ADMBase::metric[metric_timelevels], ADMBase::curv[metric_timelevels], ADMBase::lapse[lapse_timelevels], ADMBase::shift[shift_timelevels]

# Next allocate storage for all 3 gridfunction groups used in BaikalETK
STORAGE: evol_variables[3]     # Evolution variables
STORAGE: evol_variables_rhs[1] # Variables storing right-hand-sides
STORAGE: aux_variables[3]      # Diagnostics variables
STORAGE: auxevol_variables[1]  # Single-timelevel storage of variables needed for evolutions.

# The following scheduler is based on Lean/LeanBSSNMoL/schedule.ccl

schedule BaikalETK_Banner at STARTUP
{
  LANG: C
  OPTIONS: meta
} "Output ASCII art banner"

schedule BaikalETK_RegisterSlicing at STARTUP after BaikalETK_Banner
{
  LANG: C
  OPTIONS: meta
} "Register 3+1 slicing condition"

schedule BaikalETK_Symmetry_registration at BASEGRID
{
  LANG: C
  OPTIONS: Global
} "Register symmetries, the CartGrid3D way."

schedule BaikalETK_zero_rhss at BASEGRID after BaikalETK_Symmetry_registration
{
  LANG: C
} "Idea from Lean: set all rhs functions to zero to prevent spurious nans"

schedule BaikalETK_ADM_to_BSSN at CCTK_INITIAL after ADMBase_PostInitial
{
  LANG: C
  OPTIONS: Local
  SYNC: evol_variables
} "Convert initial data into BSSN variables"

schedule GROUP ApplyBCs as BaikalETK_ApplyBCs at CCTK_INITIAL after BaikalETK_ADM_to_BSSN
{
} "Apply boundary conditions"


# MoL: registration

schedule BaikalETK_MoL_registration in MoL_Register
{
  LANG: C
  OPTIONS: META
} "Register variables for MoL"


# MoL: compute RHSs, etc
""")
            if add_stress_energy_source_terms == True:
                file.write("""
schedule driver_BSSN_T4UU in MoL_CalcRHS as BaikalETK_T4UU before BaikalETK_BSSN_to_ADM
{
  LANG: C
} "MoL: Compute T4UU, needed for BSSN RHSs."

schedule BaikalETK_BSSN_to_ADM in MoL_CalcRHS after BaikalETK_T4UU before BaikalETK_Ricci
{
  LANG: C
} "Perform BSSN-to-ADM conversion. Needed for HydroBase coupling."
""")
            file.write("""
schedule driver_pt1_BSSN_Ricci in MoL_CalcRHS as BaikalETK_Ricci before BaikalETK_RHS
{
  LANG: C
} "MoL: Compute Ricci tensor"

schedule driver_pt2_BSSN_RHSs in MoL_CalcRHS as BaikalETK_RHS after BaikalETK_Ricci
{
  LANG: C
} "MoL: Evaluate BSSN RHSs"

schedule BaikalETK_NewRad in MoL_CalcRHS after BaikalETK_RHS
{
  LANG: C
} "NewRad boundary conditions, scheduled right after RHS eval."

schedule enforce_detgammabar_constraint in MoL_PostStep before BC_Update
{
  LANG: C
} "Enforce detgammabar = detgammahat (= 1 in Cartesian)"

schedule BaikalETK_BoundaryConditions_evolved_gfs in MoL_PostStep
{
  LANG: C
  OPTIONS: LEVEL
  SYNC: evol_variables
} "Apply boundary conditions and perform AMR+interprocessor synchronization"

schedule GROUP ApplyBCs as BaikalETK_ApplyBCs in MoL_PostStep after BaikalETK_BoundaryConditions_evolved_gfs
{
} "Group for applying boundary conditions"


# Next update ADM quantities

schedule BaikalETK_BSSN_to_ADM in MoL_PostStep after BaikalETK_ApplyBCs before ADMBase_SetADMVars
{
  LANG: C
  OPTIONS: Local
} "Perform BSSN-to-ADM conversion. Useful for diagnostics."

# Compute Hamiltonian & momentum constraints
""")
            if add_stress_energy_source_terms == True:
                file.write("""
schedule driver_BSSN_T4UU in MoL_PseudoEvolution before BaikalETK_BSSN_constraints
{
  LANG: C
  OPTIONS: Local
} "MoL_PseudoEvolution: Compute T4UU, needed for BSSN constraints"
""")
            file.write("""

schedule BaikalETK_BSSN_constraints in MoL_PseudoEvolution
{
  LANG: C
  OPTIONS: Local
} "Compute BSSN (Hamiltonian and momentum) constraints"

schedule BaikalETK_BoundaryConditions_aux_gfs in MoL_PseudoEvolution after BaikalETK_BSSN_constraints
{
  LANG: C
  OPTIONS: LOCAL # Needed so that cctk_nghostzones[0] (the number of boundary points) is defined.
                 #  In other words, don't use LEVEL mode here, or the number of boundary points
                 #  filled may not match the actual number of ghost zones. Weird, huh?
  SYNC: aux_variables
} "Enforce symmetry BCs in constraint computation"

""")
            if add_stress_energy_source_terms == True:
                file.write("""
schedule BaikalETK_BSSN_to_ADM in MoL_PseudoEvolution after BaikalETK_BoundaryConditions_aux_gfs
{
  LANG: C
  OPTIONS: Local
} "Perform BSSN-to-ADM conversion in MoL_PseudoEvolution. Needed for proper HydroBase integration."
""")
            file.write("""
schedule GROUP ApplyBCs as BaikalETK_auxgfs_ApplyBCs in MoL_PseudoEvolution after BaikalETK_BoundaryConditions_aux_gfs
{
} "Apply boundary conditions"
""")

    # Step 4: C driver functions for ETK registration & NRPy+-generated kernels
    make_code_defn_list = []

    def append_to_make_code_defn_list(filename):
        if filename not in make_code_defn_list:
            make_code_defn_list.append(filename)
        return os.path.join(outdir, filename)

    with open(append_to_make_code_defn_list("RegisterSlicing.c"), "w") as file:
        file.write("""    
#include "cctk.h"

#include "Slicing.h"

int BaikalETK_RegisterSlicing (void)
{
  Einstein_RegisterSlicing ("BaikalETK");
  return 0;
}""")
    # First the ETK banner code, proudly showing the NRPy+ banner
    import NRPy_logo as logo

    with open(append_to_make_code_defn_list("Banner.c"),"w") as file:
        file.write("""
#include <stdio.h>

void BaikalETK_Banner() 
{
    """)
        logostr = logo.print_logo(print_to_stdout=False)
        file.write("printf(\"BaikalETK: another Einstein Toolkit thorn generated by\\n\");\n")
        for line in logostr.splitlines():
            file.write("    printf(\""+line+"\\n\");\n")
        file.write("}\n")

    # Next register symmetries
    with open(append_to_make_code_defn_list("Symmetry_registration_oldCartGrid3D.c"), "w") as file:
        file.write("""
#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"
#include "Symmetry.h"

void BaikalETK_Symmetry_registration(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  // Stores gridfunction parity across x=0, y=0, and z=0 planes, respectively
  int sym[3];

  // Next register parities for each gridfunction based on its name 
  //    (to ensure this algorithm is robust, gridfunctions with integers
  //     in their base names are forbidden in NRPy+).
""")
        full_gfs_list = []
        full_gfs_list.extend(evol_gfs_list)
        full_gfs_list.extend(auxevol_gfs_list)
        full_gfs_list.extend(aux_gfs_list)
        for gf in full_gfs_list:
            file.write("""
  // Default to scalar symmetry:
  sym[0] = 1; sym[1] = 1; sym[2] = 1;
  // Now modify sym[0], sym[1], and/or sym[2] as needed 
  //    to account for gridfunction parity across 
  //    x=0, y=0, and/or z=0 planes, respectively
""")
            # If gridfunction name does not end in a digit, by NRPy+ syntax, it must be a scalar
            if gf[len(gf) - 1].isdigit() == False:
                pass  # Scalar = default
            elif len(gf) > 2:
                # Rank-1 indexed expression (e.g., vector)
                if gf[len(gf) - 2].isdigit() == False:
                    if int(gf[-1]) > 2:
                        print("Error: Found invalid gridfunction name: " + gf)
                        sys.exit(1)
                    symidx = gf[-1]
                    file.write("  sym[" + symidx + "] = -1;\n")
                # Rank-2 indexed expression
                elif gf[len(gf) - 2].isdigit() == True:
                    if len(gf) > 3 and gf[len(gf) - 3].isdigit() == True:
                        print(
                            "Error: Found a Rank-3 or above gridfunction: " + gf + ", which is at the moment unsupported.")
                        print("It should be easy to support this if desired.")
                        sys.exit(1)
                    symidx0 = gf[-2]
                    file.write("  sym[" + symidx0 + "] *= -1;\n")
                    symidx1 = gf[-1]
                    file.write("  sym[" + symidx1 + "] *= -1;\n")
            else:
                print(
                    "Don't know how you got this far with a gridfunction named " + gf + ", but I'll take no more of this nonsense.")
                print("   Please follow best-practices and rename your gridfunction to be more descriptive")
                sys.exit(1)
            file.write("  SetCartSymVN(cctkGH, sym, \"BaikalETK::" + gf + "\");\n")
        file.write("}\n")
    # Next register symmetries
    with open(append_to_make_code_defn_list("zero_rhss.c"), "w") as file:
        file.write("""
#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"
#include "Symmetry.h"

void BaikalETK_zero_rhss(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
""")
        set_rhss_to_zero = ""
        for gf in rhs_list:
            set_rhss_to_zero += gf + "[CCTK_GFINDEX3D(cctkGH,i0,i1,i2)] = 0.0;\n"

        file.write(lp.loop(["i2", "i1", "i0"], ["0", "0", "0"],
                           ["cctk_lsh[2]", "cctk_lsh[1]", "cctk_lsh[0]"],
                           ["1", "1", "1"],
                           ["#pragma omp parallel for", "", "", ], "", set_rhss_to_zero))
        file.write("}\n")
    # Next registration with the Method of Lines thorn
    with open(append_to_make_code_defn_list("MoL_registration.c"), "w") as file:
        file.write("""
//--------------------------------------------------------------------------
// Register with the Method of Lines time stepper
// (MoL thorn, found in arrangements/CactusBase/MoL)
// MoL documentation located in arrangements/CactusBase/MoL/doc
//--------------------------------------------------------------------------
#include <stdio.h>

#include "cctk.h"
#include "cctk_Parameters.h"
#include "cctk_Arguments.h"

#include "Symmetry.h"

void BaikalETK_MoL_registration(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  CCTK_INT ierr = 0, group, rhs;

  // Register evolution & RHS gridfunction groups with MoL, so it knows

  group = CCTK_GroupIndex("BaikalETK::evol_variables");
  rhs = CCTK_GroupIndex("BaikalETK::evol_variables_rhs");
  ierr += MoLRegisterEvolvedGroup(group, rhs);

  if (ierr) CCTK_ERROR("Problems registering with MoL");
}
""")
    # Next register with the boundary conditions thorns.
    # PART 1: Set BC type to "none" for all variables
    # Since we choose NewRad boundary conditions, we must register all
    #   gridfunctions to have boundary type "none". This is because
    #   NewRad is seen by the rest of the Toolkit as a modification to the
    #   RHSs.

    # This code is based on Kranc's McLachlan/ML_BSSN/src/Boundaries.cc code.
    with open(append_to_make_code_defn_list("BoundaryConditions.c"), "w") as file:
        file.write("""
#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"
#include "cctk_Faces.h"
#include "util_Table.h"
#include "Symmetry.h"

// Set `none` boundary conditions on BSSN RHSs, as these are set via NewRad.
void BaikalETK_BoundaryConditions_evolved_gfs(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  CCTK_INT ierr CCTK_ATTRIBUTE_UNUSED = 0;
""")
        for gf in evol_gfs_list:
            file.write("""
  ierr = Boundary_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, -1, "BaikalETK::""" + gf + """", "none");
  if (ierr < 0) CCTK_ERROR("Failed to register BC for BaikalETK::""" + gf + """!");
""")
        file.write("""
}

// Set `flat` boundary conditions on BSSN constraints, similar to what Lean does.
void BaikalETK_BoundaryConditions_aux_gfs(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  CCTK_INT ierr CCTK_ATTRIBUTE_UNUSED = 0;

""")
        for gf in aux_gfs_list:
            file.write("""
  ierr = Boundary_SelectVarForBC(cctkGH, CCTK_ALL_FACES, cctk_nghostzones[0], -1, "BaikalETK::""" + gf + """", "flat");
  if (ierr < 0) CCTK_ERROR("Failed to register BC for BaikalETK::""" + gf + """!");
""")
        file.write("}\n")

    # PART 2: Set C code for calling NewRad BCs
    #   As explained in lean_public/LeanBSSNMoL/src/calc_bssn_rhs.F90,
    #   the function NewRad_Apply takes the following arguments:
    #   NewRad_Apply(cctkGH, var, rhs, var0, v0, radpower),
    #     which implement the boundary condition:
    #       var  =  var_at_infinite_r + u(r-var_char_speed*t)/r^var_radpower
    #  Obviously for var_radpower>0, var_at_infinite_r is the value of
    #    the variable at r->infinity. var_char_speed is the propagation
    #    speed at the outer boundary, and var_radpower is the radial
    #    falloff rate.

    with open(append_to_make_code_defn_list("BoundaryCondition_NewRad.c"), "w") as file:
        file.write("""
#include <math.h>

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

void BaikalETK_NewRad(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

""")
        for gf in evol_gfs_list:
            var_at_infinite_r = "0.0"
            var_char_speed = "1.0"
            var_radpower = "1.0"

            if gf == "alpha":
                var_at_infinite_r = "1.0"
                if LapseCondition == "OnePlusLog":
                    var_char_speed = "sqrt(2.0)"
                else:
                    pass  # 1.0 (default) is fine
            if "aDD" in gf or "trK" in gf:  # consistent with Lean code.
                var_radpower = "2.0"

            file.write(
                "  NewRad_Apply(cctkGH, " + gf + ", " + gf.replace("GF", "") + "_rhsGF, " + var_at_infinite_r + ", " +
                var_char_speed + ", " + var_radpower + ");\n")
        file.write("}\n")

    # Step 4.b: BSSN to/from ADM conversions
    # Step 4.b.i: ADM to BSSN (i.e., BSSN in terms of ADM variables)
    # First we convert from ADM to BSSN, as is required to convert initial data
    #    (given using) ADM quantities, to the BSSN evolved variables
    import BSSN.ADM_Numerical_Spherical_or_Cartesian_to_BSSNCurvilinear as atob

    IDhDD, IDaDD, IDtrK, IDvetU, IDbetU, IDalpha, IDcf, IDlambdaU = \
        atob.Convert_Spherical_or_Cartesian_ADM_to_BSSN_curvilinear("Cartesian", "DoNotOutputADMInputFunction", outdir)

    alphaSphorCart = gri.register_gridfunctions("AUXEVOL", "alphaSphorCart")
    betaSphorCartU = ixp.register_gridfunctions_for_single_rank1("AUXEVOL", "betaSphorCartU")
    BSphorCartU = ixp.register_gridfunctions_for_single_rank1("AUXEVOL", "BSphorCartU")
    gammaSphorCartDD = ixp.register_gridfunctions_for_single_rank2("AUXEVOL", "gammaSphorCartDD", "sym01")
    KSphorCartDD = ixp.register_gridfunctions_for_single_rank2("AUXEVOL", "KSphorCartDD", "sym01")

    # Step : Output ADM to BSSN conversion.
    with open(append_to_make_code_defn_list("ADM_to_BSSN.c"), "w") as file:
        file.write("""
#include <math.h>

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

void BaikalETK_ADM_to_BSSN(CCTK_ARGUMENTS) {
    DECLARE_CCTK_ARGUMENTS;
    DECLARE_CCTK_PARAMETERS;

    CCTK_REAL *alphaSphorCartGF = alp;
""")
        # It's ugly if we output code in the following ordering, so we'll first
        #   output to a string and then sort the string to beautify the code a bit.
        outstr = []
        for i in range(DIM):
            outstr.append("    CCTK_REAL *betaSphorCartU" + str(i) + "GF = beta" + chr(ord('x') + i) + ";\n")
            outstr.append("    CCTK_REAL *BSphorCartU" + str(i) + "GF = dtbeta" + chr(ord('x') + i) + ";\n")
            for j in range(i, DIM):
                outstr.append("    CCTK_REAL *gammaSphorCartDD" + str(i) + str(j) + "GF = g" + chr(ord('x') + i) + chr(
                    ord('x') + j) + ";\n")
                outstr.append("    CCTK_REAL *KSphorCartDD" + str(i) + str(j) + "GF = k" + chr(ord('x') + i) + chr(
                    ord('x') + j) + ";\n")
        outstr.sort()
        for line in outstr:
            file.write(line)
        file.write("""
    const CCTK_REAL invdx0 = 1.0/CCTK_DELTA_SPACE(0);
    const CCTK_REAL invdx1 = 1.0/CCTK_DELTA_SPACE(1);
    const CCTK_REAL invdx2 = 1.0/CCTK_DELTA_SPACE(2);
""")

        all_but_lambdaU_expressions = [
            lhrh(lhs=gri.gfaccess("in_gfs", "hDD00"), rhs=IDhDD[0][0]),
            lhrh(lhs=gri.gfaccess("in_gfs", "hDD01"), rhs=IDhDD[0][1]),
            lhrh(lhs=gri.gfaccess("in_gfs", "hDD02"), rhs=IDhDD[0][2]),
            lhrh(lhs=gri.gfaccess("in_gfs", "hDD11"), rhs=IDhDD[1][1]),
            lhrh(lhs=gri.gfaccess("in_gfs", "hDD12"), rhs=IDhDD[1][2]),
            lhrh(lhs=gri.gfaccess("in_gfs", "hDD22"), rhs=IDhDD[2][2]),
            lhrh(lhs=gri.gfaccess("in_gfs", "aDD00"), rhs=IDaDD[0][0]),
            lhrh(lhs=gri.gfaccess("in_gfs", "aDD01"), rhs=IDaDD[0][1]),
            lhrh(lhs=gri.gfaccess("in_gfs", "aDD02"), rhs=IDaDD[0][2]),
            lhrh(lhs=gri.gfaccess("in_gfs", "aDD11"), rhs=IDaDD[1][1]),
            lhrh(lhs=gri.gfaccess("in_gfs", "aDD12"), rhs=IDaDD[1][2]),
            lhrh(lhs=gri.gfaccess("in_gfs", "aDD22"), rhs=IDaDD[2][2]),
            lhrh(lhs=gri.gfaccess("in_gfs", "trK"), rhs=IDtrK),
            lhrh(lhs=gri.gfaccess("in_gfs", "vetU0"), rhs=IDvetU[0]),
            lhrh(lhs=gri.gfaccess("in_gfs", "vetU1"), rhs=IDvetU[1]),
            lhrh(lhs=gri.gfaccess("in_gfs", "vetU2"), rhs=IDvetU[2]),
            lhrh(lhs=gri.gfaccess("in_gfs", "betU0"), rhs=IDbetU[0]),
            lhrh(lhs=gri.gfaccess("in_gfs", "betU1"), rhs=IDbetU[1]),
            lhrh(lhs=gri.gfaccess("in_gfs", "betU2"), rhs=IDbetU[2]),
            lhrh(lhs=gri.gfaccess("in_gfs", "alpha"), rhs=IDalpha),
            lhrh(lhs=gri.gfaccess("in_gfs", "cf"), rhs=IDcf)]

        outCparams = "preindent=1,outCfileaccess=a,outCverbose=False,includebraces=False"
        all_but_lambdaU_outC = fin.FD_outputC("returnstring", all_but_lambdaU_expressions, outCparams)
        file.write(lp.loop(["i2", "i1", "i0"], ["0", "0", "0"], ["cctk_lsh[2]", "cctk_lsh[1]", "cctk_lsh[0]"],
                           ["1", "1", "1"], ["#pragma omp parallel for", "", ""], "", all_but_lambdaU_outC))

        outCparams = "preindent=1,outCfileaccess=a,outCverbose=False,includebraces=False"
        lambdaU_expressions = [lhrh(lhs=gri.gfaccess("in_gfs", "lambdaU0"), rhs=IDlambdaU[0]),
                               lhrh(lhs=gri.gfaccess("in_gfs", "lambdaU1"), rhs=IDlambdaU[1]),
                               lhrh(lhs=gri.gfaccess("in_gfs", "lambdaU2"), rhs=IDlambdaU[2])]
        lambdaU_expressions_FDout = fin.FD_outputC("returnstring", lambdaU_expressions, outCparams)

        file.write(lp.loop(["i2", "i1", "i0"], ["cctk_nghostzones[2]", "cctk_nghostzones[1]", "cctk_nghostzones[0]"],
                           ["cctk_lsh[2]-cctk_nghostzones[2]", "cctk_lsh[1]-cctk_nghostzones[1]",
                            "cctk_lsh[0]-cctk_nghostzones[0]"],
                           ["1", "1", "1"], ["#pragma omp parallel for", "", ""], "", lambdaU_expressions_FDout))

        file.write("""
    ExtrapolateGammas(cctkGH,lambdaU0GF);
    ExtrapolateGammas(cctkGH,lambdaU1GF);
    ExtrapolateGammas(cctkGH,lambdaU2GF);
}
""")

    # Step 4.b.ii: BSSN to ADM (i.e., ADM in terms of BSSN variables)
    import BSSN.ADM_in_terms_of_BSSN as btoa

    btoa.ADM_in_terms_of_BSSN()
    Bq.BSSN_basic_tensors()  # Gives us betaU & BU

    with open(append_to_make_code_defn_list("BSSN_to_ADM.c"), "w") as file:
        file.write("""
#include <math.h>

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

void BaikalETK_BSSN_to_ADM(CCTK_ARGUMENTS) {
    DECLARE_CCTK_ARGUMENTS;
    DECLARE_CCTK_PARAMETERS;

""")
        btoa_lhrh = []
        for i in range(DIM):
            for j in range(i, DIM):
                btoa_lhrh.append(
                    lhrh(lhs="g" + chr(ord('x') + i) + chr(ord('x') + j) + "[CCTK_GFINDEX3D(cctkGH,i0,i1,i2)]",
                         rhs=btoa.gammaDD[i][j]))
        for i in range(DIM):
            for j in range(i, DIM):
                btoa_lhrh.append(
                    lhrh(lhs="k" + chr(ord('x') + i) + chr(ord('x') + j) + "[CCTK_GFINDEX3D(cctkGH,i0,i1,i2)]",
                         rhs=btoa.KDD[i][j]))
        btoa_lhrh.append(lhrh(lhs="alp[CCTK_GFINDEX3D(cctkGH,i0,i1,i2)]", rhs=Bq.alpha))

        for i in range(DIM):
            btoa_lhrh.append(lhrh(lhs="beta" + chr(ord('x') + i) + "[CCTK_GFINDEX3D(cctkGH,i0,i1,i2)]",
                                  rhs=Bq.betaU[i]))

        for i in range(DIM):
            btoa_lhrh.append(lhrh(lhs="dtbeta" + chr(ord('x') + i) + "[CCTK_GFINDEX3D(cctkGH,i0,i1,i2)]",
                                  rhs=Bq.BU[i]))

        outCparams = "preindent=1,outCfileaccess=a,outCverbose=False,includebraces=False"
        bssn_to_adm_Ccode = fin.FD_outputC("returnstring", btoa_lhrh, outCparams)
        file.write(lp.loop(["i2", "i1", "i0"], ["0", "0", "0"], ["cctk_lsh[2]", "cctk_lsh[1]", "cctk_lsh[0]"],
                           ["1", "1", "1"], ["#pragma omp parallel for", "", ""], "", bssn_to_adm_Ccode))

        file.write("}\n")

    # Step 4.c: C driver functions for evaluating BSSN RHSs & related quantities
    # Step 4.c.i: Evaluate Ricci tensor
    with open(append_to_make_code_defn_list("driver_pt1_BSSN_Ricci.c"), "w") as file:
        file.write("""
#include <math.h>

#include "SIMD/SIMD_intrinsics.h"

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

void driver_pt1_BSSN_Ricci(CCTK_ARGUMENTS) {
    DECLARE_CCTK_ARGUMENTS;

    const CCTK_REAL NOSIMDinvdx0 = 1.0/CCTK_DELTA_SPACE(0);
    const REAL_SIMD_ARRAY invdx0 = ConstSIMD(NOSIMDinvdx0);
    const CCTK_REAL NOSIMDinvdx1 = 1.0/CCTK_DELTA_SPACE(1);
    const REAL_SIMD_ARRAY invdx1 = ConstSIMD(NOSIMDinvdx1);
    const CCTK_REAL NOSIMDinvdx2 = 1.0/CCTK_DELTA_SPACE(2);
    const REAL_SIMD_ARRAY invdx2 = ConstSIMD(NOSIMDinvdx2);
#include "BSSN_Ricci.h"
}\n""")

    # Step 4.c.ii: Evaluate BSSN RHSs, using Ricci tensor as input
    def SIMD_declare_C_params():
        SIMD_declare_C_params_str = ""
        for i in range(len(par.glb_Cparams_list)):
            # keep_param is a boolean indicating whether we should accept or reject
            #    the parameter. singleparstring will contain the string indicating
            #    the variable type.
            keep_param, singleparstring = keep_param__return_type(par.glb_Cparams_list[i])

            if (keep_param) and ("CCTK_REAL" in singleparstring):
                parname = par.glb_Cparams_list[i].parname
                SIMD_declare_C_params_str += "    const " + singleparstring + "*NOSIMD" + parname + \
                                             " = CCTK_ParameterGet(\"" + parname + "\",\"BaikalETK\",NULL);\n"
                SIMD_declare_C_params_str += "    const REAL_SIMD_ARRAY " + parname + " = ConstSIMD(*NOSIMD" + parname + ");\n"
        return SIMD_declare_C_params_str


    with open(append_to_make_code_defn_list("driver_pt2_BSSN_RHSs.c"), "w") as file:
        file.write("""
#include <math.h>

#include "SIMD/SIMD_intrinsics.h"

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

//void BSSN_RHSs()

void driver_pt2_BSSN_RHSs(CCTK_ARGUMENTS) {
    DECLARE_CCTK_ARGUMENTS;

    const CCTK_REAL NOSIMDinvdx0 = 1.0/CCTK_DELTA_SPACE(0);
    const REAL_SIMD_ARRAY invdx0 = ConstSIMD(NOSIMDinvdx0);
    const CCTK_REAL NOSIMDinvdx1 = 1.0/CCTK_DELTA_SPACE(1);
    const REAL_SIMD_ARRAY invdx1 = ConstSIMD(NOSIMDinvdx1);
    const CCTK_REAL NOSIMDinvdx2 = 1.0/CCTK_DELTA_SPACE(2);
    const REAL_SIMD_ARRAY invdx2 = ConstSIMD(NOSIMDinvdx2);
""" + SIMD_declare_C_params() + """
#include "BSSN_RHSs.h"
}\n""")

    # Step 4.d: Enforcing conformal 3-metric constraint det(gammabar) = det(gammahat),
    #           where det(gammahat) = 1 in Cartesian coordinates
    with open(append_to_make_code_defn_list("enforce_detgammabar_constraint.c"), "w") as file:
        file.write("""
#include <math.h>

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

void enforce_detgammabar_constraint(CCTK_ARGUMENTS) {
    DECLARE_CCTK_ARGUMENTS;
    DECLARE_CCTK_PARAMETERS;

#include "enforce_detgammabar_constraint.h"
}\n""")

    # Step 4.e: Diagnostics: Computing Hamiltonian & momentum constraints
    with open(append_to_make_code_defn_list("BSSN_constraints.c"), "w") as file:
        file.write("""
#include <math.h>

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

void BaikalETK_BSSN_constraints(CCTK_ARGUMENTS) {
    DECLARE_CCTK_ARGUMENTS;
    DECLARE_CCTK_PARAMETERS;

    const CCTK_REAL invdx0 = 1.0/CCTK_DELTA_SPACE(0);
    const CCTK_REAL invdx1 = 1.0/CCTK_DELTA_SPACE(1);
    const CCTK_REAL invdx2 = 1.0/CCTK_DELTA_SPACE(2);

#include "BSSN_constraints.h"
}\n""")
    # Step 4.f: driver_BSSN_T4UU(): Compute T^{mu nu} from TmunuBase's T_{mu nu}
    if add_stress_energy_source_terms == True:
        # Declare T4DD as a set of gridfunctions. These won't
        #    actually appear in interface.ccl, as interface.ccl
        #    was set above. Thus before calling the code output
        #    by FD_outputC(), we'll have to set pointers
        #    to the actual gridfunctions they reference.
        #    (In this case the eTab's.)
        T4DD = ixp.register_gridfunctions_for_single_rank2("AUXEVOL", "T4DD", "sym01", DIM=4)
        import BSSN.ADMBSSN_tofrom_4metric as AB4m

        AB4m.g4UU_ito_BSSN_or_ADM("BSSN")

        T4UUraised = ixp.zerorank2(DIM=4)
        for mu in range(4):
            for nu in range(4):
                for delta in range(4):
                    for gamma in range(4):
                        T4UUraised[mu][nu] += AB4m.g4UU[mu][delta] * AB4m.g4UU[nu][gamma] * T4DD[delta][gamma]

        T4UU_expressions = [
            lhrh(lhs=gri.gfaccess("in_gfs", "T4UU00"), rhs=T4UUraised[0][0]),
            lhrh(lhs=gri.gfaccess("in_gfs", "T4UU01"), rhs=T4UUraised[0][1]),
            lhrh(lhs=gri.gfaccess("in_gfs", "T4UU02"), rhs=T4UUraised[0][2]),
            lhrh(lhs=gri.gfaccess("in_gfs", "T4UU03"), rhs=T4UUraised[0][3]),
            lhrh(lhs=gri.gfaccess("in_gfs", "T4UU11"), rhs=T4UUraised[1][1]),
            lhrh(lhs=gri.gfaccess("in_gfs", "T4UU12"), rhs=T4UUraised[1][2]),
            lhrh(lhs=gri.gfaccess("in_gfs", "T4UU13"), rhs=T4UUraised[1][3]),
            lhrh(lhs=gri.gfaccess("in_gfs", "T4UU22"), rhs=T4UUraised[2][2]),
            lhrh(lhs=gri.gfaccess("in_gfs", "T4UU23"), rhs=T4UUraised[2][3]),
            lhrh(lhs=gri.gfaccess("in_gfs", "T4UU33"), rhs=T4UUraised[3][3])]

        outCparams = "outCverbose=False,includebraces=False,preindent=2,SIMD_enable=True"
        T4UUstr = fin.FD_outputC("returnstring", T4UU_expressions, outCparams)
        T4UUstr_loop = lp.loop(["i2", "i1", "i0"], ["0", "0", "0"], ["cctk_lsh[2]", "cctk_lsh[1]", "cctk_lsh[0]"],
                               ["1", "1", "SIMD_width"], ["#pragma omp parallel for", "", ""], "", T4UUstr)

        with open(append_to_make_code_defn_list("driver_BSSN_T4UU.c"), "w") as file:
            file.write("""
#include <math.h>

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

#include "SIMD/SIMD_intrinsics.h"

void driver_BSSN_T4UU(CCTK_ARGUMENTS) {
    DECLARE_CCTK_ARGUMENTS;
    DECLARE_CCTK_PARAMETERS;

    const CCTK_REAL *restrict T4DD00GF = eTtt;
    const CCTK_REAL *restrict T4DD01GF = eTtx;
    const CCTK_REAL *restrict T4DD02GF = eTty;
    const CCTK_REAL *restrict T4DD03GF = eTtz;
    const CCTK_REAL *restrict T4DD11GF = eTxx;
    const CCTK_REAL *restrict T4DD12GF = eTxy;
    const CCTK_REAL *restrict T4DD13GF = eTxz;
    const CCTK_REAL *restrict T4DD22GF = eTyy;
    const CCTK_REAL *restrict T4DD23GF = eTyz;
    const CCTK_REAL *restrict T4DD33GF = eTzz;
""" + T4UUstr_loop + """
}\n""")

    # Step 4.g: make.code.defn: List of all C driver functions needed to compile BaikalETK
    with open(os.path.join(outdir, "make.code.defn"), "w") as file:
        file.write("""
# Main make.code.defn file for thorn BaikalETK

# Source files in this directory
SRCS =""")
        filestring = ""
        for i in range(len(make_code_defn_list)):
            filestring += "      " + make_code_defn_list[i]
            if i != len(make_code_defn_list) - 1:
                filestring += " \\\n"
            else:
                filestring += "\n"
        file.write(filestring)
