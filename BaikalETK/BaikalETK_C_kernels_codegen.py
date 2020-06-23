
# Step 1: Import needed core NRPy+ modules
from outputC import lhrh,outC_function_dict  # NRPy+: Core C code output module
import finite_difference as fin  # NRPy+: Finite difference C code generation module
import NRPy_param_funcs as par   # NRPy+: Parameter interface
import grid as gri               # NRPy+: Functions having to do with numerical grids
import loop as lp                # NRPy+: Generate C code loops
import indexedexp as ixp         # NRPy+: Symbolic indexed expression (e.g., tensors, vectors, etc.) support
import reference_metric as rfm   # NRPy+: Reference metric support
import os, sys                   # Standard Python modules for multiplatform OS-level functions
import time                      # Standard Python module; useful for benchmarking
import BSSN.BSSN_RHSs as rhs
import BSSN.BSSN_gauge_RHSs as gaugerhs

def BaikalETK_C_kernels_codegen_onepart(params=
        "WhichPart=BSSN_RHSs,ThornName=Baikal,FD_order=4,enable_stress_energy_source_terms=True"):
    # Set default parameters
    WhichPart = "BSSN_RHSs"
    ThornName = "Baikal"
    FD_order = 4
    enable_stress_energy_source_terms = True
    LapseCondition  = "OnePlusLog" # Set the standard 1+log lapse condition
    # Set the standard, second-order advecting-shift, Gamma-driving shift condition:
    ShiftCondition  = "GammaDriving2ndOrder_NoCovariant"
    # Default Kreiss-Oliger dissipation strength
    default_KO_strength = 0.1

    import re
    if params != "":
        params2 = re.sub("^,","",params)
        params = params2.strip()
        splitstring = re.split("=|,", params)

        if len(splitstring) % 2 != 0:
            print("outputC: Invalid params string: "+params)
            sys.exit(1)

        parnm = []
        value = []
        for i in range(int(len(splitstring)/2)):
            parnm.append(splitstring[2*i])
            value.append(splitstring[2*i+1])

        for i in range(len(parnm)):
            parnm.append(splitstring[2*i])
            value.append(splitstring[2*i+1])

        for i in range(len(parnm)):
            if parnm[i] == "WhichPart":
                WhichPart = value[i]
            elif parnm[i] == "ThornName":
                ThornName = value[i]
            elif parnm[i] == "FD_order":
                FD_order = int(value[i])
            elif parnm[i] == "enable_stress_energy_source_terms":
                enable_stress_energy_source_terms = False
                if value[i] == "True":
                    enable_stress_energy_source_terms = True
            elif parnm[i] == "LapseCondition":
                LapseCondition = value[i]
            elif parnm[i] == "ShiftCondition":
                ShiftCondition = value[i]
            elif parnm[i] == "default_KO_strength":
                default_KO_strength = float(value[i])
            else:
                print("BaikalETK Error: Could not parse input param: "+parnm[i])
                sys.exit(1)

    # Set output directory for C kernels
    outdir = os.path.join(ThornName, "src")  # Main C code output directory

    # Set spatial dimension (must be 3 for BSSN)
    par.set_parval_from_str("grid::DIM", 3)

    # Step 2: Set some core parameters, including CoordSystem MoL timestepping algorithm,
    #                                 FD order, floating point precision, and CFL factor:
    # Choices are: Spherical, SinhSpherical, SinhSphericalv2, Cylindrical, SinhCylindrical,
    #              SymTP, SinhSymTP
    # NOTE: Only CoordSystem == Cartesian makes sense here; new
    #       boundary conditions are needed within the ETK for
    #       Spherical, etc. coordinates.
    CoordSystem     = "Cartesian"

    par.set_parval_from_str("reference_metric::CoordSystem", CoordSystem)
    rfm.reference_metric()  # Create ReU, ReDD needed for rescaling B-L initial data, generating BSSN RHSs, etc.

    # Set the gridfunction memory access type to ETK-like, so that finite_difference
    #    knows how to read and write gridfunctions from/to memory.
    par.set_parval_from_str("grid::GridFuncMemAccess", "ETK")

    par.set_parval_from_str("BSSN.BSSN_gauge_RHSs::ShiftEvolutionOption", ShiftCondition)
    par.set_parval_from_str("BSSN.BSSN_gauge_RHSs::LapseEvolutionOption", LapseCondition)

    T4UU = None
    if enable_stress_energy_source_terms == True:
        registered_already = False
        for i in range(len(gri.glb_gridfcs_list)):
            if gri.glb_gridfcs_list[i].name == "T4UU00":
                registered_already = True
        if not registered_already:
            T4UU = ixp.register_gridfunctions_for_single_rank2("AUXEVOL", "T4UU", "sym01", DIM=4)
        else:
            T4UU = ixp.declarerank2("T4UU", "sym01", DIM=4)

    # Register the BSSN constraints (Hamiltonian & momentum constraints) as gridfunctions.
    registered_already = False
    for i in range(len(gri.glb_gridfcs_list)):
        if gri.glb_gridfcs_list[i].name == "H":
            registered_already = True
    if not registered_already:
        # We ignore return values for register_gridfunctions...() calls below
        #    as they are unused.
        gri.register_gridfunctions("AUX", "H")
        ixp.register_gridfunctions_for_single_rank1("AUX", "MU")

    def BSSN_RHSs__generate_symbolic_expressions():
        ######################################
        # START: GENERATE SYMBOLIC EXPRESSIONS
        print("Generating symbolic expressions for BSSN RHSs...")
        start = time.time()
        # Enable rfm_precompute infrastructure, which results in
        #   BSSN RHSs that are free of transcendental functions,
        #   even in curvilinear coordinates, so long as
        #   ConformalFactor is set to "W" (default).
        par.set_parval_from_str("reference_metric::enable_rfm_precompute","True")
        par.set_parval_from_str("reference_metric::rfm_precompute_Ccode_outdir",os.path.join(outdir,"rfm_files/"))

        # Evaluate BSSN + BSSN gauge RHSs with rfm_precompute enabled:
        import BSSN.BSSN_quantities as Bq
        par.set_parval_from_str("BSSN.BSSN_quantities::LeaveRicciSymbolic","True")

        rhs.BSSN_RHSs()

        if T4UU != None:
            import BSSN.BSSN_stress_energy_source_terms as Bsest
            Bsest.BSSN_source_terms_for_BSSN_RHSs(T4UU)
            rhs.trK_rhs += Bsest.sourceterm_trK_rhs
            for i in range(3):
                # Needed for Gamma-driving shift RHSs:
                rhs.Lambdabar_rhsU[i] += Bsest.sourceterm_Lambdabar_rhsU[i]
                # Needed for BSSN RHSs:
                rhs.lambda_rhsU[i]    += Bsest.sourceterm_lambda_rhsU[i]
                for j in range(3):
                    rhs.a_rhsDD[i][j] += Bsest.sourceterm_a_rhsDD[i][j]

        gaugerhs.BSSN_gauge_RHSs()

        # Add Kreiss-Oliger dissipation to the BSSN RHSs:
        thismodule = "KO_Dissipation"
        diss_strength = par.Cparameters("REAL", thismodule, "diss_strength", default_KO_strength)

        alpha_dKOD = ixp.declarerank1("alpha_dKOD")
        cf_dKOD    = ixp.declarerank1("cf_dKOD")
        trK_dKOD   = ixp.declarerank1("trK_dKOD")
        betU_dKOD    = ixp.declarerank2("betU_dKOD","nosym")
        vetU_dKOD    = ixp.declarerank2("vetU_dKOD","nosym")
        lambdaU_dKOD = ixp.declarerank2("lambdaU_dKOD","nosym")
        aDD_dKOD = ixp.declarerank3("aDD_dKOD","sym01")
        hDD_dKOD = ixp.declarerank3("hDD_dKOD","sym01")
        for k in range(3):
            gaugerhs.alpha_rhs += diss_strength*alpha_dKOD[k]*rfm.ReU[k] # ReU[k] = 1/scalefactor_orthog_funcform[k]
            rhs.cf_rhs         += diss_strength*   cf_dKOD[k]*rfm.ReU[k] # ReU[k] = 1/scalefactor_orthog_funcform[k]
            rhs.trK_rhs        += diss_strength*  trK_dKOD[k]*rfm.ReU[k] # ReU[k] = 1/scalefactor_orthog_funcform[k]
            for i in range(3):
                if "2ndOrder" in ShiftCondition:
                    gaugerhs.bet_rhsU[i] += diss_strength*   betU_dKOD[i][k]*rfm.ReU[k] # ReU[k] = 1/scalefactor_orthog_funcform[k]
                gaugerhs.vet_rhsU[i]     += diss_strength*   vetU_dKOD[i][k]*rfm.ReU[k] # ReU[k] = 1/scalefactor_orthog_funcform[k]
                rhs.lambda_rhsU[i]       += diss_strength*lambdaU_dKOD[i][k]*rfm.ReU[k] # ReU[k] = 1/scalefactor_orthog_funcform[k]
                for j in range(3):
                    rhs.a_rhsDD[i][j] += diss_strength*aDD_dKOD[i][j][k]*rfm.ReU[k] # ReU[k] = 1/scalefactor_orthog_funcform[k]
                    rhs.h_rhsDD[i][j] += diss_strength*hDD_dKOD[i][j][k]*rfm.ReU[k] # ReU[k] = 1/scalefactor_orthog_funcform[k]

        # We use betaU as our upwinding control vector:
        Bq.BSSN_basic_tensors()
        betaU = Bq.betaU

        # Now that we are finished with all the rfm hatted
        #           quantities in generic precomputed functional
        #           form, let's restore them to their closed-
        #           form expressions.
        par.set_parval_from_str("reference_metric::enable_rfm_precompute","False") # Reset to False to disable rfm_precompute.
        rfm.ref_metric__hatted_quantities()
        par.set_parval_from_str("BSSN.BSSN_quantities::LeaveRicciSymbolic","False")
        end = time.time()
        print("(BENCH) Finished BSSN RHS symbolic expressions in "+str(end-start)+" seconds.")
        # END: GENERATE SYMBOLIC EXPRESSIONS
        ######################################

        BSSN_RHSs_SymbExpressions = [lhrh(lhs=gri.gfaccess("rhs_gfs","aDD00"),   rhs=rhs.a_rhsDD[0][0]),
                                     lhrh(lhs=gri.gfaccess("rhs_gfs","aDD01"),   rhs=rhs.a_rhsDD[0][1]),
                                     lhrh(lhs=gri.gfaccess("rhs_gfs","aDD02"),   rhs=rhs.a_rhsDD[0][2]),
                                     lhrh(lhs=gri.gfaccess("rhs_gfs","aDD11"),   rhs=rhs.a_rhsDD[1][1]),
                                     lhrh(lhs=gri.gfaccess("rhs_gfs","aDD12"),   rhs=rhs.a_rhsDD[1][2]),
                                     lhrh(lhs=gri.gfaccess("rhs_gfs","aDD22"),   rhs=rhs.a_rhsDD[2][2]),
                                     lhrh(lhs=gri.gfaccess("rhs_gfs","alpha"),   rhs=gaugerhs.alpha_rhs),
                                     lhrh(lhs=gri.gfaccess("rhs_gfs","betU0"),   rhs=gaugerhs.bet_rhsU[0]),
                                     lhrh(lhs=gri.gfaccess("rhs_gfs","betU1"),   rhs=gaugerhs.bet_rhsU[1]),
                                     lhrh(lhs=gri.gfaccess("rhs_gfs","betU2"),   rhs=gaugerhs.bet_rhsU[2]),
                                     lhrh(lhs=gri.gfaccess("rhs_gfs","cf"),      rhs=rhs.cf_rhs),
                                     lhrh(lhs=gri.gfaccess("rhs_gfs","hDD00"),   rhs=rhs.h_rhsDD[0][0]),
                                     lhrh(lhs=gri.gfaccess("rhs_gfs","hDD01")   ,rhs=rhs.h_rhsDD[0][1]),
                                     lhrh(lhs=gri.gfaccess("rhs_gfs","hDD02"),   rhs=rhs.h_rhsDD[0][2]),
                                     lhrh(lhs=gri.gfaccess("rhs_gfs","hDD11"),   rhs=rhs.h_rhsDD[1][1]),
                                     lhrh(lhs=gri.gfaccess("rhs_gfs","hDD12"),   rhs=rhs.h_rhsDD[1][2]),
                                     lhrh(lhs=gri.gfaccess("rhs_gfs","hDD22"),   rhs=rhs.h_rhsDD[2][2]),
                                     lhrh(lhs=gri.gfaccess("rhs_gfs","lambdaU0"),rhs=rhs.lambda_rhsU[0]),
                                     lhrh(lhs=gri.gfaccess("rhs_gfs","lambdaU1"),rhs=rhs.lambda_rhsU[1]),
                                     lhrh(lhs=gri.gfaccess("rhs_gfs","lambdaU2"),rhs=rhs.lambda_rhsU[2]),
                                     lhrh(lhs=gri.gfaccess("rhs_gfs","trK"),     rhs=rhs.trK_rhs),
                                     lhrh(lhs=gri.gfaccess("rhs_gfs","vetU0"),   rhs=gaugerhs.vet_rhsU[0]),
                                     lhrh(lhs=gri.gfaccess("rhs_gfs","vetU1"),   rhs=gaugerhs.vet_rhsU[1]),
                                     lhrh(lhs=gri.gfaccess("rhs_gfs","vetU2"),   rhs=gaugerhs.vet_rhsU[2]) ]

        return [betaU,BSSN_RHSs_SymbExpressions]

    def Ricci__generate_symbolic_expressions():
        ######################################
        # START: GENERATE SYMBOLIC EXPRESSIONS
        print("Generating symbolic expressions for Ricci tensor...")
        start = time.time()
        # Enable rfm_precompute infrastructure, which results in
        #   BSSN RHSs that are free of transcendental functions,
        #   even in curvilinear coordinates, so long as
        #   ConformalFactor is set to "W" (default).
        par.set_parval_from_str("reference_metric::enable_rfm_precompute","True")
        par.set_parval_from_str("reference_metric::rfm_precompute_Ccode_outdir",os.path.join(outdir,"rfm_files/"))

        # Evaluate BSSN + BSSN gauge RHSs with rfm_precompute enabled:
        import BSSN.BSSN_quantities as Bq

        # Next compute Ricci tensor
        # par.set_parval_from_str("BSSN.BSSN_quantities::LeaveRicciSymbolic","False")
        RbarDD_already_registered = False
        for i in range(len(gri.glb_gridfcs_list)):
            if "RbarDD00" in gri.glb_gridfcs_list[i].name:
                RbarDD_already_registered = True
        if not RbarDD_already_registered:
            # We ignore the return value of ixp.register_gridfunctions_for_single_rank2() below
            #    as it is unused.
            ixp.register_gridfunctions_for_single_rank2("AUXEVOL","RbarDD","sym01")
        rhs.BSSN_RHSs()
        Bq.RicciBar__gammabarDD_dHatD__DGammaUDD__DGammaU()

        # Now that we are finished with all the rfm hatted
        #           quantities in generic precomputed functional
        #           form, let's restore them to their closed-
        #           form expressions.
        par.set_parval_from_str("reference_metric::enable_rfm_precompute","False") # Reset to False to disable rfm_precompute.
        rfm.ref_metric__hatted_quantities()
        end = time.time()
        print("(BENCH) Finished Ricci symbolic expressions in "+str(end-start)+" seconds.")
        # END: GENERATE SYMBOLIC EXPRESSIONS
        ######################################

        Ricci_SymbExpressions = [lhrh(lhs=gri.gfaccess("auxevol_gfs","RbarDD00"),rhs=Bq.RbarDD[0][0]),
                                 lhrh(lhs=gri.gfaccess("auxevol_gfs","RbarDD01"),rhs=Bq.RbarDD[0][1]),
                                 lhrh(lhs=gri.gfaccess("auxevol_gfs","RbarDD02"),rhs=Bq.RbarDD[0][2]),
                                 lhrh(lhs=gri.gfaccess("auxevol_gfs","RbarDD11"),rhs=Bq.RbarDD[1][1]),
                                 lhrh(lhs=gri.gfaccess("auxevol_gfs","RbarDD12"),rhs=Bq.RbarDD[1][2]),
                                 lhrh(lhs=gri.gfaccess("auxevol_gfs","RbarDD22"),rhs=Bq.RbarDD[2][2])]

        return [Ricci_SymbExpressions]

    def BSSN_RHSs__generate_Ccode(all_RHSs_Ricci_exprs_list):
        betaU                     = all_RHSs_Ricci_exprs_list[0]
        BSSN_RHSs_SymbExpressions = all_RHSs_Ricci_exprs_list[1]

        print("Generating C code for BSSN RHSs (FD_order="+str(FD_order)+",Tmunu="+str(enable_stress_energy_source_terms)+") in "+par.parval_from_str("reference_metric::CoordSystem")+" coordinates.")
        start = time.time()

        # Store original finite-differencing order:
        FD_order_orig = par.parval_from_str("finite_difference::FD_CENTDERIVS_ORDER")
        # Set new finite-differencing order:
        par.set_parval_from_str("finite_difference::FD_CENTDERIVS_ORDER", FD_order)

        BSSN_RHSs_string = fin.FD_outputC("returnstring",BSSN_RHSs_SymbExpressions,
                                          params="outCverbose=False,SIMD_enable=True,GoldenKernelsEnable=True",
                                          upwindcontrolvec=betaU)

        filename = "BSSN_RHSs_enable_Tmunu_"+str(enable_stress_energy_source_terms)+"_FD_order_"+str(FD_order)+".h"
        with open(os.path.join(outdir,filename), "w") as file:
            file.write(lp.loop(["i2","i1","i0"],["cctk_nghostzones[2]","cctk_nghostzones[1]","cctk_nghostzones[0]"],
           ["cctk_lsh[2]-cctk_nghostzones[2]","cctk_lsh[1]-cctk_nghostzones[1]","cctk_lsh[0]-cctk_nghostzones[0]"],
                               ["1","1","SIMD_width"],
                                ["#pragma omp parallel for",
                                 "#include \"rfm_files/rfm_struct__SIMD_outer_read2.h\"",
                                 r"""    #include "rfm_files/rfm_struct__SIMD_outer_read1.h"
    #if (defined __INTEL_COMPILER && __INTEL_COMPILER_BUILD_DATE >= 20180804)
        #pragma ivdep         // Forces Intel compiler (if Intel compiler used) to ignore certain SIMD vector dependencies
        #pragma vector always // Forces Intel compiler (if Intel compiler used) to vectorize
    #endif"""],"",
                    "#include \"rfm_files/rfm_struct__SIMD_inner_read0.h\"\n"+BSSN_RHSs_string))

        # Restore original finite-differencing order:
        par.set_parval_from_str("finite_difference::FD_CENTDERIVS_ORDER", FD_order_orig)

        end = time.time()
        print("(BENCH) Finished BSSN_RHS C codegen (FD_order="+str(FD_order)+",Tmunu="+str(enable_stress_energy_source_terms)+") in " + str(end - start) + " seconds.")

    def Ricci__generate_Ccode(Ricci_exprs_list):
        Ricci_SymbExpressions     = Ricci_exprs_list[0]

        print("Generating C code for Ricci tensor (FD_order="+str(FD_order)+") in "+par.parval_from_str("reference_metric::CoordSystem")+" coordinates.")
        start = time.time()

        # Store original finite-differencing order:
        FD_order_orig = par.parval_from_str("finite_difference::FD_CENTDERIVS_ORDER")
        # Set new finite-differencing order:
        par.set_parval_from_str("finite_difference::FD_CENTDERIVS_ORDER", FD_order)

        Ricci_string = fin.FD_outputC("returnstring", Ricci_SymbExpressions,
                                       params="outCverbose=False,SIMD_enable=True,GoldenKernelsEnable=True")

        filename = "BSSN_Ricci_FD_order_"+str(FD_order)+".h"
        with open(os.path.join(outdir, filename), "w") as file:
            file.write(lp.loop(["i2","i1","i0"],["cctk_nghostzones[2]","cctk_nghostzones[1]","cctk_nghostzones[0]"],
           ["cctk_lsh[2]-cctk_nghostzones[2]","cctk_lsh[1]-cctk_nghostzones[1]","cctk_lsh[0]-cctk_nghostzones[0]"],
                               ["1","1","SIMD_width"],
                                ["#pragma omp parallel for",
                                 "#include \"rfm_files/rfm_struct__SIMD_outer_read2.h\"",
                                 r"""    #include "rfm_files/rfm_struct__SIMD_outer_read1.h"
    #if (defined __INTEL_COMPILER && __INTEL_COMPILER_BUILD_DATE >= 20180804)
        #pragma ivdep         // Forces Intel compiler (if Intel compiler used) to ignore certain SIMD vector dependencies
        #pragma vector always // Forces Intel compiler (if Intel compiler used) to vectorize
    #endif"""],"",
                    "#include \"rfm_files/rfm_struct__SIMD_inner_read0.h\"\n"+Ricci_string))

        # Restore original finite-differencing order:
        par.set_parval_from_str("finite_difference::FD_CENTDERIVS_ORDER", FD_order_orig)

        end = time.time()
        print("(BENCH) Finished Ricci C codegen (FD_order="+str(FD_order)+") in " + str(end - start) + " seconds.")

    def BSSN_constraints__generate_symbolic_expressions_and_C_code():
        ######################################
        # START: GENERATE SYMBOLIC EXPRESSIONS
        # Define the Hamiltonian constraint and output the optimized C code.
        import BSSN.BSSN_constraints as bssncon

        bssncon.BSSN_constraints(add_T4UUmunu_source_terms=False)
        if T4UU != None:
            import BSSN.BSSN_stress_energy_source_terms as Bsest
            Bsest.BSSN_source_terms_for_BSSN_RHSs(T4UU)
            Bsest.BSSN_source_terms_for_BSSN_constraints(T4UU)
            bssncon.H += Bsest.sourceterm_H
            for i in range(3):
                bssncon.MU[i] += Bsest.sourceterm_MU[i]
        # END: GENERATE SYMBOLIC EXPRESSIONS
        ######################################

        # Store original finite-differencing order:
        FD_order_orig = par.parval_from_str("finite_difference::FD_CENTDERIVS_ORDER")
        # Set new finite-differencing order:
        par.set_parval_from_str("finite_difference::FD_CENTDERIVS_ORDER", FD_order)

        start = time.time()
        print("Generating optimized C code for Ham. & mom. constraints. May take a while, depending on CoordSystem.")
        Ham_mom_string = fin.FD_outputC("returnstring",
                                        [lhrh(lhs=gri.gfaccess("aux_gfs", "H"),   rhs=bssncon.H),
                                         lhrh(lhs=gri.gfaccess("aux_gfs", "MU0"), rhs=bssncon.MU[0]),
                                         lhrh(lhs=gri.gfaccess("aux_gfs", "MU1"), rhs=bssncon.MU[1]),
                                         lhrh(lhs=gri.gfaccess("aux_gfs", "MU2"), rhs=bssncon.MU[2])],
                                        params="outCverbose=False")

        with open(os.path.join(outdir,"BSSN_constraints_enable_Tmunu_"+str(enable_stress_energy_source_terms)+"_FD_order_"+str(FD_order)+".h"), "w") as file:
            file.write(lp.loop(["i2","i1","i0"],["cctk_nghostzones[2]","cctk_nghostzones[1]","cctk_nghostzones[0]"],
           ["cctk_lsh[2]-cctk_nghostzones[2]","cctk_lsh[1]-cctk_nghostzones[1]","cctk_lsh[0]-cctk_nghostzones[0]"],
                               ["1","1","1"],["#pragma omp parallel for","",""], "", Ham_mom_string))

        # Restore original finite-differencing order:
        par.set_parval_from_str("finite_difference::FD_CENTDERIVS_ORDER", FD_order_orig)

        end = time.time()
        print("(BENCH) Finished Hamiltonian & momentum constraint C codegen (FD_order="+str(FD_order)+",Tmunu="+str(enable_stress_energy_source_terms)+") in " + str(end - start) + " seconds.")

    def enforce_detgammabar_eq_detgammahat__generate_symbolic_expressions_and_C_code():
        ######################################
        # START: GENERATE SYMBOLIC EXPRESSIONS

        # Enable rfm_precompute infrastructure, which results in
        #   BSSN RHSs that are free of transcendental functions,
        #   even in curvilinear coordinates, so long as
        #   ConformalFactor is set to "W" (default).
        par.set_parval_from_str("reference_metric::enable_rfm_precompute","True")
        par.set_parval_from_str("reference_metric::rfm_precompute_Ccode_outdir",os.path.join(outdir,"rfm_files/"))

        import BSSN.Enforce_Detgammabar_Constraint as EGC
        enforce_detg_constraint_symb_expressions = EGC.Enforce_Detgammabar_Constraint_symb_expressions()

        # Now that we are finished with all the rfm hatted
        #           quantities in generic precomputed functional
        #           form, let's restore them to their closed-
        #           form expressions.
        par.set_parval_from_str("reference_metric::enable_rfm_precompute","False") # Reset to False to disable rfm_precompute.
        rfm.ref_metric__hatted_quantities()
        # END: GENERATE SYMBOLIC EXPRESSIONS
        ######################################

        start = time.time()
        print("Generating optimized C code (FD_order="+str(FD_order)+") for gamma constraint. May take a while, depending on CoordSystem.")
        enforce_gammadet_string = fin.FD_outputC("returnstring", enforce_detg_constraint_symb_expressions,
                                                 params="outCverbose=False,preindent=0,includebraces=False")

        with open(os.path.join(outdir,"enforcedetgammabar_constraint.h"), "w") as file:
            file.write(lp.loop(["i2","i1","i0"],["0", "0", "0"],
                               ["cctk_lsh[2]","cctk_lsh[1]","cctk_lsh[0]"],
                               ["1","1","1"],
                                ["#pragma omp parallel for",
                                     "#include \"rfm_files/rfm_struct__read2.h\"",
                                     "#include \"rfm_files/rfm_struct__read1.h\""],"",
                                     "#include \"rfm_files/rfm_struct__read0.h\"\n"+enforce_gammadet_string))
        end = time.time()
        print("(BENCH) Finished gamma constraint C codegen (FD_order="+str(FD_order)+") in " + str(end - start) + " seconds.")

    if WhichPart=="BSSN_RHSs":
        BSSN_RHSs__generate_Ccode(BSSN_RHSs__generate_symbolic_expressions())
    elif WhichPart=="Ricci":
        Ricci__generate_Ccode(Ricci__generate_symbolic_expressions())
    elif WhichPart=="BSSN_constraints":
        BSSN_constraints__generate_symbolic_expressions_and_C_code()
    elif WhichPart=="detgammabar_constraint":
        enforce_detgammabar_eq_detgammahat__generate_symbolic_expressions_and_C_code()
    else:
        print("Error: WhichPart = "+WhichPart+" is not recognized.")
        sys.exit(1)

    # Store all NRPy+ environment variables to an output string so NRPy+ environment from within this subprocess can be easily restored
    import pickle  # Standard Python module for converting arbitrary data structures to a uniform format.
    # https://www.pythonforthelab.com/blog/storing-binary-data-and-serializing/
    outstr = []
    outstr.append(pickle.dumps(len(gri.glb_gridfcs_list)))
    for lst in gri.glb_gridfcs_list:
        outstr.append(pickle.dumps(lst.gftype))
        outstr.append(pickle.dumps(lst.name))
        outstr.append(pickle.dumps(lst.rank))
        outstr.append(pickle.dumps(lst.DIM))

    outstr.append(pickle.dumps(len(par.glb_params_list)))
    for lst in par.glb_params_list:
        outstr.append(pickle.dumps(lst.type))
        outstr.append(pickle.dumps(lst.module))
        outstr.append(pickle.dumps(lst.parname))
        outstr.append(pickle.dumps(lst.defaultval))

    outstr.append(pickle.dumps(len(par.glb_Cparams_list)))
    for lst in par.glb_Cparams_list:
        outstr.append(pickle.dumps(lst.type))
        outstr.append(pickle.dumps(lst.module))
        outstr.append(pickle.dumps(lst.parname))
        outstr.append(pickle.dumps(lst.defaultval))

    outstr.append(pickle.dumps(len(outC_function_dict)))
    for Cfuncname, Cfunc in outC_function_dict.items():
        outstr.append(pickle.dumps(Cfuncname))
        outstr.append(pickle.dumps(Cfunc))

    return outstr
