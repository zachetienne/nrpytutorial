# Step 0: Add NRPy's directory to the path
# https://stackoverflow.com/questions/16780014/import-file-from-parent-directory
import os,sys
nrpy_dir_path = os.path.join("..")
if nrpy_dir_path not in sys.path:
    sys.path.append(nrpy_dir_path)

from outputC import *            # NRPy+: Core C code output module
import finite_difference as fin  # NRPy+: Finite difference C code generation module
import NRPy_param_funcs as par   # NRPy+: Parameter interface
import grid as gri               # NRPy+: Functions having to do with numerical grids
import loop as lp                # NRPy+: Generate C code loops
import indexedexp as ixp         # NRPy+: Symbolic indexed expression (e.g., tensors, vectors, etc.) support
import reference_metric as rfm   # NRPy+: Reference metric support
import cmdline_helper as cmd     # NRPy+: Multi-platform Python command-line interface
import GRHD.equations as GRHD    # NRPy+: Generate general relativistic hydrodynamics equations
import GRFFE.equations as GRFFE  # NRPy+: Generate general relativisitic force-free electrodynamics equations
import GiRaFFE_NRPy.GiRaFFE_NRPy_Metric_Face_Values as FCVAL
import GiRaFFE_NRPy.GiRaFFE_NRPy_PPM as PPM
import GiRaFFE_NRPy.Afield_flux as Af
import GiRaFFE_NRPy.Stilde_flux as Sf
import GiRaFFE_NRPy.GiRaFFE_NRPy_BCs as BC
import GiRaFFE_NRPy.GiRaFFE_NRPy_A2B as A2B
import GiRaFFE_NRPy.GiRaFFE_NRPy_C2P_P2C as C2P_P2C

thismodule = "GiRaFFE_NRPy_Main_Driver"

CoordSystem = "Cartesian"

par.set_parval_from_str("reference_metric::CoordSystem",CoordSystem)
rfm.reference_metric() # Create ReU, ReDD needed for rescaling B-L initial data, generating BSSN RHSs, etc.

# Default Kreiss-Oliger dissipation strength
default_KO_strength = 0.1
diss_strength = par.Cparameters("REAL", thismodule, "diss_strength", default_KO_strength)

def GiRaFFE_NRPy_Main_Driver_generate_all(out_dir):
    cmd.mkdir(out_dir)
    
    gammaDD = ixp.register_gridfunctions_for_single_rank2("AUXEVOL","gammaDD","sym01",DIM=3)
    betaU = ixp.register_gridfunctions_for_single_rank1("AUXEVOL","betaU",DIM=3)
    alpha = gri.register_gridfunctions("AUXEVOL","alpha")
    AD = ixp.register_gridfunctions_for_single_rank1("EVOL","AD")
    BU = ixp.register_gridfunctions_for_single_rank1("AUXEVOL","BU")
    ValenciavU = ixp.register_gridfunctions_for_single_rank1("AUXEVOL","ValenciavU")
    psi6Phi = gri.register_gridfunctions("EVOL","psi6Phi")
    StildeD = ixp.register_gridfunctions_for_single_rank1("EVOL","StildeD")

    PhievolParenU = ixp.register_gridfunctions_for_single_rank1("AUXEVOL","PhievolParenU",DIM=3)
    AevolParen = gri.register_gridfunctions("AUXEVOL","AevolParen")

    GRHD.compute_sqrtgammaDET(gammaDD)
    GRFFE.compute_AD_source_term_parenthetical_for_FD(GRHD.sqrtgammaDET,betaU,alpha,psi6Phi,AD)
    GRFFE.compute_psi6Phi_rhs_parenthetical(gammaDD,GRHD.sqrtgammaDET,betaU,alpha,AD,psi6Phi)

    parens_to_print = [\
                       lhrh(lhs=gri.gfaccess("auxevol_gfs","AevolParen"),rhs=GRFFE.AevolParen),\
                       lhrh(lhs=gri.gfaccess("auxevol_gfs","PhievolParenU0"),rhs=GRFFE.PhievolParenU[0]),\
                       lhrh(lhs=gri.gfaccess("auxevol_gfs","PhievolParenU1"),rhs=GRFFE.PhievolParenU[1]),\
                       lhrh(lhs=gri.gfaccess("auxevol_gfs","PhievolParenU2"),rhs=GRFFE.PhievolParenU[2]),\
                      ]

    subdir = "RHSs"
    cmd.mkdir(os.path.join(out_dir, subdir))
    desc = "Calculate quantities to be finite-differenced for the GRFFE RHSs"
    name = "calculate_parentheticals_for_RHSs"
    outCfunction(
        outfile  = os.path.join(out_dir,subdir,name+".h"), desc=desc, name=name,
        params   ="const paramstruct *restrict params,const REAL *restrict in_gfs,REAL *restrict auxevol_gfs",
        body     = fin.FD_outputC("returnstring",parens_to_print,params="outCverbose=False").replace("IDX4","IDX4S"),
        loopopts ="AllPoints",
        rel_path_for_Cparams=os.path.join("../"))

    xi_damping = par.Cparameters("REAL",thismodule,"xi_damping",0.1)
    GRFFE.compute_psi6Phi_rhs_damping_term(alpha,psi6Phi,xi_damping)

    AevolParen_dD = ixp.declarerank1("AevolParen_dD",DIM=3)
    PhievolParenU_dD = ixp.declarerank2("PhievolParenU_dD","nosym",DIM=3)

    A_rhsD = ixp.zerorank1()
    psi6Phi_rhs = GRFFE.psi6Phi_damping

    for i in range(3):
        A_rhsD[i] += -AevolParen_dD[i]
        psi6Phi_rhs += -PhievolParenU_dD[i][i]

    # Add Kreiss-Oliger dissipation to the GRFFE RHSs:
#     psi6Phi_dKOD = ixp.declarerank1("psi6Phi_dKOD")
#     AD_dKOD    = ixp.declarerank2("AD_dKOD","nosym")
#     for i in range(3):
#         psi6Phi_rhs += diss_strength*psi6Phi_dKOD[i]*rfm.ReU[i] # ReU[i] = 1/scalefactor_orthog_funcform[i]
#         for j in range(3):
#             A_rhsD[j] += diss_strength*AD_dKOD[j][i]*rfm.ReU[i] # ReU[i] = 1/scalefactor_orthog_funcform[i]

    RHSs_to_print = [\
                     lhrh(lhs=gri.gfaccess("rhs_gfs","AD0"),rhs=A_rhsD[0]),\
                     lhrh(lhs=gri.gfaccess("rhs_gfs","AD1"),rhs=A_rhsD[1]),\
                     lhrh(lhs=gri.gfaccess("rhs_gfs","AD2"),rhs=A_rhsD[2]),\
                     lhrh(lhs=gri.gfaccess("rhs_gfs","psi6Phi"),rhs=psi6Phi_rhs),\
                    ]

    desc = "Calculate AD gauge term and psi6Phi RHSs"
    name = "calculate_AD_gauge_psi6Phi_RHSs"
    source_Ccode = outCfunction(
        outfile  = "returnstring", desc=desc, name=name,
        params   ="const paramstruct *params,const REAL *in_gfs,const REAL *auxevol_gfs,REAL *rhs_gfs",
        body     = fin.FD_outputC("returnstring",RHSs_to_print,params="outCverbose=False").replace("IDX4","IDX4S"),
        loopopts ="InteriorPoints",
        rel_path_for_Cparams=os.path.join("../")).replace("=NGHOSTS","=NGHOSTS_A2B").replace("NGHOSTS+Nxx0","Nxx_plus_2NGHOSTS0-NGHOSTS_A2B").replace("NGHOSTS+Nxx1","Nxx_plus_2NGHOSTS1-NGHOSTS_A2B").replace("NGHOSTS+Nxx2","Nxx_plus_2NGHOSTS2-NGHOSTS_A2B")
    # Note the above .replace() functions. These serve to expand the loop range into the ghostzones, since 
    # the second-order FD needs fewer than some other algorithms we use do. 
    with open(os.path.join(out_dir,subdir,name+".h"),"w") as file:
        file.write(source_Ccode)
    
    # Declare all the Cparameters we will need
    metricderivDDD = ixp.declarerank3("metricderivDDD","sym01",DIM=3)
    shiftderivUD = ixp.declarerank2("shiftderivUD","nosym",DIM=3)
    lapsederivD = ixp.declarerank1("lapsederivD",DIM=3)

    general_access = """const REAL gammaDD00 = auxevol_gfs[IDX4S(GAMMADD00GF,i0,i1,i2)];
const REAL gammaDD01 = auxevol_gfs[IDX4S(GAMMADD01GF,i0,i1,i2)];
const REAL gammaDD02 = auxevol_gfs[IDX4S(GAMMADD02GF,i0,i1,i2)];
const REAL gammaDD11 = auxevol_gfs[IDX4S(GAMMADD11GF,i0,i1,i2)];
const REAL gammaDD12 = auxevol_gfs[IDX4S(GAMMADD12GF,i0,i1,i2)];
const REAL gammaDD22 = auxevol_gfs[IDX4S(GAMMADD22GF,i0,i1,i2)];
const REAL betaU0 = auxevol_gfs[IDX4S(BETAU0GF,i0,i1,i2)];
const REAL betaU1 = auxevol_gfs[IDX4S(BETAU1GF,i0,i1,i2)];
const REAL betaU2 = auxevol_gfs[IDX4S(BETAU2GF,i0,i1,i2)];
const REAL alpha = auxevol_gfs[IDX4S(ALPHAGF,i0,i1,i2)];
const REAL ValenciavU0 = auxevol_gfs[IDX4S(VALENCIAVU0GF,i0,i1,i2)];
const REAL ValenciavU1 = auxevol_gfs[IDX4S(VALENCIAVU1GF,i0,i1,i2)];
const REAL ValenciavU2 = auxevol_gfs[IDX4S(VALENCIAVU2GF,i0,i1,i2)];
const REAL BU0 = auxevol_gfs[IDX4S(BU0GF,i0,i1,i2)];
const REAL BU1 = auxevol_gfs[IDX4S(BU1GF,i0,i1,i2)];
const REAL BU2 = auxevol_gfs[IDX4S(BU2GF,i0,i1,i2)];
"""
    metric_deriv_access = ixp.zerorank1(DIM=3)
    metric_deriv_access[0] = """const REAL metricderivDDD000 = (auxevol_gfs[IDX4S(GAMMA_FACEDD00GF,i0+1,i1,i2)]-auxevol_gfs[IDX4S(GAMMA_FACEDD00GF,i0,i1,i2)])/dxx0;
const REAL metricderivDDD010 = (auxevol_gfs[IDX4S(GAMMA_FACEDD01GF,i0+1,i1,i2)]-auxevol_gfs[IDX4S(GAMMA_FACEDD01GF,i0,i1,i2)])/dxx0;
const REAL metricderivDDD020 = (auxevol_gfs[IDX4S(GAMMA_FACEDD02GF,i0+1,i1,i2)]-auxevol_gfs[IDX4S(GAMMA_FACEDD02GF,i0,i1,i2)])/dxx0;
const REAL metricderivDDD110 = (auxevol_gfs[IDX4S(GAMMA_FACEDD11GF,i0+1,i1,i2)]-auxevol_gfs[IDX4S(GAMMA_FACEDD11GF,i0,i1,i2)])/dxx0;
const REAL metricderivDDD120 = (auxevol_gfs[IDX4S(GAMMA_FACEDD12GF,i0+1,i1,i2)]-auxevol_gfs[IDX4S(GAMMA_FACEDD12GF,i0,i1,i2)])/dxx0;
const REAL metricderivDDD220 = (auxevol_gfs[IDX4S(GAMMA_FACEDD22GF,i0+1,i1,i2)]-auxevol_gfs[IDX4S(GAMMA_FACEDD22GF,i0,i1,i2)])/dxx0;
const REAL shiftderivUD00 = (auxevol_gfs[IDX4S(BETA_FACEU0GF,i0+1,i1,i2)]-auxevol_gfs[IDX4S(BETA_FACEU0GF,i0,i1,i2)])/dxx0;
const REAL shiftderivUD10 = (auxevol_gfs[IDX4S(BETA_FACEU1GF,i0+1,i1,i2)]-auxevol_gfs[IDX4S(BETA_FACEU1GF,i0,i1,i2)])/dxx0;
const REAL shiftderivUD20 = (auxevol_gfs[IDX4S(BETA_FACEU2GF,i0+1,i1,i2)]-auxevol_gfs[IDX4S(BETA_FACEU2GF,i0,i1,i2)])/dxx0;
const REAL lapsederivD0 = (auxevol_gfs[IDX4S(ALPHA_FACEGF,i0+1,i1,i2)]-auxevol_gfs[IDX4S(ALPHA_FACEGF,i0,i1,i2)])/dxx0;
REAL Stilde_rhsD0;
"""
    metric_deriv_access[1] = """const REAL metricderivDDD001 = (auxevol_gfs[IDX4S(GAMMA_FACEDD00GF,i0,i1+1,i2)]-auxevol_gfs[IDX4S(GAMMA_FACEDD00GF,i0,i1,i2)])/dxx1;
const REAL metricderivDDD011 = (auxevol_gfs[IDX4S(GAMMA_FACEDD01GF,i0,i1+1,i2)]-auxevol_gfs[IDX4S(GAMMA_FACEDD01GF,i0,i1,i2)])/dxx1;
const REAL metricderivDDD021 = (auxevol_gfs[IDX4S(GAMMA_FACEDD02GF,i0,i1+1,i2)]-auxevol_gfs[IDX4S(GAMMA_FACEDD02GF,i0,i1,i2)])/dxx1;
const REAL metricderivDDD111 = (auxevol_gfs[IDX4S(GAMMA_FACEDD11GF,i0,i1+1,i2)]-auxevol_gfs[IDX4S(GAMMA_FACEDD11GF,i0,i1,i2)])/dxx1;
const REAL metricderivDDD121 = (auxevol_gfs[IDX4S(GAMMA_FACEDD12GF,i0,i1+1,i2)]-auxevol_gfs[IDX4S(GAMMA_FACEDD12GF,i0,i1,i2)])/dxx1;
const REAL metricderivDDD221 = (auxevol_gfs[IDX4S(GAMMA_FACEDD22GF,i0,i1+1,i2)]-auxevol_gfs[IDX4S(GAMMA_FACEDD22GF,i0,i1,i2)])/dxx1;
const REAL shiftderivUD01 = (auxevol_gfs[IDX4S(BETA_FACEU0GF,i0,i1+1,i2)]-auxevol_gfs[IDX4S(BETA_FACEU0GF,i0,i1,i2)])/dxx1;
const REAL shiftderivUD11 = (auxevol_gfs[IDX4S(BETA_FACEU1GF,i0,i1+1,i2)]-auxevol_gfs[IDX4S(BETA_FACEU1GF,i0,i1,i2)])/dxx1;
const REAL shiftderivUD21 = (auxevol_gfs[IDX4S(BETA_FACEU2GF,i0,i1+1,i2)]-auxevol_gfs[IDX4S(BETA_FACEU2GF,i0,i1,i2)])/dxx1;
const REAL lapsederivD1 = (auxevol_gfs[IDX4S(ALPHA_FACEGF,i0,i1+1,i2)]-auxevol_gfs[IDX4S(ALPHA_FACEGF,i0,i1,i2)])/dxx1;
REAL Stilde_rhsD1;
"""
    metric_deriv_access[2] = """const REAL metricderivDDD002 = (auxevol_gfs[IDX4S(GAMMA_FACEDD00GF,i0,i1,i2+1)]-auxevol_gfs[IDX4S(GAMMA_FACEDD00GF,i0,i1,i2)])/dxx2;
const REAL metricderivDDD012 = (auxevol_gfs[IDX4S(GAMMA_FACEDD01GF,i0,i1,i2+1)]-auxevol_gfs[IDX4S(GAMMA_FACEDD01GF,i0,i1,i2)])/dxx2;
const REAL metricderivDDD022 = (auxevol_gfs[IDX4S(GAMMA_FACEDD02GF,i0,i1,i2+1)]-auxevol_gfs[IDX4S(GAMMA_FACEDD02GF,i0,i1,i2)])/dxx2;
const REAL metricderivDDD112 = (auxevol_gfs[IDX4S(GAMMA_FACEDD11GF,i0,i1,i2+1)]-auxevol_gfs[IDX4S(GAMMA_FACEDD11GF,i0,i1,i2)])/dxx2;
const REAL metricderivDDD122 = (auxevol_gfs[IDX4S(GAMMA_FACEDD12GF,i0,i1,i2+1)]-auxevol_gfs[IDX4S(GAMMA_FACEDD12GF,i0,i1,i2)])/dxx2;
const REAL metricderivDDD222 = (auxevol_gfs[IDX4S(GAMMA_FACEDD22GF,i0,i1,i2+1)]-auxevol_gfs[IDX4S(GAMMA_FACEDD22GF,i0,i1,i2)])/dxx2;
const REAL shiftderivUD02 = (auxevol_gfs[IDX4S(BETA_FACEU0GF,i0,i1,i2+1)]-auxevol_gfs[IDX4S(BETA_FACEU0GF,i0,i1,i2)])/dxx2;
const REAL shiftderivUD12 = (auxevol_gfs[IDX4S(BETA_FACEU1GF,i0,i1,i2+1)]-auxevol_gfs[IDX4S(BETA_FACEU1GF,i0,i1,i2)])/dxx2;
const REAL shiftderivUD22 = (auxevol_gfs[IDX4S(BETA_FACEU2GF,i0,i1,i2+1)]-auxevol_gfs[IDX4S(BETA_FACEU2GF,i0,i1,i2)])/dxx2;
const REAL lapsederivD2 = (auxevol_gfs[IDX4S(ALPHA_FACEGF,i0,i1,i2+1)]-auxevol_gfs[IDX4S(ALPHA_FACEGF,i0,i1,i2)])/dxx2;
REAL Stilde_rhsD2;
"""
    write_final_quantity = ixp.zerorank1(DIM=3)
    write_final_quantity[0] = """rhs_gfs[IDX4S(STILDED0GF,i0,i1,i2)] += Stilde_rhsD0;
"""
    write_final_quantity[1] = """rhs_gfs[IDX4S(STILDED1GF,i0,i1,i2)] += Stilde_rhsD1;
"""
    write_final_quantity[2] = """rhs_gfs[IDX4S(STILDED2GF,i0,i1,i2)] += Stilde_rhsD2;
"""

    # Declare this symbol:
    sqrt4pi = par.Cparameters("REAL",thismodule,"sqrt4pi","sqrt(4.0*M_PI)")

    # We need to rerun a few of these functions with the reset lists to make sure these functions
    # don't cheat by using analytic expressions
    GRHD.u4U_in_terms_of_ValenciavU__rescale_ValenciavU_by_applying_speed_limit(alpha, betaU, gammaDD, ValenciavU)
    GRFFE.compute_smallb4U(gammaDD, betaU, alpha, GRHD.u4U_ito_ValenciavU, BU, sqrt4pi)
    GRFFE.compute_smallbsquared(gammaDD, betaU, alpha, GRFFE.smallb4U)
    GRFFE.compute_TEM4UU(gammaDD,betaU,alpha, GRFFE.smallb4U, GRFFE.smallbsquared,GRHD.u4U_ito_ValenciavU)
    GRHD.compute_g4DD_zerotimederiv_dD(gammaDD,betaU,alpha, metricderivDDD,shiftderivUD,lapsederivD)
    GRHD.compute_S_tilde_source_termD(alpha, GRHD.sqrtgammaDET,GRHD.g4DD_zerotimederiv_dD, GRFFE.TEM4UU)
    for i in range(3):
        desc = "Adds the source term to StildeD"+str(i)+"."
        name = "calculate_StildeD"+str(i)+"_source_term"
        outCfunction(
            outfile  = os.path.join(out_dir,subdir,name+".h"), desc=desc, name=name,
            params   ="const paramstruct *params,const REAL *auxevol_gfs, REAL *rhs_gfs",
            body     = general_access \
                      +metric_deriv_access[i]\
                      +outputC(GRHD.S_tilde_source_termD[i],"Stilde_rhsD"+str(i),"returnstring",params="outCverbose=False").replace("IDX4","IDX4S")\
                      +write_final_quantity[i],
            loopopts ="InteriorPoints",
            rel_path_for_Cparams=os.path.join("../"))

    subdir = "FCVAL"
    cmd.mkdir(os.path.join(out_dir, subdir))
    FCVAL.GiRaFFE_NRPy_FCVAL(os.path.join(out_dir,subdir))
    
    subdir = "PPM"
    cmd.mkdir(os.path.join(out_dir, subdir))
    PPM.GiRaFFE_NRPy_PPM(os.path.join(out_dir,subdir))
    

    # We will pass values of the gridfunction on the cell faces into the function. This requires us
    # to declare them as C parameters in NRPy+. We will denote this with the _face infix/suffix.
    alpha_face = gri.register_gridfunctions("AUXEVOL","alpha_face")
    gamma_faceDD = ixp.register_gridfunctions_for_single_rank2("AUXEVOL","gamma_faceDD","sym01")
    beta_faceU = ixp.register_gridfunctions_for_single_rank1("AUXEVOL","beta_faceU")

    # We'll need some more gridfunctions, now, to represent the reconstructions of BU and ValenciavU
    # on the right and left faces
    Valenciav_rU = ixp.register_gridfunctions_for_single_rank1("AUXEVOL","Valenciav_rU",DIM=3)
    B_rU = ixp.register_gridfunctions_for_single_rank1("AUXEVOL","B_rU",DIM=3)
    Valenciav_lU = ixp.register_gridfunctions_for_single_rank1("AUXEVOL","Valenciav_lU",DIM=3)
    B_lU = ixp.register_gridfunctions_for_single_rank1("AUXEVOL","B_lU",DIM=3)

    subdir = "RHSs"
    Af.generate_Afield_flux_function_files(out_dir,subdir,alpha_face,gamma_faceDD,beta_faceU,\
                                           Valenciav_rU,B_rU,Valenciav_lU,B_lU,True)

    Memory_Read = """const double alpha_face = auxevol_gfs[IDX4S(ALPHA_FACEGF, i0,i1,i2)];
const double gamma_faceDD00 = auxevol_gfs[IDX4S(GAMMA_FACEDD00GF, i0,i1,i2)];
const double gamma_faceDD01 = auxevol_gfs[IDX4S(GAMMA_FACEDD01GF, i0,i1,i2)];
const double gamma_faceDD02 = auxevol_gfs[IDX4S(GAMMA_FACEDD02GF, i0,i1,i2)];
const double gamma_faceDD11 = auxevol_gfs[IDX4S(GAMMA_FACEDD11GF, i0,i1,i2)];
const double gamma_faceDD12 = auxevol_gfs[IDX4S(GAMMA_FACEDD12GF, i0,i1,i2)];
const double gamma_faceDD22 = auxevol_gfs[IDX4S(GAMMA_FACEDD22GF, i0,i1,i2)];
const double beta_faceU0 = auxevol_gfs[IDX4S(BETA_FACEU0GF, i0,i1,i2)];
const double beta_faceU1 = auxevol_gfs[IDX4S(BETA_FACEU1GF, i0,i1,i2)];
const double beta_faceU2 = auxevol_gfs[IDX4S(BETA_FACEU2GF, i0,i1,i2)];
const double Valenciav_rU0 = auxevol_gfs[IDX4S(VALENCIAV_RU0GF, i0,i1,i2)];
const double Valenciav_rU1 = auxevol_gfs[IDX4S(VALENCIAV_RU1GF, i0,i1,i2)];
const double Valenciav_rU2 = auxevol_gfs[IDX4S(VALENCIAV_RU2GF, i0,i1,i2)];
const double B_rU0 = auxevol_gfs[IDX4S(B_RU0GF, i0,i1,i2)];
const double B_rU1 = auxevol_gfs[IDX4S(B_RU1GF, i0,i1,i2)];
const double B_rU2 = auxevol_gfs[IDX4S(B_RU2GF, i0,i1,i2)];
const double Valenciav_lU0 = auxevol_gfs[IDX4S(VALENCIAV_LU0GF, i0,i1,i2)];
const double Valenciav_lU1 = auxevol_gfs[IDX4S(VALENCIAV_LU1GF, i0,i1,i2)];
const double Valenciav_lU2 = auxevol_gfs[IDX4S(VALENCIAV_LU2GF, i0,i1,i2)];
const double B_lU0 = auxevol_gfs[IDX4S(B_LU0GF, i0,i1,i2)];
const double B_lU1 = auxevol_gfs[IDX4S(B_LU1GF, i0,i1,i2)];
const double B_lU2 = auxevol_gfs[IDX4S(B_LU2GF, i0,i1,i2)];
REAL Stilde_fluxD0 = 0; REAL Stilde_fluxD1 = 0; REAL Stilde_fluxD2 = 0;
"""
    Memory_Write = """rhs_gfs[IDX4S(STILDED0GF, i0, i1, i2)] += invdx0*Stilde_fluxD0;
rhs_gfs[IDX4S(STILDED1GF, i0, i1, i2)] += invdx0*Stilde_fluxD1;
rhs_gfs[IDX4S(STILDED2GF, i0, i1, i2)] += invdx0*Stilde_fluxD2;
"""

    indices = ["i0","i1","i2"]
    indicesp1 = ["i0+1","i1+1","i2+1"]
    assignment = "+="
    assignmentp1 = "-="
    invdx = ["invdx0","invdx1","invdx2"]

    for flux_dirn in range(3):
        Sf.calculate_Stilde_flux(flux_dirn,True,alpha_face,gamma_faceDD,beta_faceU,\
                                 Valenciav_rU,B_rU,Valenciav_lU,B_lU,sqrt4pi)
        Stilde_flux_to_print = [\
                                Sf.Stilde_fluxD[0],\
                                Sf.Stilde_fluxD[1],\
                                Sf.Stilde_fluxD[2],\
                               ]
        Stilde_flux_names = [\
                             "Stilde_fluxD0",\
                             "Stilde_fluxD1",\
                             "Stilde_fluxD2",\
                            ]

        desc = "Compute the flux of all 3 components of tilde{S}_i on the right face in the " + str(flux_dirn) + "."
        name = "calculate_Stilde_flux_D" + str(flux_dirn) + "_right"
        outCfunction(
            outfile  = os.path.join(out_dir,subdir,name+".h"), desc=desc, name=name,
            params   ="const paramstruct *params,const REAL *auxevol_gfs,REAL *rhs_gfs",
            body     =  Memory_Read \
                   +outputC(Stilde_flux_to_print,Stilde_flux_names,"returnstring",params="outCverbose=False").replace("IDX4","IDX4S")\
                       +Memory_Write.replace(invdx[0],invdx[flux_dirn]),
            loopopts ="InteriorPoints",
            rel_path_for_Cparams=os.path.join("../"))

        desc = "Compute the flux of all 3 components of tilde{S}_i on the left face in the " + str(flux_dirn) + "."
        name = "calculate_Stilde_flux_D" + str(flux_dirn) + "_left"
        outCfunction(
            outfile  = os.path.join(out_dir,subdir,name+".h"), desc=desc, name=name,
            params   ="const paramstruct *params,const REAL *auxevol_gfs,REAL *rhs_gfs",
            body     =  Memory_Read.replace(indices[flux_dirn],indicesp1[flux_dirn]) \
                   +outputC(Stilde_flux_to_print,Stilde_flux_names,"returnstring",params="outCverbose=False").replace("IDX4","IDX4S")\
                       +Memory_Write.replace(invdx[0],invdx[flux_dirn]).replace(assignment,assignmentp1),
            loopopts ="InteriorPoints",
            rel_path_for_Cparams=os.path.join("../"))

    subdir = "boundary_conditions"
    cmd.mkdir(os.path.join(out_dir,subdir))
    BC.GiRaFFE_NRPy_BCs(os.path.join(out_dir,subdir))
    
    subdir = "A2B"
    cmd.mkdir(os.path.join(out_dir,subdir))
    A2B.GiRaFFE_NRPy_A2B(os.path.join(out_dir,subdir),gammaDD,AD,BU)
    
    C2P_P2C.GiRaFFE_NRPy_C2P(StildeD,BU,gammaDD,betaU,alpha)

    values_to_print = [\
                       lhrh(lhs=gri.gfaccess("in_gfs","StildeD0"),rhs=C2P_P2C.outStildeD[0]),\
                       lhrh(lhs=gri.gfaccess("in_gfs","StildeD1"),rhs=C2P_P2C.outStildeD[1]),\
                       lhrh(lhs=gri.gfaccess("in_gfs","StildeD2"),rhs=C2P_P2C.outStildeD[2]),\
                       lhrh(lhs=gri.gfaccess("auxevol_gfs","ValenciavU0"),rhs=C2P_P2C.ValenciavU[0]),\
                       lhrh(lhs=gri.gfaccess("auxevol_gfs","ValenciavU1"),rhs=C2P_P2C.ValenciavU[1]),\
                       lhrh(lhs=gri.gfaccess("auxevol_gfs","ValenciavU2"),rhs=C2P_P2C.ValenciavU[2])\
                      ]

    subdir = "C2P"
    cmd.mkdir(os.path.join(out_dir,subdir))
    desc = "Apply fixes to \tilde{S}_i and recompute the velocity to match with current sheet prescription."
    name = "GiRaFFE_NRPy_cons_to_prims"
    outCfunction(
        outfile  = os.path.join(out_dir,subdir,name+".h"), desc=desc, name=name,
        params   ="const paramstruct *params,REAL *xx[3],REAL *auxevol_gfs,REAL *in_gfs",
        body     = fin.FD_outputC("returnstring",values_to_print,params="outCverbose=False").replace("IDX4","IDX4S"),
        loopopts ="AllPoints,Read_xxs",
        rel_path_for_Cparams=os.path.join("../"))

    C2P_P2C.GiRaFFE_NRPy_P2C(gammaDD,betaU,alpha,  ValenciavU,BU, sqrt4pi)

    values_to_print = [\
                       lhrh(lhs=gri.gfaccess("in_gfs","StildeD0"),rhs=C2P_P2C.StildeD[0]),\
                       lhrh(lhs=gri.gfaccess("in_gfs","StildeD1"),rhs=C2P_P2C.StildeD[1]),\
                       lhrh(lhs=gri.gfaccess("in_gfs","StildeD2"),rhs=C2P_P2C.StildeD[2]),\
                      ]

    desc = "Recompute StildeD after current sheet fix to Valencia 3-velocity to ensure consistency between conservative & primitive variables."
    name = "GiRaFFE_NRPy_prims_to_cons"
    outCfunction(
        outfile  = os.path.join(out_dir,subdir,name+".h"), desc=desc, name=name,
        params   ="const paramstruct *params,REAL *auxevol_gfs,REAL *in_gfs",
        body     = fin.FD_outputC("returnstring",values_to_print,params="outCverbose=False").replace("IDX4","IDX4S"),
        loopopts ="AllPoints",
        rel_path_for_Cparams=os.path.join("../"))

    # Write out the main driver itself:
    with open(os.path.join(out_dir,"GiRaFFE_NRPy_Main_Driver.h"),"w") as file:
        file.write("""// Structure to track ghostzones for PPM:
typedef struct __gf_and_gz_struct__ {
  REAL *gf;
  int gz_lo[4],gz_hi[4];
} gf_and_gz_struct;
// Some additional constants needed for PPM:
const int VX=0,VY=1,VZ=2,BX=3,BY=4,BZ=5;
const int NUM_RECONSTRUCT_GFS = 6;

// Include ALL functions needed for evolution
#include "RHSs/calculate_parentheticals_for_RHSs.h"
#include "RHSs/calculate_AD_gauge_psi6Phi_RHSs.h"
#include "PPM/reconstruct_set_of_prims_PPM_GRFFE_NRPy.c"
#include "FCVAL/interpolate_metric_gfs_to_cell_faces.h"
#include "RHSs/calculate_StildeD0_source_term.h"
#include "RHSs/calculate_StildeD1_source_term.h"
#include "RHSs/calculate_StildeD2_source_term.h"
// #include "RHSs/calculate_E_field_D0_right.h"
// #include "RHSs/calculate_E_field_D0_left.h"
// #include "RHSs/calculate_E_field_D1_right.h"
// #include "RHSs/calculate_E_field_D1_left.h"
// #include "RHSs/calculate_E_field_D2_right.h"
// #include "RHSs/calculate_E_field_D2_left.h"
#include "../calculate_E_field_flat_all_in_one.h"
#include "RHSs/calculate_Stilde_flux_D0_right.h"
#include "RHSs/calculate_Stilde_flux_D0_left.h"
#include "RHSs/calculate_Stilde_flux_D1_right.h"
#include "RHSs/calculate_Stilde_flux_D1_left.h"
#include "RHSs/calculate_Stilde_flux_D2_right.h"
#include "RHSs/calculate_Stilde_flux_D2_left.h"
#include "boundary_conditions/GiRaFFE_boundary_conditions.h"
#include "A2B/driver_AtoB.h"
#include "C2P/GiRaFFE_NRPy_cons_to_prims.h"
#include "C2P/GiRaFFE_NRPy_prims_to_cons.h"

void override_BU_with_old_GiRaFFE(const paramstruct *restrict params,REAL *restrict auxevol_gfs,const int n) {
#include "set_Cparameters.h"
    char filename[100];
    sprintf(filename,"BU0_override-%08d.bin",n);
    FILE *out2D = fopen(filename, "rb");
    fread(auxevol_gfs+BU0GF*Nxx_plus_2NGHOSTS0*Nxx_plus_2NGHOSTS1*Nxx_plus_2NGHOSTS2,
          sizeof(double),Nxx_plus_2NGHOSTS0*Nxx_plus_2NGHOSTS1*Nxx_plus_2NGHOSTS2,out2D);
    fclose(out2D);
    sprintf(filename,"BU1_override-%08d.bin",n);
    out2D = fopen(filename, "rb");
    fread(auxevol_gfs+BU1GF*Nxx_plus_2NGHOSTS0*Nxx_plus_2NGHOSTS1*Nxx_plus_2NGHOSTS2,
          sizeof(double),Nxx_plus_2NGHOSTS0*Nxx_plus_2NGHOSTS1*Nxx_plus_2NGHOSTS2,out2D);
    fclose(out2D);
    sprintf(filename,"BU2_override-%08d.bin",n);
    out2D = fopen(filename, "rb");
    fread(auxevol_gfs+BU2GF*Nxx_plus_2NGHOSTS0*Nxx_plus_2NGHOSTS1*Nxx_plus_2NGHOSTS2,
          sizeof(double),Nxx_plus_2NGHOSTS0*Nxx_plus_2NGHOSTS1*Nxx_plus_2NGHOSTS2,out2D);
    fclose(out2D);
}

void GiRaFFE_NRPy_RHSs(const paramstruct *restrict params,REAL *restrict auxevol_gfs,const REAL *restrict in_gfs,REAL *restrict rhs_gfs) {
#include "set_Cparameters.h"
    // First thing's first: initialize the RHSs to zero!
#pragma omp parallel for
    for(int ii=0;ii<Nxx_plus_2NGHOSTS0*Nxx_plus_2NGHOSTS1*Nxx_plus_2NGHOSTS2*NUM_EVOL_GFS;ii++) {
        rhs_gfs[ii] = 0.0;
    }
    // Next calculate the easier source terms that don't require flux directions
    // This will also reset the RHSs for each gf at each new timestep.
    calculate_parentheticals_for_RHSs(params,in_gfs,auxevol_gfs);
    calculate_AD_gauge_psi6Phi_RHSs(params,in_gfs,auxevol_gfs,rhs_gfs);
    
    // Now, we set up a bunch of structs of pointers to properly guide the PPM algorithm.
    // They also count the number of ghostzones available.
    gf_and_gz_struct in_prims[NUM_RECONSTRUCT_GFS], out_prims_r[NUM_RECONSTRUCT_GFS], out_prims_l[NUM_RECONSTRUCT_GFS];
    int which_prims_to_reconstruct[NUM_RECONSTRUCT_GFS],num_prims_to_reconstruct;
    const int Nxxp2NG012 = Nxx_plus_2NGHOSTS0*Nxx_plus_2NGHOSTS1*Nxx_plus_2NGHOSTS2;
    
    REAL *temporary = auxevol_gfs + Nxxp2NG012*AEVOLPARENGF; //We're not using this anymore
    // This sets pointers to the portion of auxevol_gfs containing the relevant gridfunction.
    int ww=0;
    in_prims[ww].gf      = auxevol_gfs + Nxxp2NG012*VALENCIAVU0GF; 
      out_prims_r[ww].gf = auxevol_gfs + Nxxp2NG012*VALENCIAV_RU0GF; 
      out_prims_l[ww].gf = auxevol_gfs + Nxxp2NG012*VALENCIAV_LU0GF; 
    ww++;
    in_prims[ww].gf      = auxevol_gfs + Nxxp2NG012*VALENCIAVU1GF; 
      out_prims_r[ww].gf = auxevol_gfs + Nxxp2NG012*VALENCIAV_RU1GF; 
      out_prims_l[ww].gf = auxevol_gfs + Nxxp2NG012*VALENCIAV_LU1GF; 
    ww++;
    in_prims[ww].gf      = auxevol_gfs + Nxxp2NG012*VALENCIAVU2GF; 
      out_prims_r[ww].gf = auxevol_gfs + Nxxp2NG012*VALENCIAV_RU2GF; 
      out_prims_l[ww].gf = auxevol_gfs + Nxxp2NG012*VALENCIAV_LU2GF; 
    ww++;
    in_prims[ww].gf      = auxevol_gfs + Nxxp2NG012*BU0GF; 
      out_prims_r[ww].gf = auxevol_gfs + Nxxp2NG012*B_RU0GF; 
      out_prims_l[ww].gf = auxevol_gfs + Nxxp2NG012*B_LU0GF; 
    ww++;
    in_prims[ww].gf      = auxevol_gfs + Nxxp2NG012*BU1GF; 
      out_prims_r[ww].gf = auxevol_gfs + Nxxp2NG012*B_RU1GF; 
      out_prims_l[ww].gf = auxevol_gfs + Nxxp2NG012*B_LU1GF; 
    ww++;
    in_prims[ww].gf      = auxevol_gfs + Nxxp2NG012*BU2GF; 
      out_prims_r[ww].gf = auxevol_gfs + Nxxp2NG012*B_RU2GF; 
      out_prims_l[ww].gf = auxevol_gfs + Nxxp2NG012*B_LU2GF; 
    ww++;

    // Prims are defined AT ALL GRIDPOINTS, so we set the # of ghostzones to zero:
    for(int i=0;i<NUM_RECONSTRUCT_GFS;i++) for(int j=1;j<=3;j++) { in_prims[i].gz_lo[j]=0; in_prims[i].gz_hi[j]=0; }
    // Left/right variables are not yet defined, yet we set the # of gz's to zero by default:
    for(int i=0;i<NUM_RECONSTRUCT_GFS;i++) for(int j=1;j<=3;j++) { out_prims_r[i].gz_lo[j]=0; out_prims_r[i].gz_hi[j]=0; }
    for(int i=0;i<NUM_RECONSTRUCT_GFS;i++) for(int j=1;j<=3;j++) { out_prims_l[i].gz_lo[j]=0; out_prims_l[i].gz_hi[j]=0; }

    ww=0;
    which_prims_to_reconstruct[ww]=VX; ww++;
    which_prims_to_reconstruct[ww]=VY; ww++;
    which_prims_to_reconstruct[ww]=VZ; ww++;
    which_prims_to_reconstruct[ww]=BX; ww++;
    which_prims_to_reconstruct[ww]=BY; ww++;
    which_prims_to_reconstruct[ww]=BZ; ww++;
    num_prims_to_reconstruct=ww;

    // In each direction, perform the PPM reconstruction procedure.
    // Then, add the fluxes to the RHS as appropriate.
    int count;
    for(int flux_dirn=0;flux_dirn<3;flux_dirn++) {
        // In each direction, interpolate the metric gfs (gamma,beta,alpha) to cell faces.
        interpolate_metric_gfs_to_cell_faces(params,auxevol_gfs,flux_dirn+1);
        // Then, reconstruct the primitive variables on the cell faces.
        // This function is housed in the file: "reconstruct_set_of_prims_PPM_GRFFE_NRPy.c"
        reconstruct_set_of_prims_PPM_GRFFE_NRPy(params, auxevol_gfs, flux_dirn+1, num_prims_to_reconstruct,                                                          
                                                which_prims_to_reconstruct, in_prims, out_prims_r, out_prims_l, temporary);
        // For example, if flux_dirn==0, then at gamma_faceDD00(i,j,k) represents gamma_{xx}
        // at (i-1/2,j,k), Valenciav_lU0(i,j,k) is the x-component of the velocity at (i-1/2-epsilon,j,k),
        // and Valenciav_rU0(i,j,k) is the x-component of the velocity at (i-1/2+epsilon,j,k).
        
        if(flux_dirn==0) {
            // Next, we calculate the source term for StildeD. Again, this also resets the rhs_gfs array at
            // each new timestep.
            calculate_StildeD0_source_term(params,auxevol_gfs,rhs_gfs);
            // Now, compute the electric field on each face of a cell and add it to the RHSs as appropriate
            //calculate_E_field_D0_right(params,auxevol_gfs,rhs_gfs);
            //calculate_E_field_D0_left(params,auxevol_gfs,rhs_gfs);
            // Finally, we calculate the flux of StildeD and add the appropriate finite-differences 
            // to the RHSs.
            calculate_Stilde_flux_D0_right(params,auxevol_gfs,rhs_gfs);
            calculate_Stilde_flux_D0_left(params,auxevol_gfs,rhs_gfs);
        }
        else if(flux_dirn==1) {
            calculate_StildeD1_source_term(params,auxevol_gfs,rhs_gfs);
            //calculate_E_field_D1_right(params,auxevol_gfs,rhs_gfs);
            //calculate_E_field_D1_left(params,auxevol_gfs,rhs_gfs);
            calculate_Stilde_flux_D1_right(params,auxevol_gfs,rhs_gfs);
            calculate_Stilde_flux_D1_left(params,auxevol_gfs,rhs_gfs);
        }
        else {
            calculate_StildeD2_source_term(params,auxevol_gfs,rhs_gfs);
            //calculate_E_field_D2_right(params,auxevol_gfs,rhs_gfs);
            //calculate_E_field_D2_left(params,auxevol_gfs,rhs_gfs);
            calculate_Stilde_flux_D2_right(params,auxevol_gfs,rhs_gfs);
            calculate_Stilde_flux_D2_left(params,auxevol_gfs,rhs_gfs);
        }
        for(int count=0;count<=1;count++) {
            // This function is written to be general, using notation that matches the forward permutation added to AD2,
            // i.e., [F_HLL^x(B^y)]_z corresponding to flux_dirn=0, count=1. By cyclically permuting with flux_dirn, we 
            // get contributions to the other components, and by incrementing count, we get the backward permutations:
            // flux_dirn | count | [F_HLL^i(B^j)]_k OR A_k += vi*Bj - vj*Bi
            //      0    |    0  |  0 2 1
            //      0    |    1  |  0 1 2
            //      1    |    0  |  1 0 2
            //      1    |    1  |  1 2 0
            //      2    |    0  |  2 1 0
            //      2    |    1  |  2 0 1
            calculate_E_field_flat_all_in_one(params,
              &auxevol_gfs[IDX4ptS(VALENCIAV_RU0GF+(flux_dirn)%3, 0)],&auxevol_gfs[IDX4ptS(VALENCIAV_RU0GF+(flux_dirn-count+2)%3, 0)],
              &auxevol_gfs[IDX4ptS(VALENCIAV_LU0GF+(flux_dirn)%3, 0)],&auxevol_gfs[IDX4ptS(VALENCIAV_LU0GF+(flux_dirn-count+2)%3, 0)],
              &auxevol_gfs[IDX4ptS(B_RU0GF        +(flux_dirn)%3, 0)],&auxevol_gfs[IDX4ptS(B_RU0GF        +(flux_dirn-count+2)%3, 0)],
              &auxevol_gfs[IDX4ptS(B_LU0GF        +(flux_dirn)%3, 0)],&auxevol_gfs[IDX4ptS(B_LU0GF        +(flux_dirn-count+2)%3, 0)],
              &rhs_gfs[IDX4ptS(AD0GF+(flux_dirn+1+count)%3,0)], 2.0*((REAL)count)-1.0, flux_dirn);
// Let's suppose flux_dirn = 0. Then we will need to update Ay (count=0) and Az (count=1):
//     flux_dirn=count=0 -> AD0GF+(flux_dirn+1+count)%3 = AD0GF + (0+1+0)%3=AD1GF <- Updating Ay!
//        (flux_dirn)%3 = (0)%3 = 0               Vx
//        (flux_dirn-count+2)%3 = (0-0+2)%3 = 2   Vz .  Inputs Vx, Vz -> SIGN = -1 ; 2.0*((REAL)count)-1.0=-1 check!
// Let's suppose flux_dirn = 0. Then we will need to update Ay (count=0) and Az (count=1):
//     flux_dirn=0,count=1 -> AD0GF+(flux_dirn+1+count)%3 = AD0GF + (0+1+1)%3=AD2GF <- Updating Az!
//        (flux_dirn)%3 = (0)%3 = 0               Vx
//        (flux_dirn-count+2)%3 = (0-1+2)%3 = 1   Vy .  Inputs Vx, Vy -> SIGN = +1 ; 2.0*((REAL)count)-1.0=2-1=+1 check!
            // SIGN = -1.0 if count=0, 1.0 if count=1
            // This is necessary because 
            // -E_z(x_i,y_j,z_k) &= 0.25 ( [F_HLL^x(B^y)]_z(i+1/2,j,k)+[F_HLL^x(B^y)]_z(i-1/2,j,k)
            //                            -[F_HLL^y(B^x)]_z(i,j+1/2,k)-[F_HLL^y(B^x)]_z(i,j-1/2,k) )
            // Note the negative signs on the reversed permuation terms!

        }

    }
}

void GiRaFFE_NRPy_post_step(const paramstruct *restrict params,REAL *xx[3],REAL *restrict auxevol_gfs,REAL *restrict evol_gfs,const int n) {
    // First, apply BCs to AD and psi6Phi. Then calculate BU from AD
    apply_bcs_potential(params,evol_gfs);
    driver_A_to_B(params,evol_gfs,auxevol_gfs);
    //override_BU_with_old_GiRaFFE(params,auxevol_gfs,n);
    // Apply fixes to StildeD, then recompute the velocity at the new timestep. 
    // Apply the current sheet prescription to the velocities
    GiRaFFE_NRPy_cons_to_prims(params,xx,auxevol_gfs,evol_gfs);
    // Then, recompute StildeD to be consistent with the new velocities
    //GiRaFFE_NRPy_prims_to_cons(params,auxevol_gfs,evol_gfs);
    // Finally, apply outflow boundary conditions to the velocities.
    apply_bcs_velocity(params,auxevol_gfs);
}
""")
