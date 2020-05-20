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
# import GiRaFFE_NRPy.Afield_flux as Af
import GiRaFFE_NRPy.Stilde_flux as Sf
import GiRaFFE_NRPy.GiRaFFE_NRPy_BCs as BC
# import GiRaFFE_NRPy.GiRaFFE_NRPy_A2B as A2B
import GiRaFFE_NRPy.GiRaFFE_NRPy_C2P_P2C as C2P_P2C

thismodule = "GiRaFFE_NRPy_Main_Driver"

CoordSystem = "Cartesian"

par.set_parval_from_str("reference_metric::CoordSystem",CoordSystem)
rfm.reference_metric() # Create ReU, ReDD needed for rescaling B-L initial data, generating BSSN RHSs, etc.

xi_damping = par.Cparameters("REAL",thismodule,"xi_damping",0.1)

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
    BstaggerU = ixp.register_gridfunctions_for_single_rank1("AUXEVOL","BstaggerU")
    ValenciavU = ixp.register_gridfunctions_for_single_rank1("AUXEVOL","ValenciavU")
    psi6Phi = gri.register_gridfunctions("EVOL","psi6Phi")
    StildeD = ixp.register_gridfunctions_for_single_rank1("EVOL","StildeD")

    PhievolParenU = ixp.register_gridfunctions_for_single_rank1("AUXEVOL","PhievolParenU",DIM=3)
    AevolParen = gri.register_gridfunctions("AUXEVOL","AevolParen")

    GRHD.compute_sqrtgammaDET(gammaDD)
    
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
    subdir = "RHSs"
    cmd.mkdir(os.path.join(out_dir,subdir))
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
    
    Valenciav_rrU = ixp.register_gridfunctions_for_single_rank1("AUXEVOL","Valenciav_rrU",DIM=3)
    Valenciav_rlU = ixp.register_gridfunctions_for_single_rank1("AUXEVOL","Valenciav_rlU",DIM=3)
    Valenciav_lrU = ixp.register_gridfunctions_for_single_rank1("AUXEVOL","Valenciav_lrU",DIM=3)
    Valenciav_llU = ixp.register_gridfunctions_for_single_rank1("AUXEVOL","Valenciav_llU",DIM=3)
    Bstagger_rU = ixp.register_gridfunctions_for_single_rank1("AUXEVOL","Bstagger_rU",DIM=3)
    Bstagger_lU = ixp.register_gridfunctions_for_single_rank1("AUXEVOL","Bstagger_lU",DIM=3)


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

    subdir = "RHSs"
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
static const int VX=0,VY=1,VZ=2,
  BX_CENTER=3,BY_CENTER=4,BZ_CENTER=5,BX_STAGGER=6,BY_STAGGER=7,BZ_STAGGER=8,
  VXR=9,VYR=10,VZR=11,VXL=12,VYL=13,VZL=14;  //<-- Be _sure_ to define MAXNUMVARS appropriately!
const int NUM_RECONSTRUCT_GFS = 15;

// Include ALL functions needed for evolution
#include "PPM/reconstruct_set_of_prims_PPM_GRFFE_NRPy.c"
#include "FCVAL/interpolate_metric_gfs_to_cell_faces.h"
#include "RHSs/calculate_StildeD0_source_term.h"
#include "RHSs/calculate_StildeD1_source_term.h"
#include "RHSs/calculate_StildeD2_source_term.h"
#include "../Lorenz_psi6phi_rhs__add_gauge_terms_to_A_i_rhs.h"
#include "../A_i_rhs_no_gauge_terms.h"
#include "../compute_B_and_Bstagger_from_A.h"
#include "RHSs/calculate_Stilde_flux_D0_right.h"
#include "RHSs/calculate_Stilde_flux_D0_left.h"
#include "RHSs/calculate_Stilde_flux_D1_right.h"
#include "RHSs/calculate_Stilde_flux_D1_left.h"
#include "RHSs/calculate_Stilde_flux_D2_right.h"
#include "RHSs/calculate_Stilde_flux_D2_left.h"
#include "boundary_conditions/GiRaFFE_boundary_conditions.h"
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
    in_prims[ww].gf      = auxevol_gfs + Nxxp2NG012*BstaggerU0GF; 
      out_prims_r[ww].gf = auxevol_gfs + Nxxp2NG012*Bstagger_RU0GF; 
      out_prims_l[ww].gf = auxevol_gfs + Nxxp2NG012*Bstagger_LU0GF; 
    ww++;
    in_prims[ww].gf      = auxevol_gfs + Nxxp2NG012*BstaggerU1GF; 
      out_prims_r[ww].gf = auxevol_gfs + Nxxp2NG012*Bstagger_RU1GF; 
      out_prims_l[ww].gf = auxevol_gfs + Nxxp2NG012*Bstagger_LU1GF; 
    ww++;
    in_prims[ww].gf      = auxevol_gfs + Nxxp2NG012*BstaggerU2GF; 
      out_prims_r[ww].gf = auxevol_gfs + Nxxp2NG012*Bstagger_RU2GF; 
      out_prims_l[ww].gf = auxevol_gfs + Nxxp2NG012*Bstagger_LU2GF; 
    ww++;
    in_prims[ww].gf      = auxevol_gfs + Nxxp2NG012*VALENCIAV_RU0GF; 
      out_prims_r[ww].gf = auxevol_gfs + Nxxp2NG012*VALENCIAV_RRU0GF; 
      out_prims_l[ww].gf = auxevol_gfs + Nxxp2NG012*VALENCIAV_RLU0GF; 
    ww++;
    in_prims[ww].gf      = auxevol_gfs + Nxxp2NG012*VALENCIAV_RU1GF; 
      out_prims_r[ww].gf = auxevol_gfs + Nxxp2NG012*VALENCIAV_RRU1GF; 
      out_prims_l[ww].gf = auxevol_gfs + Nxxp2NG012*VALENCIAV_RLU1GF; 
    ww++;
    in_prims[ww].gf      = auxevol_gfs + Nxxp2NG012*VALENCIAV_RU2GF; 
      out_prims_r[ww].gf = auxevol_gfs + Nxxp2NG012*VALENCIAV_RRU2GF; 
      out_prims_l[ww].gf = auxevol_gfs + Nxxp2NG012*VALENCIAV_RLU2GF; 
    ww++;
    in_prims[ww].gf      = auxevol_gfs + Nxxp2NG012*VALENCIAV_LU0GF; 
      out_prims_r[ww].gf = auxevol_gfs + Nxxp2NG012*VALENCIAV_LRU0GF; 
      out_prims_l[ww].gf = auxevol_gfs + Nxxp2NG012*VALENCIAV_LLU0GF; 
    ww++;
    in_prims[ww].gf      = auxevol_gfs + Nxxp2NG012*VALENCIAV_LU1GF; 
      out_prims_r[ww].gf = auxevol_gfs + Nxxp2NG012*VALENCIAV_LRU1GF; 
      out_prims_l[ww].gf = auxevol_gfs + Nxxp2NG012*VALENCIAV_LLU1GF; 
    ww++;
    in_prims[ww].gf      = auxevol_gfs + Nxxp2NG012*VALENCIAV_LU2GF; 
      out_prims_r[ww].gf = auxevol_gfs + Nxxp2NG012*VALENCIAV_LRU2GF; 
      out_prims_l[ww].gf = auxevol_gfs + Nxxp2NG012*VALENCIAV_LLU2GF; 
    ww++;

    // Prims are defined AT ALL GRIDPOINTS, so we set the # of ghostzones to zero:
    for(int i=0;i<NUM_RECONSTRUCT_GFS;i++) for(int j=1;j<=3;j++) { in_prims[i].gz_lo[j]=0; in_prims[i].gz_hi[j]=0; }
    // Left/right variables are not yet defined, yet we set the # of gz's to zero by default:
    for(int i=0;i<NUM_RECONSTRUCT_GFS;i++) for(int j=1;j<=3;j++) { out_prims_r[i].gz_lo[j]=0; out_prims_r[i].gz_hi[j]=0; }
    for(int i=0;i<NUM_RECONSTRUCT_GFS;i++) for(int j=1;j<=3;j++) { out_prims_l[i].gz_lo[j]=0; out_prims_l[i].gz_hi[j]=0; }

  int flux_dirn;
  flux_dirn=0;
  // ftilde = 0 in GRFFE, since P=rho=0.

  /* There are two stories going on here:
   * 1) Computation of \partial_x on RHS of \partial_t {mhd_st_{x,y,z}}, 
   *    via PPM reconstruction onto (i-1/2,j,k), so that 
   *    \partial_y F = [ F(i+1/2,j,k) - F(i-1/2,j,k) ] / dx
   * 2) Computation of \partial_t A_i, where A_i are *staggered* gridfunctions,
   *    where A_x is defined at (i,j+1/2,k+1/2), A_y at (i+1/2,j,k+1/2), etc.
   *    Ai_rhs = \partial_t A_i = \epsilon_{ijk} \psi^{6} (v^j B^k - v^j B^k),
   *    where \epsilon_{ijk} is the flat-space antisymmetric operator.
   * 2A) Az_rhs is defined at (i+1/2,j+1/2,k), and it depends on {Bx,By,vx,vy},
   *     so the trick is to reconstruct {Bx,By,vx,vy} cleverly to get to these
   *     staggered points. For example:
   * 2Aa) vx and vy are at (i,j,k), and we reconstruct them to (i-1/2,j,k) below. After
   *      this, we'll reconstruct again in the y-dir'n to get {vx,vy} at (i-1/2,j-1/2,k)
   * 2Ab) By_stagger is at (i,j+1/2,k), and we reconstruct below to (i-1/2,j+1/2,k). */
  ww=0;
  which_prims_to_reconstruct[ww]=VX;        ww++;
  which_prims_to_reconstruct[ww]=VY;        ww++;
  which_prims_to_reconstruct[ww]=VZ;        ww++;
  //which_prims_to_reconstruct[ww]=BX_CENTER; ww++;
  which_prims_to_reconstruct[ww]=BY_CENTER; ww++;
  which_prims_to_reconstruct[ww]=BZ_CENTER; ww++;
  which_prims_to_reconstruct[ww]=BY_STAGGER;ww++;
  num_prims_to_reconstruct=ww;
  // This function is housed in the file: "reconstruct_set_of_prims_PPM_GRFFE.C"
  reconstruct_set_of_prims_PPM_GRFFE_NRPy(params, auxevol_gfs, flux_dirn+1, num_prims_to_reconstruct,                                                          
                                          which_prims_to_reconstruct, in_prims, out_prims_r, out_prims_l, temporary);
  //Right and left face values of BI_CENTER are used in GRFFE__S_i__flux computation (first to compute b^a).
  //   Instead of reconstructing, we simply set B^x face values to be consistent with BX_STAGGER.
#pragma omp parallel for
  for(int k=0;k<Nxx_plus_2NGHOSTS2;k++) for(int j=0;j<Nxx_plus_2NGHOSTS1;j++) for(int i=0;i<Nxx_plus_2NGHOSTS0;i++) {
        const int index=IDX3S(i,j,k), indexim1=IDX3S(i-1+(i==0),j,k); /* indexim1=0 when i=0 */
        out_prims_r[BX_CENTER].gf[index]=out_prims_l[BX_CENTER].gf[index]=in_prims[BX_STAGGER].gf[indexim1]; }
  // Then add fluxes to RHS for hydro variables {vx,vy,vz}:
  // This function is housed in the file: "add_fluxes_and_source_terms_to_hydro_rhss.C"
  calculate_StildeD0_source_term(params,auxevol_gfs,rhs_gfs);
  calculate_Stilde_flux_D0_right(params,auxevol_gfs,rhs_gfs);
  calculate_Stilde_flux_D0_left(params,auxevol_gfs,rhs_gfs);

  // Note that we have already reconstructed vx and vy along the x-direction, 
  //   at (i-1/2,j,k). That result is stored in v{x,y}{r,l}.  Bx_stagger data
  //   are defined at (i+1/2,j,k).
  // Next goal: reconstruct Bx, vx and vy at (i+1/2,j+1/2,k).
  flux_dirn=1;
  // ftilde = 0 in GRFFE, since P=rho=0.

  // in_prims[{VXR,VXL,VYR,VYL}].gz_{lo,hi} ghostzones are set to all zeros, which
  //    is incorrect. We fix this below.
  // [Note that this is a cheap operation, copying only 8 integers and a pointer.]
  in_prims[VXR]=out_prims_r[VX];
  in_prims[VXL]=out_prims_l[VX];
  in_prims[VYR]=out_prims_r[VY];
  in_prims[VYL]=out_prims_l[VY];

  /* There are two stories going on here: 
   * 1) Computation of \partial_y on RHS of \partial_t {mhd_st_{x,y,z}}, 
   *    via PPM reconstruction onto (i,j-1/2,k), so that 
   *    \partial_y F = [ F(i,j+1/2,k) - F(i,j-1/2,k) ] / dy
   * 2) Computation of \partial_t A_i, where A_i are *staggered* gridfunctions,
   *    where A_x is defined at (i,j+1/2,k+1/2), A_y at (i+1/2,j,k+1/2), etc.
   *    Ai_rhs = \partial_t A_i = \epsilon_{ijk} \psi^{6} (v^j B^k - v^j B^k),
   *    where \epsilon_{ijk} is the flat-space antisymmetric operator.
   * 2A) Az_rhs is defined at (i+1/2,j+1/2,k), and it depends on {Bx,By,vx,vy},
   *     so the trick is to reconstruct {Bx,By,vx,vy} cleverly to get to these
   *     staggered points. For example:
   * 2Aa) VXR = [right-face of vx reconstructed along x-direction above] is at (i-1/2,j,k),
   *      and we reconstruct it to (i-1/2,j-1/2,k) below. Similarly for {VXL,VYR,VYL}
   * 2Ab) Bx_stagger is at (i+1/2,j,k), and we reconstruct to (i+1/2,j-1/2,k) below
   * 2Ac) By_stagger is at (i-1/2,j+1/2,k) already for Az_rhs, from the previous step.
   * 2B) Ax_rhs is defined at (i,j+1/2,k+1/2), and it depends on {By,Bz,vy,vz}.
   *     Again the trick is to reconstruct these onto these staggered points.
   * 2Ba) Bz_stagger is at (i,j,k+1/2), and we reconstruct to (i,j-1/2,k+1/2) below */
  ww=0;
  // NOTE! The order of variable reconstruction is important here,
  //   as we don't want to overwrite {vxr,vxl,vyr,vyl}!
  which_prims_to_reconstruct[ww]=VXR;       ww++;
  which_prims_to_reconstruct[ww]=VYR;       ww++;
  which_prims_to_reconstruct[ww]=VXL;       ww++;
  which_prims_to_reconstruct[ww]=VYL;       ww++;
  num_prims_to_reconstruct=ww;
  // This function is housed in the file: "reconstruct_set_of_prims_PPM_GRFFE.C"
  reconstruct_set_of_prims_PPM_GRFFE_NRPy(params, auxevol_gfs, flux_dirn+1, num_prims_to_reconstruct,                                                          
                                          which_prims_to_reconstruct, in_prims, out_prims_r, out_prims_l, temporary);
  ww=0;
  // Reconstruct other primitives last!
  which_prims_to_reconstruct[ww]=VX;        ww++;
  which_prims_to_reconstruct[ww]=VY;        ww++;
  which_prims_to_reconstruct[ww]=VZ;        ww++;
  which_prims_to_reconstruct[ww]=BX_CENTER; ww++;
  //which_prims_to_reconstruct[ww]=BY_CENTER; ww++;
  which_prims_to_reconstruct[ww]=BZ_CENTER; ww++;
  which_prims_to_reconstruct[ww]=BX_STAGGER;ww++;
  which_prims_to_reconstruct[ww]=BZ_STAGGER;ww++;
  num_prims_to_reconstruct=ww;
  // This function is housed in the file: "reconstruct_set_of_prims_PPM_GRFFE.C"
  reconstruct_set_of_prims_PPM_GRFFE_NRPy(params, auxevol_gfs, flux_dirn+1, num_prims_to_reconstruct,                                                          
                                          which_prims_to_reconstruct, in_prims, out_prims_r, out_prims_l, temporary);
  //Right and left face values of BI_CENTER are used in GRFFE__S_i__flux computation (first to compute b^a).
  //   Instead of reconstructing, we simply set B^y face values to be consistent with BY_STAGGER.
#pragma omp parallel for
  for(int k=0;k<Nxx_plus_2NGHOSTS2;k++) for(int j=0;j<Nxx_plus_2NGHOSTS1;j++) for(int i=0;i<Nxx_plus_2NGHOSTS0;i++) {
        const int index=IDX3S(i,j,k), indexjm1=IDX3S(i,j-1+(j==0),k); /* indexjm1=0 when j=0 */
        out_prims_r[BY_CENTER].gf[index]=out_prims_l[BY_CENTER].gf[index]=in_prims[BY_STAGGER].gf[indexjm1]; }
  // Then add fluxes to RHS for hydro variables {vx,vy,vz}:
  // This function is housed in the file: "add_fluxes_and_source_terms_to_hydro_rhss.C"
  calculate_StildeD1_source_term(params,auxevol_gfs,rhs_gfs);
  calculate_Stilde_flux_D1_right(params,auxevol_gfs,rhs_gfs);
  calculate_Stilde_flux_D1_left(params,auxevol_gfs,rhs_gfs);

  /*****************************************
   * COMPUTING RHS OF A_z, BOOKKEEPING NOTE:
   * We want to compute 
   * \partial_t A_z - [gauge terms] = \psi^{6} (v^x B^y - v^y B^x).
   * A_z is defined at (i+1/2,j+1/2,k).
   * ==========================
   * Where defined  | Variables
   * (i-1/2,j-1/2,k)| {vxrr,vxrl,vxlr,vxll,vyrr,vyrl,vylr,vyll}
   * (i+1/2,j-1/2,k)| {Bx_stagger_r,Bx_stagger_l} (see Table 1 in arXiv:1007.2848)
   * (i-1/2,j+1/2,k)| {By_stagger_r,By_stagger_l} (see Table 1 in arXiv:1007.2848)
   * (i,j,k)        | {phi}
   * ==========================
   ******************************************/
  // Interpolates to i+1/2
#define IPH(METRICm1,METRICp0,METRICp1,METRICp2) (-0.0625*((METRICm1) + (METRICp2)) + 0.5625*((METRICp0) + (METRICp1)))
  // Next compute phi at (i+1/2,j+1/2,k):
#pragma omp parallel for
  for(int k=0;k<Nxx_plus_2NGHOSTS2;k++) for(int j=1;j<Nxx_plus_2NGHOSTS1-2;j++) for(int i=1;i<Nxx_plus_2NGHOSTS0-2;i++) {
        temporary[IDX3S(i,j,k)]= 
          IPH(IPH(phi_bssn[IDX3S(i-1,j-1,k)],phi_bssn[IDX3S(i,j-1,k)],phi_bssn[IDX3S(i+1,j-1,k)],phi_bssn[IDX3S(i+2,j-1,k)]),
              IPH(phi_bssn[IDX3S(i-1,j  ,k)],phi_bssn[IDX3S(i,j  ,k)],phi_bssn[IDX3S(i+1,j  ,k)],phi_bssn[IDX3S(i+2,j  ,k)]),
              IPH(phi_bssn[IDX3S(i-1,j+1,k)],phi_bssn[IDX3S(i,j+1,k)],phi_bssn[IDX3S(i+1,j+1,k)],phi_bssn[IDX3S(i+2,j+1,k)]),
              IPH(phi_bssn[IDX3S(i-1,j+2,k)],phi_bssn[IDX3S(i,j+2,k)],phi_bssn[IDX3S(i+1,j+2,k)],phi_bssn[IDX3S(i+2,j+2,k)]));
      }

  int A_directionz=3;
  A_i_rhs_no_gauge_terms(A_directionz,params,out_prims_r,out_prims_l,temporary,cmax_x,cmin_x,cmax_y,cmin_y, Az_rhs);

  // in_prims[{VYR,VYL,VZR,VZL}].gz_{lo,hi} ghostzones are not correct, so we fix 
  //    this below.
  // [Note that this is a cheap operation, copying only 8 integers and a pointer.]
  in_prims[VYR]=out_prims_r[VY];
  in_prims[VYL]=out_prims_l[VY];
  in_prims[VZR]=out_prims_r[VZ];
  in_prims[VZL]=out_prims_l[VZ];

  flux_dirn=2;
  // ftilde = 0 in GRFFE, since P=rho=0.

  /* There are two stories going on here: 
   * 1) Single reconstruction to (i,j,k-1/2) for {vx,vy,vz,Bx,By,Bz} to compute
   *    z-dir'n advection terms in \partial_t {mhd_st_{x,y,z}} at (i,j,k)
   * 2) Multiple reconstructions for *staggered* gridfunctions A_i:
   *    \partial_t A_i = \epsilon_{ijk} \psi^{6} (v^j B^k - v^j B^k),
   *    where \epsilon_{ijk} is the flat-space antisymmetric operator.
   * 2A) Ax_rhs is defined at (i,j+1/2,k+1/2), depends on v{y,z} and B{y,z}
   * 2Aa) v{y,z}{r,l} are at (i,j-1/2,k), so we reconstruct here to (i,j-1/2,k-1/2)
   * 2Ab) Bz_stagger{r,l} are at (i,j-1/2,k+1/2) already.
   * 2Ac) By_stagger is at (i,j+1/2,k), and below we reconstruct its value at (i,j+1/2,k-1/2)
   * 2B) Ay_rhs is defined at (i+1/2,j,k+1/2), depends on v{z,x} and B{z,x}.
   * 2Ba) v{x,z} are reconstructed to (i,j,k-1/2). Later we'll reconstruct again to (i-1/2,j,k-1/2).
   * 2Bb) Bz_stagger is at (i,j,k+1/2). Later we will reconstruct to (i-1/2,j,k+1/2).
   * 2Bc) Bx_stagger is at (i+1/2,j,k), and below we reconstruct its value at (i+1/2,j,k-1/2)
   */
  ww=0;
  // NOTE! The order of variable reconstruction is important here,
  //   as we don't want to overwrite {vxr,vxl,vyr,vyl}!
  which_prims_to_reconstruct[ww]=VYR;       ww++;
  which_prims_to_reconstruct[ww]=VZR;       ww++;
  which_prims_to_reconstruct[ww]=VYL;       ww++;
  which_prims_to_reconstruct[ww]=VZL;       ww++;
  num_prims_to_reconstruct=ww;
  // This function is housed in the file: "reconstruct_set_of_prims_PPM_GRFFE.C"
  reconstruct_set_of_prims_PPM_GRFFE_NRPy(params, auxevol_gfs, flux_dirn+1, num_prims_to_reconstruct,                                                          
                                          which_prims_to_reconstruct, in_prims, out_prims_r, out_prims_l, temporary);
  // Reconstruct other primitives last!
  ww=0;
  which_prims_to_reconstruct[ww]=VX;        ww++;
  which_prims_to_reconstruct[ww]=VY;        ww++;
  which_prims_to_reconstruct[ww]=VZ;        ww++;
  which_prims_to_reconstruct[ww]=BX_CENTER; ww++;
  which_prims_to_reconstruct[ww]=BY_CENTER; ww++;
  //which_prims_to_reconstruct[ww]=BZ_CENTER; ww++;
  which_prims_to_reconstruct[ww]=BX_STAGGER; ww++;
  which_prims_to_reconstruct[ww]=BY_STAGGER; ww++;
  num_prims_to_reconstruct=ww;
  // This function is housed in the file: "reconstruct_set_of_prims_PPM_GRFFE.C"
  reconstruct_set_of_prims_PPM_GRFFE_NRPy(params, auxevol_gfs, flux_dirn+1, num_prims_to_reconstruct,                                                          
                                          which_prims_to_reconstruct, in_prims, out_prims_r, out_prims_l, temporary);
  //Right and left face values of BI_CENTER are used in GRFFE__S_i__flux computation (first to compute b^a).
  //   Instead of reconstructing, we simply set B^z face values to be consistent with BZ_STAGGER.
#pragma omp parallel for
  for(int k=0;k<Nxx_plus_2NGHOSTS2;k++) for(int j=0;j<Nxx_plus_2NGHOSTS1;j++) for(int i=0;i<Nxx_plus_2NGHOSTS0;i++) {
        const int index=IDX3S(i,j,k), indexkm1=IDX3S(i,j,k-1+(k==0)); /* indexkm1=0 when k=0 */
        out_prims_r[BZ_CENTER].gf[index]=out_prims_l[BZ_CENTER].gf[index]=in_prims[BZ_STAGGER].gf[indexkm1]; }
  // Then add fluxes to RHS for hydro variables {vx,vy,vz}:
  // This function is housed in the file: "add_fluxes_and_source_terms_to_hydro_rhss.C"
  calculate_StildeD2_source_term(params,auxevol_gfs,rhs_gfs);
  calculate_Stilde_flux_D2_right(params,auxevol_gfs,rhs_gfs);
  calculate_Stilde_flux_D2_left(params,auxevol_gfs,rhs_gfs);

  // in_prims[{VYR,VYL,VZR,VZL}].gz_{lo,hi} ghostzones are not set correcty.
  //    We fix this below.
  // [Note that this is a cheap operation, copying only 8 integers and a pointer.]
  in_prims[VXR]=out_prims_r[VX]; 
  in_prims[VZR]=out_prims_r[VZ]; 
  in_prims[VXL]=out_prims_l[VX];
  in_prims[VZL]=out_prims_l[VZ];
  // FIXME: lines above seem to be inconsistent with lines below.... Possible bug, not major enough to affect evolutions though.
  in_prims[VZR].gz_lo[1]=in_prims[VZR].gz_hi[1]=0;
  in_prims[VXR].gz_lo[1]=in_prims[VXR].gz_hi[1]=0;
  in_prims[VZL].gz_lo[1]=in_prims[VZL].gz_hi[1]=0;
  in_prims[VXL].gz_lo[1]=in_prims[VXL].gz_hi[1]=0;
  /*****************************************
   * COMPUTING RHS OF A_x, BOOKKEEPING NOTE:
   * We want to compute 
   * \partial_t A_x - [gauge terms] = \psi^{6} (v^y B^z - v^z B^y).
   * A_x is defined at (i,j+1/2,k+1/2).
   * ==========================
   * Where defined  | Variables
   * (i,j-1/2,k-1/2)| {vyrr,vyrl,vylr,vyll,vzrr,vzrl,vzlr,vzll}
   * (i,j+1/2,k-1/2)| {By_stagger_r,By_stagger_l} (see Table 1 in arXiv:1007.2848)
   * (i,j-1/2,k+1/2)| {Bz_stagger_r,Bz_stagger_l} (see Table 1 in arXiv:1007.2848)
   * (i,j,k)        | {phi}
   * ==========================
   ******************************************/
  // Next compute phi at (i,j+1/2,k+1/2):
#pragma omp parallel for
  for(int k=1;k<Nxx_plus_2NGHOSTS2-2;k++) for(int j=1;j<Nxx_plus_2NGHOSTS1-2;j++) for(int i=0;i<Nxx_plus_2NGHOSTS0;i++) {
        temporary[IDX3S(i,j,k)]= 
          IPH(IPH(phi_bssn[IDX3S(i,j-1,k-1)],phi_bssn[IDX3S(i,j,k-1)],phi_bssn[IDX3S(i,j+1,k-1)],phi_bssn[IDX3S(i,j+2,k-1)]),
              IPH(phi_bssn[IDX3S(i,j-1,k  )],phi_bssn[IDX3S(i,j,k  )],phi_bssn[IDX3S(i,j+1,k  )],phi_bssn[IDX3S(i,j+2,k  )]),
              IPH(phi_bssn[IDX3S(i,j-1,k+1)],phi_bssn[IDX3S(i,j,k+1)],phi_bssn[IDX3S(i,j+1,k+1)],phi_bssn[IDX3S(i,j+2,k+1)]),
              IPH(phi_bssn[IDX3S(i,j-1,k+2)],phi_bssn[IDX3S(i,j,k+2)],phi_bssn[IDX3S(i,j+1,k+2)],phi_bssn[IDX3S(i,j+2,k+2)]));
      }

  int A_directionx=1;
  A_i_rhs_no_gauge_terms(A_directionx,params,out_prims_r,out_prims_l,temporary,cmax_y,cmin_y,cmax_z,cmin_z, Ax_rhs);

  // We reprise flux_dirn=1 to finish up computations of Ai_rhs's!
  flux_dirn=0;
  // ftilde = 0 in GRFFE, since P=rho=0.

  ww=0;
  // NOTE! The order of variable reconstruction is important here,
  //   as we don't want to overwrite {vxr,vxl,vyr,vyl}!
  which_prims_to_reconstruct[ww]=VXR;       ww++;
  which_prims_to_reconstruct[ww]=VZR;       ww++;
  which_prims_to_reconstruct[ww]=VXL;       ww++;
  which_prims_to_reconstruct[ww]=VZL;       ww++;
  which_prims_to_reconstruct[ww]=BZ_STAGGER;ww++;
  num_prims_to_reconstruct=ww;
  // This function is housed in the file: "reconstruct_set_of_prims_PPM_GRFFE.C"
  reconstruct_set_of_prims_PPM_GRFFE_NRPy(params, auxevol_gfs, flux_dirn+1, num_prims_to_reconstruct,                                                          
                                          which_prims_to_reconstruct, in_prims, out_prims_r, out_prims_l, temporary);

  /*****************************************
   * COMPUTING RHS OF A_y, BOOKKEEPING NOTE:
   * We want to compute 
   * \partial_t A_y - [gauge terms] = \psi^{6} (v^z B^x - v^x B^z).
   * A_y is defined at (i+1/2,j,k+1/2).
   * ==========================
   * Where defined  | Variables
   * (i-1/2,j,k-1/2)| {vyrr,vyrl,vylr,vyll,vzrr,vzrl,vzlr,vzll}
   * (i+1/2,j,k-1/2)| {By_stagger_r,By_stagger_l} (see Table 1 in arXiv:1007.2848)
   * (i-1/2,j,k+1/2)| {Bz_stagger_r,Bz_stagger_l} (see Table 1 in arXiv:1007.2848)
   * (i,j,k)        | {phi}
   * ==========================
   ******************************************/
  // Next compute phi at (i+1/2,j,k+1/2):
#pragma omp parallel for
  for(int k=1;k<Nxx_plus_2NGHOSTS2-2;k++) for(int j=0;j<Nxx_plus_2NGHOSTS1;j++) for(int i=1;i<Nxx_plus_2NGHOSTS0-2;i++) {
        temporary[IDX3S(i,j,k)]= 
          IPH(IPH(phi_bssn[IDX3S(i-1,j,k-1)],phi_bssn[IDX3S(i,j,k-1)],phi_bssn[IDX3S(i+1,j,k-1)],phi_bssn[IDX3S(i+2,j,k-1)]),
              IPH(phi_bssn[IDX3S(i-1,j,k  )],phi_bssn[IDX3S(i,j,k  )],phi_bssn[IDX3S(i+1,j,k  )],phi_bssn[IDX3S(i+2,j,k  )]),
              IPH(phi_bssn[IDX3S(i-1,j,k+1)],phi_bssn[IDX3S(i,j,k+1)],phi_bssn[IDX3S(i+1,j,k+1)],phi_bssn[IDX3S(i+2,j,k+1)]),
              IPH(phi_bssn[IDX3S(i-1,j,k+2)],phi_bssn[IDX3S(i,j,k+2)],phi_bssn[IDX3S(i+1,j,k+2)],phi_bssn[IDX3S(i+2,j,k+2)]));
      }
  
  int A_directiony=2;
  A_i_rhs_no_gauge_terms(A_directiony,params,out_prims_r,out_prims_l,temporary,cmax_z,cmin_z,cmax_x,cmin_x, Ay_rhs);

  // Next compute psi6phi_rhs, and add gauge terms to A_i_rhs terms!
  //   Note that in the following function, we don't bother with reconstruction, instead interpolating.
  // We need A^i, but only have A_i. So we add gtupij to the list of input variables.
  REAL *interp_vars[MAXNUMINTERP];
  ww=0;
  interp_vars[ww]=auxevol_gfs+Nxxp2NG012*BETAU0GF;   ww++;
  interp_vars[ww]=auxevol_gfs+Nxxp2NG012*BETAU1GF;   ww++;
  interp_vars[ww]=auxevol_gfs+Nxxp2NG012*BETAU2GF;   ww++;
  interp_vars[ww]=auxevol_gfs+Nxxp2NG012*GAMMADD00GF;  ww++;
  interp_vars[ww]=auxevol_gfs+Nxxp2NG012*GAMMADD01GF;  ww++;
  interp_vars[ww]=auxevol_gfs+Nxxp2NG012*GAMMADD02GF;  ww++;
  interp_vars[ww]=auxevol_gfs+Nxxp2NG012*GAMMADD11GF;  ww++;
  interp_vars[ww]=auxevol_gfs+Nxxp2NG012*GAMMADD12GF;  ww++;
  interp_vars[ww]=auxevol_gfs+Nxxp2NG012*GAMMADD22GF;  ww++;
  interp_vars[ww]=temporary;ww++;
  interp_vars[ww]=auxevol_gfs_Nxxp2NG012*ALPHAGF;   ww++;
  interp_vars[ww]=evol_gfs+Nxxp2NG012*AD0GF;      ww++;
  interp_vars[ww]=evol_gfs+Nxxp2NG012*AD1GF;      ww++;
  interp_vars[ww]=evol_gfs+Nxxp2NG012*AD2GF;      ww++;
  const int max_num_interp_variables=ww;
//   if(max_num_interp_variables>MAXNUMINTERP) {CCTK_VError(VERR_DEF_PARAMS,"Error: Didn't allocate enough space for interp_vars[]."); }
  // We are FINISHED with v{x,y,z}{r,l} and P{r,l} so we use these 8 gridfunctions' worth of space as temp storage.
  Lorenz_psi6phi_rhs__add_gauge_terms_to_A_i_rhs(params,interp_vars,psi6phi, 
                                                 auxevol_gfs+Nxxp2NG012*VALENCIAV_R0GF,  // WARNING: 
                                                 auxevol_gfs+Nxxp2NG012*VALENCIAV_R1GF,  // ALL VARIABLES
                                                 auxevol_gfs+Nxxp2NG012*VALENCIAV_R2GF,  // ON THESE LINES
                                                 auxevol_gfs+Nxxp2NG012*VALENCIAV_L0GF,  // ARE OVERWRITTEN
                                                 auxevol_gfs+Nxxp2NG012*VALENCIAV_L1GF,  // FOR TEMP STORAGE
                                                 auxevol_gfs+Nxxp2NG012*VALENCIAV_L2GF,  // .
                                                 auxevol_gfs+Nxxp2NG012*VALENCIAV_RR0GF, // .
                                                 auxevol_gfs+Nxxp2NG012*VALENCIAV_RL0GF, // .
                                                 rhs_gfs+Nxxp2NG012*PSI6PHIGF,
                                                 rhs_gfs+Nxxp2NG012*AD0GF,
                                                 rhs_gfs+Nxxp2NG012*AD1GF,
                                                 rhs_gfs+Nxxp2NG012*AD2GF);
#pragma omp parallel for
  for(int k=0;k<Nxx_plus_2NGHOSTS2;k++) for(int j=0;j<Nxx_plus_2NGHOSTS1;j++) for(int i=0;i<Nxx_plus_2NGHOSTS0;i++) {
        const int index=IDX3S(i,j,k);
        if(r[index]<min_radius_inside_of_which_conserv_to_prims_FFE_and_FFE_evolution_is_DISABLED) {
          st_x_rhs[index]=0.0;
          st_y_rhs[index]=0.0;
          st_z_rhs[index]=0.0;

          psi6phi_rhs[index] = 0.0;
          Ax_rhs[index] = 0.0;
          Ay_rhs[index] = 0.0;
          Az_rhs[index] = 0.0;
        }
      }
}

void GiRaFFE_NRPy_post_step(const paramstruct *restrict params,REAL *xx[3],REAL *restrict auxevol_gfs,REAL *restrict evol_gfs,const int n) {
    // First, apply BCs to AD and psi6Phi. Then calculate BU from AD
    apply_bcs_potential(params,evol_gfs);
GiRaFFE_compute_B_and_Bstagger_from_A(const paramstruct *params,
                                      auxevol_gfs+Nxxp2NG012*GAMMADD00GF,
                                      auxevol_gfs+Nxxp2NG012*GAMMADD01GF,
                                      auxevol_gfs+Nxxp2NG012*GAMMADD02GF,
                                      auxevol_gfs+Nxxp2NG012*GAMMADD11GF,
                                      auxevol_gfs+Nxxp2NG012*GAMMADD12GF,
                                      auxevol_gfs+Nxxp2NG012*GAMMADD22GF,
                                      auxevol_gfs + Nxxp2NG012*AEVOLPARENGF, /* Temporary storage */ 
                                      evol_gfs+Nxxp2NG012*AD0GF,
                                      evol_gfs+Nxxp2NG012*AD1GF,
                                      evol_gfs+Nxxp2NG012*AD2GF,
                                      auxevol_gfs+Nxxp2NG012*BU0GF,
                                      auxevol_gfs+Nxxp2NG012*BU1GF,
                                      auxevol_gfs+Nxxp2NG012*BU2GF,
                                      auxevol_gfs+Nxxp2NG012*BSTAGGERU0GF,
                                      auxevol_gfs+Nxxp2NG012*BSTAGGERU1GF,
                                      auxevol_gfs+Nxxp2NG012*BSTAGGERU2GF);
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
