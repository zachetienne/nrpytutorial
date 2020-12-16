# Step 1: The StildeD RHS *source* term
import os
from outputC import outputC, outCfunction # NRPy+: Core C code output module
import indexedexp as ixp         # NRPy+: Symbolic indexed expression (e.g., tensors, vectors, etc.) support
import GRHD.equations as GRHD    # NRPy+: Generate general relativistic hydrodynamics equations
import GRFFE.equations as GRFFE  # NRPy+: Generate general relativistic force-free electrodynamics equations

thismodule = __name__

def generate_memory_access_code(gammaDD,betaU,alpha):
    # There are several pieces of C code that we will write ourselves because we need to do things
    # a little bit outside of what NRPy+ is built for.
    # First, we will write general memory access. We will read in values from memory at a given point
    # for each quantity we care about.
    global general_access
    general_access = ""
    for var in ["GAMMADD00", "GAMMADD01", "GAMMADD02",
                "GAMMADD11", "GAMMADD12", "GAMMADD22",
                "BETAU0", "BETAU1", "BETAU2","ALPHA",
                "BU0","BU1","BU2",
                "VALENCIAVU0","VALENCIAVU1","VALENCIAVU2"]:
        lhsvar = var.lower().replace("dd","DD").replace("u","U").replace("bU","BU").replace("valencia","Valencia")
        # e.g.,
        # const REAL gammaDD00dD0 = auxevol_gfs[IDX4S(GAMMA_FACEDD00GF,i0,i1,i2)];
        general_access += "const REAL "+lhsvar+" = auxevol_gfs[IDX4S("+var+"GF,i0,i1,i2)];\n"

    # This quick function returns a nearby point for memory access. We need this because derivatives are not local operations.
    def idxp1(dirn):
        if dirn==0:
            return "i0+1,i1,i2"
        if dirn==1:
            return "i0,i1+1,i2"
        if dirn==2:
            return "i0,i1,i2+1"

    # Next we evaluate needed derivatives of the metric, based on their values at cell faces
    global metric_deriv_access
    metric_deriv_access = []
#     for dirn in range(3):
#         metric_deriv_access.append("")
#         for var in ["GAMMA_FACEDDdD00", "GAMMA_FACEDDdD01", "GAMMA_FACEDDdD02",
#                     "GAMMA_FACEDDdD11", "GAMMA_FACEDDdD12", "GAMMA_FACEDDdD22",
#                     "BETA_FACEUdD0", "BETA_FACEUdD1", "BETA_FACEUdD2","ALPHA_FACEdD"]:
#             lhsvar = var.lower().replace("dddd","DDdD").replace("udd","UdD").replace("dd","dD").replace("u","U").replace("_face","")
#             rhsvar = var.replace("dD","")
#             # e.g.,
#             # const REAL gammaDDdD000 = (auxevol_gfs[IDX4S(GAMMA_FACEDD00GF,i0+1,i1,i2)]-auxevol_gfs[IDX4S(GAMMA_FACEDD00GF,i0,i1,i2)])/dxx0;
#             metric_deriv_access[dirn] += "const REAL "+lhsvar+str(dirn)+" = (auxevol_gfs[IDX4S("+rhsvar+"GF,"+idxp1(dirn)+")]-auxevol_gfs[IDX4S("+rhsvar+"GF,i0,i1,i2)])/dxx"+str(dirn)+";\n"
#         metric_deriv_access[dirn] += "REAL Stilde_rhsD"+str(dirn)+";\n"
    # For this workaround, instead of taking the derivative of the metric components and then building the
    # four-metric, we build the four-metric and then take derivatives. Do this at i and i+1
    for dirn in range(3):
        metric_deriv_access.append("")
        for var in ["GAMMA_FACEDD00", "GAMMA_FACEDD01", "GAMMA_FACEDD02",
                    "GAMMA_FACEDD11", "GAMMA_FACEDD12", "GAMMA_FACEDD22",
                    "BETA_FACEU0", "BETA_FACEU1", "BETA_FACEU2","ALPHA_FACE"]:
            lhsvar = var.lower().replace("dd","DD").replace("u","U")
            rhsvar = var
            # e.g.,
            # const REAL gammaDD00 = auxevol_gfs[IDX4S(GAMMA_FACEDD00GF,i0,i1,i2)];
            metric_deriv_access[dirn] += "const REAL "+lhsvar+" = auxevol_gfs[IDX4S("+rhsvar+"GF,i0,i1,i2)];\n"
        # Read in at the next grid point
        for var in ["GAMMA_FACEDD00", "GAMMA_FACEDD01", "GAMMA_FACEDD02",
                    "GAMMA_FACEDD11", "GAMMA_FACEDD12", "GAMMA_FACEDD22",
                    "BETA_FACEU0", "BETA_FACEU1", "BETA_FACEU2","ALPHA_FACE"]:
            lhsvar = var.lower().replace("dd","DD").replace("u","U").replace("_face","_facep1")
            rhsvar = var
            # e.g.,
            # const REAL gammaDD00 = auxevol_gfs[IDX4S(GAMMA_FACEDD00GF,i0+1,i1,i2)];
            metric_deriv_access[dirn] += "const REAL "+lhsvar+" = auxevol_gfs[IDX4S("+rhsvar+"GF,"+idxp1(dirn)+")];\n"
        metric_deriv_access[dirn] += "REAL Stilde_rhsD"+str(dirn)+";\n"
    import BSSN.ADMBSSN_tofrom_4metric as AB4m
    AB4m.g4DD_ito_BSSN_or_ADM("ADM",gammaDD,betaU,alpha)
    four_metric_vars = [
                        AB4m.g4DD[0][0],
                        AB4m.g4DD[0][1],
                        AB4m.g4DD[0][2],
                        AB4m.g4DD[0][3],
                        AB4m.g4DD[1][1],
                        AB4m.g4DD[1][2],
                        AB4m.g4DD[1][3],
                        AB4m.g4DD[2][2],
                        AB4m.g4DD[2][3],
                        AB4m.g4DD[3][3]
                       ]
    four_metric_names = [
                         "g4DD00",
                         "g4DD01",
                         "g4DD02",
                         "g4DD03",
                         "g4DD11",
                         "g4DD12",
                         "g4DD13",
                         "g4DD22",
                         "g4DD23",
                         "g4DD33"
                        ]
    global four_metric_C, four_metric_Cp1
    four_metric_C = outputC(four_metric_vars,four_metric_names,"returnstring",params="outCverbose=False,CSE_sorting=none")
    for ii in range(len(four_metric_names)):
        four_metric_names[ii] += "p1"
    four_metric_Cp1 = outputC(four_metric_vars,four_metric_names,"returnstring",params="outCverbose=False,CSE_sorting=none")
    four_metric_C = four_metric_C.replace("gamma","gamma_face").replace("beta","beta_face").replace("alpha","alpha_face").replace("{","").replace("}","").replace("g4","const REAL g4").replace("tmp_","tmp_deriv")
    four_metric_Cp1 = four_metric_Cp1.replace("gamma","gamma_facep1").replace("beta","beta_facep1").replace("alpha","alpha_facep1").replace("{","").replace("}","").replace("g4","const REAL g4").replace("tmp_","tmp_derivp")

    global four_metric_deriv
    four_metric_deriv = []
    for dirn in range(3):
        four_metric_deriv.append("")
        for var in ["g4DDdD00", "g4DDdD01", "g4DDdD02", "g4DDdD03", "g4DDdD11",
                    "g4DDdD12", "g4DDdD13", "g4DDdD22", "g4DDdD23", "g4DDdD33"]:
            lhsvar   = var + str(dirn+1)
            rhsvar   = var.replace("dD","")
            rhsvarp1 = rhsvar + "p1"
            # e.g.,
            # const REAL g44DDdD000 = (g4DD00p1 - g4DD00)/dxx0;
            four_metric_deriv[dirn] += "const REAL "+lhsvar+" = ("+rhsvarp1+" - "+rhsvar+")/dxx"+str(dirn)+";\n"

    # This creates the C code that writes to the Stilde_rhs direction specified.
    global write_final_quantity
    write_final_quantity = []
    for dirn in range(3):
        write_final_quantity.append("")
        write_final_quantity[dirn] += "rhs_gfs[IDX4S(STILDED"+str(dirn)+"GF,i0,i1,i2)] += Stilde_rhsD"+str(dirn)+";"

def write_out_functions_for_StildeD_source_term(outdir,outCparams,gammaDD,betaU,alpha,ValenciavU,BU,sqrt4pi):
    generate_memory_access_code(gammaDD,betaU,alpha)
    # First, we declare some dummy tensors that we will use for the codegen.
    gammaDDdD  = ixp.declarerank3("gammaDDdD","sym01",DIM=3)
    betaUdD = ixp.declarerank2("betaUdD","nosym",DIM=3)
    alphadD = ixp.declarerank1("alphadD",DIM=3)
    g4DDdD = ixp.declarerank3("g4DDdD","sym01",DIM=4)

    # We need to rerun a few of these functions with the reset lists to make sure these functions
    # don't cheat by using analytic expressions
    GRHD.compute_sqrtgammaDET(gammaDD)
    GRHD.u4U_in_terms_of_ValenciavU__rescale_ValenciavU_by_applying_speed_limit(alpha, betaU, gammaDD, ValenciavU)
    GRFFE.compute_smallb4U(gammaDD, betaU, alpha, GRHD.u4U_ito_ValenciavU, BU, sqrt4pi)
    GRFFE.compute_smallbsquared(gammaDD, betaU, alpha, GRFFE.smallb4U)
    GRFFE.compute_TEM4UU(gammaDD,betaU,alpha, GRFFE.smallb4U, GRFFE.smallbsquared,GRHD.u4U_ito_ValenciavU)
#     GRHD.compute_g4DD_zerotimederiv_dD(gammaDD,betaU,alpha, gammaDDdD,betaUdD,alphadD)
    GRHD.compute_S_tilde_source_termD(alpha, GRHD.sqrtgammaDET,g4DDdD, GRFFE.TEM4UU)
    for i in range(3):
        desc = "Adds the source term to StildeD"+str(i)+"."
        name = "calculate_StildeD"+str(i)+"_source_term"
        outCfunction(
            outfile  = os.path.join(outdir,name+".h"), desc=desc, name=name,
            params   ="const paramstruct *params,const REAL *auxevol_gfs, REAL *rhs_gfs",
            body     = general_access \
                      +metric_deriv_access[i]\
                      +four_metric_C\
                      +four_metric_Cp1\
                      +four_metric_deriv[i]\
                      +outputC(GRHD.S_tilde_source_termD[i],"Stilde_rhsD"+str(i),"returnstring",params=outCparams).replace("IDX4","IDX4S")\
                      +write_final_quantity[i],
            loopopts ="InteriorPoints",
            rel_path_for_Cparams=os.path.join("../"))
