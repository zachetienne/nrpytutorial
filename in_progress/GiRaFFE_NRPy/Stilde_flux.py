from outputC import outCfunction, outputC # NRPy+: Core C code output module
import sympy as sp               # SymPy: The Python computer algebra package upon which NRPy+ depends
import NRPy_param_funcs as par   # NRPy+: Parameter interface
import grid as gri               # NRPy+: Functions having to do with numerical grids
import indexedexp as ixp         # NRPy+: Symbolic indexed expression (e.g., tensors, vectors, etc.) support
import os, sys           # Standard Python modules for multiplatform OS-level functions
import GiRaFFE_NRPy.GiRaFFE_NRPy_Characteristic_Speeds as chsp # GRFFE: the characteristic speeds

par.set_parval_from_str("grid::DIM", 3)
DIM = par.parval_from_str("grid::DIM")

thismodule = __name__

# We'll rewrite this assuming that we've passed the entire reconstructed
# gridfunctions. You could also do this with only one point, but then you'd
# need to declare everything as a Cparam in NRPy+
nrpy_dir_path = os.path.join("..")
if nrpy_dir_path not in sys.path:
    sys.path.append(nrpy_dir_path)

# We'll rewrite this assuming that we've passed the entire reconstructed
# gridfunctions. You could also do this with only one point, but then you'd
# need to declare everything as a Cparam in NRPy+

import GRHD.equations as GRHD
import GRFFE.equations as GRFFE

def calculate_GRFFE_Tmunu_and_contractions(flux_dirn, mom_comp, gammaDD,betaU,alpha,ValenciavU,BU,sqrt4pi):
    GRHD.compute_sqrtgammaDET(gammaDD)

    GRHD.u4U_in_terms_of_ValenciavU__rescale_ValenciavU_by_applying_speed_limit(alpha, betaU, gammaDD, ValenciavU)
    GRFFE.compute_smallb4U(gammaDD, betaU, alpha, GRHD.u4U_ito_ValenciavU, BU, sqrt4pi)
    GRFFE.compute_smallbsquared(gammaDD, betaU, alpha, GRFFE.smallb4U)

    GRFFE.compute_TEM4UU(gammaDD, betaU, alpha, GRFFE.smallb4U, GRFFE.smallbsquared, GRHD.u4U_ito_ValenciavU)
    GRFFE.compute_TEM4UD(gammaDD, betaU, alpha, GRFFE.TEM4UU)

    # Compute conservative variables in terms of primitive variables
    GRHD.compute_S_tildeD(alpha, GRHD.sqrtgammaDET, GRFFE.TEM4UD)

    global U,F
    # Flux F = alpha*sqrt{gamma}*T^i_j
    F = alpha*GRHD.sqrtgammaDET*GRFFE.TEM4UD[flux_dirn+1][mom_comp+1]
    # U = alpha*sqrt{gamma}*T^0_j = Stilde_j
    U = GRHD.S_tildeD[mom_comp]


def HLLE_solver(cmax, cmin, Fr, Fl, Ur, Ul):
    # This solves the Riemann problem for the mom_comp component of the momentum
    # flux StildeD in the flux_dirn direction.

    # st_j_flux = (c_\min f_R + c_\max f_L - c_\min c_\max ( st_j_r - st_j_l )) / (c_\min + c_\max)
    return (cmin*Fr + cmax*Fl - cmin*cmax*(Ur-Ul) )/(cmax + cmin)

def calculate_Stilde_flux(flux_dirn,alpha_face,gamma_faceDD,beta_faceU,\
                          Valenciav_rU,B_rU,Valenciav_lU,B_lU,sqrt4pi):

    chsp.find_cmax_cmin(flux_dirn,gamma_faceDD,beta_faceU,alpha_face)

    global Stilde_fluxD
    Stilde_fluxD = ixp.zerorank3()
    for mom_comp in range(3):
        calculate_GRFFE_Tmunu_and_contractions(flux_dirn, mom_comp, gamma_faceDD,beta_faceU,alpha_face,\
                                               Valenciav_rU,B_rU,sqrt4pi)
        Fr = F
        Ur = U
        calculate_GRFFE_Tmunu_and_contractions(flux_dirn, mom_comp, gamma_faceDD,beta_faceU,alpha_face,\
                                               Valenciav_lU,B_lU,sqrt4pi)
        Fl = F
        Ul = U
        Stilde_fluxD[mom_comp] = HLLE_solver(chsp.cmax, chsp.cmin, Fr, Fl, Ur, Ul)

def write_C_Code_to_read_memory(write_cmax_cmin=False):
    Memory_Read = ""
    # Read in all the inputs:
    for var in ["GAMMA_FACEDD00", "GAMMA_FACEDD01", "GAMMA_FACEDD02",
                "GAMMA_FACEDD11", "GAMMA_FACEDD12", "GAMMA_FACEDD22",
                "BETA_FACEU0", "BETA_FACEU1", "BETA_FACEU2","ALPHA_FACE",
                "B_RU0","B_RU1","B_RU2","B_LU0","B_LU1","B_LU2",
                "VALENCIAV_RU0","VALENCIAV_RU1","VALENCIAV_RU2",
                "VALENCIAV_LU0","VALENCIAV_LU1","VALENCIAV_LU2"]:
        lhsvar = var.lower().replace("dd","DD").replace("u","U").replace("b_","B_").replace("valencia","Valencia")
        # e.g.,
        # const REAL gammaDD00dD0 = auxevol_gfs[IDX4S(GAMMA_FACEDD00GF,i0,i1,i2)];
        Memory_Read += "const REAL "+lhsvar+" = auxevol_gfs[IDX4S("+var+"GF,i0,i1,i2)];\n"
    # Storage for the outputs:
    for var in ["Stilde_fluxD0","Stilde_fluxD1","Stilde_fluxD2"]:
        # e.g.,
        # REAL Stilde_fluxD0 = 0;
        Memory_Read += "REAL "+var+" = 0;\n"
    if write_cmax_cmin:
        # In the staggered case, we will also want to output cmax and cmin:
        for var in ["cmax","cmin"]:
            Memory_Read += "REAL "+var+" = 0;\n"
    return Memory_Read

# This quick function returns a nearby point for memory access. We need this because derivatives are not local operations.
def idxm1(dirn):
    if dirn==0:
        return "i0-1, i1, i2"
    if dirn==1:
        return "i0, i1-1, i2"
    if dirn==2:
        return "i0, i1, i2-1"

def write_C_code_to_write_results(write_cmax_cmin=False):
    # Write the outputs:
    Memory_Write = []
    for dirn in range(3):
        for comp in range(3):
            Memory_Write.append("")
            # e.g.,
            # rhs_gfs[IDX4S(STILDED0GF, i0, i1, i2)] += invdx0*Stilde_fluxD0;
            # (Note that the invdx0 gets substituted separately! It shouldn't necessarily match Stilde!)
            Memory_Write[dirn] += "rhs_gfs[IDX4S(STILDED"+str(comp)+"GF, i0, i1, i2)] += invdx"+str(dirn)+"*Stilde_fluxD"+str(comp)+";\n"
            # e.g.,
            # rhs_gfs[IDX4S(STILDED0GF, i0-1, i1, i2)] -= invdx0*Stilde_fluxD0;
            Memory_Write[dirn] += "rhs_gfs[IDX4S(STILDED"+str(comp)+"GF, "+idxm1(dirn)+")] -= invdx"+str(dirn)+"*Stilde_fluxD"+str(comp)+";\n"
    if write_cmax_cmin:
        name_suffixes = ["_X","_Y","_Z"]
        for dirn in range(3):
            for var in ["cmax","cmin"]:
                Memory_Write[dirn] += "auxevol_gfs[IDX4S("+var.upper()+name_suffixes[dirn]+"GF, i0, i1, i2)] = "+var+";\n"
    return Memory_Write

def generate_C_code_for_Stilde_flux(out_dir,inputs_provided = False, alpha_face=None, gamma_faceDD=None, beta_faceU=None,
                                    Valenciav_rU=None, B_rU=None, Valenciav_lU=None, B_lU=None, sqrt4pi=None,
                                    outCparams = "outCverbose=False,CSE_sorting=none", write_cmax_cmin=False):
    if not inputs_provided:
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
        sqrt4pi = par.Cparameters("REAL",thismodule,"sqrt4pi","sqrt(4.0*M_PI)")

    Memory_Read = write_C_Code_to_read_memory(write_cmax_cmin)
    Memory_Write = write_C_code_to_write_results(write_cmax_cmin)
    if write_cmax_cmin:
        # In the staggered case, we will also want to output cmax and cmin
        # If we want to write cmax and cmin, we will need to be able to change auxevol_gfs:
        input_params_for_Stilde_flux = "const paramstruct *params,REAL *auxevol_gfs,REAL *rhs_gfs"
    else:
        input_params_for_Stilde_flux = "const paramstruct *params,const REAL *auxevol_gfs,REAL *rhs_gfs"

    for flux_dirn in range(3):
        calculate_Stilde_flux(flux_dirn,alpha_face,gamma_faceDD,beta_faceU,\
                              Valenciav_rU,B_rU,Valenciav_lU,B_lU,sqrt4pi)
        Stilde_flux_to_print = [
                                Stilde_fluxD[0],
                                Stilde_fluxD[1],
                                Stilde_fluxD[2]
                               ]
        Stilde_flux_names = [
                             "Stilde_fluxD0",
                             "Stilde_fluxD1",
                             "Stilde_fluxD2"
                            ]

        if write_cmax_cmin:
            Stilde_flux_to_print = Stilde_flux_to_print + [chsp.cmax,chsp.cmin]
            Stilde_flux_names = Stilde_flux_names + ["cmax","cmin"]

        desc = "Compute the flux term of all 3 components of tilde{S}_i on the right face in the " + str(flux_dirn) + "direction."
        name = "calculate_Stilde_flux_D" + str(flux_dirn)
        Ccode_function = outCfunction(
            outfile  = "returnstring", desc=desc, name=name,
            params   = input_params_for_Stilde_flux,
            body     =  Memory_Read \
                       +outputC(Stilde_flux_to_print,Stilde_flux_names,"returnstring",params=outCparams).replace("IDX4","IDX4S")\
                       +Memory_Write[flux_dirn],
            loopopts ="InteriorPoints",
            rel_path_for_Cparams=os.path.join("../")).replace("NGHOSTS+Nxx0","NGHOSTS+Nxx0+1").replace("NGHOSTS+Nxx1","NGHOSTS+Nxx1+1").replace("NGHOSTS+Nxx2","NGHOSTS+Nxx2+1")

        with open(os.path.join(out_dir,name+".h"),"w") as file:
            file.write(Ccode_function)