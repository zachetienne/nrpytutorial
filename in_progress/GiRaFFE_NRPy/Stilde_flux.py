from outputC import outCfunction, lhrh # NRPy+: Core C code output module
import sympy as sp               # SymPy: The Python computer algebra package upon which NRPy+ depends
import NRPy_param_funcs as par   # NRPy+: Parameter interface
import finite_difference as fin  # NRPy+: Finite difference C code generation module
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
        if mom_comp==0:
            global F_out,U_out,smallb_out
            F_out=F
            U_out=U
            smallb_out = GRFFE.smallbsquared
        calculate_GRFFE_Tmunu_and_contractions(flux_dirn, mom_comp, gamma_faceDD,beta_faceU,alpha_face,\
                                               Valenciav_lU,B_lU,sqrt4pi)
        Fl = F
        Ul = U
        Stilde_fluxD[mom_comp] = HLLE_solver(chsp.cmax, chsp.cmin, Fr, Fl, Ur, Ul)

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

        # We'll also need to store the results of the HLLE step between functions.
        ixp.register_gridfunctions_for_single_rank1("AUXEVOL","Stilde_flux_HLLED")

    if write_cmax_cmin:
        # In the staggered case, we will also want to output cmax and cmin
        # If we want to write cmax and cmin, we will need to be able to change auxevol_gfs:
        input_params_for_Stilde_flux = "const paramstruct *params,REAL *auxevol_gfs,REAL *rhs_gfs"
    else:
        input_params_for_Stilde_flux = "const paramstruct *params,const REAL *auxevol_gfs,REAL *rhs_gfs"

    if write_cmax_cmin:
        name_suffixes = ["_x","_y","_z"]

    for flux_dirn in range(3):
        calculate_Stilde_flux(flux_dirn,alpha_face,gamma_faceDD,beta_faceU,\
                              Valenciav_rU,B_rU,Valenciav_lU,B_lU,sqrt4pi)

        Stilde_flux_to_print = [
                                lhrh(lhs=gri.gfaccess("out_gfs","Stilde_flux_HLLED0"),rhs=Stilde_fluxD[0]),
                                lhrh(lhs=gri.gfaccess("out_gfs","Stilde_flux_HLLED1"),rhs=Stilde_fluxD[1]),
                                lhrh(lhs=gri.gfaccess("out_gfs","Stilde_flux_HLLED2"),rhs=Stilde_fluxD[2])
                               ]

        if write_cmax_cmin:
            Stilde_flux_to_print = Stilde_flux_to_print \
                                  +[
                                    lhrh(lhs=gri.gfaccess("out_gfs","cmax"+name_suffixes[flux_dirn]),rhs=chsp.cmax),
                                    lhrh(lhs=gri.gfaccess("out_gfs","cmin"+name_suffixes[flux_dirn]),rhs=chsp.cmin)
                                   ]

        desc = "Compute the flux term of all 3 components of tilde{S}_i on the left face in the " + str(flux_dirn) + "direction for all components."
        name = "calculate_Stilde_flux_D" + str(flux_dirn)
        Ccode_function = outCfunction(
            outfile  = "returnstring", desc=desc, name=name,
            params   = input_params_for_Stilde_flux,
            body     = fin.FD_outputC("returnstring",Stilde_flux_to_print,params=outCparams).replace("IDX4","IDX4S"),
            loopopts ="InteriorPoints",
            rel_path_to_Cparams=os.path.join("../")).replace("NGHOSTS+Nxx0","NGHOSTS+Nxx0+1").replace("NGHOSTS+Nxx1","NGHOSTS+Nxx1+1").replace("NGHOSTS+Nxx2","NGHOSTS+Nxx2+1")

        with open(os.path.join(out_dir,name+".h"),"w") as file:
            file.write(Ccode_function)

    pre_body = """// Notice in the loop below that we go from 3 to cctk_lsh-3 for i, j, AND k, even though
    //   we are only computing the flux in one direction. This is because in the end,
    //   we only need the rhs's from 3 to cctk_lsh-3 for i, j, and k.
    const REAL invdxi[4] = {1e100,invdx0,invdx1,invdx2};
    const REAL invdx = invdxi[flux_dirn];"""

    FD_body = """const int index = IDX3S(i0,i1,i2);
const int indexp1 = IDX3S(i0+kronecker_delta[flux_dirn][0],i1+kronecker_delta[flux_dirn][1],i2+kronecker_delta[flux_dirn][2]);

rhs_gfs[IDX4ptS(STILDED0GF,index)] += (auxevol_gfs[IDX4ptS(STILDE_FLUX_HLLED0GF,index)]     - auxevol_gfs[IDX4ptS(STILDE_FLUX_HLLED0GF,indexp1)]    ) * invdx;
rhs_gfs[IDX4ptS(STILDED1GF,index)] += (auxevol_gfs[IDX4ptS(STILDE_FLUX_HLLED1GF,index)]     - auxevol_gfs[IDX4ptS(STILDE_FLUX_HLLED1GF,indexp1)]    ) * invdx;
rhs_gfs[IDX4ptS(STILDED2GF,index)] += (auxevol_gfs[IDX4ptS(STILDE_FLUX_HLLED2GF,index)]     - auxevol_gfs[IDX4ptS(STILDE_FLUX_HLLED2GF,indexp1)]    ) * invdx;"""

    desc = "Compute the difference in the flux of StildeD on the opposite faces in flux_dirn for all components."
    name = "calculate_Stilde_rhsD"
    outCfunction(
        outfile  = os.path.join(out_dir,name+".h"), desc=desc, name=name,
        params   = "const int flux_dirn,const paramstruct *params,const REAL *auxevol_gfs,REAL *rhs_gfs",
        preloop  = pre_body,
        body     = FD_body,
        loopopts = "InteriorPoints",
        rel_path_to_Cparams=os.path.join("../")
    )
