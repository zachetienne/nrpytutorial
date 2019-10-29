from outputC import *            # NRPy+: Core C code output module
import finite_difference as fin  # NRPy+: Finite difference C code generation module
import NRPy_param_funcs as par   # NRPy+: Parameter interface
import grid as gri               # NRPy+: Functions having to do with numerical grids
import loop as lp                # NRPy+: Generate C code loops
import indexedexp as ixp         # NRPy+: Symbolic indexed expression (e.g., tensors, vectors, etc.) support
import reference_metric as rfm   # NRPy+: Reference metric support
import cmdline_helper as cmd     # NRPy+: Multi-platform Python command-line interface
import shutil, os, sys           # Standard Python modules for multiplatform OS-level functions

par.set_parval_from_str("grid::DIM", 3)
DIM = par.parval_from_str("grid::DIM")

thismodule = __name__

# This function rewrites the original version from Tutorial-u0_smallb_Poynting-Cartesian.ipynb
# without the if statement.
def compute_u0_noif(gammaDD,alpha,ValenciavU):
    # Inputs:  Metric gammaDD, lapse alpha, Valencia 3-velocity ValenciavU
    # Outputs: u^0, speed-limited ValenciavU
    
    # R = gamma_{ij} v^i v^j
    R = par.Cparameters("#define",thismodule,"TINYDOUBLE",1e-100)
    # We initialize R to TINYDOUBLE instead of 0 to avoid a 0/0 below,
    #    in the case that ValenciavU == 0 for all Valencia 3-velocities.
    # "Those tiny *doubles* make me warm all over
    #  with a feeling that I'm gonna love you till the end of time."
    #    - Adapted from Connie Francis' "Tiny Bubbles"
    for i in range(DIM):
        for j in range(DIM):
            R += gammaDD[i][j]*ValenciavU[i]*ValenciavU[j]

    # The default value isn't terribly important here, since we'll overwrite in the main C code
    GAMMA_SPEED_LIMIT = par.Cparameters("REAL",thismodule,"GAMMA_SPEED_LIMIT",10.0) # Default value based on
                                                                                    # IllinoisGRMHD.
                                                                                    # GiRaFFE default = 2000.0
    # Rmax = 1 - 1/Gamma_{\max}^2
    Rmax = 1 - 1/(GAMMA_SPEED_LIMIT*GAMMA_SPEED_LIMIT)
    # Now, we set Rmax = min(Rmax,R):
    # If Rmax>R, then Rmax = 0.5*(Rmax+R-Rmax+R) = R
    # If R>Rmax, then Rmax = 0.5*(Rmax+R+Rmax-R) = Rmax
    # If R==TINYDOUBLE, then Rmax = 0.5*(Rmax-Rmax) = 0, since, e.g., 10 +/- 1e-100 = 10 exactly in double precision
    Rmax =  sp.Rational(1,2)*(Rmax+R-sp.Abs(Rmax-R))

    # With our rescaled Rmax, v^i = sqrt{Rmax/R} v^i
    global rescaledValenciavU,rescaledu0

    rescaledValenciavU = ixp.zerorank1()
    for i in range(DIM):
        # If R == TINYDOUBLE, then Rmax=0,R=1e-100, and we get Rmax/R=0.
        #   If your velocities are of order 1e-100 and this is physically
        #   meaningful, there must be something wrong with your unit conversion.
        rescaledValenciavU[i] = ValenciavU[i]*sp.sqrt(Rmax/R)

    # We stick with the "rescaled" version since we redefine Rmax as detailed above
    # u^0 = 1/(alpha-sqrt(1-R_{\max}))
    rescaledu0 = 1/(alpha*sp.sqrt(1-Rmax))


# We'll write this as a function so that we can calculate the expressions on-demand for any choice of i
def find_cp_cm(lapse,shifti,gupii):
    # Inputs:  u0,vi,lapse,shift,gammadet,gupii
    # Outputs: cplus,cminus 
    
    # a = 1/(alpha^2)
    a = 1/(lapse*lapse)
    # b = 2 beta^i / alpha^2
    b = 2 * shifti /(lapse*lapse)
    # c = -g^{ii} + (beta^i)^2 / alpha^2
    c = - gupii + shifti*shifti/(lapse*lapse)
    
    # Now, we are free to solve the quadratic equation as usual. We take care to avoid passing a
    # negative value to the sqrt function.
    detm = b*b - 4*a*c
    detm = sp.sqrt(sp.Rational(1,2)*(detm + sp.Abs(detm)))
    global cplus,cminus
    cplus  = sp.Rational(1,2)*(-b/a + detm/a)
    cminus = sp.Rational(1,2)*(-b/a - detm/a)

# We'll write this as a function, and call it within HLLE_solver, below.
def find_cmax_cmin(flux_dirn,gamma_faceUU,beta_faceU,alpha_face,gammadet_face):
    # Inputs:  flux direction flux_dirn, Inverse metric gamma_faceUU, shift beta_faceU,
    #          lapse alpha_face, metric determinant gammadet_face
    # Outputs: maximum and minimum characteristic speeds cmax and cmin
    # First, we need to find the characteristic speeds on each face
    find_cp_cm(alpha_face,beta_faceU[flux_dirn],gamma_faceUU[flux_dirn][flux_dirn])
    cpr = cplus
    cmr = cminus
    find_cp_cm(alpha_face,beta_faceU[flux_dirn],gamma_faceUU[flux_dirn][flux_dirn])
    cpl = cplus
    cml = cminus
    
    # The following algorithms have been verified with random floats:
    
    global cmax,cmin
    # Now, we need to set cmax to the larger of cpr,cpl, and 0
    cmax = sp.Rational(1,2)*(cpr+cpl+sp.Abs(cpr-cpl))
    cmax = sp.Rational(1,2)*(cmax+sp.Abs(cmax))
    
    # And then, set cmin to the smaller of cmr,cml, and 0
    cmin =  sp.Rational(1,2)*(cmr+cml-sp.Abs(cmr-cml))
    cmin = -sp.Rational(1,2)*(cmin-sp.Abs(cmin))

# We'll rewrite this assuming that we've passed the entire reconstructed
# gridfunctions. You could also do this with only one point, but then you'd 
# need to declare everything as a Cparam in NRPy+

def HLLE_solver(flux_dirn, mom_comp, alpha_face, gammadet_face, gamma_faceDD, gamma_faceUU, beta_faceU, Valenciav_rU, B_rU, Valenciav_lU, B_lU): 
    # This solves the Riemann problem for the mom_comp component of the momentum
    # flux StildeD in the flux_dirn direction.

    alpha_sqrt_gamma = alpha_face*sp.sqrt(gammadet_face)

    # We already did something like this in GiRaFFE_Higher_Order_v2; we'll 
    # paste that code here, but using the face variables. This assumes that
    # we've already interpolated the metric quantities on the faces.
    import u0_smallb_Poynting__Cartesian.u0_smallb_Poynting__Cartesian as u0b

    # Using the function coded above, we find the speed-limited v^i and u^0 on both 
    # left and right faces.
    compute_u0_noif(gamma_faceDD,alpha_face,Valenciav_rU)
    Valenciav_rU = rescaledValenciavU
    u4upperZero_r = rescaledu0
    compute_u0_noif(gamma_faceDD,alpha_face,Valenciav_lU)
    Valenciav_lU = rescaledValenciavU
    u4upperZero_l = rescaledu0

    # Note the Kronecker delta symbol in the above; let's create it real quick:
    # We'll zero-offset everything here, unlike in the original GiRaFFE.
    kronecker_delta = ixp.zerorank2()
    for i in range(DIM): 
        kronecker_delta[i][i] = 1

    # After importing the relevant module, we'll run the function contained therein 
    # with the interpolated metric face values and the reconstructed velocities
    # and magnetic fields. Then, we'll simply read in the expressions that that 
    # function set, renaming them as convenient and necessary. In particular, we 
    # add _r and _l tags and replace u0 with u4upperZero, which is necessary to avoid
    # NRPy+ interpreting u0 as something it's not.
    # LEFT
    u0b.compute_u0_smallb_Poynting__Cartesian(gamma_faceDD,beta_faceU,alpha_face,Valenciav_lU,B_lU)
    u_lD = ixp.zerorank1()
    u_lU = ixp.zerorank1()

    for i in range(DIM):
        u_lD[i] = u0b.uD[i].subs(u0b.u0,u4upperZero_l)
        u_lU[i] = u0b.uU[i].subs(u0b.u0,u4upperZero_l)

    smallb4_lU = ixp.zerorank1(DIM=4)
    smallb4_lD = ixp.zerorank1(DIM=4)
    for mu in range(4):
        smallb4_lU[mu] = u0b.smallb4U[mu].subs(u0b.u0,u4upperZero_l)
        smallb4_lD[mu] = u0b.smallb4D[mu].subs(u0b.u0,u4upperZero_l)

    smallb2_l = u0b.smallb2etk.subs(u0b.u0,u4upperZero_l)

    # We must repeat this process on both faces.
    # RIGHT
    u0b.compute_u0_smallb_Poynting__Cartesian(gamma_faceDD,beta_faceU,alpha_face,Valenciav_rU,B_rU)
    u_rD = ixp.zerorank1()
    u_rU = ixp.zerorank1()

    for i in range(DIM):
        u_rD[i] = u0b.uD[i].subs(u0b.u0,u4upperZero_r)
        u_rU[i] = u0b.uU[i].subs(u0b.u0,u4upperZero_r)

    smallb4_rU = ixp.zerorank1(DIM=4)
    smallb4_rD = ixp.zerorank1(DIM=4)
    for mu in range(4):
        smallb4_rU[mu] = u0b.smallb4U[mu].subs(u0b.u0,u4upperZero_r)
        smallb4_rD[mu] = u0b.smallb4D[mu].subs(u0b.u0,u4upperZero_r)

    smallb2_r = u0b.smallb2etk.subs(u0b.u0,u4upperZero_r)

    
    # There's a caveat here to consider with the indices here: remember the general rule
    # that a latin index must be incremented if it's in a 4-vector. 
    # f = alpha sqrt{gamma} ( b^2 u^m u_j + (1/2) b^2 delta^m_j - b^m b_j )
    Fr = alpha_sqrt_gamma * (smallb2_r * u_rU[flux_dirn] * u_rD[mom_comp] 
                             + smallb2_r * kronecker_delta[flux_dirn][mom_comp] / 2.0
                             - smallb4_rU[flux_dirn+1] * smallb4_rD[mom_comp+1])
    
    # f = alpha sqrt{gamma} ( b^2 u^m u_j + (1/2) b^2 delta^m_j - b^m b_j )
    Fl = alpha_sqrt_gamma * (smallb2_l * u_lU[flux_dirn] * u_lD[mom_comp] 
                             + smallb2_l * kronecker_delta[flux_dirn][mom_comp] / 2.0
                             - smallb4_lU[flux_dirn+1] * smallb4_lD[mom_comp+1])
    
    # st_j = alpha sqrt{gamma} ( b^2 u^0 u_j - b^0 b_j )
    st_j_r = alpha_sqrt_gamma * (smallb2_r*u4upperZero_r*u_rD[mom_comp] - smallb4_rU[0]*smallb4_rD[mom_comp+1])
    st_j_l = alpha_sqrt_gamma * (smallb2_l*u4upperZero_l*u_lD[mom_comp] - smallb4_lU[0]*smallb4_lD[mom_comp+1])
    
    find_cmax_cmin(flux_dirn,gamma_faceUU,beta_faceU,alpha_face,gammadet_face)
    
    # st_j_flux = (c_\min f_R + c_\max f_L - c_\min c_\max ( st_j_r - st_j_l )) / (c_\min + c_\max)
    return (cmin*Fr + cmax*Fl - cmin*cmax*(st_j_r-st_j_l) )/(cmax + cmin)

global calculate_Stilde_flux
# Zach says: Patrick -- check on whether gamma_faceDD is an exact (within roundoff) inverse of gamma_faceUU. Otherwise make it so!
def calculate_Stilde_flux(flux_dirn, inputs_provided=False, alpha_face=None, gammadet_face=None, gamma_faceDD=None, gamma_faceUU=None, beta_faceU=None, Valenciav_rU=None, B_rU=None, Valenciav_lU=None, B_lU=None):

    if not inputs_provided:
        # We will pass values of the gridfunction on the cell faces into the function. This requires us
        # to declare them as C parameters in NRPy+. We will denote this with the _face infix/suffix.
        alpha_face,gammadet_face = gri.register_gridfunctions("AUXEVOL",["alpha_face","gammadet_face"])
        gamma_faceDD = ixp.register_gridfunctions_for_single_rank2("AUXEVOL","gamma_faceDD","sym01")
        gamma_faceUU = ixp.register_gridfunctions_for_single_rank2("AUXEVOL","gamma_faceUU","sym01")
        beta_faceU = ixp.register_gridfunctions_for_single_rank1("AUXEVOL","beta_faceU")

        # We'll need some more gridfunctions, now, to represent the reconstructions of BU and ValenciavU
        # on the right and left faces
        Valenciav_rU = ixp.register_gridfunctions_for_single_rank1("AUXEVOL","Valenciav_rU",DIM=3)
        B_rU = ixp.register_gridfunctions_for_single_rank1("AUXEVOL","B_rU",DIM=3)
        Valenciav_lU = ixp.register_gridfunctions_for_single_rank1("AUXEVOL","Valenciav_lU",DIM=3)
        B_lU = ixp.register_gridfunctions_for_single_rank1("AUXEVOL","B_lU",DIM=3)



    # We'll also compute some powers of the conformal factor psi:
    # Since psi^6 = sqrt(gammadet) by definition,
    # psi = sp.sqrt(gammadet_face)**(1.0/6.0)
    # psim4 = psi**(-4.0)

    global Stilde_fluxD
    Stilde_fluxD = ixp.zerorank3()

    for mom_comp in range(DIM):
        Stilde_fluxD[mom_comp] = HLLE_solver(flux_dirn,mom_comp,alpha_face, gammadet_face, gamma_faceDD, gamma_faceUU, beta_faceU, Valenciav_rU, B_rU, Valenciav_lU, B_lU)
