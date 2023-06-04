from traceback import print_exc
import grid
import os
import subprocess
from sympy import sympify
from sympy import Eq, Ne, Piecewise, Rational
from sympy import pi
from sympy.functions.elementary.miscellaneous import cbrt, Max, sqrt

import NRPy_param_funcs as par
import indexedexp as ixp
from cactusthorn import CactusThorn, loop
from outputC import lhrh, outCparams

import safewrite
import sympy as sp
from sugar import *

def main():
    # safewrite.nochange = True
    safewrite.verbose = True

    outCparams.CSE_enable = "True"

    # "CCTK_REAL" or "CCTK_REALVEC"
    par.set_parval_from_str("PRECISION", "CCTK_REALVEC")

    # Current options are Carpet and CarpetX
    grid.ET_driver = "CarpetX"

    thorn = CactusThorn("CarpetXNRPy", "Z4cNRPy")

    fd_order = thorn.declare_param('fd_order', default=4, vmin=2, vmax=8, doc="Finite differencing order")

    fd_order = 4
    par.set_parval_from_str("finite_difference::FD_CENTDERIVS_ORDER", fd_order)

    centering='VVV'

    # EVOL: evolved grid functions (possibly multiple time levels)
    # AUXEVOL: needed for evolution, can be freed after evaluating rhs (1 time level)
    # AUX: e.g. RHS (1 time level)
    # TILE_TMP: not actually a grid function, only temporary tile
    # SCALAR_TMP: not actually a grid function, only a double

    # Generic helpers

    def namefun(symbol, index, shape, prefix):
        symbol = prefix
        result = [sp.Symbol(symbol + ''.join(ixnam(n) for n in index + [i]))
                if symbol else sp.sympify(0) for i in range(shape[0])]
        return result

    set_Theta_zero = thorn.declare_param('set_Theta_zero', default=False, doc="set_Theta_zero")
    kappa1 = thorn.declare_param('kappa1', default=0.02, doc="kappa1")
    kappa2 = thorn.declare_param('kappa2', default=0.0, doc="kappa2")
    f_mu_L = thorn.declare_param('f_mu_L', default=2.0, doc="f_mu_L")
    f_mu_S = thorn.declare_param('f_mu_S', default=1.0, doc="f_mu_S")
    eta = thorn.declare_param('eta', default=2.0, doc="eta")
    alphaG_floor = thorn.declare_param('alphaG_floor', default=1.0e-10, vmin=0.0, doc="alphaG_floor")
    chi_floor = thorn.declare_param('chi_floor', default=1.0e-10, vmin=0.0, doc="chi_floor")
    eps_diss = thorn.declare_param('eps_diss', default=0.2, vmin=0.0, vmax=1.0, doc="eps_diss")
    zero = thorn.declare_param('zero', default=0.0, vmin=0.0, vmax=0.0, doc="zero")


    decl_indexes()
    set_coords("x","y","z")

    # Need to make a better way of setting
    # this kind of thing up
    one = sp.IndexedBase("one")
    definitions["oneU"] = [1]*3

    gfparams(gf_type="EXTERNAL",symmetries="sym01",centering=centering,external_module="ADMBase",namefun=name_xyz)
    gflatex(r"g_{i j} k_{i j} dtk_{i j} alp dtalp beta^i dtbeta^i dtdtbeta^i")
    gfdecl("dt2alp",[])
    gfparams(gf_type="EXTERNAL",symmetries="sym01",centering=centering,external_module="TmunuBase",namefun=name_xyz)
    gfdecl("eTtt",[],"eTt",[la],"eT",[la,lb])

    gfparams(gf_type="EVOL",symmetries="sym01",centering=centering)
    gfdecl("chi",[],"gammatilde",[la,lb],"Khat",[],"Atilde",[la,lb],"Gammatilde",[ua],"Theta","alphaG",[],"betaG",[ua])

    gfparams(gf_type="AUX",symmetries="sym01",centering=centering)
    gfdecl("rhs_chi",[],"rhs_gammatilde",[la,lb],"rhs_Khat",[],"rhs_Atilde",[la,lb], \
        "rhs_Gammatilde",[ua],"rhs_Theta","rhs_alphaG",[], "rhs_betaG", [ua],"ZetatildeC",[ua],\
        "HC",[],"MtildeC",[ua],"allC",[])

    # Derivatives as calculated by NRPy

    deriv_decl(chi, ("_d",[li]), ("_dup",[li]), ("_ddn",[li]), ("_dKO",[li]), ("_d",[li,lj]))
    deriv_decl(gammatilde[la,lb], ("_d",[lc]), ("_d",[lc,ld]), ("_dKO",[lc]))
    deriv_decl(Khat, ("_d",[la]), ("_dKO",[la]))
    deriv_decl(Atilde[la,lb], ("_d",[la]), ("_dKO",[la]))
    deriv_decl(Theta, ("_d",[la]), ("_dKO",[la]))
    deriv_decl(Gammatilde[ua], ("_d",[lb]), ("_dKO",[lb]))

    deriv_decl(alphaG, ("_d",[la]), ("_d",[la,lb]), ("_dKO",[la]))
    deriv_decl(betaG[uc], ("_d",[la]), ("_d",[la,lb]), ("_dKO",[la]))

    # Our tile-local derivatives
    gfparams(gf_type="TILE_TMP",symmetries="sym01",centering=centering)
    gfdecl("dchi",[la],"ddchi",[la,lb],"dbetaGchi","disschi",[],"dgammatilde",[la,lb,lc])

    gfparams(gf_type="TILE_TMP",symmetries="sym01_sym23",centering=centering)
    gfdecl("ddgammatilde",[la,lb,lc,ld])

    gfparams(gf_type="TILE_TMP",symmetries="sym01",centering=centering)
    gfdecl("dissgammatilde",[la,lb],"dKhat",[la],"dissKhat",[],"dAtilde",[la,lb,lc],"dissAtilde",[la,lb],\
        "dTheta",[la],"dissTheta",[])

    gfparams(gf_type="TILE_TMP",symmetries="",centering=centering)
    gfdecl("dGammatilde",[ua,lb],"dissGammatilde",[ua],"dalphaG",[la])

    gfparams(gf_type="TILE_TMP",symmetries="sym01",centering=centering)
    gfdecl("ddalphaG",[la,lb],"dissalphaG",[])

    gfparams(gf_type="TILE_TMP",symmetries="",centering=centering)
    gfdecl("dbetaG",[ua,lb])

    gfparams(gf_type="TILE_TMP",symmetries="sym12",centering=centering)
    gfdecl("ddbetaG",[ua,lb,lc],"dissbetaG",[ua])

    # Expressions

    gUU_expr, detg_expr = ixp.symm_matrix_inverter3x3(gDD)
    gammatildeUU_expr, detgammatilde_expr = ixp.symm_matrix_inverter3x3(gammatildeDD)

    # Initial conditions

    gfparams(gf_type="SCALAR_TMP",symmetries="sym01",centering=centering)
    gfdecl("g","gammatilde",[ua,ub],"dgammatilde",[ua,ub,lc],"Theta_val","trK","chi_val",[])

    def Initial():
        initial1_eqns = [
            [lhrh(lhs=chi_val, rhs=1 / cbrt(detg_expr))],
            geneqns(lhs=g[ua,ub],values=gUU_expr),
            [lhrh(lhs=Theta_val, rhs=sympify(0))],
            geneqns(lhs=trK, rhs=g[ua,ub]*k[la,lb]),
            [lhrh(lhs=chi, rhs=chi_val)],
            geneqns(r"\tilde{\gamma}_{i j} = \text{chi_val} g_{i j}"),
            [lhrh(lhs=Theta, rhs=Theta_val)],
            [lhrh(lhs=Khat, rhs=trK - 2 * Theta_val)],
            geneqns(r'\tilde{A}_{i j} = \text{chi_val} (k_{i j} - \text{trK}/3 g_{i j})'),
            [lhrh(lhs=alphaG, rhs=alp)],
            geneqns(lhs=betaG[ui], rhs=r"\beta^i"),
        ]

        thorn.add_func(
            "Z4cNRPy_Initial1",
            body=initial1_eqns,
            where='everywhere',
            schedule_bin="Z4cNRPy_InitialGroup",
            doc="Convert ADM to Z4c variables, part 1",
            centering=centering)

        initial2_eqns = [
            geneqns(lhs=dgammatilde[li,lj,lk], values=gammatildeDD_dD),
            loop,
            geneqns(lhs=gammatilde[ua,ub],values=gammatildeUU_expr),
            geneqns(lhs=dgammatilde[ui,uj,lk], rhs=gammatilde[ui,ux] * gammatilde[uj,uy] * dgammatilde[lx,ly,lk]),
            geneqns(lhs=Gammatilde[ui], rhs=dgammatilde[ui,ux,lx])
        ]

        thorn.add_func(
            "Z4cNRPy_Initial2",
            body=initial2_eqns,
            where='interior',
            schedule_bin="Z4cNRPy_InitialGroup AFTER Z4cNRPy_Initial1",
            doc="Convert ADM to Z4c variables, part 2",
            sync="GammatildeUGF",
            centering=centering)
    Initial()

    # Enforce constraints

    gfparams(gf_type="SCALAR_TMP",centering=centering)
    gfdecl("detgammatilde","trAtilde",[])

    def Enforce():
        enforce_eqns = [
            # Enforce floors
            geneqns(lhs=gammatilde[ua,ub],values=gammatildeUU_expr),
            [lhrh(lhs=detgammatilde, rhs=detgammatilde_expr)],
            geneqns(lhs=trAtilde, rhs=gammatilde[ux,uy]*Atilde[lx,ly]),
            [lhrh(lhs=chi, rhs=Max(chi_floor, chi))],
            [lhrh(lhs=alphaG, rhs=Max(alphaG_floor, alphaG))],
            # Enforce algebraic constraints; see arXiv:1212.2901 [gr-qc]
            geneqns(lhs=gammatilde[li,lj], rhs=(1 / cbrt(detgammatilde)) * gammatilde[li,lj]),
            geneqns(lhs=Atilde[li,lj], rhs=(Atilde[li,lj] - trAtilde / 3 * gammatilde[li,lj]))
        ]

        thorn.add_func(
            "Z4cNRPy_Enforce",
            body=enforce_eqns,
            where='interior',
            schedule_bin="Z4cNRPy_PostStepGroup",
            doc="Enforce constraints and synchronize",
            # We need to sync the whole state vector here; nothing has synced it yet
            sync="chiGF gammatildeDDGF KhatGF AtildeDDGF ThetaGF GammatildeUGF alphaGGF betaGUGF",
            centering=centering
            )
    Enforce()

    # Calculate ADM variables
    def ADM():
        adm_eqns = [
            geneqns(r"g_{i j} = 1/\chi \tilde{\gamma}_{i j}"),
            geneqns(r"k_{i j} = 1/\chi (\tilde{A}_{i j} + (\hat{K} + 2 \Theta) / 3 \tilde{\gamma}_{i j})"),
            [lhrh(lhs=alp, rhs=alphaG)],
            [lhrh(lhs=dtalp, rhs=-alphaG * f_mu_L * Khat)],
            geneqns(r"\beta^i = \text{betaG}^i"),
            geneqns(lhs=dtbeta[ui], rhs=f_mu_S * Gammatilde[ui] - eta * betaG[ui]),

            # Second time derivatives
            # TODO: Use correct expressions
            geneqns(lhs=dtk[li,lj], rhs=zero * Atilde[li,lj]),
            geneqns(lhs=dt2alp, rhs=sympify(0)),
            geneqns(lhs=dtdtbeta[ui], rhs=zero * betaG[ui]),
        ]

        thorn.add_func(
            "Z4cNRPy_ADM",
            body=adm_eqns,
            # where='interior',
            where='everywhere',
            schedule_bin="Z4cNRPy_PostStepGroup AFTER Z4cNRPy_Enforce",
            doc="Calculate ADM variables",
            # sync="ADMBase::metric ADMBase::curv ADMBase::lapse ADMBase::dtlapse ADMBase::shift ADMBase::dtshift",
            centering=centering)
    ADM()

    # Calculate RHS

    gfparams(gf_type="SCALAR_TMP",symmetries="sym01",centering=centering)
    gfdecl("dg",[la,lb,lc])

    gfparams(gf_type="SCALAR_TMP",symmetries="sym12",centering=centering)
    gfdecl("Christoffel",[la,lb,lc],"Christoffel",[ua,lb,lc])

    gfparams(gf_type="SCALAR_TMP",symmetries="sym01",centering=centering)
    gfdecl("DDalphaG",[la,lb])

    gfparams(gf_type="SCALAR_TMP",symmetries="",centering=centering)
    gfdecl("DbetaG", [ua,lb], "rho", [], "Si",[la],"Sij",[la,lb],"trS",[],"Atilde",[ua,ld])

    gfparams(gf_type="SCALAR_TMP",symmetries="sym01",centering=centering)
    gfdecl("Atilde",[ua,ub],"R",[la,lb],"Rchi",[la,lb],"Rtilde",[la,lb],"alphaRicciTmunu",[la,lb])
    gfparams(gf_type="SCALAR_TMP",symmetries="",centering=centering)
    gfdecl("Christoffeltilde",[ua,lb,uc],
           "Christoffeltilde",[la,lb,uc],
           "GammatildeChristoffel",[ua],
           "dAtilde",[ua,ub,lc])
    gfdecl("DDchi",[la,lb],"ZetatildeCval",[ua],"HCval",[],"MtildeCval",[ua])
    gfparams(gf_type="SCALAR_TMP",symmetries="sym12",centering=centering)
    gfdecl("Christoffeltilde",[la,lb,lc],
           "Christoffeltilde",[ua,lb,lc],
           "Rsc",[],
           "dAtilde",[ua,lb,lc])

    def RHS():
        rhs_eqns = [
            # Derivatives
            geneqns(lhs=dchi[li], values=chi_dD),
            geneqns(lhs=ddchi[li,lj], values=chi_dDD),
            loop,
            geneqns(lhs=dgammatilde[li,lj,lk], values=gammatildeDD_dD),
            loop,
            geneqns(lhs=ddgammatilde[li,lj,lk,ll], values=gammatildeDD_dDD),
            loop,
            geneqns(lhs=dGammatilde[ui,lj], values=GammatildeU_dD),
            geneqns(lhs=dissGammatilde[ui], rhs=Gammatilde1_dKO1[ui,lx]*one[ux]),
            geneqns(lhs=dKhat[li], rhs=Khat_d1[li]),
            geneqns(lhs=dAtilde[li,lj,lk], rhs=Atilde2_d1[li,lj,lk],loop=True),
            geneqns(lhs=dTheta[li], rhs=Theta_d1[li],loop=True),
            geneqns(lhs=dalphaG[li], rhs=alphaG_d1[li],loop=True),
            geneqns(lhs=ddalphaG[li,lj], rhs=alphaG_d2[li,lj],loop=True),
            geneqns(lhs=dbetaG[ui,lj], rhs=betaG1_d1[ui,lj],loop=True),
            geneqns(lhs=ddbetaG[ui,lj,lk], rhs=betaG1_d2[ui,lj,lk],loop=True),

            # RHS
            geneqns(lhs=gammatilde[ui,uj], values=gammatildeUU_expr),
            geneqns(lhs=g[ui,uj], rhs=chi*gammatilde[ui,uj]),
            geneqns(lhs=dg[li,lj,lk],
                    rhs=(- 1 / chi**2 * dchi[lk] * gammatilde[li,lj]
                         + 1 / chi * dgammatilde[li,lj,lk])),
            geneqns(lhs=Christoffel[li,lj,lk],
                  rhs=Rational(1,2) * (dg[li,lj,lk] + dg[li,lk,lj] - dg[lj,lk,li])),
            geneqns(lhs=Christoffel[ui,lj,lk], rhs=g[ui,ux]*Christoffel[lx,lj,lk]),
            geneqns(DbetaG[ui,lj], rhs=dbetaG[ui,lj] + Christoffel[ui,lx,lj]* betaG[ux]),
            # arXiv:1212.2901 [gr-qc], (1)
            geneqns(lhs=rhs_chi, rhs=Rational(2,3)*chi*(alphaG*(Khat+2*Theta)-DbetaG[ux,lx])),

            # arXiv:1212.2901 [gr-qc], (2)
            geneqns(lhs=rhs_gammatilde[li,lj],
                    rhs=(-2 * alphaG * Atilde[li,lj]
                         + (betaG[ux] * dgammatilde[li,lj,lx])
                         + (gammatilde[lx,li] * dbetaG[ux,lj] + gammatilde[lx,lj] * dbetaG[ux,li])
                         - Rational(2,3) * (gammatilde[li,lj] * dbetaG[ux,lx])
                         # + eps_diss * dissgammatildeDD[i][j]
                         )),

            geneqns(lhs=DDalphaG[li,lj], rhs=ddalphaG[li,lj] - Christoffel[ux,li,lj] * dalphaG[lx]),
            geneqns(lhs=Atilde[ui,lj], rhs=gammatilde[ui,ux] * Atilde[lx,lj]),
            geneqns(lhs=rho, rhs=1 / alphaG**2 * (+ eTtt - 2 * betaG[ux] * eTt[lx]
                                                  + betaG[ux] * betaG[uy] * eT[lx,ly])),
            geneqns(lhs=Si[li], rhs=-1 / alphaG * (eTt[li] - betaG[ux] * eT[li,lx])),
            geneqns(lhs=Sij[li,lj], rhs=eT[li,lj]),
            geneqns(lhs=trS, rhs=g[ux,uy] * Sij[lx,ly]),
            # arXiv:1212.2901 [gr-qc], (3)
            geneqns(lhs=rhs_Khat,
                    rhs=(- (g[ux,uy] * DDalphaG[lx,ly])
                         + alphaG * (+ (Atilde[ux,ly] * Atilde[uy,lx])
                                     + Rational(1,3) * (Khat + 2 * Theta)**2)
                         + 4 * pi * alphaG * (trS + rho)
                         + alphaG * kappa1 * (1 - kappa2) * Theta
                         + (betaG[ux] * dKhat[lx])
                         # + eps_diss * dissKhat
                         )),
            geneqns(lhs=Christoffeltilde[li,lj,lk],
                    rhs=Rational(1,2) * (dgammatilde[li,lj,lk] + dgammatilde[li,lk,lj] - dgammatilde[lj,lk,li])),
            geneqns(lhs=Christoffeltilde[ui,lj,lk], rhs=gammatilde[ui,ux]*Christoffeltilde[lx,lj,lk]),
            geneqns(lhs=Christoffeltilde[li,lj,uk], rhs=gammatilde[uk,ux]*Christoffeltilde[li,lj,lx]),
            geneqns(lhs=Christoffeltilde[ui,lj,uk], rhs=Christoffeltilde[ui,lj,lx]*gammatilde[ux,uk]),
            geneqns(lhs=GammatildeChristoffel[ui], rhs=Christoffeltilde[ui,lx,ux]),
            geneqns(lhs=DDchi[li,lj], rhs=ddchi[li,lj] - Christoffeltilde[ux,li,lj] * dchi[lx]),
            # arXiv:1212.2901 [gr-qc], (8)
            geneqns(lhs=Rchi[li,lj],
                    rhs=(+ Rational(1,2) / chi * DDchi[li,lj]
                         + Rational(1,2) / chi * gammatilde[li,lj] * gammatilde[ux,uy] * DDchi[lx,ly]
                         - Rational(1,4) / chi**2 * dchi[li] * dchi[lj]
                         - Rational(3,4) / chi**2 * (gammatilde[li,lj] *
                                                     gammatilde[ux,uy] * dchi[lx] * dchi[ly]))),
            # arXiv:1212.2901 [gr-qc], (9)
            geneqns(lhs=Rtilde[li,lj],
                    rhs=(- Rational(1,2) * gammatilde[ux,uy] * ddgammatilde[li,lj,lx,ly]
                         + Rational(1,2) * (+ gammatilde[lx,li] * dGammatilde[ux,lj]
                                            + gammatilde[lx,lj] * dGammatilde[ux,li])
                         + Rational(1,2) * (+ GammatildeChristoffel[ux] * Christoffeltilde[li,lj,lx]
                                            + GammatildeChristoffel[ux] * Christoffeltilde[lj,li,lx])
                         + (+ Christoffeltilde[ux,ly,li] * Christoffeltilde[lj,lx,uy]
                            + Christoffeltilde[ux,ly,lj] * Christoffeltilde[li,lx,uy]
                            + Christoffeltilde[ux,li,ly] * Christoffeltilde[lx,lj,uy]))),
            # arXiv:1212.2901 [gr-qc], (7)
            geneqns(lhs=R[li,lj], rhs=Rchi[li,lj]+Rtilde[li,lj]),
            geneqns(lhs=Rsc, rhs=g[ux,uy]*R[lx,ly]),
            geneqns(lhs=alphaRicciTmunu[li,lj],
                  rhs=- DDalphaG[li,lj] + alphaG * (R[li,lj] - 8 * pi * Sij[li,lj])),
            # arXiv:1212.2901 [gr-qc], (4)
            geneqns(lhs=rhs_Atilde[li,lj],
                    rhs=(+ chi * (alphaRicciTmunu[li,lj]
                                  - (Rational(1,3) * gammatilde[li,lj] *
                                     gammatilde[ux,uy] * alphaRicciTmunu[lx,ly]))
                         + alphaG * (+ (Khat + 2 * Theta) * Atilde[li,lj]
                                     - 2 * Atilde[lx,li] * Atilde[ux,lj])
                         + betaG[ux] * dAtilde[li,lj,lx]
                         + Atilde[lx,li] * dbetaG[ux,lj] + Atilde[lx,lj] * dbetaG[ux,li]
                         - Rational(2,3) * Atilde[li,lj] * dbetaG[ux,lx]
                         # + eps_diss * dissAtildeDD[i][j]
                         )),

            geneqns(lhs=Atilde[ui,uj],
                    rhs=gammatilde[ui,ux] * Atilde[uj,lx]),
            # arXiv:1212.2901 [gr-qc], (5)
            geneqns(lhs=rhs_Gammatilde[ui],
                    rhs=(- 2 * Atilde[ui,ux] * dalphaG[lx]
                         + 2 * alphaG *
                         (+ Christoffeltilde[ui,lx,ly] * Atilde[ux,uy]
                          - Rational(3,2) / chi * Atilde[ui,ux] * dchi[lx]
                          - Rational(1,3) * gammatilde[ui,ux] * (2 * dKhat[lx] + dTheta[lx])
                          - 8 * pi * gammatilde[ui,ux] * Si[lx])
                         + gammatilde[ux,uy] * ddbetaG[ui,lx,ly]
                         + Rational(1,3) * gammatilde[ui,ux] * ddbetaG[uy,lx,ly]
                         + betaG[ux] * dGammatilde[ui,lx]
                         - GammatildeChristoffel[ux] * dbetaG[ui,lx]
                         + Rational(2,3) * GammatildeChristoffel[ui] * dbetaG[ux,lx]
                         - 2 * alphaG * kappa1 * Gammatilde[ui] - GammatildeChristoffel[ui]
                         # + eps_diss * dissGammatildeU[i]
                         )),

            # arXiv:1212.2901 [gr-qc], (6)
            geneqns(lhs=rhs_Theta,
                    rhs=Piecewise((0, set_Theta_zero),
                                  ((+ Rational(1,2) * alphaG * (+ Rsc
                                                                - Atilde[ux,uy] * Atilde[lx,ly]
                                                                + Rational(2,3) * (Khat + 2 * Theta)**2)
                                    - alphaG * (+ 8 * pi * rho
                                                + kappa1 * (2 + kappa2) * Theta)
                                    + betaG[ux] * dTheta[lx]
                                    # + eps_diss * dissTheta
                                    ), True))),

            [lhrh(lhs=rhs_alphaG,
                  rhs=(- alphaG * f_mu_L * Khat
                       # + eps_diss * dissalphaG
                       ))],

            geneqns(lhs=rhs_betaG[ui], rhs=(+ f_mu_S * Gammatilde[ui] - eta * betaG[ui])),
            loop,

            geneqns(lhs=rhs_chi, rhs=rhs_chi + eps_diss * chi_dKO1[lx] * one[ux]),
            loop,
            geneqns(lhs=rhs_gammatilde[li,lj], rhs=rhs_gammatilde[li,lj] + eps_diss * gammatilde2_dKO1[li,lj,lx]*one[ux]),
            loop,
            geneqns(lhs=rhs_Gammatilde[ui], rhs=rhs_Gammatilde[ui] + eps_diss * Gammatilde1_dKO1[ui,lx]*one[ux]),
            loop,
            geneqns(lhs=rhs_Khat, rhs=rhs_Khat + eps_diss *  Khat_dKO1[lx] * one[ux]),
            loop,
            geneqns(lhs=rhs_Atilde[li,lj], rhs=rhs_Atilde[li,lj] + eps_diss * Atilde2_dKO1[li,lj,lx] * one[ux]),
            loop,
            geneqns(lhs=rhs_Theta, rhs=rhs_Theta + eps_diss * Theta_dKO1[lx]*one[ux]),
            loop,
            geneqns(lhs=rhs_alphaG, rhs=rhs_alphaG + eps_diss * alphaG_dKO1[lx]*one[ux]),
            loop,
            geneqns(lhs=rhs_betaG[ui], rhs=rhs_betaG[ui] + eps_diss * betaG1_dKO1[ui,lx]*one[ux]),
            loop,
        ]

        thorn.add_func(
            "Z4cNRPy_RHS",
            body=rhs_eqns,
            where='interior',
            schedule_bin="ODESolvers_RHS",
            doc="Calculate RHS",
            sync="rhs_chiGF rhs_gammatildeDDGF rhs_KhatGF rhs_AtildeDDGF rhs_GammatildeUGF rhs_ThetaGF rhs_alphaGGF rhs_betaGUGF",
            centering=centering)
    RHS()

    # Evaluate constraints

    def Constraints():
        constraints_eqns = [
            # Derivatives
            geneqns(lhs=dchi[li], rhs=chi_d1[li]),
            geneqns(lhs=ddchi[li,lj], rhs=chi_d2[li,lj]),
            loop,
            geneqns(lhs=dgammatilde[li,lj,lk], rhs=gammatilde2_d1[li,lj,lk]),
            geneqns(lhs=ddgammatilde[li,lj,lk,ll], rhs=gammatilde2_d2[li,lj,lk,ll]),
            loop,
            geneqns(lhs=dGammatilde[ui,lj], rhs=Gammatilde1_d1[ui,lj]),
            loop,
            geneqns(lhs=dAtilde[li,lj,lk], rhs=Atilde2_d1[li,lj,lk]),
            loop,
            geneqns(lhs=dKhat[li], rhs=Khat_d1[li]),
            loop,
            geneqns(lhs=dTheta[li], rhs=Theta_d1[li]),
            loop,

            # Constraints

            # ZetatildeC, eqn. (13)
            geneqns(lhs=gammatilde[ui,uj], values=gammatildeUU_expr),
            geneqns(lhs=Christoffeltilde[li,lj,lk],
                    rhs=Rational(1,2) * (dgammatilde[li,lj,lk] + dgammatilde[li,lk,lj] - dgammatilde[lj,lk,li])),
            geneqns(lhs=Christoffeltilde[ui,lj,lk],
                    rhs=gammatilde[ui,ux] * Christoffeltilde[lx,lj,lk]),
            geneqns(lhs=Christoffeltilde[ui,lj,uk],
                    rhs=Christoffeltilde[ui,lj,lx] * gammatilde[ux,uk]),
            geneqns(lhs=GammatildeChristoffel[ui], rhs=Christoffeltilde[ui,lx,ux]),
            geneqns(lhs=ZetatildeCval[ui],
                    rhs=(Gammatilde[ui] - GammatildeChristoffel[ui]) / 2),
            geneqns(lhs=ZetatildeC[ui], rhs=(ZetatildeCval[ui])),

            # HC, eqn. (14)
            geneqns(lhs=g[ui,uj], rhs=chi*gammatilde[ui,uj]),
            geneqns(lhs=dg[li,lj,lk],
                    rhs=(- 1 / chi**2 * dchi[lk] * gammatilde[li,lj]
                         + 1 / chi * dgammatilde[li,lj,lk])),
            geneqns(lhs=Christoffel[li,lj,lk],
                  rhs=Rational(1,2) * (dg[li,lj,lk] + dg[li,lk,lj] - dg[lj,lk,li])),
            geneqns(lhs=Christoffel[ui,lj,lk], rhs=g[ui,ux]*Christoffel[lx,lj,lk]),
            geneqns(Atilde[ui,lj], rhs=gammatilde[ui,ux]*Atilde[lx,lj]),
            geneqns(lhs=rho,
                    rhs=1 / alphaG**2 * (+ eTtt
                                         - 2 * betaG[ux] * eTt[lx]
                                         + betaG[ux] * betaG[uy] * eT[lx,ly])),
            geneqns(lhs=DDchi[li,lj],
                  rhs=ddchi[li,lj] - Christoffeltilde[ux,li,lj] * dchi[lx]),
            # arXiv:1212.2901 [gr-qc], (8)
            geneqns(lhs=Rchi[li,lj],
                    rhs=(+ Rational(1,2) / chi * DDchi[li,lj]
                         + Rational(1,2) / chi * gammatilde[li,lj] * gammatilde[ux,uy] * DDchi[lx,ly]
                         - Rational(1,4) / chi**2 * dchi[li] * dchi[lj]
                         - Rational(3,4) / chi**2 * (gammatilde[li,lj] *
                                                     gammatilde[ux,uy] * dchi[lx] * dchi[ly]))),
            geneqns(lhs=Christoffeltilde[li,lj,uk], rhs=gammatilde[uk,ux]*Christoffeltilde[li,lj,lx]),
            # arXiv:1212.2901 [gr-qc], (9)
            geneqns(lhs=Rtilde[li,lj],
                    rhs=(- Rational(1,2) * gammatilde[ux,uy] * ddgammatilde[li,lj,lx,ly]
                         + Rational(1,2) * (+ gammatilde[lx,li] * dGammatilde[ux,lj]
                                            + gammatilde[lx,lj] * dGammatilde[ux,li])
                         + Rational(1,2) * (+ GammatildeChristoffel[ux] * Christoffeltilde[li,lj,lx]
                                            + GammatildeChristoffel[ux] * Christoffeltilde[lj,li,lx])
                         + (+ Christoffeltilde[ux,ly,li] * Christoffeltilde[lj,lx,uy]
                            + Christoffeltilde[ux,ly,lj] * Christoffeltilde[li,lx,uy]
                            + Christoffeltilde[ux,li,ly] * Christoffeltilde[lx,lj,uy]))),
            # arXiv:1212.2901 [gr-qc], (7)
            geneqns(lhs=R[li,lj], rhs=Rchi[li,lj] + Rtilde[li,lj]),
            geneqns(lhs=Rsc,
                    rhs=g[ux,uy] * R[lx,ly]),
            # arXiv:1212.2901 [gr-qc], (14)
            geneqns(lhs=HCval,
                    rhs=(+ Rsc
                         + Atilde[ux,ly] * Atilde[uy,lx]
                         - Rational(2,3) * (Khat + 2 * Theta)**2
                         - 16 * pi * rho)),
            [lhrh(lhs=HC,
                  rhs=HCval)],

            # MtildeC, eqn. (15)
            geneqns(lhs=Atilde[ui,uj], rhs=gammatilde[ui,ux]*Atilde[uj,lx]),
            geneqns(lhs=dgammatilde[ui,uj,lk],
                    rhs=-gammatilde[ui,ux] * gammatilde[uj,uy] * dgammatilde[lx,ly,lk]),
            geneqns(lhs=dAtilde[ui,lj,lk],
                    rhs=(+ dgammatilde[ui,ux,lk] * Atilde[lx,lj]
                         + gammatilde[ui,ux] * dAtilde[lx,lj,lk])),
            geneqns(lhs=dAtilde[ui,uj,lk],
                    rhs=(+ dgammatilde[ui,ux,lk] * Atilde[uj,lx]
                         + gammatilde[ui,ux] * dAtilde[uj,lx,lk])),
            geneqns(lhs=Si[li],
                    rhs=-1 / alphaG * (eTt[li] - betaG[ux] * eT[li,lx])),
            # arXiv:1212.2901 [gr-qc], (15)
            geneqns(lhs=MtildeCval[ui],
                    rhs=(+ dAtilde[ui,ux,lx]
                         + Christoffeltilde[ui,lx,ly] * Atilde[ux,uy]
                         - Rational(2,3) * gammatilde[ui,ux] * (dKhat[lx] + 2 * dTheta[lx])
                         - Rational(2,3) * Atilde[ui,ux] * dchi[lx] / chi
                         - 8 * pi * gammatilde[ui,ux] * Si[lx])),
            geneqns(lhs=MtildeC[ui], rhs=MtildeCval[ui]),

            # allC, arXiv:1111.2177, eqn. (73)
            geneqns(lhs=allC,
                    rhs=sqrt(Max(0, (+ HCval**2
                                     + gammatilde[lx,ly] * MtildeCval[ux] * MtildeCval[uy]
                                     + Theta**2
                                     + 2 * gammatilde[lx,ly] * ZetatildeCval[ux] * ZetatildeCval[uy])))),
        ]

        thorn.add_func(
            "Z4cNRPy_Constraints",
            body=constraints_eqns,
            where='interior',
            schedule_bin="Z4c_AnalysisGroup",
            doc="Evaluate Constraints",
            sync="ZetatildeCUGF HCGF MtildeCUGF allCGF",
            centering=centering)
    Constraints()

    # Generate thorn

    assert "CACTUS_HOME" in os.environ, "Please set the CACTUS_HOME variable to point to your Cactus installation"
    cactus_home = os.environ["CACTUS_HOME"]
    cactus_sim = os.environ.get("CACTUS_SIM", "sim")
    cactus_thornlist = os.environ.get("CACTUS_THORNLIST", None)

    schedule_raw = """
    # Define schedule groups

    SCHEDULE GROUP Z4cNRPy_InitialGroup AT initial AFTER ADMBase_PostInitial
    {
    } "Convert ADM to Z4c variables"

    SCHEDULE GROUP Z4cNRPy_PostStepGroup IN ODESolvers_PostStep BEFORE ADMBase_SetADMVars
    {
    } "Post-process Z4c variables"

    SCHEDULE GROUP Z4c_AnalysisGroup AT analysis
    {
    } "Analyse Z4c variables"
    """

    thorn.generate(cactus_home, cactus_config=cactus_sim, cactus_thornlist=cactus_thornlist, schedule_raw=schedule_raw)

    # Use correct name for ADMBase variables
    subprocess.run(["perl", "-pi", "-e", "s/dtdtalp/dt2alp/g,s/dtdtbeta/dt2beta/g",
                    cactus_home + "/arrangements/CarpetXNRPy/Z4cNRPy/schedule.ccl",
                    cactus_home + "/arrangements/CarpetXNRPy/Z4cNRPy/src/Z4cNRPy_ADM.cc",
                    ])

    # Correct schedule `reads` statements for `Enforce`
    subprocess.run(["patch", "-p1"],
                   cwd=cactus_home + "/arrangements/CarpetXNRPy",
                   input=bytes("""
    --- a/Z4cNRPy/schedule.ccl
    +++ b/Z4cNRPy/schedule.ccl
    @@ -107,20 +107,20 @@

     SCHEDULE Z4cNRPy_Enforce IN Z4cNRPy_PostStepGroup {
         LANG: C
    -    READS: Z4cNRPy::AtildeDD00GF(everywhere)
    -    READS: Z4cNRPy::AtildeDD01GF(everywhere)
    -    READS: Z4cNRPy::AtildeDD02GF(everywhere)
    -    READS: Z4cNRPy::AtildeDD11GF(everywhere)
    -    READS: Z4cNRPy::AtildeDD12GF(everywhere)
    -    READS: Z4cNRPy::AtildeDD22GF(everywhere)
    -    READS: Z4cNRPy::alphaGGF(everywhere)
    -    READS: Z4cNRPy::chiGF(everywhere)
    -    READS: Z4cNRPy::gammatildeDD00GF(everywhere)
    -    READS: Z4cNRPy::gammatildeDD01GF(everywhere)
    -    READS: Z4cNRPy::gammatildeDD02GF(everywhere)
    -    READS: Z4cNRPy::gammatildeDD11GF(everywhere)
    -    READS: Z4cNRPy::gammatildeDD12GF(everywhere)
    -    READS: Z4cNRPy::gammatildeDD22GF(everywhere)
    +    READS: Z4cNRPy::AtildeDD00GF(interior)
    +    READS: Z4cNRPy::AtildeDD01GF(interior)
    +    READS: Z4cNRPy::AtildeDD02GF(interior)
    +    READS: Z4cNRPy::AtildeDD11GF(interior)
    +    READS: Z4cNRPy::AtildeDD12GF(interior)
    +    READS: Z4cNRPy::AtildeDD22GF(interior)
    +    READS: Z4cNRPy::alphaGGF(interior)
    +    READS: Z4cNRPy::chiGF(interior)
    +    READS: Z4cNRPy::gammatildeDD00GF(interior)
    +    READS: Z4cNRPy::gammatildeDD01GF(interior)
    +    READS: Z4cNRPy::gammatildeDD02GF(interior)
    +    READS: Z4cNRPy::gammatildeDD11GF(interior)
    +    READS: Z4cNRPy::gammatildeDD12GF(interior)
    +    READS: Z4cNRPy::gammatildeDD22GF(interior)
         WRITES: Z4cNRPy::AtildeDD00GF(interior)
         WRITES: Z4cNRPy::AtildeDD01GF(interior)
         WRITES: Z4cNRPy::AtildeDD02GF(interior)
    """, "utf-8"))

if __name__ == "__main__":
    main()
