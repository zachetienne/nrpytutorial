from __future__ import annotations

import MoLtimestepping.C_Code_Generation as MoL
from MoLtimestepping.RK_Butcher_Table_Dictionary import Butcher_dict
from nrpytutorial.ccode_rendering import CprototypeFactory, CparameterFactory

from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from nrpytutorial.dirtree import Dirtree
    from nrpytutorial.integration_methods import IntegrationMethod

def _add_home(dirtree: Dirtree):
    dirtree.add_new_dir('MoLtimestepping', 'codes', 'moltimestepping')

def _create_home(dirtree: Dirtree):
    dirtree.create_node('moltimestepping')

def _setup_home(dirtree: Dirtree) -> None:
    _add_home(dirtree)
    _create_home(dirtree)

def _build_factory() -> CprototypeFactory:
    cparmater_factory = CparameterFactory()
    cprototype_factory = CprototypeFactory(cparmater_factory)
    return cprototype_factory

def _create_rhs_string(cprototype_factory: CprototypeFactory):
    ricci_eval = cprototype_factory.create_ricci_eval()
    rhs_eval = cprototype_factory.create_rhs_eval()
    rhs_string = [ricci_eval, rhs_eval]
    rhs_string = "\n".join([p.as_string() for p in rhs_string])
    return rhs_string

def _create_post_rhs_string(cprototype_factory: CprototypeFactory):
    apply_bcs_curvilinear = cprototype_factory.create_apply_bcs_curvilinear()
    enforce_detgammahat_constraint = cprototype_factory.create_enforce_detgammahat_constraint()
    post_rhs_string = [apply_bcs_curvilinear, enforce_detgammahat_constraint]
    post_rhs_string = "\n".join([p.as_string() for p in post_rhs_string])
    return post_rhs_string

def setup_rktimestepping(dirtree: Dirtree, integration_method: IntegrationMethod):
    _setup_home(dirtree)

    rk_order = integration_method.variant.attrs.order
    rk_method = f"{integration_method.variant}"
    
    cprototype_factory = _build_factory()
    rhs_string = _create_rhs_string(cprototype_factory)
    post_rhs_string = _create_post_rhs_string(cprototype_factory)
    
    home = dirtree.get('moltimestepping')
    MoL.MoL_C_Code_Generation(
        rk_method,
        RHS_string=rhs_string,
        post_RHS_string=post_rhs_string,
        outdir=str(home),
        )
    return
