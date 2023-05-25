from __future__ import annotations

import grid as gri
import indexedexp as ixp
import finite_difference as fin
import reference_metric as rfm
import NRPy_param_funcs as par
import BSSN.BSSN_RHSs as rhs
import BSSN.BSSN_gauge_RHSs as gaugerhs
import BSSN.BSSN_quantities as Bq
import ScalarField.ScalarField_RHSs as sfrhs
import ScalarField.ScalarField_Tmunu as sfTmunu
import BSSN.BSSN_stress_energy_source_terms as Bsest
import BSSN.Enforce_Detgammahat_Constraint as EGC
import BSSN.BSSN_constraints as bssncon
from outputC import lhrh,outputC,outCfunction


from typing import TYPE_CHECKING
if TYPE_CHECKING:
    from nrpytutorial.dirtree import Dirtree
    from nrpytutorial.numerical_grids import NumericalGrid
    from nrpytutorial.foliation_conditions import FoliationCondition


def _prepare_rfm_infra_home(dirtree: Dirtree):
    dirtree.add_new_dir('rfm_files', 'codes', 'rfm_files')
    dirtree.create_node('rfm_files')


def _evaluate_gauge_bssn_rhs():
    rhs.BSSN_RHSs()


def _evaluate_scalar_field_rhs():
    sfrhs.ScalarField_RHSs()


def _compute_contravariant_stress_enrgy_tensor(
        numerical_grid: NumericalGrid) -> StressEnergyTensorDto:
    sfTmunu.ScalarField_Tmunu()
    stress_energy_tensor = sfTmunu.T4UU

    Bsest.BSSN_source_terms_for_BSSN_RHSs(stress_energy_tensor)
    rhs.trK_rhs += Bsest.sourceterm_trK_rhs
    for i in range(numerical_grid.spatial_dimensions):
        # Needed for Gamma-driving shift RHSs:
        rhs.Lambdabar_rhsU[i] += Bsest.sourceterm_Lambdabar_rhsU[i]
        # Needed for BSSN RHSs:
        rhs.lambda_rhsU[i] += Bsest.sourceterm_lambda_rhsU[i]
        for j in range(numerical_grid.spatial_dimensions):
            rhs.a_rhsDD[i][j] += Bsest.sourceterm_a_rhsDD[i][j]

    gaugerhs.BSSN_gauge_RHSs()
    # We use betaU as our upwinding control vector:
    Bq.BSSN_basic_tensors()
    betaU = Bq.betaU
    enforce_detg_constraint_symb_expressions = EGC.Enforce_Detgammahat_Constraint_symb_expressions()

    return stress_energy_tensor, betaU


def _compute_ricci_tensor():
    Bq.RicciBar__gammabarDD_dHatD__DGammaUDD__DGammaU()


def _register_hamiltonian_as_gridfunction():
    H = gri.register_gridfunctions("AUX", "H")


def _define_hamiltonian_constraint_and_output_ccode(
        stress_energy_tensor: StressEnergyTensorDto):
    bssncon.BSSN_constraints(add_T4UUmunu_source_terms=False)
    Bsest.BSSN_source_terms_for_BSSN_constraints(stress_energy_tensor)
    bssncon.H += Bsest.sourceterm_H


def _add_kreiss_oliger_dissipation(numerical_grid: NumericalGrid,
                                   foliation_condition: FoliationCondition
                                   ):
    diss_strength = par.Cparameters(
        "REAL", "ScalarFieldCollapse", ["diss_strength"], 0.1)

    alpha_dKOD = ixp.declarerank1("alpha_dKOD")
    cf_dKOD = ixp.declarerank1("cf_dKOD")
    trK_dKOD = ixp.declarerank1("trK_dKOD")
    sf_dKOD = ixp.declarerank1("sf_dKOD")
    sfM_dKOD = ixp.declarerank1("sfM_dKOD")
    betU_dKOD = ixp.declarerank2("betU_dKOD", "nosym")
    vetU_dKOD = ixp.declarerank2("vetU_dKOD", "nosym")
    lambdaU_dKOD = ixp.declarerank2("lambdaU_dKOD", "nosym")
    aDD_dKOD = ixp.declarerank3("aDD_dKOD", "sym01")
    hDD_dKOD = ixp.declarerank3("hDD_dKOD", "sym01")

    for k in range(numerical_grid.spatial_dimensions):
        gaugerhs.alpha_rhs += diss_strength * alpha_dKOD[k] * rfm.ReU[k]
        rhs.cf_rhs += diss_strength * cf_dKOD[k] * rfm.ReU[k]
        rhs.trK_rhs += diss_strength * trK_dKOD[k] * rfm.ReU[k]
        sfrhs.sf_rhs += diss_strength * sf_dKOD[k] * rfm.ReU[k]
        sfrhs.sfM_rhs += diss_strength * sfM_dKOD[k] * rfm.ReU[k]
        for i in range(numerical_grid.spatial_dimensions):
            if "2ndOrder" in foliation_condition.shift.as_string():
                gaugerhs.bet_rhsU[i] += diss_strength * betU_dKOD[i][k] * rfm.ReU[k]
            gaugerhs.vet_rhsU[i] += diss_strength * vetU_dKOD[i][k] * rfm.ReU[k]
            rhs.lambda_rhsU[i] += diss_strength * lambdaU_dKOD[i][k] * rfm.ReU[k]
            for j in range(numerical_grid.spatial_dimensions):
                rhs.a_rhsDD[i][j] += diss_strength * aDD_dKOD[i][j][k] * rfm.ReU[k]
                rhs.h_rhsDD[i][j] += diss_strength * hDD_dKOD[i][j][k] * rfm.ReU[k]


def _evaluate_rhs_terms(
        numerical_grid: NumericalGrid,
        foliation_condition: FoliationCondition):
    par.set_parval_from_str("BSSN.BSSN_quantities::LeaveRicciSymbolic", "True")
    _evaluate_gauge_bssn_rhs()
    _evaluate_scalar_field_rhs()
    stress_energy_tensor, beta_u = _compute_contravariant_stress_enrgy_tensor(
        numerical_grid)
    par.set_parval_from_str(
        "BSSN.BSSN_quantities::LeaveRicciSymbolic", "False")
    _compute_ricci_tensor()
    _register_hamiltonian_as_gridfunction()
    _define_hamiltonian_constraint_and_output_ccode(stress_energy_tensor)
    _add_kreiss_oliger_dissipation(numerical_grid, foliation_condition)

    return beta_u


def calculate_bssn_vacuum_terms(
        dirtree: Dirtree,
        numerical_grid: NumericalGrid,
        foliation_condition: FoliationCondition) -> BetauDto:
    _prepare_rfm_infra_home(dirtree)

    par.set_parval_from_str("reference_metric::enable_rfm_precompute", "True")

    par.set_parval_from_str("reference_metric::rfm_precompute_Ccode_outdir",
                            str(dirtree.get('rfm_files')))
    beta_u = _evaluate_rhs_terms(numerical_grid, foliation_condition)

    par.set_parval_from_str("reference_metric::enable_rfm_precompute","False")
    rfm.ref_metric__hatted_quantities()

    return beta_u

def generate_c_code_bssn_rhs_terms(dirtree: Dirtree,
                                   numerical_grid: NumericalGrid,
                                   beta_u: BetauDto):
    # Construct the left-hand sides and right-hand-side expressions for all BSSN RHSs
    lhs_names = [        "alpha",       "cf",       "trK",         "sf",         "sfM"   ]
    rhs_exprs = [gaugerhs.alpha_rhs, rhs.cf_rhs, rhs.trK_rhs, sfrhs.sf_rhs, sfrhs.sfM_rhs]
    for i in range(numerical_grid.spatial_dimensions):
        lhs_names.append(        "betU"+str(i))
        rhs_exprs.append(gaugerhs.bet_rhsU[i])
        lhs_names.append(   "lambdaU"+str(i))
        rhs_exprs.append(rhs.lambda_rhsU[i])
        lhs_names.append(        "vetU"+str(i))
        rhs_exprs.append(gaugerhs.vet_rhsU[i])
        for j in range(i,numerical_grid.spatial_dimensions):
            lhs_names.append(   "aDD"+str(i)+str(j))
            rhs_exprs.append(rhs.a_rhsDD[i][j])
            lhs_names.append(   "hDD"+str(i)+str(j))
            rhs_exprs.append(rhs.h_rhsDD[i][j])

    # Sort the lhss list alphabetically, and rhss to match.
    #   This ensures the RHSs are evaluated in the same order
    #   they're allocated in memory:
    lhs_names,rhs_exprs = [list(x) for x in zip(*sorted(zip(lhs_names,rhs_exprs), key=lambda pair: pair[0]))]

    # Declare the list of lhrh's
    BSSN_evol_rhss = []
    for var in range(len(lhs_names)):
        BSSN_evol_rhss.append(lhrh(lhs=gri.gfaccess("rhs_gfs",lhs_names[var]),rhs=rhs_exprs[var]))

    # Set up the C function for the BSSN RHSs
    desc="Evaluate the BSSN RHSs"
    name="rhs_eval"

    dirtree.add_new_file(f"{name}.h", 'codes', f'{name}_header')
    dirtree.create_node(f'{name}_header')
    outCfunction(
            outfile  = str(dirtree.get(f'{name}_header')), desc=desc, name=name,
            params   = """rfm_struct *restrict rfmstruct,const paramstruct *restrict params,
            const REAL *restrict auxevol_gfs,const REAL *restrict in_gfs,REAL *restrict rhs_gfs""",
            body     = fin.FD_outputC("returnstring",BSSN_evol_rhss, params="outCverbose=False,enable_SIMD=True",
                                      upwindcontrolvec=beta_u),
            loopopts = "InteriorPoints,enable_SIMD,enable_rfm_precompute")
