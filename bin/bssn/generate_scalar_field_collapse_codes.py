from __future__ import annotations

from nrpytutorial.dirtree import Dirtree
from nrpytutorial.timestepping_codes import setup_rktimestepping
from nrpytutorial.integration_methods import make_integration_method, IntegrationMethodVariant
from nrpytutorial.ccode_rendering import CfloatPrecision
from nrpytutorial.coordinate_systems import build_spherical_coord_system
from nrpytutorial.numerical_grids import NumericalGrid
from nrpytutorial.pulse import build_gaussian_pulse
from nrpytutorial.foliation_conditions import FoliationCondition, LapseCondition, ShiftCondition
from nrpytutorial.bssn_vacuum_terms import calculate_bssn_vacuum_terms
from nrpytutorial.bssn_vacuum_terms import generate_c_code_bssn_rhs_terms
from nrpytutorial.bssn_vacuum_terms import generate_c_code_ricci
from nrpytutorial.bssn_vacuum_terms import generate_c_code_hamiltonian
from nrpytutorial.bssn_vacuum_terms import generate_c_code_gammadet

import reference_metric as rfm
import deprecated_reference_metric as evil_rfm
import ScalarField.ScalarField_InitialData as sfid
import BSSN.ADM_Numerical_Spherical_or_Cartesian_to_BSSNCurvilinear as AtoBnum
import NRPy_param_funcs as par
import deprecated_NRPy_param_funcs as evil_par
import CurviBoundaryConditions.CurviBoundaryConditions as cbcs
import cmdline_helper as cmd

import os
from pathlib import Path
from dataclasses import dataclass
import time
from string import Template

from typing import TYPE_CHECKING
if TYPE_CHECKING:
    from nrpytutorial.coordinate_systems import CoordSystem
    from nrpytutorial.pulse import PulseParameters


HOME = 'build/bssn/scalar_field_collapse'
INTEGRATION_METHOD = IntegrationMethodVariant.RK4
CFL_FACTOR = 0.5
FD_ORDER = 4
CFLOAT_TYPE = CfloatPrecision.DOUBLE
DOMAIN_SIZE = 32
SPATIAL_DIMENSIONS = 3
GRID_POINTS_COUNT = 30000
SYMMETRY_AXES = '12'
SIMD_HEADERS = 'SIMD/SIMD_intrinsics.h'
PULSE_AMPLITUDE = 0.4
PULSE_CENTER = 0
PULSE_WIDTH = 1
LAPSE_CONDITOIN = LapseCondition.ONE_PLUS_LOG
SHIFT_CONDITION = ShiftCondition.GAMMA_DRIVING
SINH_WIDTH = 0.2
SINHV2_CONST_DR = 0.05
SYMTP_BSCALE = 0.5

def _add_skeleton(dirtree: Dirtree):
    dirtree.add_new_dir('codes', 'root')
    dirtree.add_new_dir('build', 'root')
    dirtree.add_new_dir('bin', 'root')
    dirtree.add_new_dir('output', 'root')


def _create_skeleton(dirtree: Dirtree):
    dirtree.create()


def setup_dirtree(home) -> Dirtree:
    dirtree = Dirtree(home)
    _add_skeleton(dirtree)
    _create_skeleton(dirtree)
    return dirtree


def _create_numerical_grid(coord_system: CoordSystem):
    numerical_grid_builder = NumericalGrid.get_builder()
    numerical_grid_builder.set_spatial_dimensions(SPATIAL_DIMENSIONS)
    numerical_grid_builder.set_coord_system(coord_system)
    numerical_grid_builder.set_points_count(GRID_POINTS_COUNT)
    numerical_grid_builder.set_infinity(coord_system.domain_size * 1.1)
    return numerical_grid_builder.build()


def create_simd_headers(src: Path, dst: Dirtree):
    dst.add_new_dir('SIMD', 'codes', 'simds')
    dst.create_node('simds')

    dst.add_new_file('SIMD_intrinsics.h', 'simds', 'simd_header')
    dst.create_node('simd_header')

    header = dst.get('simd_header')
    header.path.write_text(src.read_text())
    return


def make_find_timestep_header(dirtree: Dirtree):
    dirtree.add_new_file('find_timestep.h', 'codes', 'find_timestep_header')
    dirtree.create_node('find_timestep_header')
    evil_rfm.out_timestep_func_to_file(
        str(dirtree.get('find_timestep_header')))


def make_scalar_field_init_data(dirtree: Dirtree, pulse: PulseParameters,
                                numerical_grid: NumericalGrid):
    dirtree.add_new_file("SFID.txt", 'output', 'sfidtxt')
    dirtree.create_node('sfidtxt')

    outputfilename = dirtree.get('sfidtxt')

    sfid.ScalarField_InitialData(str(outputfilename),
                                 pulse.profile.as_string(),
                                 pulse.amplitude,
                                 pulse.center,
                                 pulse.width,
                                 numerical_grid.points_count.radial,
                                 numerical_grid.infinity.radial)

    ccodesdir = dirtree.get('codes')
    ccodesdir = str(ccodesdir)
    sfid.NRPy_param_funcs_register_C_functions_and_NRPy_basic_defines(
        Ccodesdir=ccodesdir)

def make_foliation_conditions(lapse, shift) -> FoliationCondition:
    foliation_condition_builder = FoliationCondition.get_builder()
    foliation_condition_builder.set_lapse_condition(lapse)
    foliation_condition_builder.set_shift_condition(shift)
    foliation_condition = foliation_condition_builder.build()
    return foliation_condition

if __name__ == '__main__':
    dirtree = setup_dirtree(HOME)

    integration_method = make_integration_method(
        INTEGRATION_METHOD, CFL_FACTOR, FD_ORDER, CFLOAT_TYPE)

    setup_rktimestepping(dirtree, integration_method)

    coord_system = build_spherical_coord_system(DOMAIN_SIZE, SYMMETRY_AXES)
    numerical_grid = _create_numerical_grid(coord_system)

    create_simd_headers(Path(SIMD_HEADERS), dirtree)

    make_find_timestep_header(dirtree)

    pulse = build_gaussian_pulse(PULSE_AMPLITUDE, PULSE_CENTER, PULSE_WIDTH)
    make_scalar_field_init_data(dirtree, pulse, numerical_grid)

    ccodesdir = str(dirtree.get('codes'))
    AtoBnum.Convert_Spherical_or_Cartesian_ADM_to_BSSN_curvilinear(
            coord_system.kind.as_string(),
            "ID_scalarfield_ADM_quantities",
            Ccodesdir=ccodesdir,
            loopopts=""
            )

    foliation_condition = make_foliation_conditions(LAPSE_CONDITOIN, SHIFT_CONDITION)

    beta_u, detg_constraint = calculate_bssn_vacuum_terms(dirtree, numerical_grid, foliation_condition)

    generate_c_code_bssn_rhs_terms(dirtree, numerical_grid, beta_u)
    generate_c_code_ricci(dirtree)
    generate_c_code_hamiltonian(dirtree)
    generate_c_code_gammadet(dirtree, detg_constraint)

    dirtree.add_new_file('free_parameters.h', 'codes', 'free_parameters_header')
    dirtree.create_node('free_parameters_header')
    rfm.out_default_free_parameters_for_rfm(str(dirtree.get('free_parameters_header')),
                                            numerical_grid.coord_system.domain_size,
                                            SINH_WIDTH,
                                            SINHV2_CONST_DR,
                                            SYMTP_BSCALE)

    evil_rfm.set_Nxx_dxx_invdx_params__and__xx_h(str(dirtree.get('codes')))

    dirtree.add_new_file("xx_to_Cart.h", 'codes', 'xx_to_cart_header')
    dirtree.create_node('xx_to_cart_header')
    evil_rfm.xx_to_Cart_h("xx_to_Cart",
                          "./set_Cparameters.h",
                          str(dirtree.get('xx_to_cart_header'))
                          )
    evil_par.generate_Cparameters_Ccodes(str(dirtree.get('codes')))

    dirtree.add_new_dir("boundary_conditions", 'codes', 'boundary_conditions')
    dirtree.create_node('boundary_conditions')
    cbcs.Set_up_CurviBoundaryConditions(
            str(dirtree.get('boundary_conditions')),
            Cparamspath=os.path.join("../"),
            path_prefix='/home/dimitar/projects/nrpytutorial/Deprecated'
            )

    playground_header = 'ScalarFieldCollapse_Playground_REAL__NGHOSTS__CFL_FACTOR.h'
    dirtree.add_new_file(playground_header, 'codes', 'playground_header')
    dirtree.create_node('playground_header')

    playground_header_src = Path('data/templates/ScalarFieldCollapse_Playground_REAL__NGHOSTS__CFL_FACTOR.h')
    playground_header_template = Template(playground_header_src.read_text())
    dirtree.get('playground_header').path.write_text(
            playground_header_template.substitute({
                'nghosts': int(integration_method.fd_order/2)+1,
                'cfloat_precision': integration_method.cfloat_type.as_string(),
                'cfl_factor': integration_method.cfl_factor
                })
            )

    dirtree.add_new_file('ScalarField_interp.h', 'codes', 'scalar_field_interp_header')
    dirtree.create_node('scalar_field_interp_header')
    dirtree.get('scalar_field_interp_header').path.write_text(
            Path('in_progress-NRPyCritCol/ScalarField/ScalarField_interp.h').read_text()
            )

    dirtree.add_new_file('ScalarFieldCollapse_Playground.c', 'codes', 'playground')
    dirtree.create_node('playground')
    dirtree.get('playground').path.write_text(
            Path('data/templates/ScalarFieldCollapse_Playground.c').read_text()
            )

    dirtree.add_new_file('ScalarFieldCollapse_Playground', 'bin', 'mainexec')
    cmd.C_compile( str(dirtree.get('playground')),
                   str(dirtree.get('mainexec')),
                  compile_mode="optimized"
                  )
