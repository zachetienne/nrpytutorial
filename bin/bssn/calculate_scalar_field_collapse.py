from __future__ import annotations

from nrpytutorial.dirtree import Dirtree
from nrpytutorial.timestepping_codes import setup_rktimestepping
from nrpytutorial.integration_methods import make_integration_method, IntegrationMethodVariant
from nrpytutorial.ccode_rendering import CfloatPrecision
from nrpytutorial.coordinate_systems import build_spherical_coord_system
from nrpytutorial.numerical_grids import NumericalGrid

from typing import TYPE_CHECKING
if TYPE_CHECKING:
    from nrpytutorial.coordinate_systems import CoordSystem


HOME = 'build/bssn/scalar_field_collapse'
INTEGRATION_METHOD = IntegrationMethodVariant.RK4
CFL_FACTOR = 0.5
FD_ORDER = 4
CFLOAT_TYPE = CfloatPrecision.DOUBLE
DOMAIN_SIZE = 32
SPATIAL_DIMENSIONS = 3
GRID_POINTS_COUNT = 30000


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


if __name__ == '__main__':
    dirtree = setup_dirtree(HOME)

    integration_method = make_integration_method(
        INTEGRATION_METHOD, CFL_FACTOR, FD_ORDER, CFLOAT_TYPE)

    setup_rktimestepping(dirtree, integration_method)

    coord_system = build_spherical_coord_system(DOMAIN_SIZE)
    numerical_grid = _create_numerical_grid(coord_system)
