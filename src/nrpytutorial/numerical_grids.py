from __future__ import annotations

import NRPy_param_funcs as par
from dataclasses import dataclass
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from nrpytutorial.coordinate_systems import CoordSystem


@dataclass
class NumericalGrid:
    spatial_dimensions: int
    coord_system: CoordSystem
    points_count: CoordNameReference
    infinity: CoordNameReference

    @staticmethod
    def get_builder():
        return NumericalGridBuilder()

@dataclass
class NumericalGridBuilder:

    def set_spatial_dimensions(self, spatial_dimension: int):
        par.set_parval_from_str("grid::DIM", spatial_dimension)
        self._spatial_dimensions = spatial_dimension
        return self

    def set_coord_system(self, coord_system: CoordSystem):
        self._coord_system = coord_system
        return self

    def set_points_count(self, *args):
        assign_points_count = self._coord_system.get_coord_name_reference()
        self._points_count = assign_points_count(*args)
        return self

    def set_infinity(self, *args):
        assign_infinity = self._coord_system.get_coord_name_reference()
        self._infinity = assign_infinity(*args)
        return self

    def build(self):
        return NumericalGrid(self._spatial_dimensions,
                             self._coord_system,
                             self._points_count,
                             self._infinity)
