import NRPy_param_funcs as par
from dataclasses import dataclass

@dataclass
class NumericalGrid:
    spatial_dimensions: int

    @staticmethod
    def get_builder():
        return NumericalGridBuilder()

@dataclass
class NumericalGridBuilder:

    def set_spatial_dimensions(self, spatial_dimension: int):
        par.set_parval_from_str("grid::DIM", spatial_dimension)
        self._spatial_dimensions = spatial_dimension
        return self

    def build(self):
        return NumericalGrid(self._spatial_dimensions)
