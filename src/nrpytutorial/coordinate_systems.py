import enum
from dataclasses import dataclass
from collections import namedtuple
from typing import Any
from functools import partial


@enum.unique
class CoordSystemKind(enum.Enum):
    SPHERICAL = 'Spherical'
    SINHSPHERICAL = "SinhSpherical"
    SINHSPHERICALV2 = 'SinhSphericalv2'
    CYLINDRICAL = "Cylindrical"
    SINHCYLINDRICAL = "SinhCylindrical"
    SYMTP = "SymTP"
    SINHSYMTP = "SinhSymTP"

    def as_string(self):
        return self.value

def get_params_builder(kind: CoordSystemKind):
    if kind is CoordSystemKind.SPHERICAL:
        return namedtuple('SphericalParams', '')
    raise NotImplementedError(
"""
REVIEW comments to match coordsystemkinds
sinh_width sets the default value for:
  * SinhSpherical's params.SINHW
  * SinhCylindrical's params.SINHW{RHO,Z}
  * SinhSymTP's params.SINHWAA
sinh_width      = 0.2 # If Sinh* coordinates chosen

sinhv2_const_dr sets the default value for:
  * SinhSphericalv2's params.const_dr
  * SinhCylindricalv2's params.const_d{rho,z}
sinhv2_const_dr = 0.05# If Sinh*v2 coordinates chosen

SymTP_bScale sets the default value for:
  * SinhSymTP's params.bScale
SymTP_bScale    = 0.5 # If SymTP chosen
""")
@dataclass
class CoordSystem:
    kind: CoordSystemKind
    domain_size: int
    params: Any

    @staticmethod
    def get_builder():
        return CoordSystemBuilder()

@dataclass
class CoordSystemBuilder:

    def set_kind(self, kind: CoordSystemKind):
        self._kind = kind

    def set_domain_size(self, size: int):
        self._domain_size = size

    def set_free_params(self, *args):
        params_builder = get_params_builder(self._kind)
        self._params = params_builder(*args)

    def build(self):
        return CoordSystem(self._kind, self._domain_size, self._params)

def build_coord_system(kind: CoordSystemKind, domain_size: int, *args) -> CoordSystem:
    builder = CoordSystemBuilder()
    builder.set_kind(kind)
    builder.set_domain_size(domain_size)
    builder.set_free_params()
    return builder.build()

build_spherical_coord_system = partial(build_coord_system, CoordSystemKind.SPHERICAL)
