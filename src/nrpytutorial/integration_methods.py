from __future__ import annotations

import enum
from MoLtimestepping.RK_Butcher_Table_Dictionary import Butcher_dict
from collections import namedtuple
from dataclasses import dataclass
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from nrpytutorial.ccode_rendering import CfloatPrecision


RkMethodAttrs = namedtuple('RkMethodAttrs', 'coeffs order')


class IntegrationMethodVariant(enum.Enum):
    EULER = 'Euler'
    RK2HEUN = "RK2 Heun"
    RK2MP = "RK2 MP"
    RK2RALSTON = "RK2 Ralston"
    RK3 = "RK3"
    RK3HEUN = "RK3 Heun"
    RK3RALSTON = "RK3 Ralston"
    SSPRK3 = "SSPRK3"
    RK4 = "RK4"
    DP5 = "DP5"
    DP5ALT = "DP5alt"
    CK5 = "CK5"
    DP6 = "DP6"
    L6 = "L6"
    DP8 = "DP8"

    def __str__(self):
        return self.value

    @property
    def attrs(self):
        coeffs, order = Butcher_dict[self.value]
        return RkMethodAttrs(coeffs, order)

@dataclass
class IntegrationMethod:
    variant: IntegrationMethodVariant
    cfl_factor: float
    fd_order: float
    cfloat_type: CfloatPrecision

def make_integration_method(variant: IntegrationMethodVariant,
                            cfl_factor: float,
                            fd_order: float,
                            cfloat_type: CfloatPrecision
                            ):
    if cfl_factor == 0 or not cfl_factor % 2:
        raise ValueError(f"cfl_factor need be positive even number")

    return IntegrationMethod(variant, cfl_factor, fd_order, cfloat_type)
