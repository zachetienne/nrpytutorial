# needed because global state transfer
import BSSN.BSSN_gauge_RHSs as gaugerhs

import enum
from dataclasses import dataclass
import NRPy_param_funcs as par

class LapseCondition(enum.Enum):
    ONE_PLUS_LOG = "OnePlusLog"

    def as_string(self):
        return self.value


class ShiftCondition(enum.Enum):
    GAMMA_DRIVING = "GammaDriving2ndOrder_Covariant"

    def as_string(self):
        return self.value

@dataclass
class FoliationCondition:
    lapse: LapseCondition
    shift: ShiftCondition

    @staticmethod
    def get_builder():
        return FoliationConditionBuilder()

@dataclass
class FoliationConditionBuilder:

    def set_lapse_condition(self, condition: LapseCondition):
        par.set_parval_from_str("BSSN.BSSN_gauge_RHSs::LapseEvolutionOption",
                                condition.as_string())
        self._lapse = condition
        return self

    def set_shift_condition(self, condition: ShiftCondition):
        par.set_parval_from_str("BSSN.BSSN_gauge_RHSs::ShiftEvolutionOption",
                                condition.as_string()
                                )
        self._shift = condition
        return self

    def build(self):
        return FoliationCondition(self._lapse, self._shift)
