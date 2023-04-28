from dataclasses import dataclass
import enum
from typing import List


@dataclass
class Cparameter:
    name: str
    is_value: bool = True

    def as_string(self):
        return self.name if self.is_value else f"&{self.name}"


@dataclass
class Cprototype:
    name: str
    parameters: List[Cparameter]

    def as_string(self):
        parameters = (a.as_string() for a in self.parameters)
        parameters = ", ".join(parameters)
        return f"{self.name}({parameters});"

    def add_parameter(self, param: Cparameter):
        self.parameters.append(param)


class CparameterFactory:
    def create_rfmstruct(self):
        return Cparameter('rfmstruct', is_value=False)

    def create_params(self):
        return Cparameter('params', is_value=False)

    def create_rk_input_gfs(self):
        return Cparameter('RK_INPUT_GFS')

    def create_rk_output_gfs(self):
        return Cparameter('RK_OUTPUT_GFS')

    def create_auxevol_gfs(self):
        return Cparameter('auxevol_gfs')

    def create_bcstruct(self):
        return Cparameter('bcstruct', is_value=False)

    def create_num_evol_gfs(self):
        return Cparameter("NUM_EVOL_GFS")

    def create_evol_gf_parity(self):
        return Cparameter('evol_gf_parity')


@dataclass
class CprototypeFactory:
    parameter_factory: CparameterFactory

    def _create_parameters(self, targets):
        parameters = (f'create_{p}' for p in targets)
        parameters = (getattr(self.parameter_factory, p) for p in parameters)
        parameters = [p() for p in parameters]
        return parameters

    def create_ricci_eval(self):
        parameters = ['rfmstruct', 'params', 'rk_input_gfs', 'auxevol_gfs']
        parameters = self._create_parameters(parameters)
        return Cprototype('Ricci_eval', parameters)

    def create_rhs_eval(self):
        parameters = [
            'rfmstruct',
            'params',
            'auxevol_gfs',
            'rk_input_gfs',
            'rk_output_gfs']
        parameters = self._create_parameters(parameters)
        return Cprototype('rhs_eval', parameters)

    def create_apply_bcs_curvilinear(self):
        parameters = [
            'params',
            'bcstruct',
            'num_evol_gfs',
            'evol_gf_parity',
            'rk_output_gfs']
        parameters = self._create_parameters(parameters)
        return Cprototype('apply_bcs_curvilinear', parameters)

    def create_enforce_detgammahat_constraint(self):
        parameters = ['rfmstruct', 'params', 'rk_output_gfs']
        parameters = self._create_parameters(parameters)
        return Cprototype('enforce_detgammahat_constraint', parameters)


class CfloatPrecision(enum.Enum):
    FLOAT = 'float'
    DOUBLE = 'double'
    LONG_DOUBLE = 'long double'

    def as_string(self):
        return self.value
