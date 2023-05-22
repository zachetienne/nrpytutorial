from nrpytutorial.dirtree import Dirtree
from nrpytutorial.timestepping_codes import setup_rktimestepping
from nrpytutorial.integration_methods import make_integration_method, IntegrationMethodVariant
from nrpytutorial.ccode_rendering import CfloatPrecision


HOME = 'build/bssn/scalar_field_collapse'
INTEGRATION_METHOD = IntegrationMethodVariant.RK4
CFL_FACTOR = 0.5
FD_ORDER = 4
CFLOAT_TYPE = CfloatPrecision.DOUBLE


def setup_dirtree(home) -> Dirtree:
    dirtree = Dirtree(home)
    dirtree.add_new_dir('codes', 'root')
    dirtree.add_new_dir('build', 'root')
    dirtree.add_new_dir('bin', 'root')
    dirtree.add_new_dir('output', 'root')
    dirtree.create()
    return dirtree


if __name__ == '__main__':
    dirtree = setup_dirtree(HOME)
    integration_method = make_integration_method(INTEGRATION_METHOD, CFL_FACTOR,
                                                 FD_ORDER, CFLOAT_TYPE)
    setup_rktimestepping(dirtree, integration_method)
