# Necessary imports for unit testing framework
import unittest
import logging
from run_test import run_test
from functions_and_globals import functions_and_globals
from RepeatedTimer import RepeatedTimer

# TODO: Change level based on desired amount of output.
# ERROR -> Outputs minimal information -- only when there's an error
# INFO -> Outputs when starting and finishing a module, as well as everything in ERROR
# DEBUG -> Displays all pairs of values being compared, as well as everything in INFO
# NOTSET -> Displays symbolic dictionary for all modules, as well as everything in DEBUG
logging.basicConfig(level=logging.INFO)


# # https://stackoverflow.com/questions/3393612/run-certain-code-every-n-seconds/13151299
# # Creates a threaded timer object that prints to the console every 5 minutes
def print_fun(msg):
    logging.info(msg)


Timer = RepeatedTimer(300, print_fun, "\nPrinting every 5 minutes to prevent timeouts.\n")


# Python unittest class
class TestBSSNGlobals(unittest.TestCase):

    @classmethod
    def tearDownClass(cls):
        Timer.stop()

    # Testing globals for ADM in terms of BSSN module
    def test_ADM_Globals(self):
        # TODO: Import modules to be tested
        # Note: Even though it says the modules are unused, these imports are vital for run_test to work properly.
        # Their information gets passed into run_test through locals()
        import BSSN.ADM_in_terms_of_BSSN as ADM_in_terms_of_BSSN

        # TODO: Create lists of globals to calculate
        ADM_in_terms_of_BSSN_global_list = ['gammaDD', 'gammaDDdD', 'gammaDDdDD', 'gammaUU', 'detgamma',
                                      'GammaUDD', 'KDD', 'KDDdD']

        # TODO: Create Module dictionary based on imported modules, functions to initialize the modules, and globals
        # Note that the name of the modules in mod_dict MUST have the same name as the imported module.
        # Example: If you say 'import My_Modules.Module1 as M1', then mod_dict should have the entry 'M1' as a string.
        mod_dict = {
            'ADM_in_terms_of_BSSN': functions_and_globals(['ADM_in_terms_of_BSSN()'], ADM_in_terms_of_BSSN_global_list)
        }

        # TODO: Call run_test with arguments (self, mod_dict, locals())
        run_test(self, mod_dict, locals())

    # Testing globals for BSSN constraints
    def test_Constraints_Globals(self):
        # TODO: Import modules to be tested
        # Note: Even though it says the modules are unused, these imports are vital for run_test to work properly.
        # Their information gets passed into run_test through locals()
        import BSSN.BSSN_constraints as BSSN_constraints

        # TODO: Create lists of globals to calculate
        constraints_global_list = ['H', 'MU']

        # TODO: Create Module dictionary based on imported modules, functions to initialize the modules, and globals
        # Note that the name of the modules in mod_dict MUST have the same name as the imported module.
        # Example: If you say 'import MyModules.Module1 as M1', then mod_dict should have the entry 'M1' as a string.
        mod_dict = {
            'BSSN_constraints': functions_and_globals(['BSSN_constraints()'], constraints_global_list)
        }

        # TODO: Call run_test with arguments (self, mod_dict, locals())
        run_test(self, mod_dict, locals())

    # Testing globals for BSSN exact modules
    def test_Exact_Globals(self):
        # TODO: Import modules to be tested
        # Note: Even though it says the modules are unused, these imports are vital for run_test to work properly.
        # Their information gets passed into run_test through locals()
        import BSSN.BrillLindquist as BrillLindquist
        import BSSN.ShiftedKerrSchild as ShiftedKerrSchild
        import BSSN.StaticTrumpet as StaticTrumpet
        import BSSN.UIUCBlackHole as UIUCBlackHole

        # TODO: Create lists of globals to calculate
        cart_global_list = ['alphaCart', 'betaCartU', 'BCartU', 'gammaCartDD', 'KCartDD']
        sph_global_list = ['alphaSph', 'betaSphU', 'BSphU', 'gammaSphDD', 'KSphDD']

        # TODO: Create Module dictionary based on imported modules, functions to initialize the modules, and globals
        # IMPORTANT: The name of the modules in mod_dict MUST have the same name as the imported module.
        # Example: If you say 'import MyModules.Module1 as M1', then mod_dict should have the entry 'M1',not 'Module1'.
        mod_dict = {
            'BrillLindquist': functions_and_globals(['BrillLindquist(ComputeADMGlobalsOnly = True)'], cart_global_list),

            'ShiftedKerrSchild': functions_and_globals(['ShiftedKerrSchild(ComputeADMGlobalsOnly = True)'],
                                                     sph_global_list),

            'StaticTrumpet': functions_and_globals(['StaticTrumpet(ComputeADMGlobalsOnly = True)'], sph_global_list),

            'UIUCBlackHole': functions_and_globals(['UIUCBlackHole(ComputeADMGlobalsOnly = True)'], sph_global_list)
        }

        # TODO: Call run_test with arguments (self, mod_dict, locals())
        run_test(self, mod_dict, locals())

    # Testing globals for BSSN Psi4 Globals
    def test_Psi4_Globals(self):
        # TODO: Import modules to be tested
        # Note: Even though it says the modules are unused, these imports are vital for run_test to work properly.
        # Their information gets passed into run_test through locals()
        import BSSN.Psi4 as Psi4
        import BSSN.Psi4_tetrads as Psi4Tetrads

        # TODO: Create lists of globals to calculate
        psi4_global_list = ['psi4_re_pt', 'psi4_im_pt']
        psi4_tetrads_global_list = ['l4U', 'n4U', 'mre4U', 'mim4U']

        # TODO: Create Module dictionary based on imported modules, functions to initialize the modules, and globals
        # Note that the name of the modules in mod_dict MUST have the same name as the imported module.
        # Example: If you say 'import MyModules.Module1 as M1', then mod_dict should have the entry 'M1' as a string.
        mod_dict = {
            'Psi4': functions_and_globals(['Psi4(specify_tetrad=False)'], psi4_global_list),

            'Psi4Tetrads': functions_and_globals(['Psi4_tetrads()'], psi4_tetrads_global_list)
        }

        # TODO: Call run_test with arguments (self, mod_dict, locals())
        run_test(self, mod_dict, locals())

    # Testing globals for BSSN quantities
    def test_Quantities_Globals(self):
        # TODO: Import modules to be tested
        # Note: Even though it says the modules are unused, these imports are vital for run_test to work properly.
        # Their information gets passed into run_test through locals()
        import BSSN.BSSN_quantities as BSSN_quantities

        # TODO: Modules that need to be imported to pre-initialize module and their function calls
        import reference_metric as rfm
        rfm.reference_metric()

        # TODO: Create lists of globals to calculate
        quantities_function_list = ['declare_BSSN_gridfunctions_if_not_declared_already()', 'BSSN_basic_tensors()',
                 'gammabar__inverse_and_derivs()', 'detgammabar_and_derivs()', 'AbarUU_AbarUD_trAbar_AbarDD_dD()',
                 'RicciBar__gammabarDD_dHatD__DGammaUDD__DGammaU()', 'betaU_derivs()', 'phi_and_derivs()']

        quantities_global_list = ['hDD', 'aDD', 'lambdaU', 'vetU', 'betU', 'trK', 'cf', 'alpha', 'gammabarDD', 'AbarDD',
                                'LambdabarU', 'betaU', 'BU', 'gammabarUU', 'gammabarDD_dD', 'gammabarDD_dupD',
                                'gammabarDD_dDD', 'GammabarUDD', 'detgammabar', 'detgammabar_dD', 'detgammabar_dDD',
                                'AbarUU', 'AbarUD', 'trAbar', 'AbarDD_dD', 'AbarDD_dupD', 'RbarDD', 'DGammaUDD',
                                'gammabarDD_dHatD', 'DGammaU', 'betaU_dD', 'betaU_dupD', 'betaU_dDD', 'phi_dD',
                                'phi_dupD', 'phi_dDD', 'exp_m4phi', 'phi_dBarD', 'phi_dBarDD']

        # TODO: Create Module dictionary based on imported modules, functions to initialize the modules, and globals
        # Note that the name of the modules in mod_dict MUST have the same name as the imported module.
        # Example: If you say 'import MyModules.Module1 as M1', then mod_dict should have the entry 'M1' as a string.
        mod_dict = {
            'BSSN_quantities': functions_and_globals(quantities_function_list, quantities_global_list)
        }

        # TODO: Call run_test with arguments (self, mod_dict, locals())
        run_test(self, mod_dict, locals())

    # Testing globals for BSSN RHS
    def test_RHS_Globals(self):
        # TODO: Import modules to be tested
        # Note: Even though it says the modules are unused, these imports are vital for run_test to work properly.
        # Their information gets passed into run_test through locals()
        import BSSN.BSSN_RHSs as RHS
        import BSSN.BSSN_gauge_RHSs as gaugeRHS

        # TODO: Create lists of globals to calculate
        RHS_global_list = ['cf_rhs', 'trK_rhs', 'lambda_rhsU', 'a_rhsDD', 'h_rhsDD']
        gauge_RHS_global_list = ['alpha_rhs', 'bet_rhsU', 'vet_rhsU']

        # TODO: Create Module dictionary based on imported modules, functions to initialize the modules, and globals
        # Note that the name of the modules in mod_dict MUST have the same name as the imported module.
        # Example: If you say 'import MyModules.Module1 as M1', then mod_dict should have the entry 'M1' as a string.
        mod_dict = {
            'RHS': functions_and_globals(['BSSN_RHSs()'], RHS_global_list),

            'gaugeRHS': functions_and_globals(['BSSN_gauge_RHSs()'], gauge_RHS_global_list),
        }

        # TODO: Call run_test with arguments (self, mod_dict, locals())
        run_test(self, mod_dict, locals())

    # Testing globals for BSSN T4UUmunu_vars
    def test_T4UU_Globals(self):
        # TODO: Import modules to be tested
        # Note: Even though it says the modules are unused, these imports are vital for run_test to work properly.
        # Their information gets passed into run_test through locals()
        import BSSN.BSSN_T4UUmunu_vars as BSSN_T4UUmunu_vars

        # TODO: Create lists of globals to calculate
        T4UU_global_list = ['rho', 'S', 'sD', 'sDD']

        # TODO: Create Module dictionary based on imported modules, functions to initialize the modules, and globals
        # Note that the name of the modules in mod_dict MUST have the same name as the imported module.
        # Example: If you say 'import MyModules.Module1 as M1', then mod_dict should have the entry 'M1' as a string.
        mod_dict = {
            'BSSN_T4UUmunu_vars': functions_and_globals(['define_BSSN_T4UUmunu_rescaled_source_terms()'],
                                                        T4UU_global_list)
        }

        # TODO: Call run_test with arguments (self, mod_dict, locals())
        run_test(self, mod_dict, locals())


# Necessary for unittest class to work properly
if __name__ == '__main__':
    unittest.main()
