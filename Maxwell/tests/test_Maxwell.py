from UnitTesting.create_test import create_test


def test_VacuumMaxwell_Flat_Evol_Cartesian_System_I():

    module = 'Maxwell.VacuumMaxwell_Flat_Evol_Cartesian'

    module_name = 'VacuumMaxwell_Flat_Evol_Cartesian_System_I'

    function_and_global_dict = {'VacuumMaxwellRHSs()': ['ArhsU', 'ErhsU', 'C', 'psi_rhs']}

    initialization_string_dict = {'VacuumMaxwellRHSs()': 'import NRPy_param_funcs as par\npar.initialize_param(par.glb_param("char", "Maxwell.InitialData","System_to_use","System_I"))'}

    create_test(module, module_name, function_and_global_dict, initialization_string_dict=initialization_string_dict)

def test_VacuumMaxwell_Flat_Evol_Cartesian_System_II():

    module = 'Maxwell.VacuumMaxwell_Flat_Evol_Cartesian'

    module_name = 'VacuumMaxwell_Flat_Evol_Cartesian_System_II'

    function_and_global_dict = {'VacuumMaxwellRHSs()': ['ArhsU', 'ErhsU', 'C', 'psi_rhs', 'Gamma_rhs', 'G']}

    initialization_string_dict = {'VacuumMaxwellRHSs()': 'import NRPy_param_funcs as par\npar.initialize_param(par.glb_param("char", "Maxwell.InitialData","System_to_use","System_II"))'}

    create_test(module, module_name, function_and_global_dict, initialization_string_dict=initialization_string_dict)


def test_VacuumMaxwell_Flat_Evol_Curvilinear_rescaled():

    module = 'Maxwell.VacuumMaxwell_Flat_Evol_Curvilinear_rescaled'

    module_name = 'VacuumMaxwell_Flat_Evol_Curvilinear_rescaled'

    function_and_global_dict = {'VacuumMaxwellRHSs_rescaled()': ['erhsU', 'arhsU', 'psi_rhs', 'Gamma_rhs', 'C', 'G', 'EU_Cart', 'AU_Cart']}

    initialization_string_dict = {'VacuumMaxwellRHSs_rescaled()':  'import NRPy_param_funcs as par\npar.initialize_param(par.glb_param("char", "Maxwell.InitialData","System_to_use","System_II"))'}

    create_test(module, module_name, function_and_global_dict, initialization_string_dict=initialization_string_dict)



def test_InitialData_I():

    module = 'Maxwell.InitialData'

    module_name = 'InitialData_System_I'

    function_and_global_dict = {'InitialData()': ['AidU', 'EidU', 'psi_ID']}

    initialization_string_dict = {'InitialData()': 'import NRPy_param_funcs as par\npar.initialize_param(par.glb_param("char", "Maxwell.InitialData","System_to_use","System_I"))'}

    create_test(module, module_name, function_and_global_dict, initialization_string_dict=initialization_string_dict)

def test_InitialData_II():

    module = 'Maxwell.InitialData'

    module_name = 'InitialData_System_II'

    function_and_global_dict = {'InitialData()': ['AidU', 'EidU', 'psi_ID', 'Gamma_ID']} # issue here

    initialization_string_dict = {'InitialData()': 'import NRPy_param_funcs as par\npar.set_parval_from_str("Maxwell.InitialData::System_to_use","System_II")'}

    create_test(module, module_name, function_and_global_dict, initialization_string_dict=initialization_string_dict, logging_level='DEBUG')


if __name__ == '__main__':
    import sys

    if len(sys.argv) <= 3:
        failed_functions = []
        for fun in dir():
            if fun[0:5] == 'test_':
                print('\nTesting ' + str(fun) + '...\n')
                try:
                    exec(fun + '()')
                except SystemExit:
                    failed_functions.append(fun)

        if failed_functions != []:
            import sys, os
            with open(os.path.join('UnitTesting', 'failed_tests.txt'), 'a') as file:
                for function in failed_functions:
                    file.write(sys.argv[0] + ': ' + str(function) + '\n')
            sys.exit(1)

    else:
        globals()[sys.argv[4]]()
