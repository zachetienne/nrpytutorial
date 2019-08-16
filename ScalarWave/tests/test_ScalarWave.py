from UnitTesting.create_test import create_test


def test_InitialData_PlaneWave():

    module = 'ScalarWave.InitialData_PlaneWave'

    module_name = 'InitialData_PlaneWave'

    function_and_global_dict = {'InitialData_PlaneWave()': ['uu_ID', 'vv_ID']}

    create_test(module, module_name, function_and_global_dict)


def test_ScalarWave_RHSs():

    module = 'ScalarWave.ScalarWave_RHSs'

    module_name = 'ScalarWave_RHSs'

    function_and_global_dict = {'ScalarWave_RHSs()': ['wavespeed', 'uu_rhs', 'vv_rhs']}

    create_test(module, module_name, function_and_global_dict)


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
            exit(1)

    else:
        globals()[sys.argv[4]]()
