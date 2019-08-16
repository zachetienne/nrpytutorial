from UnitTesting.create_test import create_test


def test_TOV_ADM_T4UUmunu():

    module = 'TOV.TOV_ADM_T4UUmunu'

    module_name = 'TOV_ADM_T4UUmunu'

    function_and_global_dict = {'TOV_ADM_T4UUmunu(True)': ['r', 'th', 'ph', 'gammaSphDD', 'KSphDD', 'alphaSph', 'betaSphU', 'BSphU', 'T4UU']}

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
