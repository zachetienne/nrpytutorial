from UnitTesting.create_test import create_test


def test_equations():

    module = 'GRFFE.equations'

    module_name = 'GRFFE'

    function_and_global_dict = {'generate_everything_for_UnitTesting()': ['B_notildeU','smallb4U','smallbsquared','TEM4UU','TEM4UD','S_tildeD','S_tilde_fluxUD','S_tilde_source_termD']}

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
            sys.exit(1)

    else:
        globals()[sys.argv[4]]()
