from UnitTesting.create_test import create_test


def test_your_module():

    module = 'your_module_location.your_module'

    module_name = 'name_of_your_module'

    function_and_global_dict = {'function_call': ['your global', 'your next global', 'your next next global', 'etc']}

    create_test(module, module_name, function_and_global_dict)


# Ignore this -- it's to ensure bash functionality
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
