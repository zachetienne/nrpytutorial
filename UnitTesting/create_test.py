import sys
import cmdline_helper as cmd
import logging
import os


def create_test(module, module_name, function_and_global_dict, logging_level='INFO', initialization_string_dict=None):

    # If more than 3 arguments were passed in, then the 4th argument is the logging level.
    if len(sys.argv) > 3:
        logging_level = sys.argv[3]

    # If initialization_string_dict wasn't passed in, set it to an empty dict
    if initialization_string_dict is None:
        initialization_string_dict = {}

    # Making sure all arguments have the correct type
    assert type(module_name) == str, "module_name is " + str(type(module_name)) + " -- it should be type str. Please check input for module_name."
    assert type(module) == str, "module is " + str(type(module)) + " -- it should be type str. Please check input for module in " + module_name
    assert type(function_and_global_dict) == dict, "function_and_global_dict is " + str(type(function_and_global_dict)) + " -- it should be type dict. Please check input for function_and_global_dict in " + module_name
    assert len(function_and_global_dict) != 0, "function_and_global_dict is empty -- it must contain at least one entry to be tested. Please check input for function_and_global_dict in " + module_name
    assert type(logging_level) == str, "logging_level is " + str(type(logging_level)) + " -- it should be type str. Please check input for logging_level in " + module_name
    assert type(initialization_string_dict) == dict, "initialization_string is " + str(type(initialization_string_dict)) + " -- it should be type dict. Please check input for initialization_string_dict in " + module_name
    assert all(name in function_and_global_dict for name in initialization_string_dict), "initialization_string_dict contains key not in function_and_global_dict. Please check input for initialization_string_dict in " + module_name

    logging.basicConfig(level=logging_level)

    # Assume the tests haven't failed
    failed = False

    # Loop through each function in function_and_global_dict
    for function, global_list in function_and_global_dict.items():

        initialization_string = initialization_string_dict.get(function, '')

        # Making sure all arguments for this function have the correct type
        assert type(function) == str, "function in function_and_global_dict is" + str(type(function)) + " -- it should be type str. Please check input for function_and_global_dict in " + module_name
        assert type(global_list) == list, "global_list for function " + function + 'is ' + str(type(global_list)) + " -- it should be type list. Please check input for function_and_global_dict in " + module_name
        assert len(global_list) != 0, "global_list for function " + function + "is empty -- it must contain at least one entry to be tested. Please check input for function " +function + " in " + module_name
        assert all(isinstance(glob, str) for glob in global_list), "global_list in function_and_global_dict contains at least one entry that's not a str. Please check input for function " + function + " in " + module_name

        # Create file_string with the unittest format
        file_string = '''
import unittest
import logging
from UnitTesting.setup_trusted_values_dict import setup_trusted_values_dict
from UnitTesting.RepeatedTimer import RepeatedTimer

logging.basicConfig(level=logging.{})  

def timer_output():
    logging.info(' Printing every 5 minutes to prevent Travis CI timeout.')


class TestGlobals(unittest.TestCase):

    def setUp(self):
        self.path = r'{}'
        # Create trusted_values_dict.py if it doesn't exist
        logging.debug(' Calling setup_trusted_values_dict...')
        setup_trusted_values_dict(self.path)
        logging.debug(' ...Success: setup_trusted_values_dict ran without errors.\\n')
        self.timer = RepeatedTimer(300, timer_output)
    
    
    def tearDown(self):
        self.timer.stop()
        

    def test_globals(self):

        self.module = '{}'

        self.module_name = '{}'

        self.function = '{}'
        
        self.global_list = {}
        
        self.initialization_string = """{}"""
        
        self.trusted_values_dict_name = '{}globals'
        
        try:
        
{}
{}
        # Something failed
        except AssertionError:
            pass
        # Nothing failed
        else:
            import os
            file = open(os.path.join(self.path, 'success.txt'), 'w')
            file.close()
                
                
if __name__ == '__main__':
    unittest.main()
'''

        # Create trusted_values_dict_name based on module_name and function
        trusted_values_dict_name = module_name + '__' + \
                function.replace('()', '__').replace('(', '__').replace(')', '__').replace(' ', '').replace('=', '_')

        # Copying the lines from run_test.py into our test file
        with open(os.path.join('UnitTesting', 'run_test.py'), 'r') as file:
            run_test_string = file.read().split('def run_test(self):')

        # Properly formatting run_test's imports
        run_test_imports_list = run_test_string[0].split('\n')
        final_run_test_imports_string = ''
        for line in run_test_imports_list:
            if line != '':
                final_run_test_imports_string += '            ' + line + '\n'

        # Properly formatting run_test's body
        run_test_body_list = run_test_string[1].split('\n')
        final_run_test_body_string = ''
        for line in run_test_body_list:
            final_run_test_body_string += '        ' + line + '\n'

        # Formatting file_string with inputs
        file_string = file_string.format(logging_level.upper(), sys.path[0], module, module_name, function,
                                         global_list, initialization_string, trusted_values_dict_name,
                                         final_run_test_imports_string, final_run_test_body_string)

        logging.debug(' Test file for:\nmodule:   ' + module_name + '\nfunction: ' + function + '\n' + file_string)

        # Setting full+path by joining the path of the user's test file, the trusted_values_dict_name, and 'test.py'
        full_path = os.path.join(sys.path[0], trusted_values_dict_name + 'test.py')

        # If the path has a space in it, tell the user there is a space and error out
        if ' ' in full_path:
            logging.error(' ' + full_path)
            logging.error(' It seems there is a space in one of your directories listed in the path above. '
                          'Please change this.')

            raise SystemExit

        # Otherwise, the path is correct -- open it and write file_string to it
        with open(full_path, 'w') as file:
            logging.info(' Creating file ' + full_path + '...')
            file.write(file_string)
            logging.info(' ...Success: File created.\n')

        logging.debug(' $PYTHONEXEC is: ' + str(sys.argv[1]) + '\n')

        logging.debug(' Calling Execute_input_string\n')

        # Use cmdline_helper to execute the test
        # cmd.Execute_input_string('coverage run -a' + ' ' + full_path, output=False)
        cmd.Execute_input_string(sys.argv[1] + ' ' + full_path, output=False)

        # See if success.txt was created -- if so, the test passed. Otherwise, the test failed
        try:
            success_file = os.path.join(sys.path[0], 'success.txt')
            open(success_file)
            cmd.delete_existing_files(success_file)
        except IOError:
            logging.error(' Test for function ' + function + ' in module ' + module_name +
                          ' failed! Please examine test file.\n'
                          '\n----------------------------------------------------------------------\n'
                          '----------------------------------------------------------------------\n')
            failed = True
        else:
            logging.info(' Test for function ' + function + ' in module ' + module_name +
                         ' passed! Deleting test file...')
            cmd.delete_existing_files(full_path)
            logging.info(' ...Deletion successful. Test complete.\n'
                         '\n----------------------------------------------------------------------\n'
                         '----------------------------------------------------------------------\n')

    if failed:
        raise SystemExit
