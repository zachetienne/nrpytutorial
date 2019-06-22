
# Takes in a module dictionary [mod_dict] and locals [old_locals] and returns a dictionary [result_dict] containing
# the symbolic expressions for each global specified in mod_dict

# Requires: Modules can't have the same name. Two globals for a given module can't have the same name.
#           Functions must be able to be called on their respective module. Globals must 
#           be defined in their respective modules.
# Returns: Returns [resultDict], which is a dictionary whose keys are the modules that were passed 
#          through modDict and the values are dictionaries containing the values for the specified globals
# Note: Must pass in locals() as second argument to insure that all imports that have been done are accessible
#       by the evaluateGlobals module

# Called by run_test


def evaluate_globals(mod_dict, old_locals):

    result_dict = dict()
    
    for mod,func_glob_dict in mod_dict.items():

        # Initializing string of execution
        stringexec = ''

        # Calling all functions and assigning all globals
        for function in func_glob_dict['function_list']:
            stringexec += mod + '.' + function + '\n'
        for glob in func_glob_dict['global_list']:
            stringexec += glob + '=' + mod + '.' + glob + '\n'

        # Initializing location 
        loc = {}

        # Executing string of execution with current globals and storing resulting globals in loc
        exec(stringexec, old_locals, loc)
        
        # Storing the module-variable pair [mod],[loc] into the dictionary [resultDict]
        result_dict[mod] = loc
        
    return result_dict
