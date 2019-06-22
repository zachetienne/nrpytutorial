
# functions_and_globals takes a list of functions [function_list] and a list of globals [global_list] and returns a
# dictionary with keys 'function_list' and 'global_list', and respective values [function_list] and [global_list].
# Mainly used for creating value of module dictionary.
# Throws AssertionError if either functionList or globalList is not a list,
# or if any element in either list is not a string.


def functions_and_globals(function_list, global_list):

    assert(isinstance(function_list, list))
    assert(isinstance(global_list,   list))

    for function in function_list:
        assert(isinstance(function, str))
    for glob in global_list:
        assert(isinstance(glob, str))

    return {'function_list': function_list, 'global_list': global_list}