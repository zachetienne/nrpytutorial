from trusted_values_dict import trusted_values_dict


# createTrustedGlobalsDict takes in a module dictionary [mod_dict] and a boolean list [first_times].
# For each module, if [first_time] is True, then an empty dict is returned.
# This ensures that when the dictionary is passed into [calc_error], there will be an error.
# If [first_time] is False, then a dictionary that contains every module in ModDict as keys, and each module's
# respective dictionary from trustedValuesDict as values. The naming convention for the dictionaries is as follows:
#   trustedValuesDict['(MODULE_NAME)Globals'] -- The module name with 'Globals' concatenated on the end.
#   This is consistent throughout all files.
# Throws an [AssertionError] if [mod_dict] and [first_times] contain different keys.
# Throws a [KeyError] if [first_time] is true for a module whose entry isn't in [trusted_values_dict].

# Called by run_test


def create_trusted_globals_dict(mod_dict, first_times):

    assert set(mod_dict.keys()) == set(first_times.keys())

    trusted_dict = dict()

    for mod in mod_dict:

        first_time = first_times[mod]

        if first_time:
            trusted_dict[mod] = dict()
        else:
            dict_string = mod + 'Globals'
            trusted_dict[mod] = trusted_values_dict[dict_string]

    return trusted_dict
