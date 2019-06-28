# isFirstTime takes in a module dictionary [mod_dict] and a trusted values dict [trusted_values_dict],
# and determines if it is the first time the code is being run based off the existence of trusted values in
# [trusted_values_dict] for each module in [mod_dict].
# Requires: The name of the trusted values dictionary is the same as the convention set by create_trusted_globals_dict

# Called by run_test


def is_first_time(mod_dict, trusted_values_dict):

    boolean_dict = dict()

    for mod in mod_dict:
        # Boolean stating whether or not the module has an entry with the proper name in trustedValuesDict
        boolean_dict[mod] = ((mod + 'Globals') not in trusted_values_dict)

    return boolean_dict
