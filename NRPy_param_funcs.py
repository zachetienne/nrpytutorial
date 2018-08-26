glb_params_list = []  # = where we store the parameters and default values of parameters. A list of named tuples
glb_paramsvals_list = []  # = where we store parameter values.
from collections import namedtuple
glb_param = namedtuple('glb_param', 'type module parname defaultval')
import sympy as sp

def initialize_param(input):
    if get_params_idx(input) == -1:
        glb_params_list.append(input)
        glb_paramsvals_list.append(input.defaultval)
    else:
        print("initialize_param() minor warning: Did nothing; already initialized parameter "+input.module+"::"+input.parname)

# Given the named tuple `input` and list of named tuples `params`,
#    defined according to namedtuple('param', 'type module name defaultval'),
#    where in the case of `input`, defaultval need not be set,
#    return the list index of `params` that matches `input`.
# On error returns -1
def get_params_idx(input):
    # inspired by: https://stackoverflow.com/questions/2917372/how-to-search-a-list-of-tuples-in-python:
    list = [i for i, v in enumerate(glb_params_list)
            if (input.type=="ignoretype" or input.type==v[0]) and input.module == v[1] and input.parname == v[2]]
    if list == []:
        return -1 # No match found => error out!
    else:
        if len(list) > 1:
            print("Error: Found multiple parameters matching "+str(input))
            exit(1)
        return list.pop() # pop() returns the index

def get_params_value(input):
    idx = get_params_idx(input)
    if idx < 0:
        print("Error: could not find a parameter matching:",input)
        print("Full list of modules:\n",glb_params_list)
        exit(1)
    else:
        return glb_paramsvals_list[idx]

#
def idx_from_str(varname,modname=""):
    if "::" in varname:
        splitstring = re.split('::', varname)
        modname=splitstring[0]
        varname=splitstring[1]
        
    # inspired by: https://stackoverflow.com/questions/2917372/how-to-search-a-list-of-tuples-in-python:    
    if modname == "":
        list = [i for i, v in enumerate(glb_params_list) if v[2] == varname]
    else:
        list = [i for i, v in enumerate(glb_params_list) if (v[1] == modname and v[2] == varname)]
    if list == []:
        print("Error: Could not find a parameter matching \""+varname+"\" in ",glb_params_list)
        exit(1)
    if len(list) > 1:
        print("Error: Found more than one parameter named \""+varname+"\". Use get_params_value() instead.")
        exit(1)
    return list.pop()

def parval_from_str(string):
    return glb_paramsvals_list[idx_from_str(string)]

def set_parval_from_str(string,value):
    glb_paramsvals_list[idx_from_str(string)] = value

# parse_param_string__set__params_and_paramsvars:
# Summary: This function parses a string like this
# module::variablename=value
# into its component parts: module,variablename,value
# then hunts for module,variablename in the
# params list. When it finds the matching index in the
# list "idx", it sets paramsvals[idx] = value.
# Usage comment: This function parses both parameter file
#   inputs as well as parameter inputs from the command
#   line. You should set filename = "" if reading from the
#   command line. This ensures the error messages are
#   appropriate for the context.
import re
def set_paramsvals_value(line,filename="", FindMainModuleMode=False):
    MainModuleFound = True
    if FindMainModuleMode == True:
        MainModuleFound = False
    # First remove carriage return and leading whitespace:
    stripped_line_of_text = line.strip()
    # Only process lines that do NOT start with a hash
    if not stripped_line_of_text.startswith("#"):
        # Valid lines take the form:
        #    module::variable = value # Comment here
        # Thus, using the delimiters "::", "=", and "#",
        # with the regex split command, the first three
        # items in single_param_def will be [module, variablename, value]
        single_param_def = re.split('::|=|#', stripped_line_of_text)

        # First verify that single_param_def has at least 3 parts:
        if len(single_param_def) < 3:
            if filename != "":
                print("Error: the line " + line + " in parameter file " + filename + " is not in the form")
                print("\"module::variable = value\"")
            else:
                print("Error: the command-line argument " + stripped_line_of_text + " is not in the form")
                print("\"module::variable=value\"   <-- NOTICE NO SPACES ALLOWED!")
            exit(1)

        # Next remove all leading/trailing whitespace from single_param_def
        for i in range(len(single_param_def)):
            single_param_def[i] = single_param_def[i].strip()

        if FindMainModuleMode == False:
            # Next find the parameter in the params list and set
            # the corresponding element in paramsvals to the value.
            idx = get_params_idx(glb_param("ignoretype", single_param_def[0], single_param_def[1], "ignoredefval"))
            # If parameter is not found, print useful error message, then exit:
            if idx == -1:
                if filename != "":
                    print("Error: when reading line \"" + line + "\" in parameter file \"" + filename + "\":")
                else:
                    print("Error: when parsing command-line argument \"" + stripped_line_of_text + "\":")
                print("\t\tcould not find parameter \""+ single_param_def[1] + "\" in \""+single_param_def[0]+"\" module.")
                exit(1)
            # If parameter is found at index idx, set paramsval[idx] to the value specified in the file.
            if glb_params_list[idx].defaultval != "RUNTIME":
                partype = glb_params_list[idx].type
                if partype == "bool":
                    if single_param_def[2] == "True":
                        glb_paramsvals_list[idx] = True
                    elif single_param_def[2] == "False":
                        glb_paramsvals_list[idx] = False
                    else:
                        print("Error: \"bool\" type can only take values of \"True\" or \"False\"")
                        exit(1)
                elif partype == "INT":
                    glb_paramsvals_list[idx] = int(single_param_def[2])
                elif partype is "REAL" or \
                    partype is "REALARRAY" or \
                    partype is "char" or \
                    partype is "char *":
                    glb_paramsvals_list[idx] = single_param_def[2]
                else:
                    print("Error: type \""+partype+"\" on variable \""+ glb_params_list[idx].parname +"\" is unsupported.")
                    print("Supported types include: bool, INT, REAL, REALARRAY, char, and char *")
                    exit(1)
#                    glb_paramsvals_list[idx] = single_param_def[2]
            else:
                print("Error: Tried to set the parameter "
                      + single_param_def[0] + "::" + single_param_def[1] +
                      " with default value RUNTIME")
                print("Such a parameter is defined by NRPy+, but must be "
                      "set at C code runtime. Go fix your C code parameter file!")
                exit(1)
        elif FindMainModuleMode == True and MainModuleFound == False:
            if single_param_def[0] == "NRPy" and single_param_def[1] == "MainModule":
                idx = get_params_idx(glb_param("ignoretype", single_param_def[0], single_param_def[1], "ignoredefval"))
                if idx == -1:
                    print("Critical error: NRPy::MainModule is uninitialized!")
                    exit(1)
                glb_paramsvals_list[idx] = single_param_def[2]

def Cparameters(type,module,names,assumption="Real"):
    output = []
    # if names is not a list, make it a list, to
    #      simplify the remainder of this routine.
    if not isinstance(names,list):
        names = [names]
    for i in range(len(names)):
        initialize_param(glb_param(type, module, names[i], "SetAtCRuntime"))
        if assumption == "Real":
            tmp = sp.Symbol(names[i], real=True) # Assumes all Cparameters are real.
        elif assumption == "RealPositive":
            tmp = sp.Symbol(names[i], real=True, positive=True) # Assumes all Cparameters are real and positive.
        else:
            print("Error: assumption "+str(assumption)+" not supported.")
            exit(1)
        output.append(tmp)
    if len(names) == 1:
        return output[0]
    return output

def Ccode__declare_params(filename):
    Coutput = ""
    for i in range(len(glb_params_list)):
        partype = glb_params_list[i].type
        if partype != "bool" and \
           partype != "char" and \
           partype != "INT" and \
           partype != "REAL":
            print("Error: parameter "+glb_params_list[i].module+"::"+glb_params_list[i].parname+" has unsupported type: \""
                  + glb_params_list[i].type + "\"")
            exit(1)
        if partype == "char":
            Ctype = "char *"
        else:
            Ctype = partype
        Coutput += Ctype + " " + glb_params_list[i].module + "__" + glb_params_list[i].parname
        if isinstance(glb_paramsvals_list[i], (bool,int,float)):
            Coutput += " = " + str(glb_paramsvals_list[i]).lower() + ";\n"
        elif isinstance(glb_paramsvals_list[i], (str)):
            Coutput += " = \"" + str(glb_paramsvals_list[i]).lower() + "\";\n"
        else:
            Coutput += " = " + str(glb_paramsvals_list[i]) + ";\n"
    if filename == "stdout":
        print(Coutput)
        return
    elif filename == "returnstring":
        return Coutput
    else:
        # Output to the file specified by outCfilename
        with open(filename, "w") as file:
            file.write(Coutput)
        print("Wrote to file \"" + filename + "\"")