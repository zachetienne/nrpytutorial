# As documented in the NRPy+ tutorial module
#   Tutorial-Coutput__Parameter_Interface.ipynb
#   this core NRPy+ module is used for
#   generating C code and functions.

# Authors: Zachariah B. Etienne; zachetie **at** gmail **dot* com
#          Ken Sible; ksible **at** outlook **dot* com

# Step 0.a: Define __all__, which is the complete
#           list of symbols imported when
#           "from outputC import *" is called.
__all__ = ['lhrh', 'outCparams', 'nrpyAbs', 'superfast_uniq', 'check_if_string__error_if_not',
           'outputC', 'parse_outCparams_string',
           'outC_NRPy_basic_defines_h_dict',
           'outC_function_prototype_dict', 'outC_function_dict', 'Cfunction', 'add_to_Cfunction_dict', 'outCfunction']

import loop as lp                             # NRPy+: C code loop interface
import NRPy_param_funcs as par                # NRPy+: parameter interface
from SIMD import expr_convert_to_SIMD_intrins # NRPy+: SymPy expression => SIMD intrinsics interface
from cse_helpers import cse_preprocess,cse_postprocess  # NRPy+: CSE preprocessing and postprocessing
import sympy as sp                            # SymPy: The Python computer algebra package upon which NRPy+ depends
import re, sys, os, stat                      # Standard Python: regular expressions, system, and multiplatform OS funcs
from collections import namedtuple            # Standard Python: Enable namedtuple data type
from suffixes import dosubs
from grid import find_gftype
from var_access import var_from_access

lhrh = namedtuple('lhrh', 'lhs rhs')
outCparams = namedtuple('outCparams', 'preindent includebraces declareoutputvars outCfileaccess outCverbose CSE_enable CSE_varprefix CSE_sorting CSE_preprocess enable_SIMD SIMD_find_more_subs SIMD_find_more_FMAsFMSs SIMD_debug enable_TYPE')

# Sometimes SymPy has problems evaluating complicated expressions involving absolute
#    values, resulting in hangs. So instead of using sp.Abs(), if we instead use
#    nrpyAbs, we can sidestep the internal SymPy evaluation and force the C
#    codegen to output our desired fabs().
nrpyAbs = sp.Function('nrpyAbs')
custom_functions_for_SymPy_ccode = {
    "nrpyAbs": "fabs",
    'Pow': [(lambda b, e: e == 0.5, lambda b, e: 'sqrt(%s)'     % (b)),
            (lambda b, e: e ==-0.5, lambda b, e: '(1.0/sqrt(%s))'     % (b)),
            (lambda b, e: e == sp.S.One/3, lambda b, e: 'cbrt(%s)' % (b)),
            (lambda b, e: e ==-sp.S.One/3, lambda b, e: '(1.0/cbrt(%s))' % (b)),
            (lambda b, e: e == 2, lambda b, e: '((%s)*(%s))'                % (b,b)),
            (lambda b, e: e == 3, lambda b, e: '((%s)*(%s)*(%s))'           % (b,b,b)),
            (lambda b, e: e == 4, lambda b, e: '((%s)*(%s)*(%s)*(%s))'      % (b,b,b,b)),
            (lambda b, e: e == 5, lambda b, e: '((%s)*(%s)*(%s)*(%s)*(%s))' % (b,b,b,b,b)),
            (lambda b, e: e ==-1, lambda b, e: '(1.0/(%s))'                       % (b)),
            (lambda b, e: e ==-2, lambda b, e: '(1.0/((%s)*(%s)))'                % (b,b)),
            (lambda b, e: e ==-3, lambda b, e: '(1.0/((%s)*(%s)*(%s)))'           % (b,b,b)),
            (lambda b, e: e ==-4, lambda b, e: '(1.0/((%s)*(%s)*(%s)*(%s)))'      % (b,b,b,b)),
            (lambda b, e: e ==-5, lambda b, e: '(1.0/((%s)*(%s)*(%s)*(%s)*(%s)))' % (b,b,b,b,b)),
            (lambda b, e: e !=-5, 'pow')]
##    (lambda b, e: e != 2, 'pow')]
}


# Parameter initialization is called once, within nrpy.py.
par.initialize_param(par.glb_param("char", __name__, "PRECISION", "double")) # __name__ = "outputC", this module's name.
# par.initialize_param(par.glb_param("bool", thismodule, "enable_SIMD", False))

# super fast 'uniq' function:
# f8() function from https://www.peterbe.com/plog/uniqifiers-benchmark
def superfast_uniq(seq): # Author: Dave Kirby
    # Order preserving
    seen = set()
    return [x for x in seq if x not in seen and not seen.add(x)]


def indent_Ccode(Ccode, indent="  "):
    Ccodesplit = Ccode.splitlines()
    outstring = ""
    for i in range(len(Ccodesplit)):
        if Ccodesplit[i] != "":
            if Ccodesplit[i].lstrip().startswith("#"):
                # Remove all indentation from preprocessor statements (lines that start with "#")
                outstring += Ccodesplit[i].lstrip() + '\n'
            else:
                outstring += indent + Ccodesplit[i] + '\n'
        else:
            outstring += '\n'
    return outstring.rstrip(" ")  # make sure to remove trailing whitespace!


def check_if_string__error_if_not(allegedstring, stringdesc):
    if sys.version_info[0] == 3:
        string_types = str
    else:
        string_types = basestring
    if not isinstance(allegedstring, string_types):
        print("ERROR: "+str(stringdesc)+" =="+str(allegedstring)+" not a string!")
        sys.exit(1)


def ccode_postproc(string):
    PRECISION = par.parval_from_str("PRECISION")

    # In the C math library, e.g., pow(x,y) assumes x and y are doubles, and returns a double.
    #  If x and y are floats, then for consistency should use powf(x,y) instead.
    #  Similarly, in the case of x and y being long doubles, should use powl(x,y) for consistency.
    # First we find the appropriate suffix depending on the desired precision:
    cmathsuffix = ""
    if PRECISION == "double":
        pass
    elif PRECISION in ('CCTK_REAL', 'CCTK_REALVEC'):
        pass
    elif PRECISION == "long double":
        cmathsuffix = "l"
    elif PRECISION == "float":
        cmathsuffix = "f"
    else:
        print("Error: "+__name__+"::PRECISION = \""+ PRECISION +"\" not supported")
        sys.exit(1)
    # ... then we append the above suffix to standard C math library functions:
    for func in ['pow', 'sqrt', 'cbrt', 'sin', 'cos', 'tan', 'sinh', 'cosh', 'tanh', 'exp', 'log', 'fabs', 'fmin', 'fmax']:
        string2 = re.sub(func+r'\(', func + cmathsuffix+"(", string); string = string2

    # Finally, SymPy prefers to output Rationals as long-double fractions.
    #  E.g., Rational(1,3) is output as 1.0L/3.0L.
    #  The Intel compiler vectorizer complains miserably about this,
    #  and strictly speaking it is useless when we're in double precision.
    # So here we get rid of the "L" suffix on floating point numbers:
    if PRECISION!="long double":
        string2 = re.sub(r'([0-9.]+)L/([0-9.]+)L', '(\\1 / \\2)', string); string = string2

    return string

def parse_outCparams_string(params):
    # Default values:
    preindent = ""
    includebraces = "True"
    declareoutputvars = "False"
    outCfileaccess = "w"
    outCverbose = "True"
    CSE_enable = "True"
    CSE_sorting = "canonical"
    CSE_varprefix = "tmp"
    CSE_preprocess = "False"
    enable_SIMD = "False"
    SIMD_find_more_subs = "False"
    SIMD_find_more_FMAsFMSs = "True" # Finding too many FMAs/FMSs can degrade performance; currently tuned to optimize BSSN
    SIMD_debug = "False"
    enable_TYPE = "True"

    if params != "":
        params2 = re.sub("^,","",params)
        params = params2.strip()
        split_string = re.split("=|,", params)

        if len(split_string) % 2 != 0:
            print("outputC: Invalid params string: "+params)
            sys.exit(1)

        parnm = []
        value = []
        for i in range(int(len(split_string)/2)):
            parnm.append(split_string[2*i])
            value.append(split_string[2*i+1])

        for i, parname in enumerate(parnm):
            # Clean the string
            if value[i] == "true":
                value[i] = "True"
            if value[i] == "false":
                value[i] = "False"
            if parname == "preindent":
                if not value[i].isdigit():
                    print("Error: preindent must be set to an integer (corresponding to the number of tab stops). ")
                    print(value[i]+" is not an integer.")
                    sys.exit(1)
                preindent = ""
                for _j in range(int(value[i])):  # _j is unused
                    preindent += "  "
            elif parname == "includebraces":
                includebraces = value[i]
            elif parname == "declareoutputvars":
                declareoutputvars = value[i]
            elif parname == "outCfileaccess":
                outCfileaccess = value[i]
            elif parname == "outCverbose":
                outCverbose = value[i]
            elif parname == "CSE_enable":
                CSE_enable = value[i]
            elif parname == "CSE_varprefix":
                CSE_varprefix = value[i]
            elif parname == "CSE_sorting":
                CSE_sorting = value[i]
            elif parname == "CSE_preprocess":
                CSE_preprocess = value[i]
            elif parname == "enable_SIMD":
                enable_SIMD = value[i]
            elif parname == "SIMD_find_more_subs":
                SIMD_find_more_subs = value[i]
            elif parname == "SIMD_find_more_FMAsFMSs":
                SIMD_find_more_FMAsFMSs = value[i]
            elif parname == "SIMD_debug":
                SIMD_debug = value[i]
            elif parname == "enable_TYPE":
                enable_TYPE = value[i]
            elif parname == "GoldenKernelsEnable" and value[i] == "True":
                # GoldenKernelsEnable==True enables the most optimized kernels,
                #   at the expense of ~3x longer codegen runtimes.
                CSE_preprocess          = "True"
                SIMD_find_more_subs     = "True"
                SIMD_find_more_FMAsFMSs = "True"
            elif parname == "GoldenKernelsEnable" and value[i] == "False":
                pass # Do nothing; just allow user to set GoldenKernelsEnable="False".
            else:
                print("Error: outputC parameter name \""+parname+"\" unrecognized.")
                sys.exit(1)

            # CSE preprocessing does not work with SymPy < 1.3. Error out if this is chosen.
            if CSE_preprocess == "True":
                sympy_version = sp.__version__.replace('rc', '...').replace('b', '...')
                sympy_major_version = int(sympy_version.split(".")[0])
                sympy_minor_version = int(sympy_version.split(".")[1])
                if sympy_major_version < 1 or (sympy_major_version == 1 and sympy_minor_version < 3):
                    print('Warning: SymPy version', sympy_version, 'does not support CSE preprocessing. Disabling...')
                    # print('         Please update your SymPy version, or disable CSE preprocessing/GoldenKernels.')
                    CSE_preprocess = "False"

    return outCparams(preindent,includebraces,declareoutputvars,outCfileaccess,outCverbose,
                      CSE_enable,CSE_varprefix,CSE_sorting,CSE_preprocess,
                      enable_SIMD,SIMD_find_more_subs,SIMD_find_more_FMAsFMSs,SIMD_debug,
                      enable_TYPE)

# Input: sympyexpr = a single SymPy expression *or* a list of SymPy expressions
#        output_varname_str = a single output variable name *or* a list of output
#                             variable names, one per sympyexpr.
# Output: C code, as a string.
def outputC(sympyexpr, output_varname_str, filename = "stdout", params = "", prestring = "", poststring = ""):
    outCparams = parse_outCparams_string(params)
    preindent = outCparams.preindent
    TYPE = par.parval_from_str("PRECISION")

    if outCparams.enable_TYPE == "False":
        TYPE = ""

    # Step 0: Initialize
    #  commentblock: comment block containing the input SymPy string,
    #                set only if outCverbose==True
    #  outstring:    the output C code string
    commentblock = ""
    outstring = ""

    # Step 1: If enable_SIMD==True, then check if TYPE=="double". If not, error out.
    #         Otherwise set TYPE="REAL_SIMD_ARRAY", which should be #define'd
    #         within the C code. For example for AVX-256, the C code should have
    #         #define REAL_SIMD_ARRAY __m256d
    if outCparams.enable_SIMD == "True":
        if TYPE not in ('double', ''):
            print("SIMD output currently only supports double precision or typeless. Sorry!")
            sys.exit(1)
        if TYPE == "double":
            TYPE = "REAL_SIMD_ARRAY"

    # Step 2a: Apply sanity checks when either sympyexpr or
    #          output_varname_str is a list.
    if isinstance(output_varname_str, list) and not isinstance(sympyexpr, list):
        print("Error: Provided a list of output variable names, but only one SymPy expression.")
        sys.exit(1)
    if isinstance(sympyexpr, list):
        if not isinstance(output_varname_str, list):
            print("Error: Provided a list of SymPy expressions, but no corresponding list of output variable names")
            sys.exit(1)
        elif len(output_varname_str) != len(sympyexpr):
            print("Error: Length of SymPy expressions list ("+str(len(sympyexpr))+
                  ") != Length of corresponding output variable name list ("+str(len(output_varname_str))+")")
            sys.exit(1)
    # Step 2b: If sympyexpr and output_varname_str are not lists,
    #          convert them to lists of one element each, to
    #          simplify proceeding code.
    if not isinstance(output_varname_str, list) and not isinstance(sympyexpr, list):
        output_varname_strtmp = [output_varname_str]
        output_varname_str = output_varname_strtmp
        sympyexprtmp = [sympyexpr]
        sympyexpr = sympyexprtmp
    sympyexpr = sympyexpr[:]  # pass-by-value (copy list)

    # Step 3: If outCparams.verbose = True, then output the original SymPy
    #         expression(s) in code comments prior to actual C code
    if outCparams.outCverbose == "True":
        commentblock += preindent+"/*\n"+preindent+" *  Original SymPy expression"
        if len(output_varname_str)>1:
            commentblock += "s"
        commentblock += ":\n"
        for i, varname in enumerate(output_varname_str):
            if i == 0:
                if len(output_varname_str) != 1:
                    commentblock += preindent+" *  \"["
                else:
                    commentblock += preindent+" *  \""
            else:
                commentblock += preindent+" *    "
            commentblock += varname + " = " + str(sympyexpr[i])
            if i == len(output_varname_str)-1:
                if len(output_varname_str) != 1:
                    commentblock += "]\"\n"
                else:
                    commentblock += "\"\n"
            else:
                commentblock += ",\n"
        commentblock += preindent+" */\n"

    # Step 4: Add proper indentation of C code:
    if outCparams.includebraces == "True":
        indent = outCparams.preindent+"  "
    else:
        indent = outCparams.preindent+""

    # Step 5: Should the output variable, e.g., outvar, be declared?
    #         If so, start output line with e.g., "double outvar "
    outtypestring = ""
    if outCparams.declareoutputvars == "True":
        outtypestring = indent+TYPE + " "
    else:
        outtypestring = indent

    # Step 6a: If common subexpression elimination (CSE) disabled, then
    #         just output the SymPy string in the most boring way,
    #         nearly consistent with SymPy's ccode() function,
    #         though with support for float & long double types
    #         as well.
    SIMD_RATIONAL_decls = RATIONAL_decls = ""

    if outCparams.CSE_enable == "False":
        # If CSE is disabled:
        # Synthesizing `muladd` calls seems to slow down the code
        # sympyexpr = list(map(map_synthesize_muladd, sympyexpr))
        for i in range(len(sympyexpr)):
            outstring += outtypestring + ccode_postproc(sp.ccode(dosubs(sympyexpr[i]), output_varname_str[i],
                                                                 user_functions=custom_functions_for_SymPy_ccode))+"\n"
    # Step 6b: If CSE enabled, then perform CSE using SymPy and then
    #          resulting C code.
    else:
        # If CSE is enabled:
        SIMD_const_varnms = []
        SIMD_const_values = []

        varprefix = '' if outCparams.CSE_varprefix == 'tmp' else outCparams.CSE_varprefix
        if outCparams.CSE_preprocess == "True" or outCparams.enable_SIMD == "True":
            # If CSE_preprocess == True, then perform partial factorization
            # If enable_SIMD == True, then declare _NegativeOne_ in preprocessing
            factor_negative = eval(outCparams.enable_SIMD) and eval(outCparams.SIMD_find_more_subs)
            sympyexpr, map_sym_to_rat = cse_preprocess(sympyexpr, prefix=varprefix,
                declare=eval(outCparams.enable_SIMD), negative=factor_negative, factor=eval(outCparams.CSE_preprocess))
            for v in map_sym_to_rat:
                p, q = float(map_sym_to_rat[v].p), float(map_sym_to_rat[v].q)
                if outCparams.enable_SIMD == "False":
                    RATIONAL_decls += indent + 'const double ' + str(v) + ' = '
                    # Since Integer is a subclass of Rational in SymPy, we need only check whether
                    # the denominator q = 1 to determine if a rational is an integer.
                    if q != 1: RATIONAL_decls += str(p) + '/' + str(q) + ';\n'
                    else:      RATIONAL_decls += str(p) + ';\n'

        #####
        # Prior to the introduction of the SCALAR_TMP type, NRPy+
        # did not insert sympy.Eq objects into the list of functions
        # that were passed to the CSE. SCALAR_TMP places certain
        # ordering requirements on the system which can only be
        # untangled if sympy.Eq objects are used.
        ###

        # This set of symbols does include the SCALAR_TMPs
        # But includes them as sq.Eq objects.
        sympyexpr_group1 = []
        names_group1 = []

        # This set of symbols does not include the SCALAR_TMPs
        sympyexpr_group2 = []
        names_group2 = []

        # Synthesizing `muladd` calls seems to slow down the code
        # sympyexpr = list(map(map_synthesize_muladd, sympyexpr))
        for ii in range(len(sympyexpr)):
            nm = output_varname_str[ii]
            var = var_from_access(nm)
            typ = find_gftype(var,die=False)
            if typ == "SCALAR_TMP":
                sympyexpr_group1 += [ sp.Eq(sp.sympify(var), sympyexpr[ii]) ]
                names_group1 += [nm]
            else:
                sympyexpr_group1 += [sympyexpr[ii]]
                names_group1 += [nm]
                sympyexpr_group2 += [sympyexpr[ii]]
                names_group2 += [nm]

        sympyexpr_group = sympyexpr_group1
        names_group = names_group1
        # Start processing a group
        sympy_version = sp.__version__.replace('rc', '...').replace('b', '...')
        sympy_major_version = int(sympy_version.split(".")[0])
        sympy_minor_version = int(sympy_version.split(".")[1])
        if sympy_major_version < 1 or (sympy_major_version == 1 and sympy_minor_version < 3):
            print('Warning: SymPy version', sympy_version, 'does not support CSE postprocessing.')
            CSE_results = sp.cse(sympyexpr_group, sp.numbered_symbols(outCparams.CSE_varprefix),
                                 order=outCparams.CSE_sorting)
        else:
            CSE_tmp = sp.cse(sympyexpr_group, sp.numbered_symbols(outCparams.CSE_varprefix),
                                                 order=outCparams.CSE_sorting)
            CSE_results = cse_postprocess(CSE_tmp)

        # cse_postprocess moves SCALAR_TMP out of the group
        sympyexpr_group = sympyexpr_group2
        names_group = names_group2

        for commonsubexpression in CSE_results[0]:
            FULLTYPESTRING = "const " + TYPE + " "
            if outCparams.enable_TYPE == "False":
                FULLTYPESTRING = ""

            if outCparams.enable_SIMD == "True":
                outstring += indent + FULLTYPESTRING + str(commonsubexpression[0]) + " = " + \
                             str(expr_convert_to_SIMD_intrins(commonsubexpression[1],map_sym_to_rat,varprefix,outCparams.SIMD_find_more_FMAsFMSs)) + ";\n"
            else:
                outstring += indent + FULLTYPESTRING + ccode_postproc(sp.ccode(dosubs(commonsubexpression[1]), commonsubexpression[0],
                                                                user_functions=custom_functions_for_SymPy_ccode)) + "\n"

        for i, result in enumerate(CSE_results[1]):
            if outCparams.enable_SIMD == "True":
                outstring += outtypestring + names_group[i] + " = " + \
                             str(expr_convert_to_SIMD_intrins(result,map_sym_to_rat,varprefix,outCparams.SIMD_find_more_FMAsFMSs)) + ";\n"
            else:
                result = dosubs(result)
                outstring += outtypestring+ccode_postproc(sp.ccode(result,names_group[i],
                                                                   user_functions=custom_functions_for_SymPy_ccode))+"\n"
        # Finish processing a group

        # Complication: SIMD functions require numerical constants to be stored in SIMD arrays
        # Resolution: This function extends lists "SIMD_const_varnms" and "SIMD_const_values",
        #             which store the name of each constant SIMD array (e.g., _Integer_1) and
        #             the value of each variable (e.g., 1.0).
        if outCparams.enable_SIMD == "True":
            for v in map_sym_to_rat:
                p, q = float(map_sym_to_rat[v].p), float(map_sym_to_rat[v].q)
                SIMD_const_varnms.extend([str(v)])
                if q != 1: SIMD_const_values.extend([str(p) + '/' + str(q)])
                else:      SIMD_const_values.extend([str(p)])

        # Step 6b.i: If enable_SIMD == True , and
        #            there is at least one SIMD const variable,
        #            then declare the SIMD_const_varnms and SIMD_const_values arrays
        if outCparams.enable_SIMD == "True" and len(SIMD_const_varnms) != 0:
            # Step 6a) Sort the list of definitions. Idea from:
            # https://stackoverflow.com/questions/9764298/is-it-possible-to-sort-two-listswhich-reference-each-other-in-the-exact-same-w
            SIMD_const_varnms, SIMD_const_values = \
                (list(t) for t in zip(*sorted(zip(SIMD_const_varnms, SIMD_const_values))))
            # Step 6b) Remove duplicates
            uniq_varnms = superfast_uniq(SIMD_const_varnms)
            uniq_values = superfast_uniq(SIMD_const_values)
            SIMD_const_varnms = uniq_varnms
            SIMD_const_values = uniq_values
            if len(SIMD_const_varnms) != len(SIMD_const_values):
                print("Error: SIMD constant declaration arrays SIMD_const_varnms[] and SIMD_const_values[] have inconsistent sizes!")
                sys.exit(1)

            for i in range(len(SIMD_const_varnms)):
                if outCparams.enable_TYPE == "False":
                    SIMD_RATIONAL_decls += indent + SIMD_const_varnms[i] + " = " + SIMD_const_values[i]+";"
                else:
                    SIMD_RATIONAL_decls += indent + "const double " + "tmp" + SIMD_const_varnms[i] + " = " + SIMD_const_values[i] + ";\n"
                    SIMD_RATIONAL_decls += indent + "const REAL_SIMD_ARRAY " + SIMD_const_varnms[i] + " = ConstSIMD(" + "tmp" + SIMD_const_varnms[i] + ");\n"
                SIMD_RATIONAL_decls += "\n"

    # Step 7: Construct final output string
    final_Ccode_output_str = commentblock
    # Step 7a: Output C code in indented curly brackets if
    #          outCparams.includebraces = True
    if outCparams.includebraces == "True": final_Ccode_output_str += outCparams.preindent+"{\n"
    final_Ccode_output_str += prestring + RATIONAL_decls + SIMD_RATIONAL_decls + outstring + poststring
    if outCparams.includebraces == "True": final_Ccode_output_str += outCparams.preindent+"}\n"

    # Step 8: If filename == "stdout", then output
    #         C code to standard out (useful for copy-paste or interactive
    #         mode). Otherwise output to file specified in variable name.
    if filename == "stdout":
        # Output to standard out (stdout; "the screen")
        print(final_Ccode_output_str)
    elif filename == "returnstring":
        return final_Ccode_output_str
    else:
        # Output to the file specified by the function input parameter string 'filename':
        with open(filename, outCparams.outCfileaccess) as file:
            file.write(final_Ccode_output_str)
        successstr = ""
        if outCparams.outCfileaccess == "a":
            successstr = "Appended "
        elif outCparams.outCfileaccess == "w":
            successstr = "Wrote "
        print(successstr + "to file \"" + filename + "\"")

from defines_dict import outC_NRPy_basic_defines_h_dict

outC_function_prototype_dict = {}
outC_function_dict           = {}
outC_function_outdir_dict    = {}
outC_function_element = namedtuple('outC_function_element', 'includes prefunc desc c_type name params preloop body loopopts postloop enableCparameters rel_path_to_Cparams')
outC_function_master_list = []

def Cfunction(includes=None, prefunc="", desc="", c_type="void", name=None, params=None, preloop="", body=None,
              loopopts="", postloop="", enableCparameters=True, rel_path_to_Cparams=os.path.join("./")):
    if name is None or params is None or body is None: # use "is None" instead of "==None", as the former is more correct.
        print("Cfunction() error: strings must be provided for function name, parameters, and body")
        sys.exit(1)
    func_prototype = c_type+" "+name+"("+params+")"

    include_Cparams_str = ""
    if enableCparameters:
        if "enable_SIMD" in loopopts or "SIMD_width" in body:  # If using manual SIMD looping, SIMD_width will appear in body.
            include_Cparams_str = "#include \"" + os.path.join(rel_path_to_Cparams, "set_Cparameters-SIMD.h") + "\"\n"
        else:
            include_Cparams_str = "#include \"" + os.path.join(rel_path_to_Cparams, "set_Cparameters.h") + "\"\n"

    complete_func = ""
    if includes is not None:
        if not isinstance(includes, list):
            print("Error in Cfunction(name="+name+"): includes must be set to a list of strings")
            print("e.g., includes=[\"stdio.h\",\"stdlib.h\"] ;  or None (default)")
            print("Found includes = " + str(includes))
            sys.exit(1)
        for inc in includes:
            if "<" in inc:  # for C++ style code e.g., #include <cstdio>
                complete_func += "#include " + inc + "\n"
            else:
                if "NRPy_basic_defines.h" in inc or "NRPy_function_prototypes.h" in inc or \
                        "SIMD_intrinsics.h" in inc:
                    inc = os.path.join(rel_path_to_Cparams, inc)
                complete_func += "#include \"" + inc + "\"\n"
        # complete_func += "\n"

    if prefunc != "":
        complete_func += prefunc + "\n"

    def indent_Ccode(indent, Ccode):
        Ccodesplit = Ccode.splitlines()
        outstring = ""
        for i in range(len(Ccodesplit)):
            outstring += indent + Ccodesplit[i] + '\n'
        return outstring

    if desc != "":
        complete_func += "/*\n" + indent_Ccode(" * ", desc) + " */\n"
    complete_func += func_prototype + " {\n"+include_Cparams_str+preloop+"\n"+lp.simple_loop(loopopts, body)+postloop+"}\n"

    return func_prototype+";", complete_func


def add_to_Cfunction_dict(includes=None, prefunc="", desc="", c_type="void", name=None, params=None,
                          preloop="", body=None, loopopts="", postloop="",
                          path_from_rootsrcdir_to_this_Cfunc="default", enableCparameters=True,
                          rel_path_to_Cparams=os.path.join("./")):

    namesuffix = ""

    if outC_function_element(includes, prefunc, desc, c_type, name + namesuffix, params,
                             preloop, body, loopopts, postloop,
                             enableCparameters, rel_path_to_Cparams) in outC_function_master_list:
        print("OUCH! Found " + name + namesuffix + " in outC_function_master_list.")

    outC_function_master_list.append(outC_function_element(includes, prefunc, desc, c_type, name + namesuffix, params,
                                                           preloop, body, loopopts, postloop,
                                                           enableCparameters, rel_path_to_Cparams))

    outC_function_outdir_dict[name + namesuffix] = path_from_rootsrcdir_to_this_Cfunc
    # print(name, namesuffix, path_from_rootsrcdir_to_this_Cfunc)
    # print(outC_function_outdir_dict)
    outC_function_prototype_dict[name + namesuffix], outC_function_dict[name + namesuffix] = \
        Cfunction(includes, prefunc, desc, c_type, name + namesuffix, params, preloop, body, loopopts, postloop,
                  enableCparameters, rel_path_to_Cparams)


def outCfunction(outfile="", includes=None, prefunc="", desc="",
                 c_type="void", name=None, params=None, preloop="", body=None, loopopts="", postloop="",
                 enableCparameters=True, rel_path_to_Cparams=os.path.join("./")):
    _ignoreprototype, Cfunc = Cfunction(includes, prefunc, desc, c_type, name, params, preloop, body,
                                        loopopts, postloop, enableCparameters, rel_path_to_Cparams)
    if outfile == "returnstring":
        return Cfunc
    print("Writing file:",outfile,"in dir",os.getcwd())
    with open(outfile, "w") as file:
        file.write(Cfunc)
        print("Output C function "+name+"() to file "+outfile)


def construct_Makefile_from_outC_function_dict(Ccodesrootdir, exec_name, uses_free_parameters_h=False,
                                               compiler_opt_option="fastdebug", addl_CFLAGS=None,
                                               addl_libraries=None, mkdir_Ccodesrootdir=True, use_make=True, CC="gcc",
                                               create_lib=False,  include_dirs=None):
    if not create_lib and "main" not in outC_function_dict:
        print("construct_Makefile_from_outC_function_dict() error: C codes will not compile if main() function not defined!")
        print("    Make sure that the main() function registered to outC_function_dict has name \"main\".")
        sys.exit(1)

    if not os.path.isdir(Ccodesrootdir):
        if not mkdir_Ccodesrootdir:
            print("Error (in construct_Makefile_from_outC_function_dict): Directory \"" + Ccodesrootdir + "\" does not exist.")
            sys.exit(1)
        else:
            import cmdline_helper as cmd
            cmd.mkdir(Ccodesrootdir)

    Makefile_list_of_files = []
    def add_to_Makefile(Ccodesrootdir, path_and_file):
        Makefile_list_of_files.append(path_and_file)
        return os.path.join(Ccodesrootdir, path_and_file)

    list_of_uniq_functions = []
    for item in outC_function_master_list:
        if item.name not in list_of_uniq_functions:
            # Convention: Output all C files ending in _gridN into the gridN/ subdirectory.
            if "__rfm__" in item.name:
                subdir = item.name.split("__rfm__")[-1]
                import cmdline_helper as cmd
                cmd.mkdir(os.path.join(Ccodesrootdir, subdir))
                with open(add_to_Makefile(Ccodesrootdir, os.path.join(subdir, item.name+".c")), "w") as file:
                    file.write(outC_function_dict[item.name])
            elif outC_function_outdir_dict[item.name] != "default":
                subdir = outC_function_outdir_dict[item.name]
                with open(add_to_Makefile(Ccodesrootdir, os.path.join(subdir, item.name+".c")), "w") as file:
                    file.write(outC_function_dict[item.name])
            else:
                with open(add_to_Makefile(Ccodesrootdir, os.path.join(item.name+".c")), "w") as file:
                    file.write(outC_function_dict[item.name])
            list_of_uniq_functions += [item.name]
    CFLAGS      = " -O2 -march=native -g -fopenmp -Wall -Wno-unused-variable"
    DEBUGCFLAGS = " -O2 -g -Wall -Wno-unused-variable -Wno-unknown-pragmas"  # OpenMP requires -fopenmp, and when disabling
                                                                             # -fopenmp, unknown pragma warnings appear.
                                                                             # -Wunknown-pragmas silences these warnings
    FASTCFLAGS  = "  -O2 -march=native -fopenmp -Wall -Wno-unused-variable"
    if CC == "gcc":
        CFLAGS      += " -std=gnu99"
        DEBUGCFLAGS += " -std=gnu99"
        FASTCFLAGS   += " -std=gnu99"
    CHOSEN_CFLAGS = CFLAGS
    if compiler_opt_option == "debug":
        CHOSEN_CFLAGS = DEBUGCFLAGS
    elif compiler_opt_option == "fast":
        CHOSEN_CFLAGS = FASTCFLAGS
    if addl_CFLAGS is not None:
        if not isinstance(addl_CFLAGS, list):
            print("Error: construct_Makefile_from_outC_function_dict(): addl_CFLAGS must be a list!")
            sys.exit(1)
        for FLAG in addl_CFLAGS:
            CHOSEN_CFLAGS += " "+FLAG
            CFLAGS += " " + FLAG
            DEBUGCFLAGS += " " + FLAG
            FASTCFLAGS += " " + FLAG
    obj_dependency_str = ""
    dep_list = []
    compile_list = []
    for c_file in Makefile_list_of_files:
        object_file = c_file.replace(".c", ".o")
        obj_dependency_str += " " + object_file
        addl_headers = ""
        if uses_free_parameters_h:
            if c_file == "main.c":
                addl_headers += " free_parameters.h"
        dep_list.append(object_file + ": " + c_file + addl_headers)
        compile_list.append("\t$(CC) $(CFLAGS) $(INCLUDEDIRS)  -c " + c_file + " -o " + object_file)

    linked_libraries = ""
    if addl_libraries is not None:
        if not isinstance(addl_libraries, list):
            print("Error: construct_Makefile_from_outC_function_dict(): addl_libraries must be a list!")
            sys.exit(1)
        for lib in addl_libraries:
            linked_libraries += " " + lib
    linked_libraries += " -lm"
    if "openmp" in CHOSEN_CFLAGS:
        linked_libraries += " -lgomp"

    if use_make:
        with open(os.path.join(Ccodesrootdir, "Makefile"), "w") as Makefile:
            Makefile.write("""CC     = """ + CC + """
CFLAGS = """ + CHOSEN_CFLAGS + """
#CFLAGS = """ + CFLAGS + """
#CFLAGS = """ + DEBUGCFLAGS + """
#CFLAGS = """ + FASTCFLAGS + "\n")
            include_dirs_str = ""
            if include_dirs is not None:
                if not isinstance(include_dirs, list):
                    print("Error: construct_Makefile_from_outC_function_dict(): include_dirs must be a list!")
                    sys.exit(1)
                for include_dir in include_dirs:
                    include_dirs_str += "-I" + include_dir + " "
                include_dirs_str = include_dirs_str[:-1]
            Makefile.write("INCLUDEDIRS = " + include_dirs_str + "\n")
            Makefile.write("all: " + exec_name + " " + obj_dependency_str + "\n")
            for idx, dep in enumerate(dep_list):
                Makefile.write(dep + "\n")
                Makefile.write(compile_list[idx] + "\n\n")
            Makefile.write(exec_name + ": " + obj_dependency_str + "\n")
            ## LINKER STEP:
            Makefile.write("\t$(CC) " + obj_dependency_str + " -o " + exec_name + linked_libraries + "\n")
            ## MAKE CLEAN:
            Makefile.write("\nclean:\n\trm -f *.o */*.o *~ */*~ ./#* *.txt *.dat *.avi *.png " + exec_name + "\n")
    else:
        with open(os.path.join(Ccodesrootdir, "backup_script_nomake.sh"), "w") as backup:
            for compile_line in compile_list:
                backup.write(compile_line.replace("$(CC)", CC).replace("$(CFLAGS)", CFLAGS).replace("\t", "") + "\n")
            backup.write(CC + " " + CFLAGS + " " + obj_dependency_str + " -o " + exec_name + linked_libraries + "\n")
        os.chmod(os.path.join(Ccodesrootdir, "backup_script_nomake.sh"), stat.S_IRWXU)


def construct_NRPy_basic_defines_h(Ccodesrootdir, enable_SIMD=False, supplemental_dict=None):
    if supplemental_dict is None:
        supplemental_dict = {}
    if not os.path.isdir(Ccodesrootdir):
        print("Error (in construct_NRPy_basic_defines_h): Directory \"" + Ccodesrootdir + "\" does not exist.")
        sys.exit(1)

    def output_key(filestream, key, item):
        filestream.write("\n\n//********************************************\n")
        filestream.write("// Basic definitions for module " + key + ":\n" + item)
        filestream.write("//********************************************\n")

    with open(os.path.join(Ccodesrootdir, "NRPy_basic_defines.h"), "w") as file:
        file.write("""// NRPy+ basic definitions, automatically generated from outC_NRPy_basic_defines_h_dict within outputC,
//    and populated within NRPy+ modules. DO NOT EDIT THIS FILE BY HAND.\n\n""")
        if enable_SIMD:
            file.write("// construct_NRPy_basic_defines_h(...,enable_SIMD=True) was called so we #include SIMD intrinsics:\n")
            file.write("""#include "SIMD/SIMD_intrinsics.h"\n""")
        # The ordering here is based largely on data structure dependencies. E.g., griddata_struct contains bc_struct.
        core_modules_list = ["outputC", "NRPy_param_funcs", "finite_difference", "reference_metric",
                             "CurviBoundaryConditions", "MoL", "interpolate", "grid"]
        for key in core_modules_list:
            if key in outC_NRPy_basic_defines_h_dict:
                output_key(file, key, outC_NRPy_basic_defines_h_dict[key])

        for key in outC_NRPy_basic_defines_h_dict:
            if key not in core_modules_list:
                output_key(file, key, outC_NRPy_basic_defines_h_dict[key])

        for key in supplemental_dict:
            output_key(file, key, supplemental_dict[key])


def construct_NRPy_function_prototypes_h(Ccodesrootdir):
    if not os.path.isdir(Ccodesrootdir):
        print("Error (in construct_NRPy_function_prototypes_h): Directory \"" + Ccodesrootdir + "\" does not exist.")
        sys.exit(1)
    with open(os.path.join(Ccodesrootdir, "NRPy_function_prototypes.h"), "w") as file:
        for key, item in outC_function_prototype_dict.items():
            file.write(item + "\n")


def outputC_register_C_functions_and_NRPy_basic_defines(addl_includes=None):
    # First register C functions needed by outputC

    # Then set up the dictionary entry for outputC in NRPy_basic_defines
    Nbd_str  = r"""
#include "stdio.h"
#include "stdlib.h"
#include "math.h"
#include "string.h" // "string.h Needed for strncmp, etc.
#include "stdint.h" // "stdint.h" Needed for Windows GCC 6.x compatibility, and for int8_t
"""
    if addl_includes is not None:
        if not isinstance(addl_includes, list):
            print("Error: addl_includes must be a list!")
            sys.exit(1)
        for include in addl_includes:
            Nbd_str += "#include \"" + include + "\"\n"
    Nbd_str += """

#ifndef M_PI
#define M_PI 3.141592653589793238462643383279502884L
#endif
#ifndef M_SQRT1_2
#define M_SQRT1_2 0.707106781186547524400844362104849039L
#endif

#ifndef MIN
#define MIN(A, B) ( ((A) < (B)) ? (A) : (B) )
#endif
#ifndef MAX
#define MAX(A, B) ( ((A) > (B)) ? (A) : (B) )
#endif

#ifdef __cplusplus
#define restrict __restrict__
#endif
"""
    Nbd_str += "#define REAL " + par.parval_from_str("outputC::PRECISION") + "\n"
    outC_NRPy_basic_defines_h_dict["outputC"] = Nbd_str


def type_and_parname_from_Cparam(glb_Cparam):
    c_type = glb_Cparam.type
    if "char" in glb_Cparam.type:
        c_type = "char"

    parname = glb_Cparam.parname
    if "char" in glb_Cparam.type:
        if len(glb_Cparam.type.split(" ")) == 2:
            parname += glb_Cparam.type.split()[1]
        else:
            print("NRPy_param_funcs error: char array must specify size in type field.")
            print("                        E.g., char blah[100] -> type = \"char [100]\"")
            print("                        Found this: ", glb_Cparam)
            sys.exit(1)

    return c_type, parname

# Must be done here since outputC imports NRPy_param_funcs
def NRPy_param_funcs_register_C_functions_and_NRPy_basic_defines(directory=os.path.join(".")):
    # Set up the C function for BSSN basis transformations
    includes = ["NRPy_basic_defines.h"]
    desc = "Set Cparameters to default values specified within NRPy+."
    c_type = "void"
    name = "set_Cparameters_to_default"
    params = "paramstruct *restrict params"
    body = ""

    for i in range(len(par.glb_Cparams_list)):
        Cptype, parname = type_and_parname_from_Cparam(par.glb_Cparams_list[i])
        if Cptype != "#define":
            comment = "  // " + par.glb_Cparams_list[i].module + "::" + parname
            c_output = "  params->" + parname
            defaultval = par.glb_Cparams_list[i].defaultval
            if Cptype == "char":
                chararray_name = parname.split("[")[0]
                chararray_size = parname.split("[")[1].replace("]", "")
                c_output = "  snprintf(params->" + chararray_name + \
                           ", " + chararray_size + ", " \
                           "\"" + defaultval + "\");" + comment + "\n"
            elif isinstance(par.glb_Cparams_list[i].defaultval, (bool, int, float)):
                c_output += " = " + str(defaultval).lower() + ";" + comment + "\n"
            else:
                c_output += " = " + str(par.glb_Cparams_list[i].defaultval) + ";" + comment + "\n"
            body += c_output
    add_to_Cfunction_dict(
        includes=includes,
        desc=desc,
        c_type=c_type, name=name, params=params,
        body=body,
        loopopts="",
        enableCparameters=False)


    # Step 4: Generate C code to set C parameter constants
    #         (i.e., all ints != -12345678 and REALs != 1e300);
    #         output to filename "set_Cparameters.h" if enable_SIMD==False
    #         or "set_Cparameters-SIMD.h" if enable_SIMD==True
    # Step 4.a: Output non-SIMD version, set_Cparameters.h
    def gen_set_Cparameters(pointerEnable=True):
        returnstring = ""
        for i in range(len(par.glb_Cparams_list)):
            # C parameter type, parameter name
            Cptype, Cpparname = type_and_parname_from_Cparam(par.glb_Cparams_list[i])
            # For efficiency reasons, set_Cparameters*.h does not set char arrays;
            #   access those from the params struct directly.
            if Cptype != "char":
                pointer = "->"
                if pointerEnable==False:
                    pointer = "."

                if not ((Cptype == "REAL" and par.glb_Cparams_list[i].defaultval == 1e300) or Cptype == "#define" or Cptype == "externally_defined"):
                    comment = "  // " + par.glb_Cparams_list[i].module + "::" + Cpparname
                    Coutput = "const " + Cptype + " " + Cpparname + " = " + "params" + pointer + Cpparname + ";" + comment + "\n"
                    returnstring += Coutput
        return returnstring

    # Next output header files for setting C parameters to current values within functions.
    with open(os.path.join(directory, "set_Cparameters.h"), "w") as file:
        file.write(gen_set_Cparameters(pointerEnable=True))
    with open(os.path.join(directory, "set_Cparameters-nopointer.h"), "w") as file:
        file.write(gen_set_Cparameters(pointerEnable=False))

    # Step 4.b: Output SIMD version, set_Cparameters-SIMD.h
    with open(os.path.join(directory, "set_Cparameters-SIMD.h"), "w") as file:
        for i in range(len(par.glb_Cparams_list)):
            # SIMD does not support char arrays.
            if "char" not in par.glb_Cparams_list[i].type:
                Cptype, Cpparname = type_and_parname_from_Cparam(par.glb_Cparams_list[i])
                comment = "  // " + par.glb_Cparams_list[i].module + "::" + Cpparname
                if Cptype == "REAL" and par.glb_Cparams_list[i].defaultval != 1e300:
                    c_output =  "const REAL            NOSIMD" + Cpparname + " = " + "params->" + Cpparname + ";"+comment+"\n"
                    c_output += "const REAL_SIMD_ARRAY " + Cpparname + " = ConstSIMD(NOSIMD" + Cpparname + ");"+comment+"\n"
                    file.write(c_output)
                elif par.glb_Cparams_list[i].defaultval != 1e300 and Cptype != "#define":
                    c_output = "const "+Cptype+" "+Cpparname + " = " + "params->" + Cpparname + ";"+comment+"\n"
                    file.write(c_output)

    # Set up the dictionary entry for grid in NRPy_basic_defines
    # Generate C code to declare C paramstruct;
    #         output to "declare_Cparameters_struct.h"
    #         We want the elements of this struct to be *sorted*,
    #         to ensure that the struct is consistently ordered
    #         for checkpointing purposes.
    Nbd_str = "typedef struct __paramstruct__ {\n"
    CCodelines = []
    for i in range(len(par.glb_Cparams_list)):
        if par.glb_Cparams_list[i].type != "#define":
            Cptype, Cpparname = type_and_parname_from_Cparam(par.glb_Cparams_list[i])
            comment = "  // " + par.glb_Cparams_list[i].module + "::" + Cpparname
            CCodelines.append("  " + Cptype + " " + Cpparname + ";" + comment + "\n")
    for line in sorted(CCodelines):
        Nbd_str += line
    Nbd_str += "} paramstruct;\n"

    outC_NRPy_basic_defines_h_dict["NRPy_param_funcs"] = Nbd_str

def construct_NRPy_Cfunctions(outDir):
    for item in outC_function_master_list:
        d=item._asdict()
        d["outfile"] = os.path.join(outDir, d["name"]+".cc")
        outCfunction(**d)


# SymPy helpers

def expr_map_pre(f, expr):
    expr = f(expr)
    if hasattr(expr, 'args') and len(expr.args) > 0:
        expr = expr.func(*map(lambda arg: expr_map_pre(f, arg), expr.args))
    return expr


def expr_map_post(f, expr):
    if hasattr(expr, 'args') and len(expr.args) > 0:
        expr = expr.func(*map(lambda arg: expr_map_post(f, arg), expr.args))
    expr = f(expr)
    return expr


neg = sp.Function('neg')


def synthesize_neg(expr):
    if expr.is_Mul and expr.args[0] == -1:
        return neg(sp.Mul(*expr.args[1:]))
    return expr


def replace_neg(expr):
    if expr.is_Function and expr.func == neg:
        assert len(expr.args) == 1
        return sp.Mul(-1, expr.args[0])
    return expr


negone = sp.sympify(-1)

muladd = sp.Function('muladd')
mulsub = sp.Function('mulsub')
negmuladd = sp.Function('negmuladd')
negmulsub = sp.Function('negmulsub')


def muladdop(negmul, negadd):
    if negmul:
        return negmulsub if negadd else negmuladd
    else:
        return mulsub if negadd else muladd


def synthesize_muladd(expr):
    expr1 = expr
    outer_negmul = False
    outer_negadd = False
    if expr1.is_Function and expr1.func == neg:
        expr1 = expr1.args[0]
        outer_negmul = not outer_negmul
        outer_negadd = not outer_negadd
    # TODO: break up long chains of additions into a tree
    if expr1.is_Add:
        add_args = expr.args
        for i, arg in enumerate(add_args):
            negmul = outer_negmul
            negadd = outer_negadd
            if arg.is_Function and arg.func == neg:
                arg = arg.args[0]
                negmul = not negmul
            # TODO: find all the enclosed multiplications
            if arg.is_Mul:
                mul_args = arg.args
                assert len(mul_args) >= 2
                mul_expr1 = mul_args[0]
                mul_expr2 = sp.Mul(*mul_args[1:])

                add_expr = sp.Add(*(add_args[0:i] + add_args[i + 1:]))

                if add_expr.is_Function and add_expr.func == neg:
                    add_expr = add_expr.args[0]
                    negadd = not negadd

                if mul_expr1.is_Function and mul_expr1.func == neg:
                    mul_expr1 = mul_expr1.args[0]
                    negmul = not negmul
                if mul_expr2.is_Function and mul_expr2.func == neg:
                    mul_expr2 = mul_expr2.args[0]
                    negmul = not negmul

                if negmul:
                    op = negmulsub if negadd else negmuladd
                else:
                    op = mulsub if negadd else muladd

                return op(mul_expr1, mul_expr2, synthesize_muladd(add_expr))
    return expr


def map_synthesize_muladd(expr):
    expr = expr_map_post(synthesize_neg, expr)
    expr = expr_map_post(synthesize_muladd, expr)
    expr = expr_map_post(replace_neg, expr)
    return expr
