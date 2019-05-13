import NRPy_param_funcs as par
import re
from SIMD import expr_convert_to_SIMD_intrins

from collections import namedtuple
lhrh = namedtuple('lhrh', 'lhs rhs')
outCparams = namedtuple('outCparams', 'preindent includebraces declareoutputvars outCfileaccess outCverbose CSE_enable CSE_varprefix SIMD_enable SIMD_debug enable_TYPE')

# Parameter initialization is called once, within nrpy.py.
par.initialize_param(par.glb_param("char", __name__, "PRECISION", "double")) # __name__ = "outputC", this module's name.
# par.initialize_param(par.glb_param("bool", thismodule, "SIMD_enable", False))

# super fast 'uniq' function:
# f8() function from https://www.peterbe.com/plog/uniqifiers-benchmark
def superfast_uniq(seq): # Author: Dave Kirby
    # Order preserving
    seen = set()
    return [x for x in seq if x not in seen and not seen.add(x)]

def ccode_postproc(string):
    PRECISION = par.parval_from_str("PRECISION")

    # In the C math library, e.g., pow(x,y) assumes x and y are doubles, and returns a double.
    #  If x and y are floats, then for consistency should use powf(x,y) instead.
    #  Similarly, in the case of x and y being long doubles, should use powl(x,y) for consistency.
    # First we find the appropriate suffix depending on the desired precision:
    cmathsuffix = ""
    if PRECISION == "double":
        pass
    elif PRECISION == "long double":
        cmathsuffix = "l"
    elif PRECISION == "float":
        cmathsuffix = "f"
    else:
        print("Error: "+__name__+"::PRECISION = \""+ PRECISION +"\" not supported")
        exit(1)
    # ... then we append the above suffix to standard C math library functions:
    for func in ['pow', 'sqrt', 'sin', 'cos', 'tan', 'sinh', 'cosh', 'tanh', 'exp', 'log', 'fabs']:
        string2 = re.sub(func+'\(', func + cmathsuffix+"(", string); string = string2

    # Finally, SymPy prefers to output Rationals as long-double fractions.
    #  E.g., Rational(1,3) is output as 1.0L/3.0L.
    #  The Intel compiler vectorizer complains miserably about this,
    #  and strictly speaking it is useless when we're in double precision.
    # So here we get rid of the "L" suffix on floating point numbers:
    if PRECISION!="long double":
        string2 = re.sub('([0-9.]+)L/([0-9.]+)L', '(\\1 / \\2)', string); string = string2

    return string

def parse_outCparams_string(params):
    # Default values:
    preindent = ""
    includebraces = "True"
    declareoutputvars = "False"
    outCfileaccess = "w"
    outCverbose = "True"
    CSE_enable = "True"
    CSE_varprefix = "tmp"
    SIMD_enable = "False"
    SIMD_debug = "False"
    enable_TYPE = "True"

    if params != "":
        params2 = re.sub("^,","",params)
        params = params2.strip()
        splitstring = re.split("=|,", params)

        if len(splitstring) % 2 != 0:
            print("outputC: Invalid params string: "+params)
            exit(1)

        parnm = []
        value = []
        for i in range(int(len(splitstring)/2)):
            parnm.append(splitstring[2*i])
            value.append(splitstring[2*i+1])

        for i in range(len(parnm)):
            # Clean the string
            if value[i] == "true":
                value[i] = "True"
            if value[i] == "false":
                value[i] = "False"
            if parnm[i] == "preindent":
                if not value[i].isdigit():
                    print("Error: preindent must be set to an integer (corresponding to the number of tab stops). ")
                    print(value[i]+" is not an integer.")
                    exit(1)
                preindent = ""
                for i in range(int(value[i])):
                    preindent += "   "
            elif parnm[i] == "includebraces":
                includebraces = value[i]
            elif parnm[i] == "declareoutputvars":
                declareoutputvars = value[i]
            elif parnm[i] == "outCfileaccess":
                outCfileaccess = value[i]
            elif parnm[i] == "outCverbose":
                outCverbose = value[i]
            elif parnm[i] == "CSE_enable":
                CSE_enable = value[i]
            elif parnm[i] == "CSE_varprefix":
                CSE_varprefix = value[i]
            elif parnm[i] == "SIMD_enable":
                SIMD_enable = value[i]
            elif parnm[i] == "SIMD_debug":
                SIMD_debug = value[i]
            elif parnm[i] == "enable_TYPE":
                enable_TYPE = value[i]
            else:
                print("Error: outputC parameter name \""+parnm[i]+"\" unrecognized.")
                exit(1)

    return outCparams(preindent,includebraces,declareoutputvars,outCfileaccess,outCverbose,CSE_enable,CSE_varprefix,SIMD_enable,SIMD_debug,enable_TYPE)

import sympy as sp
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

    # Step 1: If SIMD_enable==True, then check if TYPE=="double". If not, error out.
    #         Otherwise set TYPE="REAL_SIMD_ARRAY", which should be #define'd
    #         within the C code. For example for AVX-256, the C code should have
    #         #define REAL_SIMD_ARRAY __m256d
    if outCparams.SIMD_enable == "True":
        if not (TYPE == "double" or TYPE == ""):
            print("SIMD output currently only supports double precision or typeless. Sorry!")
            exit(1)
        if TYPE == "double":
            TYPE = "REAL_SIMD_ARRAY"
        else:
            TYPE = ""

    # Step 2a: Apply sanity checks when either sympyexpr or
    #          output_varname_str is a list.
    if type(output_varname_str) is list and type(sympyexpr) is not list:
        print("Error: Provided a list of output variable names, but only one SymPy expression.")
        exit(1)
    if type(sympyexpr) is list:
        if type(output_varname_str) is not list:
            print("Error: Provided a list of SymPy expressions, but no corresponding list of output variable names")
            exit(1)
        elif len(output_varname_str) != len(sympyexpr):
            print("Error: Length of SymPy expressions list ("+str(len(sympyexpr))+
                  ") != Length of corresponding output variable name list ("+str(len(output_varname_str))+")")
            exit(1)
    # Step 2b: If sympyexpr and output_varname_str are not lists,
    #          convert them to lists of one element each, to
    #          simplify proceeding code.
    if type(output_varname_str) is not list and type(sympyexpr) is not list:
        output_varname_strtmp = [output_varname_str]
        output_varname_str = output_varname_strtmp
        sympyexprtmp = [sympyexpr]
        sympyexpr = sympyexprtmp


    # Step 3: If outCparams.verbose = True, then output the original SymPy
    #         expression(s) in code comments prior to actual C code
    if outCparams.outCverbose == "True":
        commentblock += preindent+"/*\n"+preindent+" *  Original SymPy expression"
        if len(output_varname_str)>1:
            commentblock += "s"
        commentblock += ":\n"
        for i in range(len(output_varname_str)):
            if i==0:
                if len(output_varname_str)!=1:
                    commentblock += preindent+" *  \"["
                else:
                    commentblock += preindent+" *  \""
            else:
                commentblock += preindent+" *    "
            commentblock += output_varname_str[i] + " = " + str(sympyexpr[i])
            if i==len(output_varname_str)-1:
                if len(output_varname_str)!=1:
                    commentblock += "]\"\n"
                else:
                    commentblock += "\"\n"
            else:
                commentblock += ",\n"
        commentblock += preindent+" */\n"

    # Step 4: Add proper indentation of C code:
    if outCparams.includebraces == "True":
        indent = outCparams.preindent+"   "
    else:
        indent = outCparams.preindent+""

    # Step 5: Should the output variable, e.g., outvar, be declared?
    #         If so, start output line with e.g., "double outvar "
    outtypestring = ""
    if outCparams.declareoutputvars == "True":
        outtypestring = outCparams.preindent+indent+TYPE + " "
    else:
        outtypestring = outCparams.preindent+indent

    # Step 6a: If common subexpression elimination (CSE) disabled, then
    #         just output the SymPy string in the most boring way,
    #         nearly consistent with SymPy's ccode() function,
    #         though with support for float & long double types
    #         as well.
    SIMD_decls = ""

    if outCparams.CSE_enable == "False":
        # If CSE is disabled:
        for i in range(len(sympyexpr)):
            outstring += outtypestring + ccode_postproc(sp.ccode(sympyexpr[i], output_varname_str[i]))+"\n"
    # Step 6b: If CSE enabled, then perform CSE using SymPy and then
    #          resulting C code.
    else:
        # If CSE is enabled:
        SIMD_const_varnms = []
        SIMD_const_values = []

        CSE_results = sp.cse(sympyexpr, sp.numbered_symbols(outCparams.CSE_varprefix), order='canonical')
        for commonsubexpression in CSE_results[0]:
            FULLTYPESTRING = "const " + TYPE + " "
            if outCparams.enable_TYPE == "False":
                FULLTYPESTRING = ""

            if outCparams.SIMD_enable == "True":
                outstring += outCparams.preindent + indent + FULLTYPESTRING + str(commonsubexpression[0]) + " = " + \
                             str(expr_convert_to_SIMD_intrins(commonsubexpression[1],SIMD_const_varnms,SIMD_const_values,outCparams.SIMD_debug)) + ";\n"
            else:
                outstring += outCparams.preindent+indent+FULLTYPESTRING+ccode_postproc(sp.ccode(commonsubexpression[1],commonsubexpression[0]))+"\n"
        for i,result in enumerate(CSE_results[1]):
            if outCparams.SIMD_enable == "True":
                outstring += outtypestring + output_varname_str[i] + " = " + \
                             str(expr_convert_to_SIMD_intrins(result,SIMD_const_varnms,SIMD_const_values,outCparams.SIMD_debug)) + ";\n"
            else:
                outstring += outtypestring+ccode_postproc(sp.ccode(result,output_varname_str[i]))+"\n"

        # Step 6b.i: If SIMD_enable == True , and
        #            there is at least one SIMD const variable, 
        #            then declare the SIMD_const_varnms and SIMD_const_values arrays
        if outCparams.SIMD_enable == "True" and len(SIMD_const_varnms) != 0:
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
                exit(1)

            for i in range(len(SIMD_const_varnms)):
                if outCparams.enable_TYPE == "False":
                    SIMD_decls += outCparams.preindent + indent + SIMD_const_varnms[i] + " = " + SIMD_const_values[i]+";"
                else:
                    SIMD_decls += outCparams.preindent + indent + "const double " + outCparams.CSE_varprefix + SIMD_const_varnms[i] + " = " + SIMD_const_values[i] + ";\n"
                    SIMD_decls += outCparams.preindent+indent+ "const REAL_SIMD_ARRAY " + SIMD_const_varnms[i] + " = ConstSIMD("+ outCparams.CSE_varprefix + SIMD_const_varnms[i] + ");\n"
                SIMD_decls += "\n"

    # Step 7: Construct final output string
    final_Ccode_output_str = commentblock
    # Step 7a: Output C code in indented curly brackets if
    #          outCparams.includebraces = True
    if outCparams.includebraces == "True": final_Ccode_output_str += outCparams.preindent+"{\n"
    final_Ccode_output_str += prestring + SIMD_decls + outstring + poststring
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
