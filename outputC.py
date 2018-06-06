import NRPy_param_funcs as par
import re
from SIMD import expr_convert_to_SIMD_intrins

from collections import namedtuple
lhrh = namedtuple('lhrh', 'lhs rhs')

# Parameter initialization is called once, within nrpy.py.
thismodule = __name__
par.initialize_param(par.glb_param("bool", thismodule, "SIMD_enable", False))
par.initialize_param(par.glb_param("bool", thismodule, "SIMD_debug", False))
par.initialize_param(par.glb_param("char", thismodule, "PRECISION", "double"))
#par.initialize_param(par.glb_param("bool", thismodule, "CSE_enable", True))
par.initialize_param(par.glb_param("char", thismodule, "CSE_varprefix", "tmp"))
par.initialize_param(par.glb_param("char", thismodule, "outCfileaccess", "w"))
par.initialize_param(par.glb_param("bool", thismodule, "outCverbose", True))
par.initialize_param(par.glb_param("bool", thismodule, "declareoutputvars", False))
par.initialize_param(par.glb_param("bool", thismodule, "includebraces", True))

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
    for func in ['pow', 'sqrt', 'sin', 'cos', 'tan', 'sinh', 'cosh', 'tanh', 'exp', 'log']:
        string2 = re.sub(func+'\(', func + cmathsuffix+"(", string); string = string2

    # Finally, SymPy prefers to output Rationals as long-double fractions.
    #  E.g., Rational(1,3) is output as 1.0L/3.0L.
    #  The Intel compiler vectorizer complains miserably about this,
    #  and strictly speaking it is useless when we're in double precision.
    # So here we get rid of the "L" suffix on floating point numbers:
    if PRECISION!="long double":
        string2 = re.sub('([0-9.]+)L/([0-9.]+)L', '(\\1 / \\2)', string); string = string2

    return string

import sympy as sp
# Input: sympyexpr = a single SymPy expression *or* a list of SymPy expressions
#        output_varname_str = a single output variable name *or* a list of output
#                             variable names, one per sympyexpr.
# Output: C code, as a string.
def outputC(sympyexpr, output_varname_str, filename = "stdout", CSE_enable = True, prestring = "", poststring = ""):
    TYPE = par.parval_from_str("PRECISION")

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
    if par.parval_from_str("SIMD_enable") == True:
        if TYPE != "double":
            print("SIMD output currently only supports double precision. Sorry!")
            exit(1)
        TYPE = "REAL_SIMD_ARRAY"

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


    # Step 3: If outputC::verbose = True, then output the original SymPy
    #         expression(s) in code comments prior to actual C code
    if par.parval_from_str("outCverbose") == True:
        commentblock += "/*\n *  Original SymPy expression"
        if len(output_varname_str)>1:
            commentblock += "s"
        commentblock += ":\n"
        for i in range(len(output_varname_str)):
            if i==0:
                if len(output_varname_str)!=1:
                    commentblock += " *  \"["
                else:
                    commentblock += " *  \""
            else:
                commentblock += " *    "
            commentblock += output_varname_str[i] + " = " + str(sympyexpr[i])
            if i==len(output_varname_str)-1:
                if len(output_varname_str)!=1:
                    commentblock += "]\"\n"
                else:
                    commentblock += "\"\n"
            else:
                commentblock += ",\n"
        commentblock += " */\n"

    # Step 4: Add proper indentation of C code:
    if par.parval_from_str("includebraces") == True:
        indent = "   "
    else:
        indent = ""

    # Step 5: Should the output variable, e.g., outvar, be declared?
    #         If so, start output line with e.g., "double outvar "
    outtypestring = ""
    if par.parval_from_str("declareoutputvars") == True:
        outtypestring = indent+TYPE + " "
    else:
        outtypestring = indent

    # Step 6a: If common subexpression elimination (CSE) disabled, then
    #         just output the SymPy string in the most boring way,
    #         nearly consistent with SymPy's ccode() function,
    #         though with support for float & long double types
    #         as well.
    SIMD_decls = ""

    if CSE_enable == False:
        # If CSE is disabled:
        for i in range(len(sympyexpr)):
            outstring += outtypestring + ccode_postproc(sp.ccode(sympyexpr[i], output_varname_str[i]))+"\n"
    # Step 6b: If CSE enabled, then perform CSE using SymPy and then
    #          resulting C code.
    else:
        # If CSE is enabled:
        SIMD_const_varnms = []
        SIMD_const_values = []

        CSE_varprefix = par.parval_from_str("CSE_varprefix")
        CSE_results = sp.cse(sympyexpr, sp.numbered_symbols(CSE_varprefix), order='canonical')
        for commonsubexpression in CSE_results[0]:
            if par.parval_from_str("SIMD_enable") == True:
                outstring += indent + "const " + TYPE + " " + str(commonsubexpression[0]) + " = " + \
                             str(expr_convert_to_SIMD_intrins(commonsubexpression[1],SIMD_const_varnms,SIMD_const_values)) + ";\n"
            else:
                outstring += indent+"const "+TYPE+" "+ccode_postproc(sp.ccode(commonsubexpression[1],commonsubexpression[0]))+"\n"
        for i,result in enumerate(CSE_results[1]):
            if par.parval_from_str("SIMD_enable") == True:
                outstring += outtypestring + output_varname_str[i] + " = " + \
                             str(expr_convert_to_SIMD_intrins(result,SIMD_const_varnms,SIMD_const_values)) + ";\n"
            else:
                outstring += outtypestring+ccode_postproc(sp.ccode(result,output_varname_str[i]))+"\n"

        # Step 6b.i: If SIMD_enable == True, then parse the SIMD_const_varnms and SIMD_const_values
        if par.parval_from_str("SIMD_enable") == True:
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
            else:
                for i in range(len(SIMD_const_varnms)):
                    SIMD_decls += indent+"const double " + CSE_varprefix + SIMD_const_varnms[i] + " = " + SIMD_const_values[i] + ";\n"
                    SIMD_decls += indent+"const REAL_SIMD_ARRAY " + " = Set1SIMD("+ CSE_varprefix + SIMD_const_varnms[i] + ");\n"
                    # if i != len(SIMD_const_varnms)-1:
                    #     SIMD_decls += "\n"
                SIMD_decls += "\n"

    # Step 7: Construct final output string
    final_Ccode_output_str = commentblock
    # Step 7a: Output C code in indented curly brackets if
    #          outputC::includebraces = True
    if par.parval_from_str("includebraces") == True: final_Ccode_output_str += "{\n"
    final_Ccode_output_str += prestring + SIMD_decls + outstring + poststring
    if par.parval_from_str("includebraces") == True: final_Ccode_output_str += "}"

    # Step 8: If parameter outputC::outCfilename = "stdout", then output
    #         C code to standard out (useful for copy-paste or interactive
    #         mode). Otherwise output to file specified in variable name.
    if filename == "stdout":
        # Output to standard out (stdout; "the screen")
        print(final_Ccode_output_str)
    elif filename == "returnstring":
        return final_Ccode_output_str
    else:
        # Output to the file specified by outCfilename
        with open(filename, par.parval_from_str("outCfileaccess")) as file:
            file.write(final_Ccode_output_str)
        successstr = ""
        if par.parval_from_str("outCfileaccess") == "a":
            successstr = "Appended to "
        elif par.parval_from_str("outCfileaccess") == "w":
            successstr = "Wrote "
        print(successstr + "to file \"" + par.parval_from_str("outCfilename") + "\"")
