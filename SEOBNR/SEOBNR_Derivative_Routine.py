# As documented in the NRPy+ tutorial module
#   Tutorial_SEOBNR_Derivative_Routine.ipynb,
#   this module computes partial derivatives
#   of the SEOBNRv3 Hamiltonian with respect
#   to 12 dynamic variables

# Authors: Zachariah B. Etienne & Tyler Knowles
#         zachetie **at** gmail **dot* com

# Step 1.a: import all needed modules from Python/NRPy+:
import sympy as sp                # SymPy: The Python computer algebra package upon which NRPy+ depends
import sys, os                    # Standard Python modules for multiplatform OS-level functions
from outputC import superfast_uniq, lhrh      # Remove duplicate entries from a Python array; store left- and right-
                                              #   hand sides of mathematical expressions
# Step 1.b: Check for a sufficiently new version of SymPy (for validation)
# Ignore the rc's and b's for release candidates & betas.
sympy_version = sp.__version__.replace('rc', '...').replace('b', '...')
sympy_version_decimal = float(sympy_version.split(".")[0]) + float(sympy_version.split(".")[1])/10.0
if sympy_version_decimal < 1.2:
    print('Error: NRPy+ does not support SymPy < 1.2')
    sys.exit(1)

# As of April 2021, "sp.sympify("Q+1")" fails because Q is a reserved keyword.
#   This is the workaround, courtesy Ken Sible.
custom_global_dict = {}
exec('from sympy import *', custom_global_dict)
del custom_global_dict['Q']
if sympy_version_decimal >= 1.6:
    custom_parse_expr = lambda expr: sp.parse_expr(expr, global_dict=custom_global_dict)
else:
    custom_parse_expr = lambda expr: sp.sympify(expr)

# Step 1.c: Name of the directory containing the input file
inputdir = "SEOBNR"

# Supporting function to simplify derivative expressions by removing terms equal to 0
def simplify_deriv(lhss_deriv,rhss_deriv):
    # Copy expressions into another array
    lhss_deriv_simp = []
    rhss_deriv_simp = []
    for i in range(len(rhss_deriv)):
        lhss_deriv_simp.append(lhss_deriv[i])
        rhss_deriv_simp.append(rhss_deriv[i])
    # If a right-hand side is 0, substitute value 0 for the corresponding left-hand side in later terms
    for i in range(len(rhss_deriv_simp)):
        if rhss_deriv_simp[i] == 0:
            for j in range(i+1,len(rhss_deriv_simp)):
                for var in rhss_deriv_simp[j].free_symbols:
                    if str(var) == str(lhss_deriv_simp[i]):
                        rhss_deriv_simp[j] = rhss_deriv_simp[j].subs(var,0)
    zero_elements_to_remove = []
    # Create array of indices for expressions that are zero
    for i in range(len(rhss_deriv_simp)):
        if rhss_deriv_simp[i] == sp.sympify(0):
            zero_elements_to_remove.append(i)

    # When removing terms that are zero, we need to take into account their new index (after each removal)
    count = 0
    for i in range(len(zero_elements_to_remove)):
        del lhss_deriv_simp[zero_elements_to_remove[i]+count]
        del rhss_deriv_simp[zero_elements_to_remove[i]+count]
        count -= 1
    return lhss_deriv_simp,rhss_deriv_simp

# Supporing function to convert a generic partial derivative into a partial derivative with respect to a specific variable
def deriv_onevar(lhss_deriv,rhss_deriv,variable_list,index):
    # Denote each variable with prm
    variableprm_list = []
    for variable in variable_list:
        variableprm_list.append(str(variable)+"prm")

    # Copy expressions into another array
    lhss_deriv_new = []
    rhss_deriv_new = []
    for i in range(len(rhss_deriv)):
        lhss_deriv_new.append(lhss_deriv[i])
        rhss_deriv_new.append(rhss_deriv[i])
    # For each free symbol's derivative, replace it with:
    #   1, if we are differentiating with respect to the variable, or
    #   0, if we are note differentiating with respect to that variable
    for i in range(len(rhss_deriv_new)):
        for var in variableprm_list:
            if variableprm_list.index(str(var))==index:
            #if var==(variable+"prm"):
                rhss_deriv_new[i] = rhss_deriv_new[i].subs(var,1)
            else:
                rhss_deriv_new[i] = rhss_deriv_new[i].subs(var,0)
    # Simplify derivative expressions again
    lhss_deriv_simp,rhss_deriv_simp = simplify_deriv(lhss_deriv_new,rhss_deriv_new)
    return lhss_deriv_simp,rhss_deriv_simp

def symbolic_parital_derivative():
    # Step 2.a: Read in expressions as a (single) string
    with open(os.path.join(inputdir,'Hamstring.txt'), 'r') as file:
        expressions_as_lines = file.readlines()

    # Step 2.b: Create and populate the "lr" array, which separates each line into left- and right-hand sides
    #   Each entry is a string of the form lhrh(lhs='',rhs='')
    lr = []

    for i in range(len(expressions_as_lines)):
        # Ignore lines with 2 or fewer characters and those starting with #
        if len(expressions_as_lines[i]) > 2 and expressions_as_lines[i][0] != "#":
            # Split each line by its equals sign
            split_line = expressions_as_lines[i].split("=")
            # Append the line to "lr", removing spaces, "sp." prefixes, and replacing Lambda->Lamb
            #   (Lambda is a protected keyword):
            lr.append(lhrh(lhs=split_line[0].replace(" ","").replace("Lambda","Lamb"),
                           rhs=split_line[1].replace(" ","").replace("sp.","").replace("Lambda","Lamb")))

    # Step 2.c: Separate and sympify right- and left-hand sides into separate arrays
    lhss = []
    rhss = []
    for i in range(len(lr)):
        lhss.append(custom_parse_expr(lr[i].lhs))
        rhss.append(custom_parse_expr(lr[i].rhs))

    # Step 3.a: Create `input_constants` array and populate with SymPy symbols
    m1,m2,tortoise,eta,KK,k0,k1,EMgamma,d1v2,dheffSSv2 = sp.symbols('m1 m2 tortoise eta KK k0 k1 EMgamma d1v2 dheffSSv2',
                                                                    real=True)
    input_constants = [m1,m2,tortoise,eta,KK,k0,k1,EMgamma,d1v2,dheffSSv2]

    # Step 3.b: Create `dynamic_variables` array and populate with SymPy symbols
    x,y,z,px,py,pz,s1x,s1y,s1z,s2x,s2y,s2z = sp.symbols('x y z px py pz s1x s1y s1z s2x s2y s2z', real=True)
    dynamic_variables = [x,y,z,px,py,pz,s1x,s1y,s1z,s2x,s2y,s2z]

    # Step 4.a: Prepare array of "free symbols" in the right-hand side expressions
    full_symbol_list_with_dups = []
    for i in range(len(lr)):
        for variable in rhss[i].free_symbols:
            full_symbol_list_with_dups.append(variable)

    # Step 4.b: Remove duplicate free symbols
    full_symbol_list = superfast_uniq(full_symbol_list_with_dups)

    # Step 4.c: Remove input constants from symbol list
    for inputconst in input_constants:
        for symbol in full_symbol_list:
            if str(symbol) == str(inputconst):
                full_symbol_list.remove(symbol)

    # Step 5.a: Convert each left-hand side to function notation
    #   while separating and simplifying left- and right-hand sides
    xx = sp.Symbol('xx')
    func = []
    for i in range(len(lr)):
        func.append(sp.sympify(sp.Function(lr[i].lhs)(xx)))

    # Step 5.b: Mark each free variable as a function with argument xx
    full_function_list = []
    for symb in full_symbol_list:
        func = sp.sympify(sp.Function(str(symb))(xx))
        full_function_list.append(func)
        for i in range(len(rhss)):
            for var in rhss[i].free_symbols:
                if str(var) == str(symb):
                    rhss[i] = rhss[i].subs(var,func)

    # Step 6.a: Use SymPy's diff function to differentiate right-hand sides with respect to xx
    #   and append "prm" notation to left-hand sides
    lhss_deriv = []
    rhss_deriv = []
    for i in range(len(rhss)):
        lhss_deriv.append(custom_parse_expr(str(lhss[i])+"prm"))
        newrhs = custom_parse_expr(str(sp.diff(rhss[i],xx)).replace("(xx)","").replace(", xx","prm").replace("Derivative",""))
        rhss_deriv.append(newrhs)

    # Step 7.b: Call the simplication function and then copy results
    lhss_deriv_simp,rhss_deriv_simp = simplify_deriv(lhss_deriv,rhss_deriv)
    lhss_deriv = lhss_deriv_simp
    rhss_deriv = rhss_deriv_simp

    # Step 8.b: Call the derivative function and populate dictionaries with the result
    lhss_derivative = {}
    rhss_derivative = {}
    for index in range(len(dynamic_variables)):
        lhss_temp,rhss_temp = deriv_onevar(lhss_deriv,rhss_deriv,dynamic_variables,index)
        lhss_derivative[dynamic_variables[index]] = lhss_temp
        rhss_derivative[dynamic_variables[index]] = rhss_temp

    # Step 9: Output original expression and each partial derivative expression in SymPy snytax
    with open("partial_derivatives.txt", "w") as output:
        for i in range(len(lr)):
            right_side = lr[i].rhs
            right_side_in_sp = right_side.replace("sqrt(","sp.sqrt(").replace("log(","sp.log(").replace("pi",
                                                    "sp.pi").replace("sign(","sp.sign(").replace("Abs(",
                                                    "sp.Abs(").replace("Rational(","sp.Rational(")
            output.write(str(lr[i].lhs)+" = "+right_side_in_sp)
        for var in dynamic_variables:
            for i in range(len(lhss_derivative[var])):
                right_side = str(rhss_derivative[var][i])
                right_side_in_sp = right_side.replace("sqrt(","sp.sqrt(").replace("log(","sp.log(").replace("pi",
                                                        "sp.pi").replace("sign(","sp.sign(").replace("Abs(",
                                                        "sp.Abs(").replace("Rational(","sp.Rational(").replace("prm",
                                                        "prm_"+str(var))
                output.write(str(lhss_derivative[var][i]).replace("prm","prm_"+str(var))+" = "+right_side_in_sp+"\n")
