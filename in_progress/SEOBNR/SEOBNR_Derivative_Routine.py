# As documented in the NRPy+ tutorial module
#   SEOBNR_Derivative_Routine.ipynb,
#   this module computes partial derivatives of
#   an input list of expressions with respect
#   to an input list of free variables.

# Authors: Zachariah B. Etienne & Tyler Knowles
#         zachetie **at** gmail **dot* com

# Step 1.a: import all needed modules from Python/NRPy+:
import sympy as sp                # SymPy: The Python computer algebra package upon which NRPy+ depends
import os, sys                    # Standard Python modules for multiplatform OS-level functions

# Step 1.?: check system path so can use outputC; #TylerK: remove and put outputC back with other imports
nrpy_dir_path = os.path.join("..")
if nrpy_dir_path not in sys.path:
    sys.path.append(nrpy_dir_path)

from outputC import superfast_uniq, lhrh      # Remove duplicate entries from a Python array; store left- and right-
                                              #   hand sides of mathematical expressions

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

def symbolic_parital_derivative(expression_text_file,constants_text_file,variables_text_file):
    # Step 2.a: Read in expressions as a (single) string
    with open(expression_text_file, 'r') as file:
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
        lhss.append(sp.sympify(lr[i].lhs))
        rhss.append(sp.sympify(lr[i].rhs))

    # Step 3.a: Read in constants as a (single) string
    with open(constants_text_file, 'r') as file:
        constants_as_lines = file.readlines()

    # Step 3.b: Create "input_constants" array and populate with SymPy constants
    input_constants = []
    for constant in constants_as_lines:
        constant = sp.symbols(constant,real=True)
        input_constants.append(constant)

    # Step 3.c: Read in variables with which to take derivatives
    with open(variables_text_file, 'r') as file:
        variables_as_lines = file.readlines()

    #Step 3.d: Create "dynamic_variables" array and populate with SymPy symbols
    dynamic_variables = []
    for variable in variables_as_lines:
        variable = sp.symbols(variable,real=True)
        dynamic_variables.append(variable)

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
        lhss_deriv.append(sp.sympify(str(lhss[i])+"prm"))
        newrhs = sp.sympify(str(sp.diff(rhss[i],xx)).replace("(xx)","").replace(", xx","prm").replace("Derivative",""))
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
