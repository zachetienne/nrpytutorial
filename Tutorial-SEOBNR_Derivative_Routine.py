#!/usr/bin/env python
# coding: utf-8

# # Symbolic Partial Derivative Routine
# 
# ## Authors: Zach Etienne & Tyler Knowles
# 
# ## This module contains a routine for computing partial derivatives of a mathematical expression that is written as several subexpressions.
# 
# **Notebook Status:** <font color='green'><b> Validated </b></font>
# 
# **Validation Notes:** This tutorial notebook has been confirmed to be self-consistent with its corresponding NRPy+ module, as documented [below](#code_validation). Additionally, this notebook has been validated by checking that results are consistent with exact derivative expressions used in the SEOBNRv3_opt approixment of [LALSuite](https://git.ligo.org/lscsoft/lalsuite).
# 
# ### NRPy+ Source Code for this module: [SEOBNR_Derivative_Routine.py](../edit/SEOBNR/SEOBNR_Derivative_Routine.py)
# 
# ## Introduction
# $$\label{intro}$$
# 
# This notebook documents the symbolic partial derivative routine used to generate analytic derivatives of the [SEOBNRv3](https://git.ligo.org/lscsoft/lalsuite) Hamiltonian (documented [here](../Tutorial-SEOBNR_v3_Hamiltonian.ipynb)) and described in [this article](https://arxiv.org/abs/1803.06346).  In general, this notebook takes as input a file of inter-dependent mathematical expressions (in SymPy syntax), a file listing the names of values within those expressions, and a file listing all variables with which to take partial derivatives of each expression.  The output is a text file containing the original expression and those for each partial derivative computation.  The intention is to perform CSE on these expressions to create efficient partial derivative code!

# <a id='toc'></a>
# 
# # Table of Contents
# $$\label{toc}$$
# 
# This notebook is organized as follows
# 
# 1. [Step 1](#initializenrpy): Initialize core Python/NRPy+ modules
# 1. [Step 2:](#read_expressions) Read in Hamiltonian expressions from `Hamstring.txt`
# 1. [Step 3:](#list_constants) Specify constants and variables in Hamiltonian expression
# 1. [Step 4:](#list_free_symbols) Extract free symbols
# 1. [Step 5:](#convert_to_func) Convert variables to function notation; e.g., `var` goes to `var(xx)`
# 1. [Step 6:](#differentiate) Differentiate with respect to `xx`
# 1. [Step 7:](#remove_zeros) Remove derivatives (of constants) that evaluate to zero, simplifying derivative expressions
# 1. [Step 8:](#partial_derivative) Simplify derivatives with respect to a specific variable
# 1. [Step 9:](#store_results) Store partial derivatives to SymPy notebook `partial_derivatives.txt-VALIDATION.txt`
# 1. [Step 10:](#code_validation) Validate against LALSuite and trusted `SEOBNR_Derivative_Routine` NRPy+ module
# 1. [Step 11:](#latex_pdf_output) Output this notebook to $\LaTeX$-formatted PDF file

# <a id='initializenrpy'></a>
# 
# # Step 1: Initialize core Python/NRPy+ modules \[Back to [top](#toc)\]
# $$\label{initializenrpy}$$
# 
# Let's start by importing all the needed modules from Python/NRPy+ and creating the output directory (if it does not already exist):

# In[1]:


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

# Step 1.c: Name of the directory containing the input file
inputdir = "SEOBNR"


# <a id='read_expressions'></a>
# 
# # Step 2: Read in Hamiltonian expressions from `Hamstring.txt` \[Back to [top](#toc)\]
# $$\label{read_expressions}$$
# 
# We read in the expressions of which we will compute partial derivatives in a single large string before splitting the string by line (carriage return) and by "=".  Doing so allows us to manipulate the right- and left-hand sides of the expressions appropriately.  We store the left- and right-hand sides in the array `lr`, which consists of `lhrh` arrays with left-hand sides `lhs` and right-hand sides `rhs`.  Note that `Lambda` is a protected keyword in Python, so the variable $\Lambda$ in the Hamiltonian is renamed `Lamb`.

# In[2]:


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

# As of April 2021, "sp.sympify("Q+1")" fails because Q is a reserved keyword.
#   This is the workaround, courtesy Ken Sible.
custom_global_dict = {}
exec('from sympy import *', custom_global_dict)
del custom_global_dict['Q']
print(sympy_version_decimal)
if sympy_version_decimal > 1.2:
    custom_parse_expr = lambda expr: sp.parse_expr(expr, global_dict=custom_global_dict)
else:
    custom_parse_expr = lambda expr: sp.sympify(expr)

for i in range(len(lr)):
    lhss.append(custom_parse_expr(lr[i].lhs))
    rhss.append(custom_parse_expr(lr[i].rhs))


# <a id='list_constants'></a>
# 
# # Step 3: Specify constants and variables in Hamiltonian expression \[Back to [top](#toc)\]
# $$\label{list_constants}$$
# 
# We read in and declare as SymPy symbols the constant values; derivatives with respect to these variables will be set to zero.  We then read in the variables with respect to which we want to take derivatives and declare those as SymPy variables as well.

# In[3]:


# Step 3.a: Create `input_constants` array and populate with SymPy symbols
m1,m2,tortoise,eta,KK,k0,k1,EMgamma,d1v2,dheffSSv2 = sp.symbols('m1 m2 tortoise eta KK k0 k1 EMgamma d1v2 dheffSSv2',
                                                                real=True)
input_constants = [m1,m2,tortoise,eta,KK,k0,k1,EMgamma,d1v2,dheffSSv2]

# Step 3.b: Create `dynamic_variables` array and populate with SymPy symbols
x,y,z,px,py,pz,s1x,s1y,s1z,s2x,s2y,s2z = sp.symbols('x y z px py pz s1x s1y s1z s2x s2y s2z', real=True)
dynamic_variables = [x,y,z,px,py,pz,s1x,s1y,s1z,s2x,s2y,s2z]


# <a id='list_free_symbols'></a>
# 
# # Step 4: Extract free symbols \[Back to [top](#toc)\]
# $$\label{list_free_symbols}$$
# 
# By ''free symbols'' we mean the variables in the right-hand sides.  We first create a list of all such terms (using SymPy's built-in free_symbol attribute), including duplicates, and then strip the duplicates.  We then remove input constants from the symbol list.

# In[4]:


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


# <a id='convert_to_func'></a>
# 
# # Step 5: Convert variables to function notation; e.g., `var` goes to `var(xx)` \[Back to [top](#toc)\]
# $$\label{convert_to_func}$$
# 
# In order to compute the partial derivative of each right-hand side, we mark each variable (left-hand side) and each free symbol (in right-hand sides) as a function with argument $\texttt{xx}$.

# In[5]:


# Step 5.a: Convert each left-hand side to function notation
#   while separating and simplifying left- and right-hand sides
xx = sp.Symbol('xx',real=True)
func = []
for i in range(len(lr)):
    func.append(sp.sympify(sp.Function(lr[i].lhs,real=True)(xx)))

# Step 5.b: Mark each free variable as a function with argument xx
full_function_list = []
for symb in full_symbol_list:
    func = sp.sympify(sp.Function(str(symb),real=True)(xx))
    full_function_list.append(func)
    for i in range(len(rhss)):
        for var in rhss[i].free_symbols:
            if str(var) == str(symb):
                rhss[i] = rhss[i].subs(var,func)


# <a id='differentiate'></a>
# 
# # Step 6: Differentiate with respect to `xx` \[Back to [top](#toc)\]
# $$\label{differentiate}$$
# 
# Now we differentiate the right-hand expressions with respect to `xx`.  We use the SymPy $\texttt{diff}$ command, differentiating with respect to $\texttt{xx}$.  After so doing, we remove $\texttt{(xx)}$ and "Derivative" (which is output by $\texttt{diff}$), and use "prm" suffix to denote the derivative with respect to $\texttt{xx}$.

# In[6]:


# Step 6: Use SymPy's diff function to differentiate right-hand sides with respect to xx
#   and append "prm" notation to left-hand sides
lhss_deriv = []
rhss_deriv = []
for i in range(len(rhss)):
    lhss_deriv.append(custom_parse_expr(str(lhss[i])+"prm"))
    newrhs = custom_parse_expr(str(sp.diff(rhss[i],xx)).replace("(xx)","").replace(", xx","prm").replace("Derivative",""))
    rhss_deriv.append(newrhs)


# <a id='remove_zeros'></a>
# 
# # Step 7: Remove derivatives (of constants) that evaluate to zero, simplifying derivative expressions \[Back to [top](#toc)\]
# $$\label{remove_zeros}$$
# 
# We declare a function to simply the derivative expressions.  In particular, we want to remove terms equal to zero.

# In[7]:


# Step 7.a: Define derivative simplification function
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

# Step 7.b: Call the simplication function and then copy results
lhss_deriv_simp,rhss_deriv_simp = simplify_deriv(lhss_deriv,rhss_deriv)
lhss_deriv = lhss_deriv_simp
rhss_deriv = rhss_deriv_simp


# <a id='partial_derivative'></a>
# 
# # Step 8: Simplify derivatives with respect to a specific variable \[Back to [top](#toc)\]
# $$\label{partial_derivative}$$
# 
# In [Step 6](#differentiate) we took a generic derivative of each expression, assuming all variables were functions of `xx`.  We now define a function that will select a specific dynamic variable (element of `dynamic_variables`) and set the derivative of the variable to 1 and all others to 0.

# In[8]:


# Step 8.a: Define onevar derivative function
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
                rhss_deriv_new[i] = rhss_deriv_new[i].subs(var,1)
            else:
                rhss_deriv_new[i] = rhss_deriv_new[i].subs(var,0)
    # Simplify derivative expressions again
    lhss_deriv_simp,rhss_deriv_simp = simplify_deriv(lhss_deriv_new,rhss_deriv_new)
    return lhss_deriv_simp,rhss_deriv_simp

# Step 8.b: Call the derivative function and populate dictionaries with the result
lhss_derivative = {}
rhss_derivative = {}
for index in range(len(dynamic_variables)):
    lhss_temp,rhss_temp = deriv_onevar(lhss_deriv,rhss_deriv,dynamic_variables,index)
    lhss_derivative[dynamic_variables[index]] = lhss_temp
    rhss_derivative[dynamic_variables[index]] = rhss_temp


# <a id='store_results'></a>
# 
# # Step 9: Store partial derivatives to SymPy notebook `partial_derivatives.txt-VALIDATION.txt` \[Back to [top](#toc)\]
# $$\label{store_results}$$
# 
# We write the resulting derivatives in SymPy syntax.  Each partial derivative is output in its own file, in a similar format to the input expressions.

# In[9]:


# Step 9: Output original expression and each partial derivative expression in SymPy snytax
with open(os.path.join(inputdir,'partial_derivatives.txt-VALIDATION'), 'w') as output:
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


# <a id='code_validation'></a>
# 
# # Step 10: Validate against LALSuite and trusted `SEOBNR_Derivative_Routine` NRPy+ module \[Back to [top](#toc)\]
# $$\label{code_validation}$$
# 
# We validate the output of this notebook against known LALSuite values of the Hamiltonian partial derivatives and the output of the `SEOBNR_Derivative_Routine` NRPy+ module.  We note that due to cancellations in the deriavtive terms, various versions of SymPy may result in relative errors that differ as much as an order of magnitude.  Furthermore, even changing the set of input pararameters can affect the relative error by as many as two orders of magnitude.  Therefore we look for agreement with LALSuite to at least 10 significant digits.
# 
# When comparing the notebook output to that of the NRPy+ module, we compare term-by-term using SymPy to check that each right-hand side side is equivalent.

# In[10]:


# Define a function to return a set of reasonable input parameters
# This function contains three distinct sets of input parameters, and index differentiates between them
def reset_values(tort_value, index):
    # Check that a reasonable tortoise value has been passed
    if tort_value!=1 and tort_value!=2:
        print("Error: tortoise must be 1 or 2.")
        sys.exit(1)
    # Each index corresponds to a specific set of input parameters
    if index==0:#-M 13 -m 11 -X 0.1 -Y -0.2 -Z 0.3 -x -0.3 -y 0.2 -z -0.1
        values = {'m1': 1.300000000000000e+01, 'm2': 1.100000000000000e+01, 'eta': 2.482638888888889e-01,
                  'x': 1.658426645098320e+01, 'y': 3.975021008701605e-02, 'z': -1.820682538442627e-07,
                  's1x': 2.934751675254397e-02, 's1y': -5.867672205485316e-02, 's1z': 8.802097562761332e-02,
                  's2x': -6.302678133897792e-02, 's2y': 4.200490780215727e-02, 's2z': -2.100705983874398e-02,
                  'KK': 3.913980338468737e-01, 'k0': -7.447639215330089e-01, 'k1': -6.380586501824999e-01,
                  'd1v2': -7.476323019145448e+01,'dheffSSv2':2.105103187692902e+01,
                  'EMgamma': 0.577215664901532860606512090082402431}
        # Note that we transform the momentum based on the tortoise values
        if tort_value==1:
            values.update({'px': -1.517631642228534e-03, 'py': 2.693180445886167e-01, 'pz': -1.320499830947482e-04,
                          'tortoise': 1})
        else:
            values.update({'px': -1.633028076483384e-03, 'py': 2.693177679992048e-01, 'pz': -1.320499918278832e-04,
                          'tortoise': 2})
    elif index==1:#-M 25 -m 10 -X 0.1 -Y -0.0 -Z 0.1 -x -0.2 -y 0.0 -z -0.2
        values = {'m1': 2.500000000000000e+01, 'm2': 1.000000000000000e+01, 'eta': 2.040816326530612e-01,
                  'x': 1.289689003662444e+01, 'y': 5.495441315063273e-03, 'z': -1.717482806041791e-11,
                  's1x': 5.102040816179230e-02, 's1y': 9.846215537206260e-07, 's1z': 5.102040816473832e-02,
                  's2x': -1.632653061189792e-02, 's2y': -6.762952223804450e-07, 's2z': -1.632653061259188e-02,
                  'KK': 5.642540639599580e-01, 'k0': -1.063532077165767e+00, 'k1': -8.835684149841774e-01,
                  'd1v2': -8.041179092044979e+01,'dheffSSv2':1.125986130778842e+01,
                  'EMgamma': 0.577215664901532860606512090082402431}
        if tort_value==1:
            values.update({'px': -1.898773926867491e-03, 'py': 3.160984442121970e-01, 'pz': 1.171602901570564e-07,
                           'tortoise': 1})
        else:
            values.update({'px': -2.209215477700561e-03, 'py': 3.160983119312114e-01, 'pz': 1.171602905704723e-07,
                          'tortoise': 2})
    elif index==2:#-M 7 -m 5 -X 0.01 -Y -0.5 -Z 0.03 -x -0.04 -y 0.05 -z -0.06
        values = {'m1': 7.000000000000000e+00, 'm2': 5.000000000000000e+00, 'eta': 2.430555555555556e-01,
                  'x': 2.633506161699224e+01, 'y': 7.574563213724998e-02, 'z': -2.789625823248071e-08,
                      's1x': 3.417297286269225e-03, 's1y': -1.701385963191495e-01, 's1z': 1.020835932957879e-02,
                      's2x': -6.945454346305877e-03, 's2y': 8.679766617922793e-03, 's2z': -1.041665076794264e-02,
                      'KK': 4.052853693162246e-01, 'k0': -7.706473492549312e-01, 'k1': -6.587426366263742e-01,
                      'd1v2': -7.555647472993827e+01,'dheffSSv2':1.972817669753086e+01,
                      'EMgamma': 0.577215664901532860606512090082402431}
        if tort_value==1:
            values.update({'px': -7.883793607066706e-04, 'py': 2.068742709904638e-01, 'pz': -7.338789145500886e-04,
                          'tortoise': 1})
        else:
            values.update({'px': -8.039726989861640e-04, 'py': 2.068742261404732e-01, 'pz': -7.338789145335709e-04,
                           'tortoise': 2})
    else:
        # If an improper index is passed, exit
        print("Error: invalid index (only three sets of input parameters available).")
        sys.exit(1)
    # Return the input values
    return values

# Numerically evaluate right-hand sides using input values
def evaluate_expression(left_sides,right_sides,input_values):
    new_right_sides = []
    for i in range(len(right_sides)):
        term = custom_parse_expr(str(right_sides[i]).replace("(xx)",""))
        # Only look for the free variables in each expression to reduce computation time
        free_vars = term.free_symbols
        for variable in free_vars:
            term = term.subs(variable, input_values[str(variable)])
        # Evaluate each term to reduce computation time
        new_right_sides.append(sp.sympify(term.evalf()))
        # Store each subexpression in values numerically
        input_values[str(left_sides[i])] = new_right_sides[i]
    # Return the input values dictionary with all numerical right-hand added
    return input_values

# Create array of trusted LALSuite derivative values
# Note that position in the array corresponds to the index of the corresponding input values
LALSuite_validated_values = []
#-M 13 -m 11 -X 0.1 -Y -0.2 -Z 0.3 -x -0.3 -y 0.2 -z -0.1
LALSuite_validated_values.append({'Hreal': 9.928923110195770e-01,'dHreal_dx': 9.932484846748471e-04,
                                  'dHreal_dy': 2.813294366789505e-06, 'dHreal_dz': 1.926378549762488e-06,
                                  'dHreal_dpx': -3.710666135737856e-04, 'dHreal_dpy': 6.116199124763537e-02,
                                  'dHreal_dpz': -5.600910364542288e-07, 'dHreal_ds1x': -1.438467658934620e-05,
                                  'dHreal_ds1y': -1.319462868057848e-06, 'dHreal_ds1z': 7.665413183773232e-04,
                                  'dHreal_ds2x': -2.075691477548065e-05,'dHreal_ds2y': 2.456427688083135e-06,
                                  'dHreal_ds2z': 8.762835349889455e-04})
#-M 25 -m 10 -X 0.1 -Y -0.0 -Z 0.1 -x -0.2 -y 0.0 -z -0.2
LALSuite_validated_values.append({'Hreal': 9.926852598351464e-01, 'dHreal_dx': 1.397519118422771e-03,
                                  'dHreal_dy': 1.928133240540033e-06, 'dHreal_dz': -1.215449398950413e-06,
                                  'dHreal_dpx': -4.004159849919695e-04, 'dHreal_dpy': 5.701850933742150e-02,
                                  'dHreal_dpz': 4.329487960716782e-08, 'dHreal_ds1x': 2.259457049322466e-06,
                                  'dHreal_ds1y': -2.544122765762015e-09, 'dHreal_ds1z': 9.834156257814124e-04,
                                  'dHreal_ds2x': 5.185557993931246e-06,'dHreal_ds2y': 2.437768415468806e-10,
                                  'dHreal_ds2z': 2.111169766641698e-03})
#-M 7 -m 5 -X 0.01 -Y -0.5 -Z 0.03 -x -0.04 -y 0.05 -z -0.06
LALSuite_validated_values.append({'Hreal': 9.955293642650920e-01, 'dHreal_dx': 3.734697245297603e-04,
                                  'dHreal_dy': 1.105998063449349e-06, 'dHreal_dz': 5.367207414282669e-08,
                                  'dHreal_dpx': -1.848412708548443e-04, 'dHreal_dpy': 4.754239153769983e-02,
                                  'dHreal_dpz': -3.549083643069269e-08, 'dHreal_ds1x': -4.819261725948465e-07,
                                  'dHreal_ds1y': 3.333280059627902e-06, 'dHreal_ds1z': 2.201786563823208e-04,
                                  'dHreal_ds2x': -7.576810957551029e-07,'dHreal_ds2y': 6.818093508597533e-06,
                                  'dHreal_ds2z': 2.922663340179887e-04})

# Sort variables by which tortoise value we use to compute the derivatives
variables_tort2 = [x,y,z]
variables_tort1 = [px,py,pz,s1x,s1y,s1z,s2x,s2y,s2z]

# Call evaluation function
print("Computing the difference between notebook output and trusted LALSuite derivative values...")
for index in range(3):
    values = reset_values(2,index)
    values = evaluate_expression(lhss,rhss,values)
    Hreal = values['Hreal']
    trusted_Hreal = LALSuite_validated_values[index]['Hreal']
    relative_difference = (trusted_Hreal - Hreal)/trusted_Hreal
    if abs(relative_difference) > 1e-9:
        print("The relative difference in Hreal is too large: %.15e" % relative_difference)
        sys.exit(1)

    for var in variables_tort2:
        Hrealprm = evaluate_expression(lhss_derivative[var],rhss_derivative[var],values)['Hrealprm']
        trusted_Hrealprm = LALSuite_validated_values[index]['dHreal_d'+str(var)]
        relative_difference = (trusted_Hrealprm - Hrealprm)/trusted_Hrealprm
        if abs(relative_difference) > 1e-9:
            print("The relative difference in dHreal_d%s is too large: %.15e" % (var,relative_difference))
            sys.exit(1)

    values = reset_values(1,index)
    values = evaluate_expression(lhss,rhss,values)

    for var in variables_tort1:
        Hrealprm = evaluate_expression(lhss_derivative[var],rhss_derivative[var],values)['Hrealprm']
        trusted_Hrealprm = LALSuite_validated_values[index]['dHreal_d'+str(var)]
        relative_difference = (trusted_Hrealprm - Hrealprm)/trusted_Hrealprm
        if abs(relative_difference) > 1e-9:
            print("The relative difference in dHreal_d%s is too large: %.15e" % (var,relative_difference))
            sys.exit(1)
print("Test passed: the notebook agrees with LALSuite to at least 10 significant digits!")

print("Printing difference between notebook output and trusted NRPy+ module output...")
# Open the files to compare
file = 'partial_derivatives.txt'
outfile = 'partial_derivatives.txt-VALIDATION'

print("Checking file " + outfile)
with open(os.path.join(inputdir,file), "r") as file1, open(os.path.join(inputdir,outfile), "r") as file2:
    # Read the lines of each file
    file1_lines = file1.readlines()
    file2_lines = file2.readlines()
    # Compare right-hand sides of the expressions by computing the difference between them
    num_diffs = 0
    for i in range(len(file1_lines)):
        expr_new = custom_parse_expr(file1_lines[i].split("=")[1].replace("sp.",""))
        expr_validated = custom_parse_expr(file2_lines[i].split("=")[1].replace("sp.",""))
        difference = sp.simplify(expr_new - expr_validated)
        if difference != 0:
            num_diffs += 1
            print(difference)
    if num_diffs == 0:
        print("No difference. TEST PASSED!")
    else:
        print("ERROR: Disagreement found with the trusted file. See differences above.")
        sys.exit(1)


# <a id='latex_pdf_output'></a>
# 
# # Step 11: Output this notebook to $\LaTeX$-formatted PDF file \[Back to [top](#toc)\]
# $$\label{latex_pdf_output}$$
# 
# The following code cell converts this Jupyter notebook into a proper, clickable $\LaTeX$-formatted PDF file. After the cell is successfully run, the generated PDF may be found in the root NRPy+ tutorial directory, with filename
# [Tutorial-SEOBNR_Derivative_Routine.pdf](Tutorial-SEOBNR_Derivative_Routine.pdf) (Note that clicking on this link may not work; you may need to open the PDF file through another means.)

# In[11]:


import cmdline_helper as cmd      # NRPy+: Multi-platform Python command-line interface
cmd.output_Jupyter_notebook_to_LaTeXed_PDF("Tutorial-SEOBNR_Derivative_Routine")

