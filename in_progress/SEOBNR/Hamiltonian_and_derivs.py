# NRPy+ code to generate first derivatives of the SEOBNRv3 Hamiltonian from a list of numerical expressions computing
# said Hamiltonian. Originally written by Zach Etienne; edited and commented by Tyler Knowles.

from outputC import outputC,lhrh,superfast_uniq  # NRPy+: Core C code output module
import sympy as sp      # SymPy: The Python computer algebra package upon which NRPy+ depends
import sys              # Python module for multiplatform OS-related functions

# simplify_deriv() simplifies derivative expressions by removing terms equal to zero.
def simplify_deriv(lhss_deriv, rhss_deriv):
    # Create 'simp' arrays to store and manipulate derivative expressions.
    lhss_deriv_simp = []
    rhss_deriv_simp = []
    # Append terms to 'simp' arrays.
    for i in range(len(rhss_deriv)):
        lhss_deriv_simp.append(lhss_deriv[i])
        rhss_deriv_simp.append(rhss_deriv[i])
    # For each term equal to zero, loop through all expressions and replace that variable with the number zero.
    for i in range(len(rhss_deriv_simp)):
        if rhss_deriv_simp[i] == 0:
            for j in range(i + 1, len(rhss_deriv_simp)):
                for var in rhss_deriv_simp[j].free_symbols:
                    if str(var) == str(lhss_deriv_simp[i]):
                        rhss_deriv_simp[j] = rhss_deriv_simp[j].subs(var, 0)
    # Create 'zero' array to store terms to be removed from derivative expressions.
    zero_elements_to_remove = []
    # Loop over all terms and add those equal to zero to 'zero' array.
    for i in range(len(rhss_deriv_simp)):
        if rhss_deriv_simp[i] == sp.sympify(0):
            zero_elements_to_remove.append(i)
    count = 0
    # Remove from derivative list all elements of 'zero' array.
    for i in range(len(zero_elements_to_remove)):
        del lhss_deriv_simp[zero_elements_to_remove[i] + count]
        del rhss_deriv_simp[zero_elements_to_remove[i] + count]
        count -= 1
    # Return simplified derivative expressions.
    return lhss_deriv_simp, rhss_deriv_simp


# deriv_onevar() replaces variable derivatives with 1 or 0 depending on which partial derivaitve is computed.  For
# example, pass 'xprm=1' to replace each instance of 'xprm' with 1 and 'qprm' with 0 for each q in (y,z,px,py,pz,s1x,
# s1y,s1z,s2x,s2y,s2z).  This produces expressions which compute the partial derivative of the Hamiltonian with respect
# to x.
def deriv_onevar(lhss_deriv, rhss_deriv, xprm=0, yprm=0, zprm=0, pxprm=0, pyprm=0, pzprm=0, s1xprm=0, s1yprm=0,
                 s1zprm=0, s2xprm=0, s2yprm=0, s2zprm=0):
    if xprm + yprm + zprm + pxprm + pyprm + pzprm + s1xprm + s1yprm + s1zprm + s2xprm + s2yprm + s2zprm != 1:
        print("deriv_onevar() cannot take more than one derivative at a time!")
        sys.exit()

    # Create 'new' arrays to store and manipulate derivative terms.
    lhss_deriv_new = []
    rhss_deriv_new = []
    # Append derivative terms to 'new' arrays
    for i in range(len(rhss_deriv)):
        lhss_deriv_new.append(lhss_deriv[i])
        rhss_deriv_new.append(rhss_deriv[i])
    # Replace each instance of 'qprm', q in (x,y,z,px,py,pz,s1x,s1y,s1z,s2x,s2y,s2z), with either 0 or 1.
    for i in range(len(rhss_deriv_new)):
        for var in rhss_deriv_new[i].free_symbols:
            if str(var) == "xprm":
                rhss_deriv_new[i] = rhss_deriv_new[i].subs(var, xprm)
            elif str(var) == "yprm":
                rhss_deriv_new[i] = rhss_deriv_new[i].subs(var, yprm)
            elif str(var) == "zprm":
                rhss_deriv_new[i] = rhss_deriv_new[i].subs(var, zprm)
            elif str(var) == "pxprm":
                rhss_deriv_new[i] = rhss_deriv_new[i].subs(var, pxprm)
            elif str(var) == "pyprm":
                rhss_deriv_new[i] = rhss_deriv_new[i].subs(var, pyprm)
            elif str(var) == "pzprm":
                rhss_deriv_new[i] = rhss_deriv_new[i].subs(var, pzprm)
            elif str(var) == "s1xprm":
                rhss_deriv_new[i] = rhss_deriv_new[i].subs(var, s1xprm)
            elif str(var) == "s1yprm":
                rhss_deriv_new[i] = rhss_deriv_new[i].subs(var, s1yprm)
            elif str(var) == "s1zprm":
                rhss_deriv_new[i] = rhss_deriv_new[i].subs(var, s1zprm)
            elif str(var) == "s2xprm":
                rhss_deriv_new[i] = rhss_deriv_new[i].subs(var, s2xprm)
            elif str(var) == "s2yprm":
                rhss_deriv_new[i] = rhss_deriv_new[i].subs(var, s2yprm)
            elif str(var) == "s2zprm":
                rhss_deriv_new[i] = rhss_deriv_new[i].subs(var, s2zprm)
    # Simplify the derivative expressions with simplify_deriv().
    lhss_deriv_simp, rhss_deriv_simp = simplify_deriv(lhss_deriv_new, rhss_deriv_new)
    # Return simplified derivative expression.
    return lhss_deriv_simp, rhss_deriv_simp


# replace_numpy_funcs() replaces specific SymPy function names with the corresponding NumPy function names.
def replace_numpy_funcs(expression):
    return str(expression).replace("sqrt(", "sp.sqrt(").replace("Abs(", "sp.Abs(").replace("log(",
                                                                "sp.log(").replace("sign(", "sp.sign(")


# output_H_and_derivs() is the main wrapper function for computing the SEONBRv3 Hamiltonian H and the twelve first
# partial derivatives of H with respect to x, y, z, px, py, pz, s1x, s1y, s1z, s2x, s2y, s2z.
def output_H_and_derivs():
    # Open and read the file of numerical expressions (written in SymPy syntax) computing the SEOBNRv3 Hamiltonian.
    f = open("SEOBNR/Hamstring.txt", 'r')
    Hamstring = str(f.read())
    f.close()

    # Split Hamstring by carriage returns.
    Hamterms = Hamstring.splitlines()

    # Create 'lr' array to store each left-hand side and right-hand side of Hamstring as strings.
    lr = []
    # Loop over each line in Hamstring to separate the left- and right-hand sides.
    for i in range(len(Hamterms)):
        # Ignore lines with 2 or fewer characters and those starting with #
        if len(Hamterms[i]) > 2 and Hamterms[i][0] != "#":
            # Split each line by its equals sign.
            splitHamterms = Hamterms[i].split("=")
            # Append terms to the 'lr' array, removing spaces, "sp." prefixes, and replacing Lambda->Lamb (Lambda is a
            # protected keyword)
            lr.append(lhrh(lhs=splitHamterms[0].replace(" ", "").replace("Lambda", "Lamb"),
                           rhs=splitHamterms[1].replace(" ", "").replace("sp.", "").replace("Lambda", "Lamb")))
    # Declare the symbol 'xx', which we use to denote each left-hand side as a function
    xx = sp.Symbol('xx')
    # Create arrays to store simplified left- and right-hand expressions, as well as left-hand sides designated as
    # functions.
    func = []
    lhss = []
    rhss = []
    # Affix '(xx)' to each left-hand side as a function designation; separate and simplify left- and right-hand sides
    # of the numerical expressions.
    for i in range(len(lr)):
        func.append(sp.sympify(sp.Function(lr[i].lhs)(xx)))
        lhss.append(sp.sympify(lr[i].lhs))
        rhss.append(sp.sympify(lr[i].rhs))
    # Creat array for and generate a list of all the "free symbols" in the right-hand side expressions.
    full_symbol_list_with_dups = []
    for i in range(len(lr)):
        for var in rhss[i].free_symbols:
            full_symbol_list_with_dups.append(var)

    # Remove all duplicated "free symbols" from the right-hand side expressions.
    full_symbol_list = superfast_uniq(full_symbol_list_with_dups)

    # Declare input constants.
    m1, m2, eta, KK, k0, k1, d1v2, dheffSSv2 = sp.symbols("m1 m2 eta KK k0 k1 d1v2 dheffSSv2", real=True)
    tortoise, EMgamma = sp.symbols("tortoise EMgamma", real=True)
    input_constants = [m1, m2, eta, KK, k0, k1, d1v2, dheffSSv2, tortoise, EMgamma]

    # Derivatives of input constants will always be zero, so remove them from the full_symbol_list.
    for inputconst in input_constants:
        for symbol in full_symbol_list:
            if str(symbol) == str(inputconst):
                full_symbol_list.remove(symbol)

    # Add symbols to the function list and replace right-hand side terms with their function equivalent.
    full_function_list = []
    for symb in full_symbol_list:
        func = sp.sympify(sp.Function(str(symb))(xx))
        full_function_list.append(func)
        for i in range(len(rhss)):
            for var in rhss[i].free_symbols:
                if str(var) == str(symb):
                    rhss[i] = rhss[i].subs(var, func)

    # Create left- and right-hand side 'deriv' arrays
    lhss_deriv = []
    rhss_deriv = []
    # Differentiate with respect to xx, remove '(xx)', and replace xx with 'prm' notation.
    for i in range(len(rhss)):
        lhss_deriv.append(sp.sympify(str(lhss[i]) + "prm"))
        newrhs = sp.sympify(
            str(sp.diff(rhss[i], xx)).replace("(xx)", "").replace(", xx", "prm").replace("Derivative", ""))
        rhss_deriv.append(newrhs)
    # Simplify derivative expressions with simplify_deriv()
    lhss_deriv_simp, rhss_deriv_simp = simplify_deriv(lhss_deriv, rhss_deriv)
    lhss_deriv = lhss_deriv_simp
    rhss_deriv = rhss_deriv_simp
    # Generate partial derivatives with respect to each of the twelve input variables
    lhss_deriv_x, rhss_deriv_x = deriv_onevar(lhss_deriv, rhss_deriv, xprm=1, yprm=0, zprm=0, pxprm=0, pyprm=0, pzprm=0,
                                              s1xprm=0, s1yprm=0, s1zprm=0, s2xprm=0, s2yprm=0, s2zprm=0)
    lhss_deriv_y, rhss_deriv_y = deriv_onevar(lhss_deriv, rhss_deriv, xprm=0, yprm=1, zprm=0, pxprm=0, pyprm=0, pzprm=0,
                                              s1xprm=0, s1yprm=0, s1zprm=0, s2xprm=0, s2yprm=0, s2zprm=0)
    lhss_deriv_z, rhss_deriv_z = deriv_onevar(lhss_deriv, rhss_deriv, xprm=0, yprm=0, zprm=1, pxprm=0, pyprm=0, pzprm=0,
                                              s1xprm=0, s1yprm=0, s1zprm=0, s2xprm=0, s2yprm=0, s2zprm=0)
    lhss_deriv_px, rhss_deriv_px = deriv_onevar(lhss_deriv, rhss_deriv, xprm=0, yprm=0, zprm=0, pxprm=1, pyprm=0,
                                                pzprm=0, s1xprm=0, s1yprm=0, s1zprm=0, s2xprm=0, s2yprm=0, s2zprm=0)
    lhss_deriv_py, rhss_deriv_py = deriv_onevar(lhss_deriv, rhss_deriv, xprm=0, yprm=0, zprm=0, pxprm=0, pyprm=1,
                                                pzprm=0, s1xprm=0, s1yprm=0, s1zprm=0, s2xprm=0, s2yprm=0, s2zprm=0)
    lhss_deriv_pz, rhss_deriv_pz = deriv_onevar(lhss_deriv, rhss_deriv, xprm=0, yprm=0, zprm=0, pxprm=0, pyprm=0,
                                                pzprm=1, s1xprm=0, s1yprm=0, s1zprm=0, s2xprm=0, s2yprm=0, s2zprm=0)
    lhss_deriv_s1x, rhss_deriv_s1x = deriv_onevar(lhss_deriv, rhss_deriv, xprm=0, yprm=0, zprm=0, pxprm=0, pyprm=0,
                                                  pzprm=0, s1xprm=1, s1yprm=0, s1zprm=0, s2xprm=0, s2yprm=0, s2zprm=0)
    lhss_deriv_s1y, rhss_deriv_s1y = deriv_onevar(lhss_deriv, rhss_deriv, xprm=0, yprm=0, zprm=0, pxprm=0, pyprm=0,
                                                  pzprm=0, s1xprm=0, s1yprm=1, s1zprm=0, s2xprm=0, s2yprm=0, s2zprm=0)
    lhss_deriv_s1z, rhss_deriv_s1z = deriv_onevar(lhss_deriv, rhss_deriv, xprm=0, yprm=0, zprm=0, pxprm=0, pyprm=0,
                                                  pzprm=0, s1xprm=0, s1yprm=0, s1zprm=1, s2xprm=0, s2yprm=0, s2zprm=0)
    lhss_deriv_s2x, rhss_deriv_s2x = deriv_onevar(lhss_deriv, rhss_deriv, xprm=0, yprm=0, zprm=0, pxprm=0, pyprm=0,
                                                  pzprm=0, s1xprm=0, s1yprm=0, s1zprm=0, s2xprm=1, s2yprm=0, s2zprm=0)
    lhss_deriv_s2y, rhss_deriv_s2y = deriv_onevar(lhss_deriv, rhss_deriv, xprm=0, yprm=0, zprm=0, pxprm=0, pyprm=0,
                                                  pzprm=0, s1xprm=0, s1yprm=0, s1zprm=0, s2xprm=0, s2yprm=1, s2zprm=0)
    lhss_deriv_s2z, rhss_deriv_s2z = deriv_onevar(lhss_deriv, rhss_deriv, xprm=0, yprm=0, zprm=0, pxprm=0, pyprm=0,
                                                  pzprm=0, s1xprm=0, s1yprm=0, s1zprm=0, s2xprm=0, s2yprm=0, s2zprm=1)
    # Prepare to output derivative expressions in C syntax
    outstring = "/* SEOBNR Hamiltonian expression: */\n"
    outstringsp = ""
    outsplhs = []
    outsprhs = []
    for i in range(len(lr)):
        outstring += outputC(sp.sympify(lr[i].rhs), lr[i].lhs, "returnstring",
                             "outCverbose=False,includebraces=False,CSE_enable=False")
        outstringsp += lr[i].lhs + " = " + lr[i].rhs + "\n"
        outsplhs.append(sp.sympify(lr[i].lhs))
        outsprhs.append(sp.sympify(lr[i].rhs))
    outstring += "\n\n\n/* SEOBNR \partial_x H expression: */\n"
    for i in range(len(lhss_deriv_x)):
        outstring += outputC(rhss_deriv_x[i], str(lhss_deriv_x[i]), "returnstring",
                             "outCverbose=False,includebraces=False,CSE_enable=False")
        outstringsp += str(lhss_deriv_x[i]) + " = " + str(rhss_deriv_x[i]) + "\n"
        outsplhs.append(lhss_deriv_x[i])
        outsprhs.append(rhss_deriv_x[i])

    with open("SEOBNR_Playground_Pycodes/numpy_expressions.py", "w") as file:
        file.write("""from __future__ import division
import numpy as np
def compute_dHdq(m1, m2, eta, x, y, z, px, py, pz, s1x, s1y, s1z, s2x, s2y, s2z, KK, k0, k1, d1v2, dheffSSv2, tortoise, EMgamma):
""")
        for i in range(len(lr) - 1):
            file.write("    " + lr[i].lhs + " = " + str(lr[i].rhs).replace("Rational(",
                                                                           "np.true_divide(").replace("sqrt(",
                                                                                                      "np.sqrt(").replace(
                "log(", "np.log(").replace("sign(", "np.sign(").replace("Abs(", "np.abs(").replace("pi",
                                                                                                   "np.pi") + "\n")
        file.write("""    return Hreal""")
    # Playground agrees to this point... now need to figure out how to do CSE and output to different files!TylerK
    # Open and write to a file that, when called, will perform CSE on Hreal and all derivative expressions.
    with open("SEOBNR_Playground_Pycodes/sympy_expression.py", "w") as file:
        file.write("""import sympy as sp
from outputC import *
def sympy_cse():
    m1,m2,x,y,z,px,py,pz,s1x,s1y,s1z,s2x,s2y,s2z = sp.symbols("m1 m2 x y z px py pz s1x s1y s1z s2x s2y s2z",real=True)
    eta,KK,k0,k1,d1v2,dheffSSv2 = sp.symbols("eta KK k0 k1 d1v2 dheffSSv2",real=True)
    tortoise,EMgamma = sp.symbols("tortoise EMgamma",real=True)
""")
        # Convert Numpy functions in expressions to SymPy expressions.
        for i in range(len(lr)):
            file.write("    " + lr[i].lhs + " = " + str(lr[i].rhs).replace("Abs(", "sp.Abs(").replace("sqrt(",
                        "sp.sqrt(").replace("log(","sp.log(").replace("sign(", "sp.sign(").replace("Rational(",
                        "sp.Rational(").replace("pi", "sp.pi") + "\n")
        # Call the CSE routine for Hreal
        file.write("""
    CSE_results = sp.cse(Hreal, sp.numbered_symbols("Htmp"), order='canonical')
    with open("SEOBNR_Playground_Pycodes/numpy_expressions.py", "a") as file:
        for commonsubexpression in CSE_results[0]:
            file.write("    "+str(commonsubexpression[0])+" = "+str(commonsubexpression[1]).replace("sqrt(","np.sqrt(").replace("log(","np.log(").replace("sign(","np.sign(").replace("Abs(", "np.abs(").replace("pi", "np.pi")+"\\n")
        for i,result in enumerate(CSE_results[1]):
            file.write("    Hreal = "+str(result).replace("sqrt(","np.sqrt(")+"\\n")
""")
        # Declare all Hamiltonian terms as symbols so they can be used in derivative computations.
        for i in range(len(lr)):
            file.write("    " + lr[i].lhs + " = " + "sp.symbols(\"" + lr[i].lhs + "\")\n")

        # Print all terms for each of the partial derivatives.
        #TylerK: make this into a loop so there's not a buch of repeated code
        for i in range(len(lhss_deriv_x)):
            file.write("    " + str(lhss_deriv_x[i]).replace("prm", "prm_x") + " = " +
                       replace_numpy_funcs(rhss_deriv_x[i]).replace("prm", "prm_x") + "\n")
        for i in range(len(lhss_deriv_y)):
            file.write("    " + str(lhss_deriv_y[i]).replace("prm", "prm_y") + " = " +
                       replace_numpy_funcs(rhss_deriv_y[i]).replace("prm", "prm_y") + "\n")
        for i in range(len(lhss_deriv_z)):
            file.write("    " + str(lhss_deriv_z[i]).replace("prm", "prm_z") + " = " +
                       replace_numpy_funcs(rhss_deriv_z[i]).replace("prm", "prm_z") + "\n")

        for i in range(len(lhss_deriv_px)):
            file.write("    " + str(lhss_deriv_px[i]).replace("prm", "prm_px") + " = " +
                       replace_numpy_funcs(rhss_deriv_px[i]).replace("prm", "prm_px") + "\n")
        for i in range(len(lhss_deriv_py)):
            file.write("    " + str(lhss_deriv_py[i]).replace("prm", "prm_py") + " = " +
                       replace_numpy_funcs(rhss_deriv_py[i]).replace("prm", "prm_py") + "\n")
        for i in range(len(lhss_deriv_pz)):
            file.write("    " + str(lhss_deriv_pz[i]).replace("prm", "prm_pz") + " = " +
                       replace_numpy_funcs(rhss_deriv_pz[i]).replace("prm", "prm_pz") + "\n")

        for i in range(len(lhss_deriv_s1x)):
            file.write("    " + str(lhss_deriv_s1x[i]).replace("prm", "prm_s1x") + " = " +
                       replace_numpy_funcs(rhss_deriv_s1x[i]).replace("prm", "prm_s1x") + "\n")
        for i in range(len(lhss_deriv_s1y)):
            file.write("    " + str(lhss_deriv_s1y[i]).replace("prm", "prm_s1y") + " = " +
                       replace_numpy_funcs(rhss_deriv_s1y[i]).replace("prm", "prm_s1y") + "\n")
        for i in range(len(lhss_deriv_s1z)):
            file.write("    " + str(lhss_deriv_s1z[i]).replace("prm", "prm_s1z") + " = " +
                       replace_numpy_funcs(rhss_deriv_s1z[i]).replace("prm", "prm_s1z") + "\n")

        for i in range(len(lhss_deriv_s2x)):
            file.write("    " + str(lhss_deriv_s2x[i]).replace("prm", "prm_s2x") + " = " +
                       replace_numpy_funcs(rhss_deriv_s2x[i]).replace("prm", "prm_s2x") + "\n")
        for i in range(len(lhss_deriv_s2y)):
            file.write("    " + str(lhss_deriv_s2y[i]).replace("prm", "prm_s2y") + " = " +
                       replace_numpy_funcs(rhss_deriv_s2y[i]).replace("prm", "prm_s2y") + "\n")
        for i in range(len(lhss_deriv_s2z)):
            file.write("    " + str(lhss_deriv_s2z[i]).replace("prm", "prm_s2z") + " = " +
                       replace_numpy_funcs(rhss_deriv_s2z[i]).replace("prm", "prm_s2z") + "\n")
        # Perform CSE
        file.write("""
    output_list = ["Hrealprm_x","Hrealprm_y","Hrealprm_z","Hrealprm_px","Hrealprm_py","Hrealprm_pz",
        "Hrealprm_s1x","Hrealprm_s1y","Hrealprm_s1z","Hrealprm_s2x","Hrealprm_s2y","Hrealprm_s2z"]
    expression_list = [Hrealprm_x,Hrealprm_y,Hrealprm_z,Hrealprm_px,Hrealprm_py,Hrealprm_pz,
        Hrealprm_s1x,Hrealprm_s1y,Hrealprm_s1z,Hrealprm_s2x,Hrealprm_s2y,Hrealprm_s2z]
    CSE_results = sp.cse(expression_list, sp.numbered_symbols("tmp"), order='canonical')
    with open("SEOBNR_Playground_Pycodes/numpy_expressions.py", "a") as file:
        for commonsubexpression in CSE_results[0]:
            file.write("    "+str(commonsubexpression[0])+" = "+str(commonsubexpression[1]).replace("sqrt(","np.sqrt(").replace("log(","np.log(").replace("sign(","np.sign(").replace("Abs(", "np.abs(")+"\\n")
        for i,result in enumerate(CSE_results[1]):
            file.write("    "+str(output_list[i])+" = "+str(result).replace("sqrt(","np.sqrt(").replace("log(","np.log(").replace("sign(","np.sign(").replace("Abs(", "np.abs(")+"\\n")
        for i,result in enumerate(CSE_results[1]):
            if i > 0:
                file.write(","+str(output_list[i]))
            else:
                file.write("    return np.array([Hreal,"+str(output_list[i]))
        file.write("])")
""")
