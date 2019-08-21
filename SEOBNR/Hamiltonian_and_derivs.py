from outputC import *
import sympy as sp


def simplify_deriv(lhss_deriv, rhss_deriv):
    lhss_deriv_simp = []
    rhss_deriv_simp = []
    for i in range(len(rhss_deriv)):
        lhss_deriv_simp.append(lhss_deriv[i])
        rhss_deriv_simp.append(rhss_deriv[i])
    for i in range(len(rhss_deriv_simp)):
        if rhss_deriv_simp[i] == 0:
            for j in range(i+1, len(rhss_deriv_simp)):
                for var in rhss_deriv_simp[j].free_symbols:
                    if str(var) == str(lhss_deriv_simp[i]):
                        rhss_deriv_simp[j] = rhss_deriv_simp[j].subs(var, 0)
    zero_elements_to_remove = []
    for i in range(len(rhss_deriv_simp)):
        if rhss_deriv_simp[i] == sp.sympify(0):
            zero_elements_to_remove.append(i)
    count = 0
    for i in range(len(zero_elements_to_remove)):
        del lhss_deriv_simp[zero_elements_to_remove[i]+count]
        del rhss_deriv_simp[zero_elements_to_remove[i]+count]
        count -= 1
    return lhss_deriv_simp, rhss_deriv_simp


def deriv_onevar(lhss_deriv, rhss_deriv, xprm=0, yprm=0, zprm=0, pxprm=0, pyprm=0, pzprm=0, s1xprm=0, s1yprm=0,
                 s1zprm=0, s2xprm=0, s2yprm=0, s2zprm=0):
    lhss_deriv_new = []
    rhss_deriv_new = []
    for i in range(len(rhss_deriv)):
        lhss_deriv_new.append(lhss_deriv[i])
        rhss_deriv_new.append(rhss_deriv[i])
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
    lhss_deriv_simp, rhss_deriv_simp = simplify_deriv(lhss_deriv_new, rhss_deriv_new)
    return lhss_deriv_simp, rhss_deriv_simp


def replace_numpy_funcs(expression):

    return str(expression).replace("sqrt(", "sp.sqrt(").replace("log(",
                                            "sp.log(").replace("Abs(", "sp.Abs(").replace("sign(", "sp.sign(")


def output_H_and_derivs():

    f = open("Hamstring.txt", 'r')
    Hamstring = str(f.read())
    f.close()

    # Split Hamstring by carriage returns:
    Hamterms = Hamstring.splitlines()

    # Create "lr" array, which will store each left-hand side and right-hand side of Hamterms as strings.
    lr = []
    # Loop over each line in Hamstring to separate the left- and right-hand sides.
    for i in range(len(Hamterms)):
        # Ignore lines with 2 or fewer characters and those starting with #
        if len(Hamterms[i]) > 2 and Hamterms[i][0] != "#":
            # Split each line by its equals sign.
            splitHamterms = Hamterms[i].split("=")
            # Append to the "lr" array, removing spaces, "sp." prefixes, and replacing Lambda->Lamb
            #       (Lambda is a protected keyword):
            lr.append(lhrh(lhs=splitHamterms[0].replace(" ", "").replace("Lambda", "Lamb"),
                           rhs=splitHamterms[1].replace(" ", "").replace("sp.A", "a").replace("sp.",
                                                            "").replace("Lambda", "Lamb")))

    xx = sp.Symbol('xx')
    func = []
    lhss = []
    rhss = []
    # Affix '(xx)' to each left-hand side as a function designation.
    for i in range(len(lr)):
        func.append(sp.sympify(sp.Function(lr[i].lhs)(xx)))
        lhss.append(sp.sympify(lr[i].lhs))
        rhss.append(sp.sympify(lr[i].rhs))

    # Generate a list of all the "free symbols" in the RHS expressions.
    full_symbol_list_with_dups = []
    for i in range(len(lr)):
        for var in rhss[i].free_symbols:
            full_symbol_list_with_dups.append(var)

    # Remove all duplicated "free symbols" from the RHS expressions.
    full_symbol_list = superfast_uniq(full_symbol_list_with_dups)

    # Declare input constants.
    m1, m2, eta = sp.symbols("m1 m2 eta", real=True)
    c0k2, c1k2, c0k3, c1k3, c0k4, c1k4, c2k4, c0k5, c1k5, c2k5 = sp.symbols(
        "c0k2 c1k2 c0k3 c1k3 c0k4 c1k4 c2k4 c0k5 c1k5 c2k5", real=True)
    KK, k5l, b3, bb3, d1, d1v2, dheffSS, dheffSSv2 = sp.symbols("KK k5l b3 bb3 d1 d1v2 dheffSS dheffSSv2", real=True)
    tortoise = sp.symbols("tortoise", real=True)
    input_constants = [m1, m2, eta, c0k2, c1k2, c0k3, c1k3, c0k4, c1k4, c2k4, c0k5, c1k5, c2k5, KK, k5l, b3, bb3, d1,
                       d1v2, dheffSS, dheffSSv2, tortoise]

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

    # Differentiate with respect to xx, remove '(xx)', and replace xx with 'prm' notation
    lhss_deriv = []
    rhss_deriv = []
    for i in range(len(rhss)):
        lhss_deriv.append(sp.sympify(str(lhss[i]) + "prm"))
        newrhs = sp.sympify(
            str(sp.diff(rhss[i], xx)).replace("(xx)", "").replace(", xx", "prm").replace("Derivative", ""))
        rhss_deriv.append(newrhs)

    lhss_deriv_simp, rhss_deriv_simp = simplify_deriv(lhss_deriv, rhss_deriv)
    lhss_deriv = lhss_deriv_simp
    rhss_deriv = rhss_deriv_simp

    lhss_deriv_x, rhss_deriv_x = deriv_onevar(lhss_deriv, rhss_deriv,
                                              xprm=1, yprm=0, zprm=0, pxprm=0, pyprm=0, pzprm=0, s1xprm=0, s1yprm=0,
                                              s1zprm=0, s2xprm=0, s2yprm=0, s2zprm=0)
    lhss_deriv_y, rhss_deriv_y = deriv_onevar(lhss_deriv, rhss_deriv,
                                              xprm=0, yprm=1, zprm=0, pxprm=0, pyprm=0, pzprm=0, s1xprm=0, s1yprm=0,
                                              s1zprm=0, s2xprm=0, s2yprm=0, s2zprm=0)
    lhss_deriv_z, rhss_deriv_z = deriv_onevar(lhss_deriv, rhss_deriv,
                                              xprm=0, yprm=0, zprm=1, pxprm=0, pyprm=0, pzprm=0, s1xprm=0, s1yprm=0,
                                              s1zprm=0, s2xprm=0, s2yprm=0, s2zprm=0)
    lhss_deriv_px, rhss_deriv_px = deriv_onevar(lhss_deriv, rhss_deriv,
                                                xprm=0, yprm=0, zprm=0, pxprm=1, pyprm=0, pzprm=0, s1xprm=0, s1yprm=0,
                                                s1zprm=0, s2xprm=0, s2yprm=0, s2zprm=0)
    lhss_deriv_py, rhss_deriv_py = deriv_onevar(lhss_deriv, rhss_deriv,
                                                xprm=0, yprm=0, zprm=0, pxprm=0, pyprm=1, pzprm=0, s1xprm=0, s1yprm=0,
                                                s1zprm=0, s2xprm=0, s2yprm=0, s2zprm=0)
    lhss_deriv_pz, rhss_deriv_pz = deriv_onevar(lhss_deriv, rhss_deriv,
                                                xprm=0, yprm=0, zprm=0, pxprm=0, pyprm=1, pzprm=1, s1xprm=0, s1yprm=0,
                                                s1zprm=0, s2xprm=0, s2yprm=0, s2zprm=0)
    lhss_deriv_s1x, rhss_deriv_s1x = deriv_onevar(lhss_deriv, rhss_deriv,
                                                  xprm=0, yprm=0, zprm=0, pxprm=0, pyprm=0, pzprm=0, s1xprm=1, s1yprm=0,
                                                  s1zprm=0, s2xprm=0, s2yprm=0, s2zprm=0)
    lhss_deriv_s1y, rhss_deriv_s1y = deriv_onevar(lhss_deriv, rhss_deriv,
                                                  xprm=0, yprm=0, zprm=0, pxprm=0, pyprm=0, pzprm=0, s1xprm=0, s1yprm=1,
                                                  s1zprm=0, s2xprm=0, s2yprm=0, s2zprm=0)
    lhss_deriv_s1z, rhss_deriv_s1z = deriv_onevar(lhss_deriv, rhss_deriv,
                                                  xprm=0, yprm=0, zprm=0, pxprm=0, pyprm=0, pzprm=0, s1xprm=0, s1yprm=0,
                                                  s1zprm=1, s2xprm=0, s2yprm=0, s2zprm=0)
    lhss_deriv_s2x, rhss_deriv_s2x = deriv_onevar(lhss_deriv, rhss_deriv,
                                                  xprm=0, yprm=0, zprm=0, pxprm=0, pyprm=0, pzprm=0, s1xprm=0, s1yprm=0,
                                                  s1zprm=0, s2xprm=1, s2yprm=0, s2zprm=0)
    lhss_deriv_s2y, rhss_deriv_s2y = deriv_onevar(lhss_deriv, rhss_deriv,
                                                  xprm=0, yprm=0, zprm=0, pxprm=0, pyprm=0, pzprm=0, s1xprm=0, s1yprm=0,
                                                  s1zprm=0, s2xprm=0, s2yprm=1, s2zprm=0)
    lhss_deriv_s2z, rhss_deriv_s2z = deriv_onevar(lhss_deriv, rhss_deriv,
                                                  xprm=0, yprm=0, zprm=0, pxprm=0, pyprm=0, pzprm=0, s1xprm=0, s1yprm=0,
                                                  s1zprm=0, s2xprm=0, s2yprm=0, s2zprm=1)

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

    with open("numpy_expressions.py", "w") as file:
        file.write("import numpy as np\n")
        file.write("def compute_dHdq(x, y, z, px, py, pz, s1x, s1y, s1z, s2x, s2y, s2z):\\n")
        for i in range(len(lr)-1):
            file.write("\\t" + lr[i].lhs + " = " + str(lr[i].rhs).replace("Abs(", "abs(").replace("Rational(","divide(") + "\n")

    with open("sympy_expression.py", "w") as file:
        file.write("""
import sympy as sp
from outputC import *
m1,m2,x,y,z,px,py,pz,s1x,s1y,s1z,s2x,s2y,s2z = sp.symbols("m1 m2 x y z px py pz s1x s1y s1z s2x s2y s2z",real=True)
c0k2,c1k2,c0k3,c1k3,c0k4 = sp.symbols("c0k2 c1k2 c0k3 c1k3 c0k4",real=True)
c1k4,c2k4,c0k5,c1k5,c2k5 = sp.symbols("c1k4 c2k4 c0k5 c1k5 c2k5",real=True)
eta,KK,k5l,b3,bb3,d1,d1v2,dheffSS,dheffSSv2 = sp.symbols("eta KK k5l b3 bb3 d1 d1v2 dheffSS dheffSSv2",real=True)
tortoise = sp.symbols("tortoise",real=True)
""")

        for i in range(len(lr)):
            file.write(lr[i].lhs + " = " + str(lr[i].rhs).replace("sqrt(","sp.sqrt(").replace("log(",
                "sp.log(").replace("Abs(", "sp.Abs(").replace("sign(", "sp.sign(").replace("Rational(",
                "sp.Rational(") + "\n")
        file.write("""
CSE_results = sp.cse(Hreal, sp.numbered_symbols("Htmp"), order='canonical')
with open("numpy_expressions.py", "a") as file:
    for commonsubexpression in CSE_results[0]:
        file.write("\\t"+str(commonsubexpression[0])+" = "+str(commonsubexpression[1]).replace("Abs", "abs")+"\\n")
    for i,result in enumerate(CSE_results[1]):
        file.write("\\tHreal = "+str(result)+"\\n")
""")

        for i in range(len(lr)):
            file.write(lr[i].lhs + " = " + "sp.symbols(\"" + lr[i].lhs + "\")\n")

        for i in range(len(lhss_deriv_x)):
#            print(rhss_deriv_x[i])
#            replace_numpy_funcs(rhss_deriv_x[i]).replace("prm", "prm_x")
            file.write(str(lhss_deriv_x[i]).replace("prm", "prm_x") + " = " +
                       replace_numpy_funcs(rhss_deriv_x[i]).replace("prm", "prm_x") + "\n")
#        for i in range(len(lhss_deriv_x)):
#            file.write(str(lhss_deriv_x[i]).replace("prm", "prm_x") + " = " + str(rhss_deriv_x[i]).replace("sqrt(",
#                "sp.sqrt(").replace("log(", "sp.log(").replace("Abs(", "sp.Abs(").replace("sign(", "sp.sign(").replace(
#                "prm", "prm_x") + "\n")
        for i in range(len(lhss_deriv_y)):
            file.write(str(lhss_deriv_y[i]).replace("prm", "prm_y") + " = " +
                       replace_numpy_funcs(rhss_deriv_y[i]).replace("prm", "prm_y") + "\n")
        for i in range(len(lhss_deriv_z)):
            file.write(str(lhss_deriv_z[i]).replace("prm", "prm_z") + " = " +
                       replace_numpy_funcs(rhss_deriv_z[i]).replace("prm", "prm_z") + "\n")

        for i in range(len(lhss_deriv_px)):
            file.write(str(lhss_deriv_px[i]).replace("prm", "prm_px") + " = " +
                       replace_numpy_funcs(rhss_deriv_px[i]).replace("prm", "prm_px") + "\n")
        for i in range(len(lhss_deriv_py)):
            file.write(str(lhss_deriv_py[i]).replace("prm", "prm_py") + " = " +
                       replace_numpy_funcs(rhss_deriv_py[i]).replace("prm", "prm_py") + "\n")
        for i in range(len(lhss_deriv_pz)):
            file.write(str(lhss_deriv_pz[i]).replace("prm", "prm_pz") + " = " +
                       replace_numpy_funcs(rhss_deriv_pz[i]).replace("prm", "prm_pz") + "\n")

        for i in range(len(lhss_deriv_s1x)):
            file.write(str(lhss_deriv_s1x[i]).replace("prm", "prm_s1x") + " = " +
                       replace_numpy_funcs(rhss_deriv_s1x[i]).replace("prm", "prm_s1x") + "\n")
        for i in range(len(lhss_deriv_s1y)):
            file.write(str(lhss_deriv_s1y[i]).replace("prm", "prm_s1y") + " = " +
                       replace_numpy_funcs(rhss_deriv_s1y[i]).replace("prm", "prm_s1y") + "\n")
        for i in range(len(lhss_deriv_s1z)):
            file.write(str(lhss_deriv_s1z[i]).replace("prm", "prm_s1z") + " = " +
                       replace_numpy_funcs(rhss_deriv_s1z[i]).replace("prm", "prm_s1z") + "\n")

        for i in range(len(lhss_deriv_s2x)):
            file.write(str(lhss_deriv_s2x[i]).replace("prm", "prm_s2x") + " = " +
                       replace_numpy_funcs(rhss_deriv_s2x[i]).replace("prm", "prm_s2x") + "\n")
        for i in range(len(lhss_deriv_s2y)):
            file.write(str(lhss_deriv_s2y[i]).replace("prm", "prm_s2y") + " = " +
                       replace_numpy_funcs(rhss_deriv_s2y[i]).replace("prm", "prm_s2y") + "\n")
        for i in range(len(lhss_deriv_s2z)):
            file.write(str(lhss_deriv_s2z[i]).replace("prm", "prm_s2z") + " = " +
                       replace_numpy_funcs(rhss_deriv_s2z[i]).replace("prm", "prm_s2z") + "\n")
        file.write("""
output_list = ["Hrealprm_x","Hrealprm_y","Hrealprm_z","Hrealprm_px","Hrealprm_py","Hrealprm_pz",
         "Hrealprm_s1x","Hrealprm_s1y","Hrealprm_s1z","Hrealprm_s2x","Hrealprm_s2y","Hrealprm_s2z"]
expression_list = [Hrealprm_x,Hrealprm_y,Hrealprm_z,Hrealprm_px,Hrealprm_py,Hrealprm_pz,
         Hrealprm_s1x,Hrealprm_s1y,Hrealprm_s1z,Hrealprm_s2x,Hrealprm_s2y,Hrealprm_s2z]
CSE_results = sp.cse(expression_list, sp.numbered_symbols("tmp"), order='canonical')
with open("numpy_expressions.py", "a") as file:
    for commonsubexpression in CSE_results[0]:
        file.write("\\t"+str(commonsubexpression[0])+" = "+str(commonsubexpression[1]).replace("Abs", "abs")+"\\n")
    for i,result in enumerate(CSE_results[1]):
        file.write("\\t"+str(output_list[i])+" = "+str(result)+"\\n")
    file.write("\\treturn np.array([")
    for i,result in enumerate(CSE_results[1]):
        if i > 0:
            file.write(","+str(output_list[i]))
        else:
            file.write("\\treturn np.array([Hreal,"+str(output_list[i]))
    file.write("])")
""")