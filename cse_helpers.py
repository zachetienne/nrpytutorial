""" CSE Partial Factorization and Post-Processing

The following script will perform partial factorization on SymPy expressions,
which should occur before common subexpression elimination (CSE) to prevent the
identification of undesirable patterns, and perform post-processing on the
the resulting replaced/reduced expressions after the CSE procedure was applied.
"""
# Author: Ken Sible
# Email:  ksible@outlook.com

from SIMDExprTree import ExprTree
import sympy as sp

# Input:  expr_list = single SymPy expression or list of SymPy expressions
#         prefix    = string prefix for variable names (i.e. replacement symbols)
#         declare   = declare negative one symbol (i.e. _NegativeOne_)
#         negative  = perform partial factorization on negative one
#         factor    = perform partial factorization (excluding negative one)
# Output: modified SymPy expression(s) where all integers and rationals were replaced
#           with temporary placeholder variables that allow for partial factorization
def cse_preprocess(expr_list, prefix='', declare=False, negative=False, factor=True, debug=False):
    if not isinstance(expr_list, list):
        expr_list = [expr_list]
    def expand(a, n):
        if  n == 2:  return sp.Mul(a, a, evaluate=False)
        elif n > 2:  return sp.Mul(expand(a, n - 1), a, evaluate=False)
        return sp.Pow(expand(a, -n), -1, evaluate=False)
    _NegativeOne_ = sp.Symbol(prefix + '_NegativeOne_')
    map_sym_to_rat, map_rat_to_sym = {}, {}
    for i, expr in enumerate(expr_list):
        tree = ExprTree(expr)
        # Expand power function, preventing replacement of exponent argument
        for subtree in tree.preorder(tree.root):
            subexpr = subtree.expr
            if subexpr.func == sp.Pow:
                exponent = subtree.children[1].expr
                if exponent.func == sp.Integer and abs(exponent) > 1:
                    subtree.expr = expand(*subexpr.args)
                    tree.build(subtree, clear=True)
        # Search through expression tree for integers/rationals
        for subtree in tree.preorder():
            subexpr = subtree.expr
            if isinstance(subexpr, sp.Rational) and subexpr != sp.S.NegativeOne:
                # If rational < 0, factor out negative and declare positive rational
                sign = 1 if subexpr >= 0 else -1
                subexpr *= sign
                # Check whether rational was already declared, otherwise declare rational
                try: repl = map_rat_to_sym[subexpr]
                except KeyError:
                    p, q = subexpr.p, subexpr.q
                    var_name = prefix + '_Rational_' + str(p) + '_' + str(q) \
                        if q != 1 else prefix + '_Integer_' + str(p)
                    repl = sp.Symbol(var_name)
                    map_sym_to_rat[repl], map_rat_to_sym[subexpr] = subexpr, repl
                subtree.expr = repl * sign
                if sign < 0: tree.build(subtree, clear=True)
            # If declare == True, then declare symbol for -1 or extracted negative
            elif declare and subexpr == sp.S.NegativeOne:
                try: subtree.expr = map_rat_to_sym[sp.S.NegativeOne]
                except KeyError:
                    repl = _NegativeOne_
                    map_sym_to_rat[repl], map_rat_to_sym[subexpr] = subexpr, repl
                    subtree.expr = repl
        # If exponent was replaced with symbol (usually -1), then back-substitute
        for subtree in tree.preorder(tree.root):
            subexpr = subtree.expr
            if subexpr.func == sp.Pow:
                exponent = subtree.children[1].expr
                if exponent.func == sp.Symbol:
                    subtree.children[1].expr = map_sym_to_rat[exponent]
        expr = tree.reconstruct()
        # If factor == True, then perform partial factoring (excluding negative one)
        if factor == True:
            # Handle the separate case of function arguments
            for subtree in tree.preorder():
                if isinstance(subtree.expr, sp.Function):
                    for var in map_sym_to_rat:
                        if var != _NegativeOne_:
                            child = subtree.children[0]
                            child.expr = sp.collect(child.expr, var)
                            child.children.clear()
            expr = tree.reconstruct()
            # Perform partial factoring on the expression(s)
            for var in map_sym_to_rat:
                if var != _NegativeOne_:
                    expr = sp.collect(expr, var)
        # If negative == True, then perform partial factoring on negative one
        if negative == True:
            expr = sp.collect(expr, _NegativeOne_)
        # If debug == True, then back-substitute everything and check difference
        if debug == True:
            def lookup_rational(arg):
                if arg.func == sp.Symbol:
                    try: arg = map_sym_to_rat[arg]
                    except KeyError: pass
                return arg
            debug_tree = ExprTree(expr)
            for subtree in debug_tree.preorder():
                subexpr = subtree.expr
                if subexpr.func == sp.Symbol:
                    subtree.expr = lookup_rational(subexpr)
            debug_expr = tree.reconstruct()
            expr_diff  = expr - debug_expr
            if sp.simplify(expr_diff) != 0:
                raise Warning('Expression Difference: ' + str(expr_diff))
        expr_list[i] = expr
    if len(expr_list) == 1:
        expr_list = expr_list[0]
    return expr_list, map_sym_to_rat

# Input:  cse_output = output from SymPy CSE with tuple format: (list of ordered pairs that 
#            contain substituted symbols and their replaced expressions, reduced SymPy expression)
# Output: output from SymPy CSE where postprocessing, such as back-substitution of addition/product
#            of symbols, has been applied to the replaced/reduced expression(s)
def cse_postprocess(cse_output):
    replaced, reduced = cse_output
    i = 0
    while i < len(replaced):
        sym, expr = replaced[i]
        # Search through replaced expressions for negative symbols
        if (expr.func == sp.Mul and len(expr.args) == 2 and \
                any((arg.func == sp.Symbol) for arg in expr.args) and \
                any((arg == sp.S.NegativeOne or '_NegativeOne_' in str(arg)) for arg in expr.args)):
            for k in range(i + 1, len(replaced)):
                if sym in replaced[k][1].free_symbols:
                    replaced[k] = (replaced[k][0], replaced[k][1].subs(sym, expr))
            for k in range(len(reduced)):
                if sym in reduced[k].free_symbols:
                    reduced[k] = reduced[k].subs(sym, expr)
            # Remove the replaced expression from the list
            replaced.pop(i)
            if i != 0: i -= 1
        # Search through replaced expressions for addition/product of 2 or less symbols
        if ((expr.func == sp.Add or expr.func == sp.Mul) and 0 < len(expr.args) < 3 and \
                all((arg.func == sp.Symbol or arg.is_integer or arg.is_rational) for arg in expr.args)) or \
                (expr.func == sp.Pow and expr.args[0].func == sp.Symbol and expr.args[1] == 2):
            sym_count = 0 # Count the number of occurrences of the substituted symbol
            for k in range(len(replaced) - i):
                # Check if the substituted symbol appears in the replaced expressions
                if sym in replaced[i + k][1].free_symbols:
                    for arg in sp.preorder_traversal(replaced[i + k][1]):
                        if arg.func == sp.Symbol and str(arg) == str(sym):
                            sym_count += 1
            for k in range(len(reduced)):
                # Check if the substituted symbol appears in the reduced expression
                if sym in reduced[k].free_symbols:
                    for arg in sp.preorder_traversal(reduced[k]):
                        if arg.func == sp.Symbol and str(arg) == str(sym):
                            sym_count += 1
            # If the number of occurrences of the substituted symbol is 2 or less, back-substitute
            if 0 < sym_count < 3:
                for k in range(i + 1, len(replaced)):
                    if sym in replaced[k][1].free_symbols:
                        replaced[k] = (replaced[k][0], replaced[k][1].subs(sym, expr))
                for k in range(len(reduced)):
                    if sym in reduced[k].free_symbols:
                        reduced[k] = reduced[k].subs(sym, expr)
                # Remove the replaced expression from the list
                replaced.pop(i); i -= 1
        i += 1
    return replaced, reduced
