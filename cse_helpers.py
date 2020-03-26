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
import sys

# Input:  expr_list = single SymPy expression or list of SymPy expressions
#         prefix    = string prefix for variable names (i.e. rational symbols)
#         declare   = declare symbol for negative one (i.e. _NegativeOne_)
#         factor    = perform partial factorization (excluding negative)
#         negative  = include negative in partial factorization
# Output: modified SymPy expression(s) where all integers and rationals were replaced
#           with temporary placeholder variables that allow for partial factorization
def cse_preprocess(expr_list, prefix='', declare=False, factor=True, negative=False, debug=False):
    if not isinstance(expr_list, list):
        expr_list = [expr_list]
    _NegativeOne_ = sp.Symbol(prefix + '_NegativeOne_')
    map_sym_to_rat, map_rat_to_sym = {}, {}
    for i, expr in enumerate(expr_list):
        tree = ExprTree(expr)
        # Search through expression tree for rational(s)
        for subtree in tree.preorder():
            subexpr = subtree.expr
            if isinstance(subexpr, sp.Rational) and subexpr != sp.S.NegativeOne:
                # Ignore replacing exponent of power function with symbol
                if subtree.func == sp.Pow: continue
                # If rational < 0, factor out negative, leaving positive rational
                sign = 1 if subexpr >= 0 else -1
                subexpr *= sign
                # Declare unique symbol for rational on first appearance
                try: repl = map_rat_to_sym[subexpr]
                except KeyError:
                    p, q = subexpr.p, subexpr.q
                    var_name = prefix + '_Rational_' + str(p) + '_' + str(q) \
                        if q != 1 else prefix + '_Integer_' + str(p)
                    repl = sp.Symbol(var_name)
                    map_sym_to_rat[repl], map_rat_to_sym[subexpr] = subexpr, repl
                subtree.expr = repl * sign
                if sign < 0: tree.build(subtree)
            # If declare == True, then declare symbol for negative one
            elif declare == True and subexpr == sp.S.NegativeOne:
                try: subtree.expr = map_rat_to_sym[sp.S.NegativeOne]
                except KeyError:
                    map_sym_to_rat[_NegativeOne_], map_rat_to_sym[subexpr] = subexpr, _NegativeOne_
                    subtree.expr = _NegativeOne_
        expr = tree.reconstruct()
        # If factor == True, then perform partial factoring (excluding _NegativeOne_)
        if factor == True:
            # Handle the separate case of function argument(s)
            for subtree in tree.preorder():
                if isinstance(subtree.expr, sp.Function):
                    for var in map_sym_to_rat:
                        if var != _NegativeOne_:
                            arg = subtree.children[0]
                            arg.expr = sp.collect(arg.expr, var)
                            arg.children.clear()
            expr = tree.reconstruct()
            # Perform partial factoring on expression(s)
            for var in map_sym_to_rat:
                if var != _NegativeOne_:
                    expr = sp.collect(expr, var)
            tree.root.expr = expr
            tree.build(tree.root, clear=True)
        # If negative == True, then perform partial factoring on _NegativeOne_
        if negative == True:
            for subtree in tree.preorder():
                if isinstance(subtree.expr, sp.Function):
                    arg = subtree.children[0]
                    arg.expr = sp.collect(arg.expr, _NegativeOne_)
                    arg.children.clear()
            expr = sp.collect(tree.reconstruct(), _NegativeOne_)
            tree.root.expr = expr
            tree.build(tree.root, clear=True)
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
        # Replace any left-over one(s) after partial factoring
        if factor == True or negative == True:
            _One_ = sp.Symbol(prefix + '_Integer_1')
            for subtree in tree.preorder():
                if subtree.expr == sp.S.One:
                    subtree.expr = _One_
            tmp_expr = tree.reconstruct()
            if tmp_expr != expr:
                try: map_rat_to_sym[sp.S.One]
                except KeyError:
                    map_sym_to_rat[_One_], map_rat_to_sym[sp.S.One] = sp.S.One, _One_
                    subtree.expr = _One_
                expr = tmp_expr
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
