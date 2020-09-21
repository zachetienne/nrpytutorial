""" CSE Partial Factorization and Post-Processing

    The following script will perform partial factorization on SymPy expressions,
    which should occur before common subexpression elimination (CSE) to prevent the
    identification of undesirable patterns, and perform post-processing on the
    the resulting replaced/reduced expressions after the CSE procedure was applied.
"""
# Author: Ken Sible
# Email:  ksible *at* outlook *dot* com

from expr_tree import ExprTree  # NRPy+: Contains expression tree data structure class definitions and manipulation functions
import sympy as sp              # SymPy: The Python computer algebra package upon which NRPy+ depends
import sys                      # Standard Python module for multiplatform OS-level functions
from collections import OrderedDict

def cse_preprocess(expr_list, prefix='', declare=False, factor=True, negative=False, debug=False):
    """ Perform CSE Preprocessing

        :arg:    single SymPy expression or list of SymPy expressions
        :arg:    string prefix for variable names (i.e. rational symbols)
        :arg:    declare symbol for negative one (i.e. _NegativeOne_)
        :arg:    perform partial factorization (excluding negative symbol)
        :arg:    include negative symbol in partial factorization
        :arg:    back-substitute and check difference for debugging
        :return: modified SymPy expression(s) where all integers and rationals were replaced
                    with temporary placeholder variables that allow for partial factorization

        >>> from sympy.abc import x, y, z
        >>> expr = -x/12 - y/12 + z
        >>> cse_preprocess(expr)
        (_Rational_1_12*(-x - y) + z, OrderedDict([(_Rational_1_12, 1/12)]))

        >>> cse_preprocess(expr, declare=True)
        (_Rational_1_12*(_NegativeOne_*x + _NegativeOne_*y) + z, OrderedDict([(_Rational_1_12, 1/12), (_NegativeOne_, -1)]))

        >>> expr = -x/12 - y/12 + z
        >>> cse_preprocess(expr, declare=True, negative=True)
        (_NegativeOne_*_Rational_1_12*(x + y) + z, OrderedDict([(_Rational_1_12, 1/12), (_NegativeOne_, -1)]))

        >>> cse_preprocess(expr, factor=False)
        ((-_Rational_1_12)*x + (-_Rational_1_12)*y + z, OrderedDict([(_Rational_1_12, 1/12)]))

        >>> cse_preprocess(expr, prefix='FD')
        (FD_Rational_1_12*(-x - y) + z, OrderedDict([(FD_Rational_1_12, 1/12)]))

        >>> from sympy import exp
        >>> expr = exp(3*x + 3*y)
        >>> cse_preprocess(expr)
        (exp(_Integer_3*(x + y)), OrderedDict([(_Integer_3, 3)]))

        >>> from sympy import Mul
        >>> expr = Mul((-1)**3, (3*x + 3*y), evaluate=False)
        >>> cse_preprocess(expr, declare=True)
        (_Integer_3*_NegativeOne_*(x + y), OrderedDict([(_NegativeOne_, -1), (_Integer_3, 3)]))
    """
    if not isinstance(expr_list, list):
        expr_list = [expr_list]
    expr_list = expr_list[:]
    _NegativeOne_ = sp.Symbol(prefix + '_NegativeOne_')
    map_sym_to_rat, map_rat_to_sym = OrderedDict(), OrderedDict()
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
                if sign < 0: tree.build(subtree, clear=False)
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
                    arg = subtree.children[0]
                    for var in map_sym_to_rat:
                        if var != _NegativeOne_:
                            arg.expr = sp.collect(arg.expr, var)
                    tree.build(arg)
            expr = tree.reconstruct()
            # Perform partial factoring on expression(s)
            for var in map_sym_to_rat:
                if var != _NegativeOne_:
                    expr = sp.collect(expr, var)
            tree.root.expr = expr
            tree.build(tree.root)
        # If negative == True, then perform partial factoring on _NegativeOne_
        if negative == True:
            for subtree in tree.preorder():
                if isinstance(subtree.expr, sp.Function):
                    arg = subtree.children[0]
                    arg.expr = sp.collect(arg.expr, _NegativeOne_)
                    tree.build(arg)
            expr = sp.collect(tree.reconstruct(), _NegativeOne_)
            tree.root.expr = expr
            tree.build(tree.root)
        # If declare == True, then simplify (-1)^n
        if declare == True:
            _One_ = sp.Symbol(prefix + '_Integer_1')
            for subtree in tree.preorder():
                subexpr = subtree.expr
                if subexpr.func == sp.Pow:
                    base, expo = subexpr.args[0], subexpr.args[1]
                    if base == _NegativeOne_:
                        subtree.expr = _One_ if expo % 2 == 0 else _NegativeOne_
                        tree.build(subtree)
            expr = tree.reconstruct()
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

def cse_postprocess(cse_output):
    """ Perform CSE Postprocessing

        :arg:    output from SymPy CSE with tuple format: (list of ordered pairs that
                    contain substituted symbols and their replaced expressions, reduced SymPy expression)
        :return: output from SymPy CSE where postprocessing, such as back-substitution of addition/product
                    of symbols, has been applied to the replaced/reduced expression(s)

        >>> from sympy.abc import x, y
        >>> from sympy import cse, cos, sin

        >>> cse_out = cse(3 + x + cos(3 + x))
        >>> cse_postprocess(cse_out)
        ([], [x + cos(x + 3) + 3])

        >>> cse_out = cse(3 + x + y + cos(3 + x + y))
        >>> cse_postprocess(cse_out)
        ([(x0, x + y + 3)], [x0 + cos(x0)])

        >>> cse_out = cse(3*x + cos(3*x))
        >>> cse_postprocess(cse_out)
        ([], [3*x + cos(3*x)])

        >>> cse_out = cse(3*x*y + cos(3*x*y))
        >>> cse_postprocess(cse_out)
        ([(x0, 3*x*y)], [x0 + cos(x0)])

        >>> cse_out = cse(x**2 + cos(x**2))
        >>> cse_postprocess(cse_out)
        ([], [x**2 + cos(x**2)])

        >>> cse_out = cse(x**3 + cos(x**3))
        >>> cse_postprocess(cse_out)
        ([(x0, x**3)], [x0 + cos(x0)])

        >>> cse_out = cse(x*y + cos(x*y) + sin(x*y))
        >>> cse_postprocess(cse_out)
        ([(x0, x*y)], [x0 + sin(x0) + cos(x0)])

        >>> from sympy import exp, log
        >>> expr = -x + exp(-x) + log(-x)
        >>> cse_pre = cse_preprocess(expr, declare=True)
        >>> cse_out = cse(cse_pre[0])
        >>> cse_postprocess(cse_out)
        ([], [_NegativeOne_*x + exp(_NegativeOne_*x) + log(_NegativeOne_*x)])
    """
    replaced, reduced = cse_output
    replaced, reduced = replaced[:], reduced[:]
    i = 0
    while i < len(replaced):
        sym, expr = replaced[i]; args = expr.args
        # Search through replaced expressions for negative symbols
        if (expr.func == sp.Mul and len(expr.args) == 2 and any(a1.func == sp.Symbol and \
               (a2 == sp.S.NegativeOne or '_NegativeOne_' in str(a2)) for a1, a2 in [args, reversed(args)])):
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

if __name__ == "__main__":
    import doctest
    sys.exit(doctest.testmod()[0])
