from SIMDExprTree import ExprTree
import sympy as sp

# Input:  sympyexpr = a single SymPy expression or list of SymPy expressions
# Output: modified SymPy expression(s) where all integers and rationals are replaced
#           with temporary placeholder variables that allow for partial factorization
def cse_preprocess(expr_list, prefix='', ignore=False, factor=True):
    def expand(a, n):
        if  n == 2:  return sp.Mul(a, a, evaluate=False)
        elif n > 2:  return sp.Mul(expand(a, n - 1), a, evaluate=False)
        return sp.Pow(expand(a, -n), -1, evaluate=False)
    if not isinstance(expr_list, list):
        expr_list = [expr_list]
    var_map, sym_map = {}, {}
    for i, expr in enumerate(expr_list):
        tree = ExprTree(expr)
        for subtree in tree.preorder(tree.root):
            subexpr = subtree.expr
            if subexpr.func == sp.Pow:
                exponent = subtree.children[1].expr
                if exponent.func == sp.Integer and abs(exponent) > 1:
                    subtree.expr = expand(*subexpr.args)
                    tree.build(subtree, clear=True)
        for subtree in tree.preorder():
            subexpr = subtree.expr
            if isinstance(subexpr, sp.Rational) and subexpr != sp.S.NegativeOne:
                sign = 1 if subexpr >= 0 else -1
                subexpr *= sign
                try:
                    repl = sym_map[subexpr]
                except KeyError:
                    p, q = subexpr.p, subexpr.q
                    var_name = prefix + '_Rational_' + str(p) + '_' + str(q) \
                        if q != 1 else prefix + '_Integer_' + str(p)
                    repl = sp.Symbol(var_name)
                    var_map[repl], sym_map[subexpr] = subexpr, repl
                subtree.expr = repl * sign
                if ignore and sign < 0:
                    try:
                        subtree.expr *= sym_map[-1] * sign
                    except KeyError:
                        repl = sp.Symbol(prefix + '_NegativeOne_')
                        var_map[repl], sym_map[-1] = sp.S.NegativeOne, repl
                        subtree.expr *= repl * sign
            elif subexpr == sp.S.NegativeOne:
                try:
                    subtree.expr = sym_map[-1]
                except KeyError:
                    repl = sp.Symbol(prefix + '_NegativeOne_')
                    var_map[repl], sym_map[-1] = sp.S.NegativeOne, repl
                    subtree.expr = repl
        for subtree in tree.preorder(tree.root):
            subexpr = subtree.expr
            if subexpr.func == sp.Pow:
                exponent = subtree.children[1].expr
                if exponent.func == sp.Symbol:
                    subtree.children[1].expr = var_map[exponent]
        expr = tree.reconstruct()
        if factor:
            for subtree in tree.preorder():
                if isinstance(subtree.expr, sp.Function):
                    for var in var_map.keys():
                        child = subtree.children[0]
                        child.expr = sp.collect(child.expr, var)
                        child.children.clear()
            expr = tree.reconstruct()
            for var in var_map:
                expr = sp.collect(expr, var)
        expr_list[i] = expr
    if len(expr_list) == 1:
        expr_list = expr_list[0]
    return expr_list, var_map

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
        if (expr.func == sp.Mul and len(expr.args) == 2 and sp.S.NegativeOne in expr.args and \
                 any((arg.func == sp.Symbol) for arg in expr.args)):
            for k in range(i + 1, len(replaced)):
                    if sym in replaced[k][1].free_symbols:
                        replaced[k] = (replaced[k][0], replaced[k][1].subs(sym, expr))
            for k in range(len(reduced)):
                if sym in reduced[k].free_symbols:
                    reduced[k] = reduced[k].subs(sym, expr)
            # Remove the replaced expression from the list
            replaced.pop(i); i -= 1
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