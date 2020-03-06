from sympy import (Integer, Rational, Float, Function, Symbol,
    Add, Mul, Pow, Abs, S, N, sign, srepr, simplify, sympify, 
    var, sin, cos, exp, log, symbols, preorder_traversal)
from SIMDExprTree import ExprTree
import re, sys

# For debugging purposes, Part 1:
# Basic arithmetic operations
def ConstSIMD_check(a):
    return Float(a, 34)
def AbsSIMD_check(a):
    return Abs(a)
def nrpyAbsSIMD_check(a):
    return Abs(a)
def AddSIMD_check(a, b):
    return a + b
def SubSIMD_check(a, b):
    return a - b
def MulSIMD_check(a, b):
    return a * b
def FusedMulAddSIMD_check(a, b, c):
    return a*b + c
def FusedMulSubSIMD_check(a, b, c):
    return a*b - c
def DivSIMD_check(a, b):
    return a / b
def signSIMD_check(a):
    return sign(a)

# For debugging purposes, Part 2:
# Transcendental operations
def PowSIMD_check(a, b):
    return a**b
def SqrtSIMD_check(a):
    return a**(Rational(1, 2))
def CbrtSIMD_check(a):
    return a**(Rational(1, 3))
def ExpSIMD_check(a):
    return exp(a)
def LogSIMD_check(a):
    return log(a)
def SinSIMD_check(a):
    return sin(a)
def CosSIMD_check(a):
    return cos(a)

# Input: SymPy expression.
# Return value: SymPy expression containing all needed SIMD compiler intrinsics
def expr_convert_to_SIMD_intrins(expr, map_sym_to_rat, prefix="", SIMD_find_more_FMAsFMSs="False", debug="False"):
    for item in preorder_traversal(expr):
        for arg in item.args:
            if isinstance(arg, Symbol):
                var(str(arg))

    def lookup_rational(arg):
        if arg.func == Symbol:
            try: arg = map_sym_to_rat[arg]
            except KeyError: pass
        return arg

    map_rat_to_sym = {map_sym_to_rat[v]:v for v in map_sym_to_rat}

    expr_orig, tree = expr, ExprTree(expr)

    AbsSIMD  = Function("AbsSIMD")
    AddSIMD  = Function("AddSIMD")
    SubSIMD  = Function("SubSIMD")
    MulSIMD  = Function("MulSIMD")
    FusedMulAddSIMD = Function("FusedMulAddSIMD")
    FusedMulSubSIMD = Function("FusedMulSubSIMD")
    DivSIMD  = Function("DivSIMD")
    SignSIMD = Function("SignSIMD")

    PowSIMD  = Function("PowSIMD")
    SqrtSIMD = Function("SqrtSIMD")
    CbrtSIMD = Function("CbrtSIMD")
    ExpSIMD  = Function("ExpSIMD")
    LogSIMD  = Function("LogSIMD")
    SinSIMD  = Function("SinSIMD")
    CosSIMD  = Function("CosSIMD")

    # Step 1: Replace transcendental, power, and division functions with SIMD equivalents
    #         Note that due to how SymPy expresses rational numbers, the following does not
    #         affect fractional expressions of integers
    for subtree in tree.preorder():
        func = subtree.expr.func
        args = subtree.expr.args
        if   func == Abs:
            subtree.expr = AbsSIMD(args[0])
        elif func == exp:
            subtree.expr = ExpSIMD(args[0])
        elif func == log:
            subtree.expr = LogSIMD(args[0])
        elif func == sin:
            subtree.expr = SinSIMD(args[0])
        elif func == cos:
            subtree.expr = CosSIMD(args[0])
        elif func == sign:
            subtree.expr = SignSIMD(args[0])
    expr = tree.reconstruct(evaluate=True)

    # Fun little recursive function for constructing integer powers:
    def IntegerPowSIMD(a, n):
        if   n == 2:
            return MulSIMD(a, a)
        elif n > 2:
            return MulSIMD(IntegerPowSIMD(a, n - 1), a)
        elif n <= -2:
            one = Symbol(prefix + '_Integer_1')
            try: map_rat_to_sym[1]
            except KeyError:
                map_sym_to_rat[one], map_rat_to_sym[1] = S.One, one
            return DivSIMD(one, IntegerPowSIMD(a, -n))
        elif n == -1:
            one = Symbol(prefix + '_Integer_1')
            try: map_rat_to_sym[1]
            except KeyError:
                map_sym_to_rat[one], map_rat_to_sym[1] = S.One, one
            return DivSIMD(one, a)

    for subtree in tree.preorder():
        func = subtree.expr.func
        args = subtree.expr.args
        if func == Pow:
            exponent = lookup_rational(args[1])
            if   exponent == 0.5:
                subtree.expr = SqrtSIMD(args[0])
                subtree.children.pop(1) # Remove 0.5
            elif exponent == -0.5:
                subtree.expr = DivSIMD(1, SqrtSIMD(args[0]))
                tree.build(subtree, clear=True)
            elif exponent == Rational(1, 3):
                subtree.expr = CbrtSIMD(args[0])
                subtree.children.pop(1) # Remove -0.5
            elif isinstance(exponent, Integer):
                subtree.expr = IntegerPowSIMD(args[0], exponent)
                tree.build(subtree, clear=True)
            else:
                subtree.expr = PowSIMD(*args)
    expr = tree.reconstruct()
    
    # We must evaluate the expression, otherwise nested multiplications
    # will arise that conflict with the following replacements in Step 3.

    # Step 2: SIMD multiplication and addition compiler intrinsics read in
    #         only two arguments at once, where SymPy's Mul() and Add()
    #         operators can read an arbitrary number of arguments.
    #         Here, we split e.g., Mul(a, b, c, d) into
    #         MulSIMD(a, MulSIMD(b, MulSIMD(c, d))),
    #         To accomplish this easily, we construct a string
    #         'MulSIMD(A, MulSIMD(B, ...', where MulSIMD(a, b) is some user-
    #         defined function that takes in only two arguments, and then
    #         evaluate the string using the eval() function.
    # Implementation detail: If we did not perform Step 2 above, the eval
    #         function would automatically evaluate all Rational expressions
    #         as though they were input as integers: e.g., 1/2 evaluates to 0.
    #         This is undesirable, so we instead define new, temporary
    #         functions IntegerTMP and RationalTMP that are undisturbed by
    #         the eval()
    for subtree in tree.preorder():
        func = subtree.expr.func
        args = subtree.expr.args
        if (func == Mul or func == Add):
            func = MulSIMD if func == Mul else AddSIMD
            subexpr = func(*args[-2:])
            for arg in args[:-2]:
                subexpr = func(arg, subexpr, evaluate=False)
            subtree.expr = subexpr
            tree.build(subtree, clear=True)
    expr = tree.reconstruct()
    
    # Step 3: Simplification patterns:
    # Step 3.a: Replace the pattern Mul(Div(1, b), a) or Mul(a, Div(1, b)) with Div(a, b):
    for subtree in tree.preorder():
        func = subtree.expr.func
        args = subtree.expr.args
        # MulSIMD(DivSIMD(1, b), a) >> DivSIMD(a, b)
        if   func == MulSIMD and args[0].func == DivSIMD and \
                lookup_rational(args[0].args[0]) == 1:
            subtree.expr = DivSIMD(args[1], args[0].args[1])
            tree.build(subtree, clear=True)
        # MulSIMD(a, DivSIMD(1, b)) >> DivSIMD(a, b)
        elif func == MulSIMD and args[1].func == DivSIMD and \
                lookup_rational(args[1].args[0]) == 1:
            subtree.expr = DivSIMD(args[0], args[1].args[1])
            tree.build(subtree, clear=True)
    expr = tree.reconstruct()

    # Step 3.b: Subtraction intrinsics. SymPy replaces all a - b with a + (-b) = Add(a, Mul(-1, b))
    #         Here, we replace
    #         a) AddSIMD(MulSIMD(-1, b), a),
    #         b) AddSIMD(MulSIMD(b, -1), a),
    #         c) AddSIMD(a, MulSIMD(-1, b)), and
    #         d) AddSIMD(a, MulSIMD(b, -1))
    #         with SubSIMD(a, b)
    for subtree in tree.preorder():
        func = subtree.expr.func
        args = subtree.expr.args
        # AddSIMD(MulSIMD(-1, b), a) >> SubSIMD(a, b)
        if   func == AddSIMD and args[0].func == MulSIMD and \
                lookup_rational(args[0].args[0]) == -1:
            subtree.expr = SubSIMD(args[1], args[0].args[1])
            tree.build(subtree, clear=True)
        # AddSIMD(MulSIMD(b, -1), a) >> SubSIMD(a, b)
        elif func == AddSIMD and args[0].func == MulSIMD and \
                lookup_rational(args[0].args[1]) == -1:
            subtree.expr = SubSIMD(args[1], args[0].args[0])
            tree.build(subtree, clear=True)
        # AddSIMD(a, MulSIMD(-1, b)) >> SubSIMD(a, b)
        elif func == AddSIMD and args[1].func == MulSIMD and \
                lookup_rational(args[1].args[0]) == -1:
            subtree.expr = SubSIMD(args[0], args[1].args[1])
            tree.build(subtree, clear=True)
        # AddSIMD(a, MulSIMD(b, -1)) >> SubSIMD(a, b)
        elif func == AddSIMD and args[1].func == MulSIMD and \
                lookup_rational(args[1].args[1]) == -1:
            subtree.expr = SubSIMD(args[0], args[1].args[0])
            tree.build(subtree, clear=True)
    expr = tree.reconstruct()

    # Step 4: Now that all multiplication and addition functions only take two
    #         arguments, we can now easily define fused-multiply-add functions,
    #         where AddSIMD(a, MulSIMD(b, c)) = b*c + a = FusedMulAddSIMD(b, c, a),
    #         or    AddSIMD(MulSIMD(b, c), a) = b*c + a = FusedMulAddSIMD(b, c, a).
    # Fused multiply add (FMA3) is standard on Intel CPUs with the AVX2
    #         instruction set, starting with Haswell processors in 2013:
    #         https://en.wikipedia.org/wiki/Haswell_(microarchitecture)

    # Step 4.a: Find double FMA patterns first [e.g., FMA(a,b,FMA(c,d,e))]:
    #           NOTE: Double FMA simplifications do not guarantee a significant performance impact when solving BSSN equations:
    if SIMD_find_more_FMAsFMSs == "True":
        for subtree in tree.preorder():
            func = subtree.expr.func
            args = subtree.expr.args
            # a + b*c + d*e -> FMA(b,c,FMA(d,e,a))
            # AddSIMD(a, AddSIMD(MulSIMD(b,c), MulSIMD(d,e))) >> FusedMulAddSIMD(b, c, FusedMulAddSIMD(d,e,a))
            # Validate:
            # x = a + b*c + d*e
            # outputC(x,"x", params="SIMD_enable=True,SIMD_debug=True")
            if  (func == AddSIMD and args[1].func == AddSIMD and args[1].args[0].func == MulSIMD and args[1].args[1].func == MulSIMD):
                subtree.expr = FusedMulAddSIMD(                args[1].args[0].args[0], args[1].args[0].args[1],
                                               FusedMulAddSIMD(args[1].args[1].args[0], args[1].args[1].args[1],
                                                               args[0]))
                tree.build(subtree, clear=True)
            # b*c + d*e + a -> FMA(b,c,FMA(d,e,a))
            # Validate:
            # x = b*c + d*e + a
            # outputC(x,"x", params="SIMD_enable=True,SIMD_debug=True")
            # AddSIMD(AddSIMD(MulSIMD(b,c), MulSIMD(d,e)),a) >> FusedMulAddSIMD(b, c, FusedMulAddSIMD(d,e,a))
            elif func == AddSIMD and args[0].func == AddSIMD and args[0].args[0].func == MulSIMD and args[0].args[1].func == MulSIMD:
                subtree.expr = FusedMulAddSIMD(                args[0].args[0].args[0], args[0].args[0].args[1],
                                               FusedMulAddSIMD(args[0].args[1].args[0], args[0].args[1].args[1],
                                                               args[1]))
                tree.build(subtree, clear=True)
        expr = tree.reconstruct()

    # Step 4.b: Next find single FMA patterns:
    for subtree in tree.preorder():
        func = subtree.expr.func
        args = subtree.expr.args
        # AddSIMD(MulSIMD(b, c), a) >> FusedMulAddSIMD(b, c, a)
        if   func == AddSIMD and args[0].func == MulSIMD:
            subtree.expr = FusedMulAddSIMD(args[0].args[0], args[0].args[1], args[1])
            tree.build(subtree, clear=True)
        # AddSIMD(a, MulSIMD(b, c)) >> FusedMulAddSIMD(b, c, a)
        elif func == AddSIMD and args[1].func == MulSIMD:
            subtree.expr = FusedMulAddSIMD(args[1].args[0], args[1].args[1], args[0])
            tree.build(subtree, clear=True)
        # SubSIMD(MulSIMD(b, c), a) >> FusedMulSubSIMD(b, c, a)
        elif func == SubSIMD and args[0].func == MulSIMD:
            subtree.expr = FusedMulSubSIMD(args[0].args[0], args[0].args[1], args[1])
            tree.build(subtree, clear=True)
    expr = tree.reconstruct()

    # Step 4.c: Leftover double FMA patterns that are difficult to find in Step 5.a:
    #           NOTE: Double FMA simplifications do not guarantee a significant performance impact when solving BSSN equations:
    if SIMD_find_more_FMAsFMSs == "True":
        for subtree in tree.preorder():
            func = subtree.expr.func
            args = subtree.expr.args
            # (b*c - d*e) + a -> AddSIMD(a, FusedMulSubSIMD(b, c, MulSIMD(d, e))) >> FusedMulSubSIMD(b, c, FusedMulSubSIMD(d,e,a))
            # Validate:
            # x = (b*c - d*e) + a
            # outputC(x,"x", params="SIMD_enable=True,SIMD_debug=True")
            if func == AddSIMD and args[1].func == FusedMulSubSIMD and args[1].args[2].func == MulSIMD:
                subtree.expr = FusedMulSubSIMD(                args[1].args[0]        ,args[1].args[1],
                                               FusedMulSubSIMD(args[1].args[2].args[0],args[1].args[2].args[1],
                                                               args[0]))
                tree.build(subtree, clear=True)
            # b*c - (a - d*e) -> SubSIMD(FusedMulAddSIMD(b, c, MulSIMD(d, e)), a) >> FMA(b,c,FMS(d,e,a))
            # Validate:
            # x = b * c - (a - d * e)
            # outputC(x, "x", params="SIMD_enable=True,SIMD_debug=True")
            elif func == SubSIMD and args[0].func == FusedMulAddSIMD and args[0].args[2].func == MulSIMD:
                subtree.expr = FusedMulAddSIMD(args[0].args[0], args[0].args[1],
                                               FusedMulSubSIMD(args[0].args[2].args[0], args[0].args[2].args[1],
                                                               args[1]))
                tree.build(subtree, clear=True)
            # (b*c - d*e) - a -> SubSIMD(FusedMulSubSIMD(b, c, MulSIMD(d, e)), a) >> FMS(b,c,FMA(d,e,a))
            # Validate:
            # x = (b*c - d*e) - a
            # outputC(x,"x", params="SIMD_enable=True,SIMD_debug=True")
            elif func == SubSIMD and args[0].func == FusedMulSubSIMD and args[0].args[2].func == MulSIMD:
                subtree.expr = FusedMulSubSIMD(args[0].args[0], args[0].args[1],
                                               FusedMulAddSIMD(args[0].args[2].args[0], args[0].args[2].args[1],
                                                               args[1]))
                tree.build(subtree, clear=True)
        expr = tree.reconstruct()

    if debug == "True":
        expr_check = eval(str(expr).replace("SIMD", "SIMD_check"))
        expr_check = expr_check.subs(-1, Symbol('_NegativeOne_'))

        expr_diff = expr_check - expr_orig
        # The eval(str(srepr())) below normalizes the expression, 
        # fixing a cancellation issue in SymPy ~0.7.4.
        expr_diff = eval(str(srepr(expr_diff)))
        tree_diff = ExprTree(expr_diff)
        for subtree in tree_diff.preorder():
            subexpr = subtree.expr
            if subexpr.func == Float:
                if abs(subexpr - Integer(subexpr)) < 1.0e-14:
                    subtree.expr = Integer(item)
        expr_diff = tree_diff.reconstruct()

        if expr_diff != 0:
            simp_expr_diff = simplify(expr_diff)
            if simp_expr_diff != 0:
                raise Warning('Expression Difference: ' + str(simp_expr_diff))
    return(expr)
