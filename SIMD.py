""" Convert Expression to SIMD Compiler Intrinsics """
# Authors: Ken Sible & Zachariah Etienne
# Emails: ksible *at* outlook *dot** com
#         zachetie *at* gmail *dot** com

from sympy import (Integer, Rational, Float, Function, Symbol,
    Add, Mul, Pow, Abs, S, sign, srepr, simplify,
    var, sin, cos, exp, log, preorder_traversal)
from expr_tree import ExprTree
from cse_helpers import cse_preprocess

# Basic Arithmetic Operations (Debugging)
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
def NegFusedMulAddSIMD_check(a, b, c):
    return -a*b + c
def NegFusedMulSubSIMD_check(a, b, c):
    return -a*b - c
def DivSIMD_check(a, b):
    return a / b
def signSIMD_check(a):
    return sign(a)

# Transcendental Operations (Debugging)
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

def expr_convert_to_SIMD_intrins(expr, map_sym_to_rat=None, prefix="", SIMD_find_more_FMAsFMSs="True", debug="False"):
    """ Convert expression to SIMD compiler intrinsics

        :arg:    SymPy expression
        :arg:    symbol to rational dictionary
        :arg:    option to find more FMA/FMS patterns
        :arg:    back-substitute and check difference
        :return: expression containing SIMD compiler intrinsics

        >>> from sympy.abc import a, b, c, d
        >>> from cse_helpers import cse_preprocess
        >>> convert = expr_convert_to_SIMD_intrins

        >>> convert(a**2)
        MulSIMD(a, a)

        >>> convert(a**(-2))
        DivSIMD(_Integer_1, MulSIMD(a, a))

        >>> convert(a**(1/2))
        SqrtSIMD(a)

        >>> convert(a**(-1/2))
        DivSIMD(1, SqrtSIMD(a))

        >>> from sympy import Rational
        >>> convert(a**Rational(1, 3))
        CbrtSIMD(a)

        >>> convert(a**b)
        PowSIMD(a, b)

        >>> convert(a - b)
        SubSIMD(a, b)

        >>> convert(a + b - c)
        AddSIMD(b, SubSIMD(a, c))

        >>> convert(a + b + c)
        AddSIMD(a, AddSIMD(b, c))

        >>> convert(a + b + c + d)
        AddSIMD(AddSIMD(a, b), AddSIMD(c, d))

        >>> convert(a*b*c)
        MulSIMD(a, MulSIMD(b, c))

        >>> convert(a*b*c*d)
        MulSIMD(MulSIMD(a, b), MulSIMD(c, d))

        >>> convert(a/b)
        DivSIMD(a, b)

        >>> convert(a*b + c)
        FusedMulAddSIMD(a, b, c)

        >>> convert(a*b - c)
        FusedMulSubSIMD(a, b, c)

        >>> convert(-a*b + c)
        NegFusedMulAddSIMD(a, b, c)

        >>> convert(-a*b - c)
        NegFusedMulSubSIMD(a, b, c)
    """
    for item in preorder_traversal(expr):
        for arg in item.args:
            if isinstance(arg, Symbol):
                var(str(arg))

    def lookup_rational(arg):
        if arg.func == Symbol:
            try: arg = map_sym_to_rat[arg]
            except KeyError: pass
        return arg

    if map_sym_to_rat is None:
        expr, map_sym_to_rat = cse_preprocess(expr)

    map_rat_to_sym = {map_sym_to_rat[v]:v for v in map_sym_to_rat}

    expr_orig, tree = expr, ExprTree(expr)

    AbsSIMD  = Function("AbsSIMD")
    AddSIMD  = Function("AddSIMD")
    SubSIMD  = Function("SubSIMD")
    MulSIMD  = Function("MulSIMD")
    FusedMulAddSIMD = Function("FusedMulAddSIMD")
    FusedMulSubSIMD = Function("FusedMulSubSIMD")
    NegFusedMulAddSIMD = Function("NegFusedMulAddSIMD")
    NegFusedMulSubSIMD = Function("NegFusedMulSubSIMD")
    DivSIMD  = Function("DivSIMD")
    SignSIMD = Function("SignSIMD")

    PowSIMD  = Function("PowSIMD")
    SqrtSIMD = Function("SqrtSIMD")
    CbrtSIMD = Function("CbrtSIMD")
    ExpSIMD  = Function("ExpSIMD")
    LogSIMD  = Function("LogSIMD")
    SinSIMD  = Function("SinSIMD")
    CosSIMD  = Function("CosSIMD")

    # Step 1: Replace transcendental functions, power functions, and division expressions.
    #   Note: SymPy does not represent fractional integers as rationals since
    #         those are explicitly declared using the rational class, and hence
    #         the following algorithm does not affect fractional integers.
    #         SymPy: srepr(a**(-2)) = Pow(a, -2)
    #         NRPy:  srepr(a**(-2)) = DivSIMD(1, MulSIMD(a, a))
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
    tree.reconstruct()

    def IntegerPowSIMD(a, n):
        # Recursive Helper Function: Construct Integer Powers
        if   n == 2:
            return MulSIMD(a, a)
        if n > 2:
            return MulSIMD(IntegerPowSIMD(a, n - 1), a)
        if n <= -2:
            one = Symbol(prefix + '_Integer_1')
            try: map_rat_to_sym[1]
            except KeyError:
                map_sym_to_rat[one], map_rat_to_sym[1] = S.One, one
            return DivSIMD(one, IntegerPowSIMD(a, -n))
        if n == -1:
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
                subtree.children.pop(1)
            elif exponent == -0.5:
                subtree.expr = DivSIMD(1, SqrtSIMD(args[0]))
                tree.build(subtree)
            elif exponent == Rational(1, 3):
                subtree.expr = CbrtSIMD(args[0])
                subtree.children.pop(1)
            elif isinstance(exponent, Integer):
                subtree.expr = IntegerPowSIMD(args[0], exponent)
                tree.build(subtree)
            else:
                subtree.expr = PowSIMD(*args)
    tree.reconstruct()

    # Step 2: Replace subtraction expressions.
    #   Note: SymPy: srepr(a - b) = Add(a, Mul(-1, b))
    #         NRPy:  srepr(a - b) = SubSIMD(a, b)
    for subtree in tree.preorder():
        func = subtree.expr.func
        args = list(subtree.expr.args)
        if func == Add:
            try:
                # Find the first occurrence of a negative product inside the addition
                i = next(i for i, arg in enumerate(args) if arg.func == Mul and \
                        any(lookup_rational(arg) == -1 for arg in args[i].args))
                # Find the first occurrence of a negative symbol inside the product
                j = next(j for j, arg in enumerate(args[i].args) if lookup_rational(arg) == -1)
                # Find the first non-negative argument of the product
                k = next(k for k in range(len(args)) if k != i)
                # Remove the negative symbol from the product
                subargs = list(args[i].args); subargs.pop(j)
                # Build the subtraction expression for replacement
                subexpr = SubSIMD(args[k], Mul(*subargs))
                args = [arg for arg in args if arg not in (args[i], args[k])]
                if len(args) > 0:
                    subexpr = Add(subexpr, *args)
                subtree.expr = subexpr
                tree.build(subtree)
            except StopIteration: pass
    tree.reconstruct()

    # Step 3: Replace addition and multiplication expressions.
    #   Note: SIMD addition and multiplication compiler intrinsics can read
    #         only two arguments at once, whereas SymPy's Mul() and Add()
    #         operators can read an arbitrary number of arguments.
    #         SymPy: srepr(a*b*c*d) = Mul(a, b, c, d)
    #         NRPy:  srepr(a*b*c*d) = MulSIMD(MulSIMD(a, b), MulSIMD(c, d))
    for subtree in tree.preorder():
        func = subtree.expr.func
        args = subtree.expr.args
        if func in (Mul, Add):
            func = MulSIMD if func == Mul else AddSIMD
            subexpr = func(*args[-2:])
            args, N = args[:-2], len(args) - 2
            for i in range(0, N, 2):
                if N - i > 1:
                    tmpexpr = func(args[i], args[i + 1])
                    subexpr = func(tmpexpr, subexpr, evaluate=False)
                else:
                    subexpr = func(args[i], subexpr, evaluate=False)
            subtree.expr = subexpr
            tree.build(subtree)
    tree.reconstruct()

    # Step 4: Replace the pattern Mul(Div(1, b), a) or Mul(a, Div(1, b)) with Div(a, b).
    for subtree in tree.preorder():
        func = subtree.expr.func
        args = subtree.expr.args
        # MulSIMD(DivSIMD(1, b), a) >> DivSIMD(a, b)
        if   func == MulSIMD and args[0].func == DivSIMD and \
                lookup_rational(args[0].args[0]) == 1:
            subtree.expr = DivSIMD(args[1], args[0].args[1])
            tree.build(subtree)
        # MulSIMD(a, DivSIMD(1, b)) >> DivSIMD(a, b)
        elif func == MulSIMD and args[1].func == DivSIMD and \
                lookup_rational(args[1].args[0]) == 1:
            subtree.expr = DivSIMD(args[0], args[1].args[1])
            tree.build(subtree)
    tree.reconstruct()

    # Step 5: Now that all multiplication and addition functions only take two
    #         arguments, we can define fused-multiply-add functions,
    #         where AddSIMD(a, MulSIMD(b, c)) = b*c + a = FusedMulAddSIMD(b, c, a),
    #         or    AddSIMD(MulSIMD(b, c), a) = b*c + a = FusedMulAddSIMD(b, c, a).
    #   Note: Fused-multiply-add (FMA3) is standard on Intel CPUs with the AVX2
    #         instruction set, starting with Haswell processors in 2013:
    #         https://en.wikipedia.org/wiki/Haswell_(microarchitecture)

    # Step 5.a: Find double FMA patterns first [e.g. FMA(a, b, FMA(c, d, e))].
    #   Note: Double FMA simplifications do not guarantee a significant performance impact when solving BSSN equations
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
                tree.build(subtree)
            # b*c + d*e + a -> FMA(b,c,FMA(d,e,a))
            # Validate:
            # x = b*c + d*e + a
            # outputC(x,"x", params="SIMD_enable=True,SIMD_debug=True")
            # AddSIMD(AddSIMD(MulSIMD(b,c), MulSIMD(d,e)),a) >> FusedMulAddSIMD(b, c, FusedMulAddSIMD(d,e,a))
            elif func == AddSIMD and args[0].func == AddSIMD and args[0].args[0].func == MulSIMD and args[0].args[1].func == MulSIMD:
                subtree.expr = FusedMulAddSIMD(                args[0].args[0].args[0], args[0].args[0].args[1],
                                               FusedMulAddSIMD(args[0].args[1].args[0], args[0].args[1].args[1],
                                                               args[1]))
                tree.build(subtree)
        tree.reconstruct()

    # Step 5.b: Find single FMA patterns.
    for subtree in tree.preorder():
        func = subtree.expr.func
        args = subtree.expr.args
        # AddSIMD(MulSIMD(b, c), a) >> FusedMulAddSIMD(b, c, a)
        if   func == AddSIMD and args[0].func == MulSIMD:
            subtree.expr = FusedMulAddSIMD(args[0].args[0], args[0].args[1], args[1])
            tree.build(subtree)
        # AddSIMD(a, MulSIMD(b, c)) >> FusedMulAddSIMD(b, c, a)
        elif func == AddSIMD and args[1].func == MulSIMD:
            subtree.expr = FusedMulAddSIMD(args[1].args[0], args[1].args[1], args[0])
            tree.build(subtree)
        # SubSIMD(MulSIMD(b, c), a) >> FusedMulSubSIMD(b, c, a)
        elif func == SubSIMD and args[0].func == MulSIMD:
            subtree.expr = FusedMulSubSIMD(args[0].args[0], args[0].args[1], args[1])
            tree.build(subtree)
        # SubSIMD(a, MulSIMD(b, c)) >> NegativeFusedMulAddSIMD(b, c, a)
        elif func == SubSIMD and args[1].func == MulSIMD:
            subtree.expr = NegFusedMulAddSIMD(args[1].args[0], args[1].args[1], args[0])
            tree.build(subtree)
        # FMS(-1, MulSIMD(a, b), c) >> NegativeFusedMulSubSIMD(b, c, a)
        func = subtree.expr.func
        args = subtree.expr.args
        if func == FusedMulSubSIMD and args[1].func == MulSIMD and lookup_rational(args[0]) == -1:
            subtree.expr = NegFusedMulSubSIMD(args[1].args[0], args[1].args[1], args[2])
            tree.build(subtree)
    tree.reconstruct()

    # Step 5.c: Remaining double FMA patterns that previously in Step 5.a were difficult to find.
    #   Note: Double FMA simplifications do not guarantee a significant performance impact when solving BSSN equations
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
                tree.build(subtree)
            # b*c - (a - d*e) -> SubSIMD(FusedMulAddSIMD(b, c, MulSIMD(d, e)), a) >> FMA(b,c,FMS(d,e,a))
            # Validate:
            # x = b * c - (a - d * e)
            # outputC(x, "x", params="SIMD_enable=True,SIMD_debug=True")
            elif func == SubSIMD and args[0].func == FusedMulAddSIMD and args[0].args[2].func == MulSIMD:
                subtree.expr = FusedMulAddSIMD(args[0].args[0], args[0].args[1],
                                               FusedMulSubSIMD(args[0].args[2].args[0], args[0].args[2].args[1],
                                                               args[1]))
                tree.build(subtree)
            # (b*c - d*e) - a -> SubSIMD(FusedMulSubSIMD(b, c, MulSIMD(d, e)), a) >> FMS(b,c,FMA(d,e,a))
            # Validate:
            # x = (b*c - d*e) - a
            # outputC(x,"x", params="SIMD_enable=True,SIMD_debug=True")
            elif func == SubSIMD and args[0].func == FusedMulSubSIMD and args[0].args[2].func == MulSIMD:
                subtree.expr = FusedMulSubSIMD(args[0].args[0], args[0].args[1],
                                               FusedMulAddSIMD(args[0].args[2].args[0], args[0].args[2].args[1],
                                                               args[1]))
                tree.build(subtree)
        tree.reconstruct()

    # Step 5.d: NegFusedMulAddSIMD(a,b,c) = -a*b + c:
    for subtree in tree.preorder():
        func = subtree.expr.func
        args = subtree.expr.args
        # FMA(a,Mul(-1,b),c) >> NFMA(a,b,c)
        if   func == FusedMulAddSIMD and args[1].func == MulSIMD and \
             lookup_rational(args[1].args[0]) == -1:
            subtree.expr = NegFusedMulAddSIMD(args[0],args[1].args[1],args[2])
            tree.build(subtree)
        # FMA(a,Mul(b,-1),c) >> NFMA(a,b,c)
        elif func == FusedMulAddSIMD and args[1].func == MulSIMD and \
             lookup_rational(args[1].args[1]) == -1:
            subtree.expr = NegFusedMulAddSIMD(args[0],args[1].args[0],args[2])
            tree.build(subtree)
        # FMA(Mul(-1,a), b,c) >> NFMA(a,b,c)
        elif func == FusedMulAddSIMD and args[0].func == MulSIMD and \
             lookup_rational(args[0].args[0]) == -1:
            subtree.expr = NegFusedMulAddSIMD(args[0].args[1],args[1],args[2])
            tree.build(subtree)
        # FMA(Mul(a,-1), b,c) >> NFMA(a,b,c)
        elif func == FusedMulAddSIMD and args[0].func == MulSIMD and \
             lookup_rational(args[0].args[1]) == -1:
            subtree.expr = NegFusedMulAddSIMD(args[0].args[0],args[1],args[2])
            tree.build(subtree)
    tree.reconstruct()

    # Step 5.e: Replace e.g., FMA(-1,b,c) with SubSIMD(c,b) and similar patterns
    for subtree in tree.preorder():
        func = subtree.expr.func
        args = subtree.expr.args
        # FMA(-1,b,c) >> SubSIMD(c,b)
        if   func == FusedMulAddSIMD and lookup_rational(args[0]) == -1:
            subtree.expr = SubSIMD(args[2], args[1])
            tree.build(subtree)
        # FMA(a,-1,c) >> SubSIMD(c,a)
        elif func == FusedMulAddSIMD and lookup_rational(args[1]) == -1:
            subtree.expr = SubSIMD(args[2], args[0])
            tree.build(subtree)
        # FMS(a,-1,c) >> MulSIMD(-1,AddSIMD(a,c))
        elif func == FusedMulSubSIMD and lookup_rational(args[1]) == -1:
            subtree.expr = MulSIMD(args[1], AddSIMD(args[0], args[2]))
            tree.build(subtree)
        # FMS(-1,b,c) >> MulSIMD(-1,AddSIMD(b,c))
        elif func == FusedMulSubSIMD and lookup_rational(args[0]) == -1:
            subtree.expr = MulSIMD(args[0], AddSIMD(args[1], args[2]))
            tree.build(subtree)
    tree.reconstruct()

    # Step 5.f: NegFusedMulSubSIMD(a,b,c) = -a*b - c:
    for subtree in tree.preorder():
        func = subtree.expr.func
        args = subtree.expr.args
        # NFMA(a,b,Mul(-1,c)) >> NFMS(a,b,c)
        if   func == NegFusedMulAddSIMD and args[2].func == MulSIMD and \
             lookup_rational(args[2].args[0]) == -1:
            subtree.expr = NegFusedMulSubSIMD(args[0],args[1],args[2].args[1])
            tree.build(subtree)
        # NFMA(a,b,Mul(c,-1)) >> NFMS(a,b,c)
        elif func == NegFusedMulAddSIMD and args[2].func == MulSIMD and \
             lookup_rational(args[2].args[1]) == -1:
            subtree.expr = NegFusedMulSubSIMD(args[0],args[1],args[2].args[0])
            tree.build(subtree)
        # FMS(a,Mul(-1,b),c) >> NFMS(a,b,c)
        elif func == FusedMulSubSIMD and args[1].func == MulSIMD and \
             lookup_rational(args[1].args[0]) == -1:
            subtree.expr = NegFusedMulSubSIMD(args[0],args[1].args[1],args[2])
            tree.build(subtree)
        # FMS(a,Mul(b,-1),c) >> NFMS(a,b,c)
        elif func == FusedMulSubSIMD and args[1].func == MulSIMD and \
             lookup_rational(args[1].args[1]) == -1:
            subtree.expr = NegFusedMulSubSIMD(args[0],args[1].args[0],args[2])
            tree.build(subtree)
        # FMS(a,Mul([something],Mul(-1,b)),c) >> NFMS(a,Mul([something],b),c)
        elif func == FusedMulSubSIMD and args[1].func == MulSIMD and \
             args[1].args[1].func == MulSIMD and lookup_rational(args[1].args[1].args[0]) == -1:
            subtree.expr = NegFusedMulSubSIMD(args[0], MulSIMD(args[1].args[0],args[1].args[1].args[1]), args[2])
            tree.build(subtree)
        # FMS(a,Mul([something],Mul(b,-1)),c) >> NFMS(a,Mul([something],b),c)
        elif func == FusedMulSubSIMD and args[1].func == MulSIMD and \
             args[1].args[1].func == MulSIMD and lookup_rational(args[1].args[1].args[1]) == -1:
            subtree.expr = NegFusedMulSubSIMD(args[0], MulSIMD(args[1].args[0],args[1].args[1].args[0]), args[2])
            tree.build(subtree)
    tree.reconstruct()

    # Step 5.g: Find single FMA patterns again, as some new ones might be found.
    for subtree in tree.preorder():
        func = subtree.expr.func
        args = subtree.expr.args
        # AddSIMD(MulSIMD(b, c), a) >> FusedMulAddSIMD(b, c, a)
        if   func == AddSIMD and args[0].func == MulSIMD:
            subtree.expr = FusedMulAddSIMD(args[0].args[0], args[0].args[1], args[1])
            tree.build(subtree)
        # AddSIMD(a, MulSIMD(b, c)) >> FusedMulAddSIMD(b, c, a)
        elif func == AddSIMD and args[1].func == MulSIMD:
            subtree.expr = FusedMulAddSIMD(args[1].args[0], args[1].args[1], args[0])
            tree.build(subtree)
        # SubSIMD(MulSIMD(b, c), a) >> FusedMulSubSIMD(b, c, a)
        elif func == SubSIMD and args[0].func == MulSIMD:
            subtree.expr = FusedMulSubSIMD(args[0].args[0], args[0].args[1], args[1])
            tree.build(subtree)
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
                if abs(subexpr - Integer(subexpr)) < 1.0e-14*subexpr:
                    subtree.expr = Integer(subexpr)
        expr_diff = tree_diff.reconstruct()

        if expr_diff != 0:
            simp_expr_diff = simplify(expr_diff)
            if simp_expr_diff != 0:
                raise Warning('Expression Difference: ' + str(simp_expr_diff))
    return(expr)

if __name__ == "__main__":
    import doctest
    doctest.testmod()
