from sympy import Integer,Symbol,symbols,simplify,Rational,Function,srepr,sin,cos,exp,log,Abs,Add,Mul,Pow,preorder_traversal,N,Float,S,var,sympify
import NRPy_param_funcs as par
import re

# For debugging purposes, Part 1:
# Basic arithmetic operations
def ConstSIMD_check(a):
    return Float(a,34)
def AbsSIMD_check(a):
    return Abs(a)
def AddSIMD_check(a,b):
    return a+b
def SubSIMD_check(a,b):
    return a-b
def MulSIMD_check(a,b):
    return a*b
def FusedMulAddSIMD_check(a,b,c):
    return a*b + c
def FusedMulSubSIMD_check(a,b,c):
    return a*b - c
def DivSIMD_check(a,b):
    return a/b

# For debugging purposes, Part 2:
# Transcendental operations
def PowSIMD_check(a, b):
    return a**b
def SqrtSIMD_check(a):
    return a**(Rational(1,2))
def CbrtSIMD_check(a):
    return a**(Rational(1,3))
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
# Complication: SIMD functions require numerical constants to be stored in SIMD arrays
# Resolution: This function extends lists "SIMD_const_varnms" and "SIMD_const_values",
#             which store the name of each constant SIMD array (e.g., _Integer_1) and
#             the value of each variable (e.g., 1.0).
def expr_convert_to_SIMD_intrins(expr,  SIMD_const_varnms,SIMD_const_values,debug="False"):

    # Declare all variables, so we can eval them in the next (AddSIMD & MulSIMD) step
    for item in preorder_traversal(expr):
        for i in range(len(item.args)):
            if isinstance(item.args[i],Symbol):
                var(str(item.args[i]))

    expr_orig = expr

    AbsSIMD = Function("AbsSIMD")
    AddSIMD = Function("AddSIMD")
    SubSIMD = Function("SubSIMD")
    MulSIMD = Function("MulSIMD")
    FusedMulAddSIMD = Function("FusedMulAddSIMD")
    FusedMulSubSIMD = Function("FusedMulSubSIMD")
    DivSIMD = Function("DivSIMD")

    PowSIMD = Function("PowSIMD")
    SqrtSIMD = Function("SqrtSIMD")
    CbrtSIMD = Function("CbrtSIMD")
    ExpSIMD = Function("ExpSIMD")
    LogSIMD = Function("LogSIMD")
    SinSIMD = Function("SinSIMD")
    CosSIMD = Function("CosSIMD")

    # Step 1: Replace transcendental, power, and division functions with SIMD equivalents
    #         Note that due to how SymPy expresses rational numbers, the following does not
    #         affect fractional expressions of integers
    for item in preorder_traversal(expr):
        if item.func == Abs:
            expr = expr.xreplace({item: AbsSIMD(item.args[0])})
        elif item.func == exp:
            expr = expr.xreplace({item: ExpSIMD(item.args[0])})
        elif item.func == log:
            expr = expr.xreplace({item: LogSIMD(item.args[0])})
        elif item.func == sin:
            expr = expr.xreplace({item: SinSIMD(item.args[0])})
        elif item.func == cos:
            expr = expr.xreplace({item: CosSIMD(item.args[0])})

    for item in preorder_traversal(expr):
        if item.func == Pow:
            if item.args[1] == 0.5:
                expr = expr.xreplace({item: SqrtSIMD(item.args[0])})
            elif item.args[1] == -0.5:
                expr = expr.xreplace({item: DivSIMD(1,SqrtSIMD(item.args[0]))})
            elif item.args[1] == Rational(1,3):
                expr = expr.xreplace({item: CbrtSIMD(item.args[0])})
            elif item.args[1] == 2:
                expr = expr.xreplace({item: MulSIMD(item.args[0],item.args[0])})
            elif item.args[1] == -2:
                expr = expr.xreplace({item: DivSIMD(1,MulSIMD(item.args[0],item.args[0]))})
            elif item.args[1] == 3: #and item.args[0].is_Symbol:
                expr = expr.xreplace({item: MulSIMD(item.args[0],MulSIMD(item.args[0],item.args[0]))})
            elif item.args[1] == -3: # and len(item.args)==1 and item.args[0].is_Symbol:
                expr = expr.xreplace({item: DivSIMD(1,MulSIMD(item.args[0], MulSIMD(item.args[0], item.args[0])))})
            elif item.args[1] == 4 and item.args[0].is_Symbol:
                expr = expr.xreplace({item: MulSIMD(item.args[0],MulSIMD(item.args[0],MulSIMD(item.args[0],item.args[0])))})
            elif item.args[1] == -4 and item.args[0].is_Symbol:
                expr = expr.xreplace({item: DivSIMD(1,MulSIMD(item.args[0], MulSIMD(item.args[0], MulSIMD(item.args[0], item.args[0]))))})
            elif item.args[1] == 5 and item.args[0].is_Symbol:
                expr = expr.xreplace({item: MulSIMD(item.args[0],MulSIMD(item.args[0],MulSIMD(item.args[0],MulSIMD(item.args[0],item.args[0]))))})
            elif item.args[1] == -5 and item.args[0].is_Symbol:
                expr = expr.xreplace({item: DivSIMD(1,MulSIMD(item.args[0], MulSIMD(item.args[0], MulSIMD(item.args[0], MulSIMD(item.args[0], item.args[0])))))})
            elif item.args[1] == 6 and item.args[0].is_Symbol:
                expr = expr.xreplace({item: MulSIMD(item.args[0],MulSIMD(item.args[0],MulSIMD(item.args[0],MulSIMD(item.args[0],MulSIMD(item.args[0],item.args[0])))))})
            elif item.args[1] == -6 and item.args[0].is_Symbol:
                expr = expr.xreplace({item: DivSIMD(1,MulSIMD(item.args[0], MulSIMD(item.args[0], MulSIMD(item.args[0], MulSIMD(item.args[0], MulSIMD(item.args[0], item.args[0]))))))})
            elif item.args[1] == 7 and item.args[0].is_Symbol:
                expr = expr.xreplace({item: MulSIMD(item.args[0],MulSIMD(item.args[0],MulSIMD(item.args[0],MulSIMD(item.args[0],MulSIMD(item.args[0],MulSIMD(item.args[0],item.args[0]))))))})
            elif item.args[1] == -7 and item.args[0].is_Symbol:
                expr = expr.xreplace({item: DivSIMD(1,MulSIMD(item.args[0], MulSIMD(item.args[0], MulSIMD(item.args[0], MulSIMD(item.args[0], MulSIMD(item.args[0], MulSIMD(item.args[0], item.args[0])))))))})
            elif item.args[1] == -1:
                expr = expr.xreplace({item: DivSIMD(1, item.args[0])})
            else:
                expr = expr.xreplace({item: PowSIMD(item.args[0],item.args[1])})

    # Step 2: Replace all rational numbers (expressed as Rational(a,b))
    #         and integers with the new functions RationalTMP and
    #         IntegerTMP, where Rational(a,b) -> RationalTMP(a,b)
    #         and Integer(a) -> IntegerTMP(a)
    RationalTMP = Function("RationalTMP")
    IntegerTMP = Function("IntegerTMP")

    string = str(srepr(expr))
    string2 = re.sub('Integer\(([+\-0-9]+)\)',
                     "(Function('IntegerTMP')('\\1'))", string)
    expr = eval(string2)

    string = str(srepr(expr))
    string2 = re.sub('Rational\(([0-9]+), ([0-9]+)\)',
                     "(Function('RationalTMP')(('\\1'),('\\2')))", string)
    expr = eval(string2)

    # Step 3: The pattern Mul(-1,Rational(a,b)) is often seen. Replace with Rational(-a,b).
    for item in preorder_traversal(expr):
        if item.func == Mul:
            idx_has_Negative1 = -10
            for i in range(len(item.args)):
                if item.args[i].func == IntegerTMP and item.args[i].args[0] == -1:
                    idx_has_Negative1 = i
            if idx_has_Negative1 >= 0:
                for i in range(len(item.args)):
                    if item.args[i].func==RationalTMP:
                        tempitem_orig = Mul(IntegerTMP(-1),RationalTMP(item.args[i].args[0], item.args[i].args[1]))
                        tempitem_new  = RationalTMP(-item.args[i].args[0], item.args[i].args[1])
                        expr = expr.subs(tempitem_orig, tempitem_new )
                        break

    # Step 4: SIMD multiplication and addition compiler intrinsics read in
    #         only two arguments at once, where SymPy's Mul() and Add()
    #         operators can read an arbitrary number of arguments.
    #         Here, we split e.g., Mul(a,b,c,d) into
    #         MulSIMD(a,MulSIMD(b,MulSIMD(c,d))),
    #         To accomplish this easily, we construct a string
    #         'MulSIMD(A,MulSIMD(B,...', where MulSIMD(a,b) is some user-
    #         defined function that takes in only two arguments, and then
    #         evaluate the string using the eval() function.
    # Implementation detail: If we did not perform Step 2 above, the eval
    #         function would automatically evaluate all Rational expressions
    #         as though they were input as integers: e.g., 1/2 evaluates to 0.
    #         This is undesirable, so we instead define new, temporary
    #         functions IntegerTMP and RationalTMP that are undisturbed by
    #         the eval()
    for item in preorder_traversal(expr):
        if item.func == Mul:
            blahtmp = symbols('blahtmp')
            tempitem = blahtmp
            for i in range(len(item.args)-1):
                tempitem = MulSIMD(item.args[i],tempitem)
            tempitem = tempitem.xreplace({blahtmp: item.args[len(item.args)-1]})
            expr = expr.xreplace({item: tempitem})
    for item in preorder_traversal(expr):
        if item.func == Add:
            blahtmp = symbols('blahtmp')
            tempitem = blahtmp
            for i in range(len(item.args)-1):
                tempitem = AddSIMD(item.args[i],tempitem)
            tempitem = tempitem.xreplace({blahtmp: item.args[len(item.args)-1]})
            expr = expr.xreplace({item: tempitem})

    # Step 5: Simplification patterns:
    # Step 5a: Replace the pattern Mul(Div(1,b),a) or Mul(a,Div(1,b)) with Div(a,b):
    for item in preorder_traversal(expr):
        if item.func == MulSIMD:
            if item.func == MulSIMD and (item.args[0].func == DivSIMD and item.args[0].args[0].func == IntegerTMP and item.args[0].args[0].args[0] == 1):
                expr = expr.xreplace({item: DivSIMD(item.args[1],item.args[0].args[1])})
            if item.func == MulSIMD and (item.args[1].func == DivSIMD and item.args[1].args[0].func == IntegerTMP and item.args[1].args[0].args[0] == 1):
                expr = expr.xreplace({item: DivSIMD(item.args[0],item.args[1].args[1])})

    # Step 5: Subtraction intrinsics. SymPy replaces all a-b with a + (-b) = Add(a,Mul(-1,b))
    #         Here, we replace
    #         a) AddSIMD(a,MulSIMD(-1,b)),
    #         b) AddSIMD(a,MulSIMD(b,-1)),
    #         c) AddSIMD(MulSIMD(-1,b),a), and
    #         d) AddSIMD(MulSIMD(b,-1),a)
    #         with SubSIMD(a,b)
    for item in preorder_traversal(expr):
        tempitem = item
        # First match patterns a) and b): AddSIMD(a,MULSIMD(.,.)):
        if item.func == AddSIMD and item.args[1].func == MulSIMD:
            # Pattern a) AddSIMD(a,MulSIMD(-1,b)) --> SubSIMD(a,b):
            if item.args[1].args[0].func == IntegerTMP and item.args[1].args[0].args[0] == -1:
                tempitem = SubSIMD(item.args[0],item.args[1].args[1])
            # Pattern b) AddSIMD(a,MulSIMD(b,-1)) --> SubSIMD(a,b):
            elif item.args[1].args[1].func == IntegerTMP and item.args[1].args[1].args[0] == -1:
                tempitem = SubSIMD(item.args[0],item.args[1].args[0])
        # Next match patterns c) and d): AddSIMD(MulSIMD(.,.),a):
        elif item.func == AddSIMD and item.args[0].func == MulSIMD:
            # Pattern c) AddSIMD(MulSIMD(-1,b),a) --> SubSIMD(a,b):
            if item.args[0].args[0].func == IntegerTMP and item.args[0].args[0].args[0] == -1:
                tempitem = SubSIMD(item.args[1],item.args[0].args[1])
            # Pattern d) AddSIMD(MulSIMD(b,-1,a)) --> SubSIMD(a,b):
            elif item.args[0].args[1].func == IntegerTMP and item.args[0].args[1].args[0] == -1:
                tempitem = SubSIMD(item.args[1],item.args[0].args[0])
        expr = expr.subs(item,tempitem)

    # Step 6: Now that all multiplication and addition functions only take two
    #         arguments, we can now easily define fused-multiply-add functions,
    #         where AddSIMD(a,MulSIMD(b,c)) = b*c + a = FusedMulAddSIMD(b,c,a),
    #         or    AddSIMD(MulSIMD(b,c),a) = b*c + a = FusedMulAddSIMD(b,c,a).
    # Fused multiply add (FMA3) is standard on Intel CPUs with the AVX2
    #         instruction set, starting with Haswell processors in 2013:
    #         https://en.wikipedia.org/wiki/Haswell_(microarchitecture)
    for item in preorder_traversal(expr):
        tempitem = item
        # If the pattern is a*b + c, replace with FMA(a,b,c)
        if item.func == AddSIMD and item.args[0].func == MulSIMD:
            tempitem = FusedMulAddSIMD(item.args[0].args[0],item.args[0].args[1],item.args[1])
        # If the pattern is c + a*b, replace with FMA(a,b,c)
        if item.func == AddSIMD and item.args[1].func == MulSIMD:
            tempitem = FusedMulAddSIMD(item.args[1].args[0],item.args[1].args[1],item.args[0])

        # If the pattern is a*b - c, replace with FMS(a,b,c)
        if item.func == SubSIMD and item.args[0].func == MulSIMD:
            tempitem = FusedMulSubSIMD(item.args[0].args[0],item.args[0].args[1],item.args[1])

        if item != tempitem: expr = expr.subs(item, tempitem)

    # Step 7: SIMD intrinsics cannot take integers or rational numbers as arguments.
    #         Therefore we must declare all integers & rational numbers as
    #         const vector doubles (e.g., const _m256d ...). To make the code
    #         more human readable, we adopt the convention
    #         RationalTMP(1,3) = 1/3 = "Rational_1_3
    #         RationalTMP(-1,3) = -1/3 = "Rational_m1_3
    # TODO: Keep track of all integers and rationals, so repeats are not computed?

    # Step 7a: Set all variable names and corresponding values.
    for item in preorder_traversal(expr):
        if item.func == RationalTMP:
            # Set variable name
            if item.args[0]*item.args[1] < 0:
                SIMD_const_varnms.extend(["_Rational_m"+str(abs(item.args[0]))+"_"+str(abs(item.args[1]))])
            elif item.args[0] > 0 and item.args[1] > 0:
                SIMD_const_varnms.extend(["_Rational_"+str(item.args[0])+"_"+str(item.args[1])])
            else:
                # E.g., doesn't make sense to have -1/-3. SymPy should have simplified this.
                print("Found a weird Rational(a,b) expression, where a<0 and b<0. Report to SymPy devels")
                print("Specifically, found that a="+str(item.args[0])+" and b="+str(item.args[1]))
                exit(1)
            # Set variable value, to 34 digits of precision
            SIMD_const_values.extend([str(N(Float(item.args[0],34)/Float(item.args[1],34),34))])
        elif item.func == IntegerTMP:
            # Set variable name
            if item.args[0] < 0:
                SIMD_const_varnms.extend(["_Integer_m"+str(-item.args[0])])
            else:
                SIMD_const_varnms.extend(["_Integer_" + str(item.args[0])])
            # Set variable value, to 34 digits of precision
            SIMD_const_values.extend([str((Float(item.args[0],34)))])

    # Step 7b: Replace all integers and rationals with the appropriate variable names:
    for item in preorder_traversal(expr):
        tempitem = item
        if item.func == RationalTMP:
            if item.args[0]*item.args[1] < 0:
                tempitem = var("_Rational_m" + str(abs(item.args[0])) + "_" + str(abs(item.args[1])))
            elif item.args[0] > 0 and item.args[1] > 0:
                tempitem = var("_Rational_" + str(item.args[0]) + "_" + str(item.args[1]))
            else:
                # E.g., doesn't make sense to have -1/-3. SymPy should have simplified this.
                print("Found a weird Rational(a,b) expression, where a<0 and b<0. Report to SymPy devels")
                print("Specifically, found that a=" + str(item.args[0]) + " and b=" + str(item.args[1]))
                exit(1)
        elif item.func == IntegerTMP:
            if item.args[0] < 0:
                tempitem = var("_Integer_m" + str(-item.args[0]))
            else:
                tempitem = var("_Integer_" + str(item.args[0]))
        if item != tempitem: expr = expr.subs(item, tempitem)

    def lookup_name_output_idx(name, list_of_names):
        for i in range(len(list_of_names)):
            if list_of_names[i] == name:
                return i
        print("I SHOULDN'T BE HERE!",name,list_of_names)
        exit(1)

    if debug=="True":
        expr_check = expr
        if "SIMD" in str(expr):
            expr_check = eval(str(expr).replace("SIMD","SIMD_check"))

        for item in preorder_traversal(expr_check):
            tempitem = item
            if item.is_Symbol and str(item)[0]=="_":
                if str(item)[:9]=="_Integer_" or str(item)[:10]=="_Rational_":
                    tempitem = SIMD_const_values[lookup_name_output_idx(str(item), SIMD_const_varnms)]
            if item != tempitem: expr_check = expr_check.subs(item, tempitem)

        expr_diff = expr_check - expr_orig
        # Some variables do not want to cancel in SymPy ~0.7.4. The eval(str(srepr())) below normalizes the expression.
        expr_diff = eval(str(srepr(expr_diff)))
        for item in preorder_traversal(expr_diff):
            if item.func == Float:
                if abs(item - Integer(item)) < 1.0e-14:
                    expr_diff = expr_diff.xreplace({item:Integer(item)})

        # Only simplify if expr_diff != 0:
        if expr_diff != 0:
            simp_expr_diff = simplify(expr_diff)
            if simp_expr_diff != 0:
                print("Warning: found possible diff",(simp_expr_diff))
    return(expr)
