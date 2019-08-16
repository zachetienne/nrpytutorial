import sympy as sp


# Generic module that can be used for testing, specifically test_functions
# function() can be called in isolation. init_function2() must be called before function2() is called.
def function(create_gamma=False):
    global alpha, betaU

    a, b, c = sp.symbols('a b c')

    alpha = a + b + c
    betaU = [0, a**2 + 2*b**2 + c**2, sp.sqrt(a + b)]

    if create_gamma:
        global gamma

        gamma = sp.atan2(b, a)


def function2(create_gamma=False):

    global alpha2, betaU2

    alpha2 = a2 + b2 + c2

    betaU2 = [0, a2**2 + 2*b2**2 + c2**2, sp.sqrt(a2 + b2)]

    if create_gamma:
        global gamma2

        gamma2 = sp.atan2(b2, a2)


def init_function2():
    global a2, b2, c2
    a2, b2, c2 = sp.symbols('a2 b2 c2')
