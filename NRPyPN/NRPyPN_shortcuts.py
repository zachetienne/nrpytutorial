# As documented in the NRPyPN notebook
# NRPyPN_shortcuts.ipynb, this Python script
# provides useful shortcuts for inputting
# post-Newtonian expressions into SymPy/NRPy+

# Basic functions:
# dot(a,b): 3-vector dot product
# cross(a,b): 3-vector cross product
# div(a,b): a shortcut for SymPy's sp.Rational(a,b), to declare rational numbers
# num_eval(expr): Numerically evaluates NRPyPN expressions

# Author:  Zach Etienne
#          zachetie **at** gmail **dot* com

# Step 0: Add NRPy's directory to the path
# https://stackoverflow.com/questions/16780014/import-file-from-parent-directory
import os,sys                    # Standard Python modules for multiplatform OS-level functions
nrpy_dir_path = os.path.join("..")
if nrpy_dir_path not in sys.path:
    sys.path.append(nrpy_dir_path)
import sympy as sp               # SymPy: The Python computer algebra package upon which NRPy+ depends
import indexedexp as ixp         # NRPy+: Symbolic indexed expression (e.g., tensors, vectors, etc.) support

# Step 1: Declare several global variables used
#         throughout NRPyPN
m1,m2 = sp.symbols('m1 m2',real=True)
S1U = ixp.declarerank1("S1U")
S2U = ixp.declarerank1("S2U")
pU = ixp.declarerank1("pU")
nU = ixp.declarerank1("nU")

drdt = sp.symbols('drdt', real=True)
Pt, Pr = sp.symbols('Pt Pr', real=True)
# Some references use r, others use q to represent the
#   distance between the two point masses. This is rather
#   confusing since q is also used to represent the
#   mass ratio m2/m1. However, q is the canonical position
#   variable name in Hamiltonian mechanics, so both are
#   well justified. It should be obvious which is which
#   throughout NRPyPN.
r, q = sp.symbols('r q', real=True)
chi1U = ixp.declarerank1('chi1U')
chi2U = ixp.declarerank1('chi2U')

# Euler-Mascheroni gamma constant:
gamma_EulerMascheroni = sp.symbols('gamma_EulerMascheroni',real=True)

# Derived quantities used in Damour et al papers:
n12U = ixp.zerorank1()
n21U = ixp.zerorank1()
p1U = ixp.zerorank1()
p2U = ixp.zerorank1()
for i in range(3):
    n12U[i] = +nU[i]
    n21U[i] = -nU[i]
    p1U[i]     = +pU[i]
    p2U[i]     = -pU[i]

# Step 2.a: Define dot and cross product of vectors
def dot(vec1,vec2):
    vec1_dot_vec2 = sp.sympify(0)
    for i in range(3):
        vec1_dot_vec2 += vec1[i]*vec2[i]
    return vec1_dot_vec2

def cross(vec1,vec2):
    vec1_cross_vec2 = ixp.zerorank1()
    LeviCivitaSymbol = ixp.LeviCivitaSymbol_dim3_rank3()
    for i in range(3):
        for j in range(3):
            for k in range(3):
                vec1_cross_vec2[i] += LeviCivitaSymbol[i][j][k]*vec1[j]*vec2[k]
    return vec1_cross_vec2

# Step 2.b: Construct rational numbers a/b via div(a,b)
def div(a,b):
    return sp.Rational(a,b)

# Step 3: num_eval(expr), a means to numerically evaluate SymPy/NRPyPN
#         expressions
def num_eval(expr,
             qmassratio =  1.0,  # must be >= 1
             nr         = 12.0,  # Orbital separation
             nchi1x     = +0.,
             nchi1y     = +0.,
             nchi1z     = +0.,
             nchi2x     = +0.,
             nchi2y     = +0.,
             nchi2z     = +0.,
             nPt=None, ndrdt=None):

    # DERIVED QUANTITIES BELOW
    # We want m1+m2 = 1, so that
    #         m2/m1 = qmassratio
    # We find below:
    nm1   =          1/(1+qmassratio)
    nm2   = qmassratio/(1+qmassratio)
    # This way nm1+nm2 = (qmassratio+1)/(1+qmassratio) = 1 CHECK
    #      and nm2/nm1 = qmassratio                        CHECK

    nS1U0 = nchi1x*nm1**2
    nS1U1 = nchi1y*nm1**2
    nS1U2 = nchi1z*nm1**2
    nS2U0 = nchi2x*nm2**2
    nS2U1 = nchi2y*nm2**2
    nS2U2 = nchi2z*nm2**2

    if nPt != None:
        expr2 = expr.subs(Pt,nPt)
        expr  = expr2
    if ndrdt != None:
        expr2 = expr.subs(drdt,ndrdt)
        expr  = expr2
    return expr\
.subs(m1,nm1).subs(m2,nm2)\
.subs(S1U[0],nS1U0).subs(S1U[1],nS1U1).subs(S1U[2],nS1U2)\
.subs(S2U[0],nS2U0).subs(S2U[1],nS2U1).subs(S2U[2],nS2U2)\
.subs(chi1U[0],nchi1x).subs(chi1U[1],nchi1y).subs(chi1U[2],nchi1z)\
.subs(chi2U[0],nchi2x).subs(chi2U[1],nchi2y).subs(chi2U[2],nchi2z)\
.subs(r,nr).subs(q,nr).subs(sp.pi,sp.N(sp.pi)\
.subs(gamma_EulerMascheroni,0.5772156649015328606065120900824024310421))
