import sympy as sp               # SymPy: The Python computer algebra package upon which NRPy+ depends
import indexedexp as ixp         # NRPy+: Symbolic indexed expression (e.g., tensors, vectors, etc.) support

thismodule = __name__

# We'll write this as a function so that we can calculate the expressions on-demand for any choice of i
def find_cp_cm(lapse,shifti,gammaUUii):
    # Inputs:  u0,vi,lapse,shift,gammadet,gupii
    # Outputs: cplus,cminus

    # a = 1/(alpha^2)
    a = sp.sympify(1)/(lapse*lapse)
    # b = 2 beta^i / alpha^2
    b = sp.sympify(2) * shifti /(lapse*lapse)
    # c = -g^{ii} + (beta^i)^2 / alpha^2
    c = - gammaUUii + shifti*shifti/(lapse*lapse)

    # Now, we are free to solve the quadratic equation as usual. We take care to avoid passing a
    # negative value to the sqrt function.
    detm = b*b - sp.sympify(4)*a*c

    import Min_Max_and_Piecewise_Expressions as noif
    detm = sp.sqrt(noif.max_noif(sp.sympify(0),detm))
    global cplus,cminus
    cplus  = sp.Rational(1,2)*(-b/a + detm/a)
    cminus = sp.Rational(1,2)*(-b/a - detm/a)

# We'll write this as a function, and call it within HLLE_solver, below.
def find_cmax_cmin(flux_dirn,gamma_faceDD,beta_faceU,alpha_face):
    # Inputs:  flux direction flux_dirn, Inverse metric gamma_faceUU, shift beta_faceU,
    #          lapse alpha_face, metric determinant gammadet_face
    # Outputs: maximum and minimum characteristic speeds cmax and cmin
    # First, we need to find the characteristic speeds on each face
    gamma_faceUU,unusedgammaDET = ixp.generic_matrix_inverter3x3(gamma_faceDD)
    find_cp_cm(alpha_face,beta_faceU[flux_dirn],gamma_faceUU[flux_dirn][flux_dirn])
    cp = cplus
    cm = cminus

    # The following algorithms have been verified with random floats:

    global cmax,cmin
    # Now, we need to set cmax to the larger of cpr,cpl, and 0

    import Min_Max_and_Piecewise_Expressions as noif
    cmax = noif.max_noif(cp,sp.sympify(0))

    # And then, set cmin to the smaller of cmr,cml, and 0
    cmin = -noif.min_noif(cm,sp.sympify(0))
