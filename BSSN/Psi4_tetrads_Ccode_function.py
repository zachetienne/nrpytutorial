# As documented in the NRPy+ tutorial module
#   Tutorial-Psi4_tetrads_Ccode_function.ipynb,
#   this module construct the C code function
#   for generating the tetrads necessary for
#   computing \psi_4 (as well as other
#   Weyl scalars and invariants in principle)

# Author: Zachariah B. Etienne
#         (zachetie **at** gmail **dot* com)

# Step 1.a: import all needed modules from NRPy+:
import sympy as sp
from outputC import *
import NRPy_param_funcs as par
import indexedexp as ixp
import grid as gri
import finite_difference as fin
import reference_metric as rfm

def Psi4_tetrads_Ccode_function(tetrad_Ccode_filename = "BSSN/Psi4_tetrad.h",TetradChoice="QuasiKinnersley"):
    # Step 1.b: Given the chosen coordinate system, set up
    #           corresponding reference metric and needed
    #           reference metric quantities
    # The following function call sets up the reference metric
    #    and related quantities, including rescaling matrices ReDD,
    #    ReU, and hatted quantities.
    rfm.reference_metric()

    # Step 1.c: Import the tetrad module
    import BSSN.Psi4_tetrads as BP4T
    par.set_parval_from_str("BSSN.Psi4_tetrads::TetradChoice", TetradChoice)

    # Step 2: Construct the C function header and
    #         convert (xx0,xx1,xx2) to the corresponding
    #         (x,y,z), as both are needed for the tetrad
    #         expressions.
    outCparams = "preindent=1,outCfileaccess=a,outCverbose=False,includebraces=False"

    with open(tetrad_Ccode_filename, "w") as file:
        file.write("""
// Taking as input (xx0,xx1,xx2), this function outputs
//   the chosen Psi4 tetrad in the (xx0,xx1,xx2) basis
void Psi4_tetrad(REAL xx0xx1xx2[3], const REAL *in_gfs,
                 REAL n4U[4],REAL mre4U[4],REAL mim4U[4]) {
    const REAL xx0 = xx0xx1xx2[0];
    const REAL xx1 = xx0xx1xx2[1];
    const REAL xx2 = xx0xx1xx2[2];
""")
    outputC([rfm.xxCart[0], rfm.xxCart[1], rfm.xxCart[2]], ["REAL x", "REAL y", "REAL z"], tetrad_Ccode_filename,
            outCparams)

    # Step 3: Output the tetrad in the reference-metric basis.

    # Step 3.a: BP4T.Psi4_tetrads() to construct the symbolic
    #           expressions for the tetrad vectors $n^\mu$,
    #           $\Re m^\mu$, and $\Im m^\mu$, which are needed
    #           to construct $\Psi_4$.

    BP4T.Psi4_tetrads()
    Psi4_tetrad_vecs = []

    # Step 3.b: As the tetrad expressions depend on BSSN
    #           gridfunctions, we pass the expressions into
    #           fin.FD_outputC() so that the needed gridfunction
    #           values are read in from memory as appropriate.
    for i in range(4):
        Psi4_tetrad_vecs.append(lhrh(lhs="n4U[" + str(i) + "]", rhs=BP4T.n4U[i]))
        Psi4_tetrad_vecs.append(lhrh(lhs="mre4U[" + str(i) + "]", rhs=BP4T.mre4U[i]))
        Psi4_tetrad_vecs.append(lhrh(lhs="mim4U[" + str(i) + "]", rhs=BP4T.mim4U[i]))

    fin.FD_outputC(tetrad_Ccode_filename, Psi4_tetrad_vecs, outCparams)

    with open(tetrad_Ccode_filename, "a") as file:
        file.write("}\n")