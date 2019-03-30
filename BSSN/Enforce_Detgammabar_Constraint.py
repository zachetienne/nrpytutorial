# This module implements the gammabar = gammahat constraint.

# Author: Zachariah B. Etienne
#         zachetie **at** gmail **dot* com

# Step P1: import all needed modules from NRPy+:
import sympy as sp
import NRPy_param_funcs as par
from outputC import *
import indexedexp as ixp
import grid as gri
import finite_difference as fin
import reference_metric as rfm
import BSSN.BSSN_quantities as Bq
import loop as lp

def output_Enforce_Detgammabar_Constraint_Ccode():
    # Set spatial dimension (must be 3 for BSSN)
    DIM = 3

    # Then we set the coordinate system for the numerical grid
    rfm.reference_metric()  # Create ReU, ReDD needed for rescaling B-L initial data, generating BSSN RHSs, etc.

    # We will need the h_{ij} quantities defined within BSSN_RHSs
    #    below when we enforce the gammahat=gammabar constraint
    # Step 1: All barred quantities are defined in terms of BSSN rescaled gridfunctions,
    #         which we declare here in case they haven't yet been declared elsewhere.

    Bq.declare_BSSN_gridfunctions_if_not_declared_already()
    hDD = Bq.hDD
    Bq.BSSN_basic_tensors()
    gammabarDD = Bq.gammabarDD

    # First define the Kronecker delta:
    KroneckerDeltaDD = ixp.zerorank2()
    for i in range(DIM):
        KroneckerDeltaDD[i][i] = sp.sympify(1)

    # The detgammabar in BSSN_RHSs is set to detgammahat when BSSN_RHSs::detgbarOverdetghat_equals_one=True (default),
    #    so we manually compute it here:
    dummygammabarUU, detgammabar = ixp.symm_matrix_inverter3x3(gammabarDD)

    # Next apply the constraint enforcement equation above.
    hprimeDD = ixp.zerorank2()
    for i in range(DIM):
        for j in range(DIM):
            hprimeDD[i][j] = \
                (sp.Abs(rfm.detgammahat) / detgammabar) ** (sp.Rational(1, 3)) * (KroneckerDeltaDD[i][j] + hDD[i][j]) \
                - KroneckerDeltaDD[i][j]

    enforce_detg_constraint_vars = [ \
        lhrh(lhs=gri.gfaccess("in_gfs", "hDD00"), rhs=hprimeDD[0][0]),
        lhrh(lhs=gri.gfaccess("in_gfs", "hDD01"), rhs=hprimeDD[0][1]),
        lhrh(lhs=gri.gfaccess("in_gfs", "hDD02"), rhs=hprimeDD[0][2]),
        lhrh(lhs=gri.gfaccess("in_gfs", "hDD11"), rhs=hprimeDD[1][1]),
        lhrh(lhs=gri.gfaccess("in_gfs", "hDD12"), rhs=hprimeDD[1][2]),
        lhrh(lhs=gri.gfaccess("in_gfs", "hDD22"), rhs=hprimeDD[2][2])]

    enforce_gammadet_string = fin.FD_outputC("returnstring", enforce_detg_constraint_vars,
                                             params="outCverbose=False,preindent=0,includebraces=False")

    with open("BSSN/enforce_detgammabar_constraint.h", "w") as file:
        indent = "   "
        file.write(
            "void enforce_detgammabar_constraint(const int Nxx_plus_2NGHOSTS[3],REAL *xx[3], REAL *in_gfs) {\n\n")
        file.write(lp.loop(["i2", "i1", "i0"], ["0", "0", "0"],
                           ["Nxx_plus_2NGHOSTS[2]", "Nxx_plus_2NGHOSTS[1]", "Nxx_plus_2NGHOSTS[0]"],
                           ["1", "1", "1"], ["#pragma omp parallel for",
                                             "    const REAL xx2 = xx[2][i2];",
                                             "        const REAL xx1 = xx[1][i1];"], "",
                           "const REAL xx0 = xx[0][i0];\n" + enforce_gammadet_string))
        file.write("}\n")

    print("Output C implementation of det(gammabar) constraint to file BSSN/enforce_detgammabar_constraint.h")