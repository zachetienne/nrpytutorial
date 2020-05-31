# This module implements the gammabar = gammahat constraint.

# Author: Zachariah B. Etienne
#         zachetie **at** gmail **dot* com

# Step P1: import all needed modules from NRPy+:
from outputC import nrpyAbs,lhrh,check_if_string__error_if_not,outCfunction  # NRPy+: Core C code output module
import finite_difference as fin   # NRPy+: Finite difference C code generation module
import grid as gri                # NRPy+: Functions having to do with numerical grids
import indexedexp as ixp          # NRPy+: Symbolic indexed expression (e.g., tensors, vectors, etc.) support
import reference_metric as rfm    # NRPy+: Reference metric support
import sympy as sp                # SymPy: The Python computer algebra package upon which NRPy+ depends
import BSSN.BSSN_quantities as Bq # NRPy+: Computes useful BSSN quantities
import os                         # Standard Python modules for multiplatform OS-level functions

def Enforce_Detgammabar_Constraint_symb_expressions():
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
                (nrpyAbs(rfm.detgammahat) / detgammabar) ** (sp.Rational(1, 3)) * (KroneckerDeltaDD[i][j] + hDD[i][j]) \
                - KroneckerDeltaDD[i][j]

    enforce_detg_constraint_symb_expressions = [
        lhrh(lhs=gri.gfaccess("in_gfs", "hDD00"), rhs=hprimeDD[0][0]),
        lhrh(lhs=gri.gfaccess("in_gfs", "hDD01"), rhs=hprimeDD[0][1]),
        lhrh(lhs=gri.gfaccess("in_gfs", "hDD02"), rhs=hprimeDD[0][2]),
        lhrh(lhs=gri.gfaccess("in_gfs", "hDD11"), rhs=hprimeDD[1][1]),
        lhrh(lhs=gri.gfaccess("in_gfs", "hDD12"), rhs=hprimeDD[1][2]),
        lhrh(lhs=gri.gfaccess("in_gfs", "hDD22"), rhs=hprimeDD[2][2])]

    return enforce_detg_constraint_symb_expressions

def output_Enforce_Detgammabar_Constraint_Ccode(outdir="BSSN/", exprs="", Read_xxs=False):
    # Step 0: Check if outdir is string; error out if not.
    check_if_string__error_if_not(outdir,"outdir")

    desc = "Enforce det(gammabar) = det(gammahat) constraint."
    name = "enforce_detgammabar_constraint"
    params = "const rfm_struct *restrict rfmstruct,const paramstruct *restrict params, REAL *restrict in_gfs"
    loopopts = "AllPoints,Enable_rfm_precompute"
    if Read_xxs:
        params = "const paramstruct *restrict params, REAL *restrict xx[3], REAL *restrict in_gfs"
        loopopts = "AllPoints,Read_xxs"
    outCfunction(
        outfile=os.path.join(outdir, name + ".h"), desc=desc, name=name, params=params,
        body=fin.FD_outputC("returnstring", exprs,
                            params="outCverbose=False,preindent=1,includebraces=False").replace("IDX4", "IDX4S"),
        loopopts=loopopts)
