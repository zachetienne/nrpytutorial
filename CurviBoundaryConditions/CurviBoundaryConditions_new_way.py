# This module provides functions for setting up Curvilinear boundary conditions,
#     as documented in Tutorial-Start_to_Finish-Curvilinear_BCs_new_way.ipynb

# Authors: Zachariah B. Etienne
#         zachetie **at** gmail **dot* com
#          Terrence Pierre Jacques

# Step P1: Import needed NRPy+ core modules:
from outputC import outputC      # NRPy+: Core C code output module
from outputC import outC_NRPy_basic_defines_h_dict  # NRPy+: Core C code output module
from outputC import add_to_Cfunction_dict  # NRPy+: Core C code output module
from outputC import Cfunction    # NRPy+: Core C code output module
from outputC import indent_Ccode # NRPy+: Core C code output module
import NRPy_param_funcs as par   # NRPy+: Parameter interface
import sympy as sp               # SymPy: The Python computer algebra package upon which NRPy+ depends
import grid as gri               # NRPy+: Functions having to do with numerical grids
import indexedexp as ixp         # NRPy+: Symbolic indexed expression (e.g., tensors, vectors, etc.) support
import reference_metric as rfm   # NRPy+: Reference metric support
import finite_difference as fin  # NRPy+: Finite-difference module
import os, sys           # Standard Python modules for multiplatform OS-level functions
from UnitTesting.assert_equal import check_zero  # NRPy+: Checks whether an expression evaluates to zero.

_unused = par.Cparameters("int", __name__, "has_outer_boundary", 0)
_unused = par.Cparameters("int", __name__, "outer_bc_type", "EXTRAPOLATION_OUTER_BCS")

# First register basic C data structures/macros inside NRPy_basic_defines.h
def NRPy_basic_defines_CurviBC_data_structures():
    return r"""
// NRPy+ Curvilinear Boundary Conditions: Core data structures
// Documented in: Tutorial-Start_to_Finish-Curvilinear_BCs.ipynb

#define EXTRAPOLATION_OUTER_BCS 0  // used to identify/specify params.outer_bc_type
#define RADIATION_OUTER_BCS     1  // used to identify/specify params.outer_bc_type

typedef struct __innerpt_bc_struct__ {
  int dstpt;  // dstpt is the 3D grid index IDX3S(i0,i1,i2) of the inner boundary point (i0,i1,i2)
  int srcpt;  // srcpt is the 3D grid index (a la IDX3S) to which the inner boundary point maps
  int8_t parity[10];  // parity[10] is a calculation of dot products for the 10 independent parity types
} innerpt_bc_struct;

typedef struct __outerpt_bc_struct__ {
  short i0,i1,i2;  // the outer boundary point grid index (i0,i1,i2), on the 3D grid
  int8_t FACEX0,FACEX1,FACEX2;  // 1-byte integers that store
  //                               FACEX0,FACEX1,FACEX2 = +1, 0, 0 if on the i0=i0min face,
  //                               FACEX0,FACEX1,FACEX2 = -1, 0, 0 if on the i0=i0max face,
  //                               FACEX0,FACEX1,FACEX2 =  0,+1, 0 if on the i1=i2min face,
  //                               FACEX0,FACEX1,FACEX2 =  0,-1, 0 if on the i1=i1max face,
  //                               FACEX0,FACEX1,FACEX2 =  0, 0,+1 if on the i2=i2min face, or
  //                               FACEX0,FACEX1,FACEX2 =  0, 0,-1 if on the i2=i2max face,
} outerpt_bc_struct;

typedef struct __bc_info_struct__ {
  int num_inner_boundary_points;  // stores total number of inner boundary points
  int num_pure_outer_boundary_points[NGHOSTS][3];  // stores number of outer boundary points on each
  //                                                  ghostzone level and direction (update min and
  //                                                  max faces simultaneously on multiple cores)
  int bc_loop_bounds[NGHOSTS][6][6];  // stores outer boundary loop bounds. Unused after bcstruct_set_up()
} bc_info_struct;

typedef struct __bc_struct__ {
  innerpt_bc_struct *restrict inner_bc_array;  // information needed for updating each inner boundary point
  outerpt_bc_struct *restrict pure_outer_bc_array[NGHOSTS*3]; // information needed for updating each outer
  //                                                             boundary point
  bc_info_struct bc_info;  // stores number of inner and outer boundary points, needed for setting loop
  //                          bounds and parallelizing over as many boundary points as possible.
} bc_struct;
"""


# Set unit-vector dot products (=parity) for each of the 10 parity condition types
def parity_conditions_symbolic_dot_products():
    parity = ixp.zerorank1(DIM=10)
    UnitVectors_inner = ixp.zerorank2()
    xx0_inbounds,xx1_inbounds,xx2_inbounds = sp.symbols("xx0_inbounds xx1_inbounds xx2_inbounds", real=True)
    for i in range(3):
        for j in range(3):
            UnitVectors_inner[i][j] = rfm.UnitVectors[i][j].subs(rfm.xx[0],xx0_inbounds).subs(rfm.xx[1],xx1_inbounds).subs(rfm.xx[2],xx2_inbounds)
    # Type 0: scalar
    parity[0] = sp.sympify(1)
    # Type 1: i0-direction vector or one-form
    # Type 2: i1-direction vector or one-form
    # Type 3: i2-direction vector or one-form
    for i in range(3):
        for Type in range(1,4):
            parity[Type] += rfm.UnitVectors[Type-1][i]*UnitVectors_inner[Type-1][i]
    # Type 4: i0i0-direction rank-2 tensor
    # parity[4] = parity[1]*parity[1]
    # Type 5: i0i1-direction rank-2 tensor
    # Type 6: i0i2-direction rank-2 tensor
    # Type 7: i1i1-direction rank-2 tensor
    # Type 8: i1i2-direction rank-2 tensor
    # Type 9: i2i2-direction rank-2 tensor
    count = 4
    for i in range(3):
        for j in range(i,3):
            parity[count] = parity[i+1]*parity[j+1]
            count = count + 1

    lhs_strings = []
    for i in range(10):
        lhs_strings.append("REAL_parity_array["+str(i)+"]")
    outstr = """
    // NRPy+ Curvilinear Boundary Conditions: Unit vector dot products for all
    //      ten parity conditions, in given coordinate system.
    //      Needed for automatically determining sign of tensor across coordinate boundary.
    // Documented in: Tutorial-Start_to_Finish-Curvilinear_BCs.ipynb
"""
    return outstr + outputC(parity, lhs_strings, filename="returnstring", params="preindent=2")


# For example, if the gridfunction name ends with "01", then (based on the table in the
# NRPy+ Jupyter notebook corresponding to this Python module) the set_parity_types()
# function below will set the parity_type of that gridfunction to 5. We can be assured
# this is a rather robust algorithm, because gri.register_gridfunctions() in grid.py
# will throw an error if a gridfunction's base name ends in an integer. This strict
# syntax was added with the express purpose of making it easier to set parity types
# based solely on the gridfunction name.
#
# After each parity type is found, we store the parity type of each gridfunction to
# const int8_t arrays evol_gf_parity and aux_gf_parity, appended to the end of
# NRPy_basic_defines.h.
def NRPy_basic_defines_set_gridfunction_defines_with_parity_types(verbose=True):
    # First add human-readable gridfunction aliases (grid.py) to NRPy_basic_defines dictionary,
    evolved_variables_list, auxiliary_variables_list, auxevol_variables_list = gri.gridfunction_lists()

    # Step 3.b: set the parity conditions on all gridfunctions in gf_list,
    #       based on how many digits are at the end of their names
    def set_parity_types(list_of_gf_names):
        parity_type = []
        for name in list_of_gf_names:
            for gf in gri.glb_gridfcs_list:
                if gf.name == name:
                    parity_type__orig_len = len(parity_type)
                    if gf.DIM < 3 or gf.DIM > 4:
                        print("Error: Cannot currently specify parity conditions on gridfunctions with DIM<3 or >4.")
                        sys.exit(1)
                    if gf.rank == 0:
                        parity_type.append(0)
                    elif gf.rank == 1:
                        if gf.DIM == 3:
                            parity_type.append(int(gf.name[-1]) + 1)  # = 1 for e.g., beta^0; = 2 for e.g., beta^1, etc.
                        elif gf.DIM == 4:
                            parity_type.append(int(gf.name[-1]))  # = 0 for e.g., b4^0; = 1 for e.g., beta^1, etc.
                    elif gf.rank == 2:
                        if gf.DIM == 3:
                            # element of a list; a[-2] the
                            # second-to-last element, etc.
                            idx0 = gf.name[-2]
                            idx1 = gf.name[-1]
                            if idx0 == "0" and idx1 == "0":
                                parity_type.append(4)
                            elif (idx0 == "0" and idx1 == "1") or (idx0 == "1" and idx1 == "0"):
                                parity_type.append(5)
                            elif (idx0 == "0" and idx1 == "2") or (idx0 == "2" and idx1 == "0"):
                                parity_type.append(6)
                            elif idx0 == "1" and idx1 == "1":
                                parity_type.append(7)
                            elif (idx0 == "1" and idx1 == "2") or (idx0 == "2" and idx1 == "1"):
                                parity_type.append(8)
                            elif idx0 == "2" and idx1 == "2":
                                parity_type.append(9)
                        elif gf.DIM == 4:
                            idx0 = gf.name[-2]
                            idx1 = gf.name[-1]
                            # g4DD00 = g_{tt} : parity type = 0
                            # g4DD01 = g_{tx} : parity type = 1
                            # g4DD02 = g_{ty} : parity type = 2
                            # g4DD0a = g_{ta} : parity type = a
                            if idx0 == "0":
                                parity_type.append(int(idx1))
                            elif idx1 == "0":
                                parity_type.append(int(idx0))
                            if idx0 == "1" and idx1 == "1":
                                parity_type.append(4)
                            elif (idx0 == "1" and idx1 == "2") or (idx0 == "2" and idx1 == "1"):
                                parity_type.append(5)
                            elif (idx0 == "1" and idx1 == "3") or (idx0 == "3" and idx1 == "1"):
                                parity_type.append(6)
                            elif idx0 == "2" and idx1 == "2":
                                parity_type.append(7)
                            elif (idx0 == "2" and idx1 == "3") or (idx0 == "3" and idx1 == "2"):
                                parity_type.append(8)
                            elif idx0 == "3" and idx1 == "3":
                                parity_type.append(9)
                    if len(parity_type) == parity_type__orig_len:
                        print("Error: Could not figure out parity type for "+gf.gftype+" gridfunction: " + gf.name,gf.DIM,gf.name[-2],gf.name[-1],gf.rank)
                        sys.exit(1)
        if len(parity_type) != len(list_of_gf_names):
            print("Error: For some reason the length of the parity types list did not match the length of the gf list.")
            sys.exit(1)
        return parity_type

    evol_parity_type = set_parity_types(evolved_variables_list)
    aux_parity_type = set_parity_types(auxiliary_variables_list)
    auxevol_parity_type = set_parity_types(auxevol_variables_list)

    # Output all gridfunctions to Ccodesrootdir/gridfunction_defines.h
    # ... then append to the file the parity type for each gridfunction.
    outstr = """
/* PARITY TYPES FOR ALL GRIDFUNCTIONS.
 * SEE \"Tutorial-Start_to_Finish-Curvilinear_BCs.ipynb\" FOR DEFINITIONS. */
"""
    if len(evolved_variables_list) > 0:
        outstr += "static const int8_t evol_gf_parity[" + str(len(evolved_variables_list)) + "] = { "
        for i in range(len(evolved_variables_list) - 1):
            outstr += str(evol_parity_type[i]) + ", "
        outstr += str(evol_parity_type[len(evolved_variables_list) - 1]) + " };\n"

    if len(auxiliary_variables_list) > 0:
        outstr += "static const int8_t aux_gf_parity[" + str(len(auxiliary_variables_list)) + "] = { "
        for i in range(len(auxiliary_variables_list) - 1):
            outstr += str(aux_parity_type[i]) + ", "
        outstr += str(aux_parity_type[len(auxiliary_variables_list) - 1]) + " };\n"

    if len(auxevol_variables_list) > 0:
        outstr += "static const int8_t auxevol_gf_parity[" + str(len(auxevol_variables_list)) + "] = { "
        for i in range(len(auxevol_variables_list) - 1):
            outstr += str(auxevol_parity_type[i]) + ", "
        outstr += str(auxevol_parity_type[len(auxevol_variables_list) - 1]) + " };\n"

    if verbose == True:
        for i in range(len(evolved_variables_list)):
            print("Evolved gridfunction \"" + evolved_variables_list[i] + "\" has parity type " + str(
                evol_parity_type[i]) + ".")
        for i in range(len(auxiliary_variables_list)):
            print("Auxiliary gridfunction \"" + auxiliary_variables_list[i] + "\" has parity type " + str(
                aux_parity_type[i]) + ".")
        for i in range(len(auxevol_variables_list)):
            print("AuxEvol gridfunction \"" + auxevol_variables_list[i] + "\" has parity type " + str(
                auxevol_parity_type[i]) + ".")
    return outstr


# EigenCoord_set_x0x1x2_inbounds__i0i1i2_inbounds_single_pt():
#      This function is documented in desc= and body= fields below.
def Cfunction__EigenCoord_set_x0x1x2_inbounds__i0i1i2_inbounds_single_pt():
    desc = """EigenCoord_set_x0x1x2_inbounds__i0i1i2_inbounds_single_pt():
  A coordinate system's "eigencoordinate" is the simplest member
  of its family; all spherical-like coordinate systems have
  Spherical as their eigencoordinate. The same is true for
  cylindrical-like (Cylindrical is eigencoordinate),
  Cartesian-like (Cartesian is the eigencoordinate), and
  SymTP-like (SymTP is the eigencoordinate) coordinates.

  For a given gridpoint (i0,i1,i2) and corresponding coordinate
  (x0,x1,x2), this function performs the dual mapping
  (x0,x1,x2) -> (Cartx,Carty,Cartz) -> (x0,x1,x2)'
  Note that (x0,x1,x2) IS NOT ALWAYS equal to (x0,x1,x2)';
  For example consider in Spherical coordinates
  (x0,x1,x2)=(r,theta,phi)=(-0.1,pi/4,pi/4).
  This point will map to (x0,x1,x2)', in which x0>0,
  because the inversion r=sqrt(Cartx^2+Carty^2+Cartz^2)
  is always positive. In this case, (x0,x1,x2) is considered
  an *inner* boundary point, and on a cell-centered grid
  is guaranteed to map to a grid point in the grid interior;
  filling in this point requires copying data, and possibly
  multiplying by a +/- 1 if the data is from a gridfunction
  storing tensors/vectors.
"""
    c_type = "void"
    name = "EigenCoord_set_x0x1x2_inbounds__i0i1i2_inbounds_single_pt"
    params = """const paramstruct *restrict params, REAL *restrict xx[3],
                                                               const int i0, const int i1, const int i2,
                                                               REAL x0x1x2_inbounds[3], int i0i1i2_inbounds[3]"""
    body = r"""
  // This is a 3-step algorithm:
  // Step 1: (x0,x1,x2) -> (Cartx,Carty,Cartz)
  //         Find the Cartesian coordinate that (x0,x1,x2)
  //         maps to, assuming (x0,x1,x2) is the eigen-
  //         coordinate. Note that we assume (x0,x1,x2)
  //         has the same grid boundaries in both the
  //         original coordinate and the eigencoordinate.
  // Step 2: (Cartx,Carty,Cartz) -> (x0,x1,x2)'
  //         Find the interior eigencoordinate point
  //         (x0,x1,x2)' to which (Cartx,Carty,Cartz)
  //         maps, as well as the corresponding
  //         gridpoint integer index (i0,i1,i2). For
  //         cell-centered grids, (x0,x1,x2) will always
  //         overlap exactly (to roundoff error) a point
  //         on the numerical grid.
  // Step 3: Sanity check
  //         Convert x0(i0_inbounds),x1(i1_inbounds),x2(i2_inbounds) -> (Cartx,Carty,Cartz),
  //         and check that
  //         (Cartx,Carty,Cartz) == (Cartx(x0(i0)),Cartx(x1(i1)),Cartx(x2(i2)))
  //         If not, error out!
"""
    # Load up the EigenCoordinate corresponding to reference_metric::CoordSystem
    CoordSystem_orig = par.parval_from_str("reference_metric::CoordSystem")
    par.set_parval_from_str("reference_metric::CoordSystem",rfm.get_EigenCoord())
    rfm.reference_metric()

    # Step 1: Output C code for the Eigen-Coordinate mapping from xx->Cartesian':
    body += r"""
  // Step 1: Convert the (curvilinear) coordinate (x0,x1,x2) to Cartesian coordinates
  REAL xCart[3];  // where (x,y,z) is output
  {
    // xx_to_Cart for EigenCoordinate """+rfm.get_EigenCoord()+r""" (orig coord = """+CoordSystem_orig+r"""):
    REAL xx0 = xx[0][i0];
    REAL xx1 = xx[1][i1];
    REAL xx2 = xx[2][i2];
"""+ \
    outputC([rfm.xx_to_Cart[0],rfm.xx_to_Cart[1],rfm.xx_to_Cart[2]],
            ["xCart[0]","xCart[1]","xCart[2]"],
            "returnstring", params="preindent=2")+"  }\n"
    body += r"""
  REAL Cartx = xCart[0];
  REAL Carty = xCart[1];
  REAL Cartz = xCart[2];"""

    # Step 2: Output C code for the Eigen-Coordinate mapping from Cartesian->xx':
    body += r"""
  // Step 2: Find the (i0_inbounds,i1_inbounds,i2_inbounds) corresponding to the above Cartesian coordinate.
  //   If (i0_inbounds,i1_inbounds,i2_inbounds) is in a ghost zone, then it must equal (i0,i1,i2), and
  //      the point is an outer boundary point.
  //   Otherwise (i0_inbounds,i1_inbounds,i2_inbounds) is in the grid interior, and data at (i0,i1,i2)
  //      must be replaced with data at (i0_inbounds,i1_inbounds,i2_inbounds), but multiplied by the
  //      appropriate parity condition (+/- 1).
  REAL Cart_to_xx0_inbounds,Cart_to_xx1_inbounds,Cart_to_xx2_inbounds;
"""
    # Step 2.a: Sanity check: First make sure that rfm.Cart_to_xx has been set. Error out if not!
    if rfm.Cart_to_xx[0] == 0 or rfm.Cart_to_xx[1] == 0 or rfm.Cart_to_xx[2] == 0:
        print("ERROR: rfm.Cart_to_xx[], which maps Cartesian -> xx, has not been set for")
        print("       reference_metric::CoordSystem = "+par.parval_from_str("reference_metric::CoordSystem"))
        print("       Boundary conditions in curvilinear coordinates REQUiRE this be set.")
        sys.exit(1)
    # Step 2.b: Output C code for the Eigen-Coordinate mapping from Cartesian->xx:
    body += """  // Cart_to_xx for EigenCoordinate """+rfm.get_EigenCoord()+r""" (orig coord = """+CoordSystem_orig+");\n"
    body += outputC([rfm.Cart_to_xx[0],rfm.Cart_to_xx[1],rfm.Cart_to_xx[2]],
                    ["Cart_to_xx0_inbounds","Cart_to_xx1_inbounds","Cart_to_xx2_inbounds"],
                    filename="returnstring", params="preindent=1")
    body += r"""
  // Next compute xxmin[i]. By definition,
  //    xx[i][j] = xxmin[i] + ((REAL)(j-NGHOSTS) + (1.0/2.0))*dxxi;
  // -> xxmin[i] = xx[i][0] - ((REAL)(0-NGHOSTS) + (1.0/2.0))*dxxi
  const REAL xxmin[3] = {
    xx[0][0] - ((REAL)(0-NGHOSTS) + (1.0/2.0))*dxx0,
    xx[1][0] - ((REAL)(0-NGHOSTS) + (1.0/2.0))*dxx1,
    xx[2][0] - ((REAL)(0-NGHOSTS) + (1.0/2.0))*dxx2 };

  // Finally compute i{0,1,2}_inbounds (add 0.5 to account for rounding down)
  const int i0_inbounds = (int)( (Cart_to_xx0_inbounds - xxmin[0] - (1.0/2.0)*dxx0 + ((REAL)NGHOSTS)*dxx0)/dxx0 + 0.5 );
  const int i1_inbounds = (int)( (Cart_to_xx1_inbounds - xxmin[1] - (1.0/2.0)*dxx1 + ((REAL)NGHOSTS)*dxx1)/dxx1 + 0.5 );
  const int i2_inbounds = (int)( (Cart_to_xx2_inbounds - xxmin[2] - (1.0/2.0)*dxx2 + ((REAL)NGHOSTS)*dxx2)/dxx2 + 0.5 );
"""

    # Restore reference_metric::CoordSystem back to the original CoordSystem
    par.set_parval_from_str("reference_metric::CoordSystem",CoordSystem_orig)
    rfm.reference_metric()

    # Step 3:
    body += """
  // Step 3: Convert x0(i0_inbounds),x1(i1_inbounds),x2(i2_inbounds) -> (Cartx,Carty,Cartz),
  //         and check that
  //         (Cartx,Carty,Cartz) == (Cartx(x0(i0)),Cartx(x1(i1)),Cartx(x2(i2)))
  //         If not, error out!

  // Step 3.a: Compute x0(i0_inbounds),x1(i1_inbounds),x2(i2_inbounds):
  const REAL x0_inbounds = xx[0][i0_inbounds];
  const REAL x1_inbounds = xx[1][i1_inbounds];
  const REAL x2_inbounds = xx[2][i2_inbounds];

  // Step 3.b: Compute {x,y,z}Cart_from_xx, as a
  //           function of i0,i1,i2
  REAL xCart_from_xx, yCart_from_xx, zCart_from_xx;
  {
    // xx_to_Cart for Coordinate """+CoordSystem_orig+r"""):
    REAL xx0 = xx[0][i0];
    REAL xx1 = xx[1][i1];
    REAL xx2 = xx[2][i2];
"""+ \
    outputC([rfm.xx_to_Cart[0],rfm.xx_to_Cart[1],rfm.xx_to_Cart[2]],
            ["xCart_from_xx","yCart_from_xx","zCart_from_xx"],
            "returnstring", params="preindent=2,includebraces=False")
    body += r"""  }

  // Step 3.c: Compute {x,y,z}Cart_from_xx_inbounds, as a
  //           function of i0_inbounds,i1_inbounds,i2_inbounds
  REAL xCart_from_xx_inbounds, yCart_from_xx_inbounds, zCart_from_xx_inbounds;
  {
    // xx_to_Cart_inbounds for Coordinate """+CoordSystem_orig+r"""):
    REAL xx0 = xx[0][i0_inbounds];
    REAL xx1 = xx[1][i1_inbounds];
    REAL xx2 = xx[2][i2_inbounds];
"""+ \
    outputC([rfm.xx_to_Cart[0],rfm.xx_to_Cart[1],rfm.xx_to_Cart[2]],
            ["xCart_from_xx_inbounds","yCart_from_xx_inbounds","zCart_from_xx_inbounds"],
            "returnstring", params="preindent=2,includebraces=False")
    body += r"""  }

  // Step 3.d: Compare xCart_from_xx to xCart_from_xx_inbounds;
  //           they should be identical!!!
#define EPS_REL 1e-8
  const REAL norm_factor = sqrt(xCart_from_xx*xCart_from_xx + yCart_from_xx*yCart_from_xx + zCart_from_xx*zCart_from_xx) + 1e-15;
  if(fabs( (double)(xCart_from_xx - xCart_from_xx_inbounds) ) > EPS_REL * norm_factor ||
     fabs( (double)(yCart_from_xx - yCart_from_xx_inbounds) ) > EPS_REL * norm_factor ||
     fabs( (double)(zCart_from_xx - zCart_from_xx_inbounds) ) > EPS_REL * norm_factor) {
    fprintf(stderr,"Error in """+CoordSystem_orig+r""" coordinate system: Cartesian disagreement: ( %.15e %.15e %.15e ) != ( %.15e %.15e %.15e ) | xx: %e %e %e -> %e %e %e | %d %d %d\n",
            (double)xCart_from_xx,(double)yCart_from_xx,(double)zCart_from_xx,
            (double)xCart_from_xx_inbounds,(double)yCart_from_xx_inbounds,(double)zCart_from_xx_inbounds,
            xx[0][i0],xx[1][i1],xx[2][i2],
            xx[0][i0_inbounds],xx[1][i1_inbounds],xx[2][i2_inbounds],
            Nxx_plus_2NGHOSTS0,Nxx_plus_2NGHOSTS1,Nxx_plus_2NGHOSTS2);
    exit(1);
  }

  // Step 4: Set output arrays.
  x0x1x2_inbounds[0] = xx[0][i0_inbounds];
  x0x1x2_inbounds[1] = xx[1][i1_inbounds];
  x0x1x2_inbounds[2] = xx[2][i2_inbounds];
  i0i1i2_inbounds[0] = i0_inbounds;
  i0i1i2_inbounds[1] = i1_inbounds;
  i0i1i2_inbounds[2] = i2_inbounds;
"""
    __function_prototype_ignore, Cfunc = Cfunction(
        desc=desc, c_type=c_type, name=name, params=params,
        body=body)
    return Cfunc


# set_parity_for_inner_boundary_single_pt():
#      This function is documented in desc= and body= fields below.
def Cfunction__set_parity_for_inner_boundary_single_pt():
    desc = """set_parity_for_inner_boundary_single_pt():
  Given (x0,x1,x2)=(xx0,xx1,xx2) and
  (x0,x1,x2)'=(x0x1x2_inbounds[0],x0x1x2_inbounds[1],x0x1x2_inbounds[2])
  (see description of
  EigenCoord_set_x0x1x2_inbounds__i0i1i2_inbounds_single_pt()
  above for more details), here we compute the parity conditions
  for all 10 tensor types supported by NRPy+.
"""
    c_type = "void"
    name = "set_parity_for_inner_boundary_single_pt"
    params = """const paramstruct *restrict params, const REAL xx0,const REAL xx1,const REAL xx2,
                                             const REAL x0x1x2_inbounds[3], const int idx,
                                             innerpt_bc_struct *restrict innerpt_bc_arr"""
    body = r"""
  const REAL xx0_inbounds = x0x1x2_inbounds[0];
  const REAL xx1_inbounds = x0x1x2_inbounds[1];
  const REAL xx2_inbounds = x0x1x2_inbounds[2];

  REAL REAL_parity_array[10];
  {
    // Evaluate dot products needed for setting parity
    //     conditions at a given point (xx0,xx1,xx2),
    //     using C code generated by NRPy+
""" + parity_conditions_symbolic_dot_products() + r"""
  }

  // Next perform sanity check on parity array output: should be +1 or -1 to within 8 significant digits:
  for(int whichparity=0;whichparity<10;whichparity++) {
    //printf("Good? Parity %d evaluated to %e\n",whichparity,(double)REAL_parity_array[whichparity]);
    if( fabs(REAL_parity_array[whichparity]) < 1 - 1e-8 || fabs(REAL_parity_array[whichparity]) > 1 + 1e-8 ) {
      fprintf(stderr,"Error at point (%e %e %e), which maps to (%e %e %e).\n",
              xx0,xx1,xx2, xx0_inbounds,xx1_inbounds,xx2_inbounds);
      fprintf(stderr,"Parity evaluated to %e , which is not within 8 significant digits of +1 or -1.\n",
              REAL_parity_array[whichparity]);
      exit(1);
    }
    // The typecast (int8_t)REAL_parity_array[whichparity] *does not work*.
    //  Thankfully, we've already checked whether REAL_parity_array[whichparity]
    //  is within 1e-8 of +/- 1, so here we just check the sign of the
    //  REAL_parity_array to find the correct value of innerpt_bc_arr[idx].parity[parity].
    for(int parity=0;parity<10;parity++) {
      innerpt_bc_arr[idx].parity[parity] = 1;
      if(REAL_parity_array[parity] < 0) innerpt_bc_arr[idx].parity[parity] = -1;
    }
  } // END for(int whichparity=0;whichparity<10;whichparity++)
"""
    __function_prototype_ignore, Cfunc = Cfunction(
        desc=desc, c_type=c_type, name=name, params=params,
        body=body)
    return Cfunc


# bcstruct_set_up():
#      This function is documented in desc= and body= fields below.
def add_to_Cfunction_dict_bcstruct_set_up(rel_path_to_Cparams=os.path.join(".")):
    includes = [os.path.join(rel_path_to_Cparams, "NRPy_basic_defines.h"),
                os.path.join(rel_path_to_Cparams, "NRPy_function_prototypes.h")]
    prefunc  = Cfunction__EigenCoord_set_x0x1x2_inbounds__i0i1i2_inbounds_single_pt()
    prefunc += Cfunction__set_parity_for_inner_boundary_single_pt()
    desc = r"""At each coordinate point (x0,x1,x2) situated at grid index (i0,i1,i2):

Step 1: Set up inner boundary structs bcstruct->inner_bc_array[].
  Recall that at each inner boundary point we must set innerpt_bc_struct:
    typedef struct __innerpt_bc_struct__ {
      int dstpt;  // dstpt is the 3D grid index IDX3S(i0,i1,i2) of the inner boundary point (i0,i1,i2)
      int srcpt;  // srcpt is the 3D grid index (a la IDX3S) to which the inner boundary point maps
      int8_t parity[10];  // parity[10] is a calculation of dot products for the 10 independent parity types
    } innerpt_bc_struct;
  At each ghostzone (i.e., each point within NGHOSTS points from grid boundary):
    Call EigenCoord_set_x0x1x2_inbounds__i0i1i2_inbounds_single_pt().
        This function converts the curvilinear coordinate (x0,x1,x2) to the corresponding
        Cartesian coordinate (x,y,z), then finds the grid point
        (i0_inbounds,i1_inbounds,i2_inbounds) in the grid interior or outer boundary
        corresponding to this Cartesian coordinate (x,y,z).
    If (i0,i1,i2) *is not* the same as (i0_inbounds,i1_inbounds,i2_inbounds),
        then we are at an inner boundary point. We must set
        Set bcstruct->inner_bc_array for this point, which requires we specify
        both (i0_inbounds,i1_inbounds,i2_inbounds) [just found!] and parity
        conditions for this gridpoint. The latter is found & specified within the
        function set_parity_for_inner_boundary_single_pt().
    If (i0,i1,i2) *is* the same as (i0_inbounds,i1_inbounds,i2_inbounds),
        then we are at an outer boundary point. Take care of outer BCs in Step 2.
Step 2: Set up outer boundary structs bcstruct->outer_bc_array[which_gz][face][idx2d]:
  Recall that at each inner boundary point we must set outerpt_bc_struct:
    typedef struct __outerpt_bc_struct__ {
      short i0,i1,i2;  // the outer boundary point grid index (i0,i1,i2), on the 3D grid
      int8_t FACEX0,FACEX1,FACEX2;  // 1-byte integers that store
      //                               FACEX0,FACEX1,FACEX2 = +1, 0, 0 if on the i0=i0min face,
      //                               FACEX0,FACEX1,FACEX2 = -1, 0, 0 if on the i0=i0max face,
      //                               FACEX0,FACEX1,FACEX2 =  0,+1, 0 if on the i1=i2min face,
      //                               FACEX0,FACEX1,FACEX2 =  0,-1, 0 if on the i1=i1max face,
      //                               FACEX0,FACEX1,FACEX2 =  0, 0,+1 if on the i2=i2min face, or
      //                               FACEX0,FACEX1,FACEX2 =  0, 0,-1 if on the i2=i2max face,
    } outerpt_bc_struct;
  Outer boundary points are filled from the inside out, two faces at a time.
    E.g., consider a Cartesian coordinate grid that has 14 points in each direction,
    including the ghostzones, with NGHOSTS=2.
    We first fill in the lower x0 face with (i0=1,i1={2,11},i2={2,11}). We fill these
    points in first, since they will in general (at least in the case of extrapolation
    outer BCs) depend on e.g., i0=2 and i0=3 points.
    Simultaneously we can fill in the upper x0 face with (i0=12,i1={2,11},i2={2,11}),
    since these points depend only on e.g., i0=11 and i0=10 (again assuming extrap. BCs).
    Next we can fill in the lower x1 face: (i0={1,12},i1=2,i2={2,11}). Notice these
    depend on i0 min and max faces being filled. The remaining pattern goes like this:
    Upper x1 face: (i0={1,12},i1=12,i2={2,11})
    Lower x2 face: (i0={1,12},i1={1,12},i2=1)
    Upper x2 face: (i0={1,12},i1={1,12},i2=12)
    Lower x0 face: (i0=0,i1={1,12},i2={1,12})
    Upper x0 face: (i0=13,i1={1,12},i2={1,12})
    Lower x1 face: (i0={0,13},i1=0,i2={2,11})
    Upper x1 face: (i0={0,13},i1=13,i2={2,11})
    Lower x2 face: (i0={0,13},i1={0,13},i2=0)
    Upper x2 face: (i0={0,13},i1={0,13},i2=13)
  Note that we allocate a outerpt_bc_struct at *all* boundary points,
    regardless of whether the point is an outer or inner point. However
    the struct is set only at outer boundary points. This is slightly
    wasteful, but only in memory, not in CPU.
"""
    c_type = "void"
    name = "bcstruct_set_up"
    params = "const paramstruct *restrict params, REAL *restrict xx[3], bc_struct *restrict bcstruct"
    body = r"""
  ////////////////////////////////////////
  // STEP 1: SET UP INNER BOUNDARY STRUCTS
  {
    // First count the number of inner points.
    int num_inner = 0;
    LOOP_OMP("omp parallel for reduction(+:num_inner)",
             i0,0,Nxx_plus_2NGHOSTS0,  i1,0,Nxx_plus_2NGHOSTS1,  i2,0,Nxx_plus_2NGHOSTS2) {
      const int i0i1i2[3] = { i0,i1,i2 };
      if(!IS_IN_GRID_INTERIOR(i0i1i2, Nxx_plus_2NGHOSTS0,Nxx_plus_2NGHOSTS1,Nxx_plus_2NGHOSTS2, NGHOSTS)) {
        REAL x0x1x2_inbounds[3];
        int i0i1i2_inbounds[3];
        EigenCoord_set_x0x1x2_inbounds__i0i1i2_inbounds_single_pt(params, xx, i0,i1,i2, x0x1x2_inbounds,i0i1i2_inbounds);
        if(i0 == i0i1i2_inbounds[0] && i1==i0i1i2_inbounds[1] && i2==i0i1i2_inbounds[2]) {
          // this is a pure outer boundary point.
        } else {
          // this is an inner boundary point, which maps either
          //  to the grid interior or to an outer boundary point
          num_inner++;
        }
      }
    }
    // Store num_inner to bc_info:
    bcstruct->bc_info.num_inner_boundary_points = num_inner;

    // Next allocate memory for inner_boundary_points:
    bcstruct->inner_bc_array = (innerpt_bc_struct *restrict)malloc( sizeof(innerpt_bc_struct)*num_inner );
  }

  // Then set inner_bc_array:
  {
    int which_inner = 0;
    LOOP_NOOMP(i0,0,Nxx_plus_2NGHOSTS0,  i1,0,Nxx_plus_2NGHOSTS1,  i2,0,Nxx_plus_2NGHOSTS2) {
      const int i0i1i2[3] = { i0,i1,i2 };
      if(!IS_IN_GRID_INTERIOR(i0i1i2, Nxx_plus_2NGHOSTS0,Nxx_plus_2NGHOSTS1,Nxx_plus_2NGHOSTS2, NGHOSTS)) {
        REAL x0x1x2_inbounds[3];
        int i0i1i2_inbounds[3];
        EigenCoord_set_x0x1x2_inbounds__i0i1i2_inbounds_single_pt(params, xx, i0,i1,i2, x0x1x2_inbounds,i0i1i2_inbounds);
        if(i0 == i0i1i2_inbounds[0] && i1==i0i1i2_inbounds[1] && i2==i0i1i2_inbounds[2]) {
          // this is a pure outer boundary point.
        } else {
          bcstruct->inner_bc_array[which_inner].dstpt = IDX3S(i0,i1,i2);
          bcstruct->inner_bc_array[which_inner].srcpt = IDX3S(i0i1i2_inbounds[0],i0i1i2_inbounds[1],i0i1i2_inbounds[2]);
          //printf("%d / %d\n",which_inner, bc_info->num_inner_boundary_points);
          set_parity_for_inner_boundary_single_pt(params, xx[0][i0],xx[1][i1],xx[2][i2],
                                                  x0x1x2_inbounds, which_inner, bcstruct->inner_bc_array);

          which_inner++;
        }
      }
    }
  }

  ////////////////////////////////////////
  // STEP 2: SET UP OUTER BOUNDARY STRUCTS
  // First set up loop bounds for outer boundary condition updates,
  //   store to bc_info->bc_loop_bounds[which_gz][face][]. Also
  //   allocate memory for outer_bc_array[which_gz][face][]:
  int imin[3] = { NGHOSTS, NGHOSTS, NGHOSTS };
  int imax[3] = { Nxx_plus_2NGHOSTS0-NGHOSTS, Nxx_plus_2NGHOSTS1-NGHOSTS, Nxx_plus_2NGHOSTS2-NGHOSTS };
  for(int which_gz=0;which_gz<NGHOSTS;which_gz++) {
    const int x0min_face_range[6] = { imin[0]-1,imin[0], imin[1],imax[1], imin[2],imax[2] };  imin[0]--;
    const int x0max_face_range[6] = { imax[0],imax[0]+1, imin[1],imax[1], imin[2],imax[2] };  imax[0]++;
    const int x1min_face_range[6] = { imin[0],imax[0], imin[1]-1,imin[1], imin[2],imax[2] };  imin[1]--;
    const int x1max_face_range[6] = { imin[0],imax[0], imax[1],imax[1]+1, imin[2],imax[2] };  imax[1]++;
    const int x2min_face_range[6] = { imin[0],imax[0], imin[1],imax[1], imin[2]-1,imin[2] };  imin[2]--;
    const int x2max_face_range[6] = { imin[0],imax[0], imin[1],imax[1], imax[2],imax[2]+1 };  imax[2]++;

    int face=0;
    ////////////////////////
    // x0min and x0max faces: Allocate memory for outer_bc_array and set bc_loop_bounds:
    //                        Note that x0min and x0max faces have exactly the same size.
    bcstruct->pure_outer_bc_array[which_gz + NGHOSTS*(face/2)] = (outerpt_bc_struct *restrict)malloc(sizeof(outerpt_bc_struct) * 2 *
                                                                                                     ((x0min_face_range[1]-x0min_face_range[0]) *
                                                                                                      (x0min_face_range[3]-x0min_face_range[2]) *
                                                                                                      (x0min_face_range[5]-x0min_face_range[4])));
    // x0min face: Can't set bc_info->bc_loop_bounds[which_gz][face] = { i0min,i0max, ... } since it's not const :(
    for(int i=0;i<6;i++) { bcstruct->bc_info.bc_loop_bounds[which_gz][face][i] = x0min_face_range[i]; }
    face++;
    // x0max face: Set loop bounds & allocate memory for outer_bc_array:
    for(int i=0;i<6;i++) { bcstruct->bc_info.bc_loop_bounds[which_gz][face][i] = x0max_face_range[i]; }
    face++;
    ////////////////////////

    ////////////////////////
    // x1min and x1max faces: Allocate memory for outer_bc_array and set bc_loop_bounds:
    //                        Note that x1min and x1max faces have exactly the same size.
    bcstruct->pure_outer_bc_array[which_gz + NGHOSTS*(face/2)] = (outerpt_bc_struct *restrict)malloc(sizeof(outerpt_bc_struct) * 2 *
                                                                                                     ((x1min_face_range[1]-x1min_face_range[0]) *
                                                                                                      (x1min_face_range[3]-x1min_face_range[2]) *
                                                                                                      (x1min_face_range[5]-x1min_face_range[4])));
    // x1min face: Can't set bc_info->bc_loop_bounds[which_gz][face] = { i0min,i0max, ... } since it's not const :(
    for(int i=0;i<6;i++) { bcstruct->bc_info.bc_loop_bounds[which_gz][face][i] = x1min_face_range[i]; }
    face++;
    // x1max face: Set loop bounds & allocate memory for outer_bc_array:
    for(int i=0;i<6;i++) { bcstruct->bc_info.bc_loop_bounds[which_gz][face][i] = x1max_face_range[i]; }
    face++;
    ////////////////////////


    ////////////////////////
    // x2min and x2max faces: Allocate memory for outer_bc_array and set bc_loop_bounds:
    //                        Note that x2min and x2max faces have exactly the same size.
    bcstruct->pure_outer_bc_array[which_gz + NGHOSTS*(face/2)] = (outerpt_bc_struct *restrict)malloc(sizeof(outerpt_bc_struct) * 2 *
                                                                                                     ((x2min_face_range[1]-x2min_face_range[0]) *
                                                                                                      (x2min_face_range[3]-x2min_face_range[2]) *
                                                                                                      (x2min_face_range[5]-x2min_face_range[4])));
    // x2min face: Can't set bc_info->bc_loop_bounds[which_gz][face] = { i0min,i0max, ... } since it's not const :(
    for(int i=0;i<6;i++) { bcstruct->bc_info.bc_loop_bounds[which_gz][face][i] = x2min_face_range[i]; }
    face++;
    // x2max face: Set loop bounds & allocate memory for outer_bc_array:
    for(int i=0;i<6;i++) { bcstruct->bc_info.bc_loop_bounds[which_gz][face][i] = x2max_face_range[i]; }
    face++;
    ////////////////////////
  }

  for(int which_gz=0;which_gz<NGHOSTS;which_gz++) for(int dirn=0;dirn<3;dirn++) {
      int idx2d = 0;
      // LOWER FACE: dirn=0 -> x0min; dirn=1 -> x1min; dirn=2 -> x2min
      {
        const int face = dirn*2;
#define IDX2D_BCS(i0,i0min,i0max, i1,i1min,i1max ,i2,i2min,i2max)       \
        ( ((i0)-(i0min)) + ((i0max)-(i0min)) * ( ((i1)-(i1min)) + ((i1max)-(i1min)) * ((i2)-(i2min)) ) )
        const int FACEX0=(face==0) - (face==1); // +1 if face==0 (x0min) ; -1 if face==1 (x0max). Otherwise 0.
        const int FACEX1=(face==2) - (face==3); // +1 if face==2 (x1min) ; -1 if face==3 (x1max). Otherwise 0.
        const int FACEX2=(face==4) - (face==5); // +1 if face==4 (x2min) ; -1 if face==5 (x2max). Otherwise 0.
        LOOP_NOOMP(i0,bcstruct->bc_info.bc_loop_bounds[which_gz][face][0],bcstruct->bc_info.bc_loop_bounds[which_gz][face][1],
                   i1,bcstruct->bc_info.bc_loop_bounds[which_gz][face][2],bcstruct->bc_info.bc_loop_bounds[which_gz][face][3],
                   i2,bcstruct->bc_info.bc_loop_bounds[which_gz][face][4],bcstruct->bc_info.bc_loop_bounds[which_gz][face][5]) {
          REAL x0x1x2_inbounds[3];
          int i0i1i2_inbounds[3];
          EigenCoord_set_x0x1x2_inbounds__i0i1i2_inbounds_single_pt(params, xx, i0,i1,i2, x0x1x2_inbounds,i0i1i2_inbounds);
          if(i0 == i0i1i2_inbounds[0] && i1==i0i1i2_inbounds[1] && i2==i0i1i2_inbounds[2]) {
            bcstruct->pure_outer_bc_array[dirn + (3*which_gz)][idx2d].i0 = i0;
            bcstruct->pure_outer_bc_array[dirn + (3*which_gz)][idx2d].i1 = i1;
            bcstruct->pure_outer_bc_array[dirn + (3*which_gz)][idx2d].i2 = i2;
            bcstruct->pure_outer_bc_array[dirn + (3*which_gz)][idx2d].FACEX0 = FACEX0;
            bcstruct->pure_outer_bc_array[dirn + (3*which_gz)][idx2d].FACEX1 = FACEX1;
            bcstruct->pure_outer_bc_array[dirn + (3*which_gz)][idx2d].FACEX2 = FACEX2;
            idx2d++;
          }
        }
      }
      // UPPER FACE: dirn=0 -> x0max; dirn=1 -> x1max; dirn=2 -> x2max
      {
        const int face = dirn*2+1;
        const int FACEX0=(face==0) - (face==1); // +1 if face==0 ; -1 if face==1. Otherwise 0.
        const int FACEX1=(face==2) - (face==3); // +1 if face==2 ; -1 if face==3. Otherwise 0.
        const int FACEX2=(face==4) - (face==5); // +1 if face==4 ; -1 if face==5. Otherwise 0.
        LOOP_NOOMP(i0,bcstruct->bc_info.bc_loop_bounds[which_gz][face][0],bcstruct->bc_info.bc_loop_bounds[which_gz][face][1],
                   i1,bcstruct->bc_info.bc_loop_bounds[which_gz][face][2],bcstruct->bc_info.bc_loop_bounds[which_gz][face][3],
                   i2,bcstruct->bc_info.bc_loop_bounds[which_gz][face][4],bcstruct->bc_info.bc_loop_bounds[which_gz][face][5]) {
          REAL x0x1x2_inbounds[3];
          int i0i1i2_inbounds[3];
          EigenCoord_set_x0x1x2_inbounds__i0i1i2_inbounds_single_pt(params, xx, i0,i1,i2, x0x1x2_inbounds,i0i1i2_inbounds);
          if(i0 == i0i1i2_inbounds[0] && i1==i0i1i2_inbounds[1] && i2==i0i1i2_inbounds[2]) {
            bcstruct->pure_outer_bc_array[dirn + (3*which_gz)][idx2d].i0 = i0;
            bcstruct->pure_outer_bc_array[dirn + (3*which_gz)][idx2d].i1 = i1;
            bcstruct->pure_outer_bc_array[dirn + (3*which_gz)][idx2d].i2 = i2;
            bcstruct->pure_outer_bc_array[dirn + (3*which_gz)][idx2d].FACEX0 = FACEX0;
            bcstruct->pure_outer_bc_array[dirn + (3*which_gz)][idx2d].FACEX1 = FACEX1;
            bcstruct->pure_outer_bc_array[dirn + (3*which_gz)][idx2d].FACEX2 = FACEX2;
            idx2d++;
          }
        }
      }
      bcstruct->bc_info.num_pure_outer_boundary_points[which_gz][dirn] = idx2d;
    }
"""
    add_to_Cfunction_dict(
        includes=includes,
        prefunc=prefunc,
        desc=desc,
        c_type=c_type, name=name, params=params,
        body=body,
        rel_path_to_Cparams=rel_path_to_Cparams)


###############################
## apply_bcs_inner_only(): Apply inner boundary conditions.
##  Function is documented below in desc= and body=.
def add_to_Cfunction_dict_apply_bcs_inner_only(rel_path_to_Cparams=os.path.join(".")):
    includes = [os.path.join(rel_path_to_Cparams, "NRPy_basic_defines.h")]
    desc = r"""
Apply BCs to inner boundary points only,
using data stored in bcstruct->inner_bc_array.
These structs are set in bcstruct_set_up().
Inner boundary points map to either the grid
interior ("pure inner") or to pure outer
boundary points ("inner maps to outer").
"""
    c_type = "void"
    name = "apply_bcs_inner_only"
    params = "const paramstruct *restrict params, const bc_struct *restrict bcstruct, REAL *restrict gfs"
    body = r"""
  // Unpack bc_info from bcstruct
  const bc_info_struct *bc_info = &bcstruct->bc_info;

#pragma omp parallel for  // spawn threads and distribute across them
  for(int pt=0;pt<bc_info->num_inner_boundary_points;pt++) {
    for(int which_gf=0;which_gf<NUM_EVOL_GFS;which_gf++) {
      const int dstpt = bcstruct->inner_bc_array[pt].dstpt;
      const int srcpt = bcstruct->inner_bc_array[pt].srcpt;

      gfs[IDX4ptS(which_gf, dstpt)] = bcstruct->inner_bc_array[pt].parity[evol_gf_parity[which_gf]] * gfs[IDX4ptS(which_gf, srcpt)];
    } // END for(int pt=0;pt<num_inner_pts;pt++)
  } // END for(int which_gf=0;which_gf<NUM_EVOL_GFS;which_gf++)
"""
    add_to_Cfunction_dict(
        includes=includes,
        desc=desc,
        c_type=c_type, name=name, params=params,
        body=body,
        rel_path_to_Cparams=rel_path_to_Cparams)


###############################
## apply_bcs_outerextrap_and_inner(): Apply extrapolation outer boundary conditions.
##  Function is documented below in desc= and body=.
def add_to_Cfunction_dict_apply_bcs_outerextrap_and_inner(rel_path_to_Cparams=os.path.join(".")):
    includes = [os.path.join(rel_path_to_Cparams, "NRPy_basic_defines.h"),
                os.path.join(rel_path_to_Cparams, "NRPy_function_prototypes.h")]
    desc = r"""
Suppose the outer boundary point is at the i0=max(i0) face. Then we fit known data at i0-3, i0-2, and i0-1
  to the unique quadratic polynomial that passes through those points, and fill the data at
  i0 with the value implied from the polynomial.
As derived in nrpytutorial's Tutorial-Start_to_Finish-Curvilinear_BCs.ipynb,
  the coefficients must be f_{i0} = f_{i0-3} - 3 f_{i0-2} + 3 f_{i0-1}.
  To check these coefficients are correct, consider
  * f(x0 = constant. Then f_{i0} = f_{i0-3} <- CHECK!
  * f(x) = x. WOLOG suppose x0=0. Then f_{i0} = (-3dx) - 3(-2dx) + 3(-dx) = + dx(-3+6-3) = 0 <- CHECK!
  * f(x) = x^2. WOLOG suppose x0=0. Then f_{i0} = (-3dx)^2 - 3(-2dx)^2 + 3(-dx)^2 = + dx^2(9-12+3) = 0 <- CHECK!
"""
    c_type = "void"
    name = "apply_bcs_outerextrap_and_inner"
    params = "const paramstruct *restrict params, const bc_struct *restrict bcstruct, REAL *restrict gfs"
    body = r"""
  // Unpack bc_info from bcstruct
  const bc_info_struct *bc_info = &bcstruct->bc_info;

  ////////////////////////////////////////////////////////
  // STEP 1 of 2: Apply BCs to pure outer boundary points.
  //              By "pure" we mean that these points are
  //              on the outer boundary and not also on
  //              an inner boundary.
  //              Here we fill in the innermost ghost zone
  //              layer first and move outward. At each
  //              layer, we fill in +/- x0 faces first,
  //              then +/- x1 faces, finally +/- x2 faces,
  //              filling in the edges as we go.
  // Spawn N OpenMP threads, either across all cores, or according to e.g., taskset.
#pragma omp parallel
  {
    for(int which_gz=0;which_gz<NGHOSTS;which_gz++) for(int dirn=0;dirn<3;dirn++) {
        // This option results in about 1.6% slower runtime for SW curvilinear at 64x24x24 on 8-core Ryzen 9 4900HS
        //#pragma omp for collapse(2)
        //for(int which_gf=0;which_gf<NUM_EVOL_GFS;which_gf++) for(int idx2d=0;idx2d<bc_info->num_pure_outer_boundary_points[which_gz][dirn];idx2d++) {
        //  {
        // Don't spawn a thread if there are no boundary points to fill; results in a nice little speedup.
        if(bc_info->num_pure_outer_boundary_points[which_gz][dirn] > 0) {
#pragma omp for  // threads have been spawned; here we distribute across them
          for(int idx2d=0;idx2d<bc_info->num_pure_outer_boundary_points[which_gz][dirn];idx2d++) {
            const short i0 = bcstruct->pure_outer_bc_array[dirn + (3*which_gz)][idx2d].i0;
            const short i1 = bcstruct->pure_outer_bc_array[dirn + (3*which_gz)][idx2d].i1;
            const short i2 = bcstruct->pure_outer_bc_array[dirn + (3*which_gz)][idx2d].i2;
            const short FACEX0 = bcstruct->pure_outer_bc_array[dirn + (3*which_gz)][idx2d].FACEX0;
            const short FACEX1 = bcstruct->pure_outer_bc_array[dirn + (3*which_gz)][idx2d].FACEX1;
            const short FACEX2 = bcstruct->pure_outer_bc_array[dirn + (3*which_gz)][idx2d].FACEX2;
            const int idx_offset0 = IDX3S(i0,i1,i2);
            const int idx_offset1 = IDX3S(i0+1*FACEX0,i1+1*FACEX1,i2+1*FACEX2);
            const int idx_offset2 = IDX3S(i0+2*FACEX0,i1+2*FACEX1,i2+2*FACEX2);
            const int idx_offset3 = IDX3S(i0+3*FACEX0,i1+3*FACEX1,i2+3*FACEX2);
            for(int which_gf=0;which_gf<NUM_EVOL_GFS;which_gf++) {
              // *** Apply 2nd-order polynomial extrapolation BCs to all outer boundary points. ***
              gfs[IDX4ptS(which_gf, idx_offset0)] =
                +3.0*gfs[IDX4ptS(which_gf, idx_offset1)]
                -3.0*gfs[IDX4ptS(which_gf, idx_offset2)]
                +1.0*gfs[IDX4ptS(which_gf, idx_offset3)];
            }
          }
        }
      }
  }

  ///////////////////////////////////////////////////////
  // STEP 2 of 2: Apply BCs to inner boundary points.
  //              These map to either the grid interior
  //              ("pure inner") or to pure outer boundary
  //              points ("inner maps to outer"). Those
  //              that map to outer require that outer be
  //              populated first; hence this being
  //              STEP 2 OF 2.
  apply_bcs_inner_only(params, bcstruct, gfs);
"""
    add_to_Cfunction_dict(
        includes=includes,
        desc=desc,
        c_type=c_type, name=name, params=params,
        body=body,
        rel_path_to_Cparams=rel_path_to_Cparams)


###############################
## RADIATION (NewRad-like) BOUNDARY CONDITIONS.
##  Functions are fully documented in nrpytutorial's
##   Tutorial-Start_to_Finish-Curvilinear_BCs.ipynb,
##   as well as below, in desc= and body=.


# partial_r f term: generate symbolic expressions for dx^i / dr
def compute_Jacobian_and_inverseJacobian_tofrom_Spherical():
    # Step 2.a: First construct Jacobian matrix:
    Jac_dUSph_dDrfmUD = ixp.zerorank2()
    for i in range(3):
        for j in range(3):
            Jac_dUSph_dDrfmUD[i][j] = sp.diff(rfm.xxSph[i], rfm.xx[j])
    Jac_dUrfm_dDSphUD, dummyDET = ixp.generic_matrix_inverter3x3(Jac_dUSph_dDrfmUD)
    return Jac_dUSph_dDrfmUD, Jac_dUrfm_dDSphUD


# r_and_partial_xi_partial_r_derivs(): Compute r(x0,x1,x2) and dx^i / dr
def setup_Cfunction_r_and_partial_xi_partial_r_derivs():
    desc = "Compute r(xx0,xx1,xx2) and partial_r x^i."
    c_type = "static inline void"
    name = "r_and_partial_xi_partial_r_derivs"
    params = """const paramstruct *restrict params,const REAL xx0,const REAL xx1,const REAL xx2,
                                                     REAL *r,
                                                     REAL *partial_x0_partial_r,REAL *partial_x1_partial_r,REAL *partial_x2_partial_r"""
    Jac_dUSph_dDrfmUD, Jac_dUrfm_dDSphUD = compute_Jacobian_and_inverseJacobian_tofrom_Spherical()
    body = outputC([rfm.xxSph[0],
                    Jac_dUrfm_dDSphUD[0][0], # sp.simplify(expr) is too slow here for SinhCylindrical
                    Jac_dUrfm_dDSphUD[1][0], # sp.simplify(expr) is too slow here for SinhCylindrical
                    Jac_dUrfm_dDSphUD[2][0]],# sp.simplify(expr) is too slow here for SinhCylindrical
                   ["*r", "*partial_x0_partial_r", "*partial_x1_partial_r", "*partial_x2_partial_r"], filename="returnstring",
                   params="preindent=1,outCverbose=False,includebraces=False")
    rel_path_to_Cparams = os.path.join(".")

    _prototype, func = Cfunction(
        includes=[],   desc=desc, c_type=c_type, name=name, params=params,    body=body,
        rel_path_to_Cparams=rel_path_to_Cparams)
    return func


# partial_r f term: generate finite-difference coefficients
#   for \partial_i f with arbitrary upwinding:
def get_arb_offset_FD_coeffs_indices(FDORDER, offset, deriv):
    # deriv = 1 <-- 1st derivative
    Minv = fin.setup_FD_matrix__return_inverse(FDORDER+1, offset)
    indices = []
    coeffs = []
    for i in range(FDORDER+1):
        indices.append(i-int(FDORDER/2) + offset)
        coeffs.append(Minv[i, deriv])
    return coeffs, indices


# partial_r f term: FD1_arbitrary_upwind(): C function to evaluate
#   partial_i f with arbitrary upwinding
def setup_Cfunction_FD1_arbitrary_upwind(dirn, radiation_BC_FD_order=-1):
    default_FDORDER = par.parval_from_str("finite_difference::FD_CENTDERIVS_ORDER")
    if radiation_BC_FD_order == -1:
        radiation_BC_FD_order = default_FDORDER

    par.set_parval_from_str("finite_difference::FD_CENTDERIVS_ORDER", radiation_BC_FD_order)

    includes = []
    desc = "Compute 1st derivative finite-difference derivative with arbitrary upwind"
    c_type = "static inline REAL"
    name = "FD1_arbitrary_upwind_x"+str(dirn)+"_dirn"
    params = """const paramstruct *restrict params, const REAL *restrict gf,
                                                const int i0,const int i1,const int i2, const int offset"""
    body = r"""switch(offset) {
"""
    tmp_list = []
    for offset in range(0, int(radiation_BC_FD_order / 2) + 1):
        tmp_list.append(offset)
        if offset > 0:
            tmp_list.append(-offset)

    for offset in tmp_list:
        body += "case " + str(offset) + ":\n"
        body += "  return ("
        coeffs, indices = get_arb_offset_FD_coeffs_indices(radiation_BC_FD_order, offset, 1)
        for i, coeff in enumerate(coeffs):
            if coeff == 0:
                continue  # skip this iteration if coeff=0
            offset = str(indices[i])
            if i > 0:
                body += "          "
            if offset == "0":
                body += "+"+str(sp.ccode(coeff))+"*gf[IDX3S(i0,i1,i2)]\n"
            else:
                if dirn == 0:
                    body += "+"+str(sp.ccode(coeff))+"*gf[IDX3S(i0+"+offset+",i1,i2)]\n"
                elif dirn == 1:
                    body += "+"+str(sp.ccode(coeff))+"*gf[IDX3S(i0,i1+"+offset+",i2)]\n"
                elif dirn == 2:
                    body += "+"+str(sp.ccode(coeff))+"*gf[IDX3S(i0,i1,i2+"+offset+")]\n"
        body = body[:-1].replace("+-", "-") + ") * invdx"+str(dirn)+";\n"
    body += """}
return 0.0 / 0.0;  // poison output if offset computed incorrectly
"""
    rel_path_to_Cparams = os.path.join(".")

    _prototype, func = Cfunction(includes=includes,
        desc=desc,
        c_type=c_type, name=name, params=params,
        body=indent_Ccode(body, "  "),
        rel_path_to_Cparams=rel_path_to_Cparams)

    par.set_parval_from_str("finite_difference::FD_CENTDERIVS_ORDER", default_FDORDER)

    return func


# partial_r f term: Numerically evaluate partial_r f,
#   calling functions defined above.
def setup_Cfunction_compute_partial_r_f(radiation_BC_FD_order=-1):
    desc = "Compute \partial_r f"
    c_type = "static inline REAL"
    name = "compute_partial_r_f"
    params = """const paramstruct *restrict params, REAL *restrict xx[3], const REAL *restrict gfs,
                                       const int which_gf, const int dest_i0,const int dest_i1,const int dest_i2,
                                       const int FACEi0,const int FACEi1,const int FACEi2,
                                       const REAL partial_x0_partial_r, const REAL partial_x1_partial_r, const REAL partial_x2_partial_r"""
    Jac_dUSph_dDrfmUD, Jac_dUrfm_dDSphUD = compute_Jacobian_and_inverseJacobian_tofrom_Spherical()

    default_FDORDER = par.parval_from_str("finite_difference::FD_CENTDERIVS_ORDER")
    if radiation_BC_FD_order == -1:
        radiation_BC_FD_order = default_FDORDER

    body = r"""  ///////////////////////////////////////////////////////////

  // FD1_stencil_radius = radiation_BC_FD_order/2 = """ + str(int(radiation_BC_FD_order/2)) + r"""
  const int FD1_stencil_radius = """ + str(int(radiation_BC_FD_order/2)) + r""";

  const int ntot = Nxx_plus_2NGHOSTS0*Nxx_plus_2NGHOSTS1*Nxx_plus_2NGHOSTS2;

  ///////////////////////////////////////////////////////////
  // Next we'll compute partial_xi f, using a maximally-centered stencil.
  //   The {i0,i1,i2}_offset parameters set the offset of the maximally-centered
  //   stencil, such that an offset=0 implies a centered stencil.

  // CHECK: Nxx_plus_2NGHOSTS0=10; FD1_stencil_radius=2. Then Nxx_plus_2NGHOSTS0-FD1_stencil_radius-1 = 7
  //  if dest_i0 = 9, we get i0_offset=7-9=-2, so the (4th order) deriv
  //  stencil is: -4,-3,-2,-1,0

  // CHECK: if FD1_stencil_radius=2 and dest_i0 = 1, we get i0_offset = FD1_stencil_radius-dest_i0 = 1,
  //  so the (4th order) deriv stencil is: -1,0,1,2,3

  // CHECK: if FD1_stencil_radius=2 and dest_i0 = 0, we get i0_offset = FD1_stencil_radius-1 = 2,
  //  so the (4th order) deriv stencil is: 0,1,2,3,4
"""
    for i in range(3):
        si = str(i)
        if check_zero(Jac_dUrfm_dDSphUD[i][0]):
            body += "  const REAL partial_x"+si+"_f=0.0;\n"
        else:
            body += "  int i"+si+"_offset = FACEi"+si+";  // Shift stencil away from the face we're updating.\n"
            body += "  // Next adjust i"+si+"_offset so that FD stencil never goes out of bounds.\n"
            body += "  if(dest_i"+si+" < FD1_stencil_radius) i"+si+"_offset = FD1_stencil_radius-dest_i"+si+";\n"
            body += "  else if(dest_i"+si+" > (Nxx_plus_2NGHOSTS"+si+"-FD1_stencil_radius-1)) i"+si+"_offset = (Nxx_plus_2NGHOSTS"+si+"-FD1_stencil_radius-1) - dest_i"+si+";\n"
            body += "  const REAL partial_x"+si+"_f=FD1_arbitrary_upwind_x"+si+"_dirn(params,&gfs[which_gf*ntot],dest_i0,dest_i1,dest_i2,i"+si+"_offset);\n"
    body += "  return partial_x0_partial_r*partial_x0_f + partial_x1_partial_r*partial_x1_f + partial_x2_partial_r*partial_x2_f;\n"
    rel_path_to_Cparams = os.path.join(".")

    _prototype, func = Cfunction(
        includes=[],   desc=desc, c_type=c_type, name=name, params=params,    body=body,
        rel_path_to_Cparams=rel_path_to_Cparams)
    return func


# radiation_bcs(): Put it all together, for a single outer boundary point.
def setup_Cfunction_radiation_bcs(radiation_BC_FD_order=-1):
    includes = []
    prefunc = ""
    Jac_dUSph_dDrfmUD, Jac_dUrfm_dDSphUD = compute_Jacobian_and_inverseJacobian_tofrom_Spherical()
    for i in range(3):
        # Do not generate FD1_arbitrary_upwind_xj_dirn() if the symbolic expression for dxj/dr == 0!
        if not check_zero(Jac_dUrfm_dDSphUD[i][0]):
            prefunc += setup_Cfunction_FD1_arbitrary_upwind(dirn=i, radiation_BC_FD_order=radiation_BC_FD_order)
    prefunc += setup_Cfunction_r_and_partial_xi_partial_r_derivs()
    prefunc += setup_Cfunction_compute_partial_r_f(radiation_BC_FD_order=radiation_BC_FD_order)
    desc = r"""*** Apply radiation BCs to all outer boundaries. ***
"""
    c_type = "static inline REAL"
    name = "radiation_bcs"
    params = """const paramstruct *restrict params, const bc_struct *restrict bcstruct,REAL *restrict xx[3],
                                 const REAL *restrict gfs, REAL *restrict gfs_rhss,
                                 const int which_gf, const REAL gf_wavespeed, const REAL gf_f_infinity,
                                 const int dest_i0,const int dest_i1,const int dest_i2,
                                 const short FACEi0,const short FACEi1,const short FACEi2"""
    body = r"""// Nearest "interior" neighbor of this gridpoint, based on current face
const int dest_i0_int=dest_i0+1*FACEi0, dest_i1_int=dest_i1+1*FACEi1, dest_i2_int=dest_i2+1*FACEi2;
REAL r, partial_x0_partial_r,partial_x1_partial_r,partial_x2_partial_r;
REAL r_int, partial_x0_partial_r_int,partial_x1_partial_r_int,partial_x2_partial_r_int;
r_and_partial_xi_partial_r_derivs(params,xx[0][dest_i0],xx[1][dest_i1],xx[2][dest_i2],
                                  &r, &partial_x0_partial_r, &partial_x1_partial_r,  &partial_x2_partial_r);
r_and_partial_xi_partial_r_derivs(params, xx[0][dest_i0_int], xx[1][dest_i1_int], xx[2][dest_i2_int],
                                  &r_int, &partial_x0_partial_r_int, &partial_x1_partial_r_int, &partial_x2_partial_r_int);
const REAL partial_r_f     = compute_partial_r_f(params,xx,gfs, which_gf,dest_i0,    dest_i1,    dest_i2,
                                                 FACEi0,FACEi1,FACEi2,
                                                 partial_x0_partial_r    ,partial_x1_partial_r    ,partial_x2_partial_r);
const REAL partial_r_f_int = compute_partial_r_f(params,xx,gfs, which_gf,dest_i0_int,dest_i1_int,dest_i2_int,
                                                 FACEi0,FACEi1,FACEi2,
                                                 partial_x0_partial_r_int,partial_x1_partial_r_int,partial_x2_partial_r_int);

const int idx3 = IDX3S(dest_i0,dest_i1,dest_i2);
const int idx3_int = IDX3S(dest_i0_int,dest_i1_int,dest_i2_int);

const REAL partial_t_f_int = gfs_rhss[IDX4ptS(which_gf, idx3_int)];

const REAL c = gf_wavespeed;
const REAL f_infinity = gf_f_infinity;
const REAL f     = gfs[IDX4ptS(which_gf, idx3)];
const REAL f_int = gfs[IDX4ptS(which_gf, idx3_int)];
const REAL partial_t_f_int_outgoing_wave = -c * (partial_r_f_int + (f_int - f_infinity) / r_int);

const REAL k = r_int*r_int*r_int * (partial_t_f_int - partial_t_f_int_outgoing_wave);

const REAL rinv = 1.0 / r;
const REAL partial_t_f_outgoing_wave = -c * (partial_r_f + (f - f_infinity) * rinv);

return partial_t_f_outgoing_wave + k * rinv*rinv*rinv;
"""
    rel_path_to_Cparams = os.path.join(".")

    _prototype, func = Cfunction(
        includes=includes, prefunc=prefunc,
        desc=desc,
        c_type=c_type, name=name, params=params,
        body=indent_Ccode(body, "  "),
        rel_path_to_Cparams=rel_path_to_Cparams)
    return func


# apply_bcs_outerradiation_and_inner():
#   Apply radiation BCs at outer boundary points, and
#   inner boundary conditions at inner boundary points.
def add_to_Cfunction_dict_apply_bcs_outerradiation_and_inner(rel_path_to_Cparams=os.path.join("."), radiation_BC_FD_order=2):
    includes = [os.path.join(rel_path_to_Cparams, "NRPy_basic_defines.h"),
                os.path.join(rel_path_to_Cparams, "NRPy_function_prototypes.h")]
    prefunc = setup_Cfunction_radiation_bcs(radiation_BC_FD_order=radiation_BC_FD_order)
    desc = ""
    c_type = "void"
    name = "apply_bcs_outerradiation_and_inner"
    params = """const paramstruct *restrict params, const bc_struct *restrict bcstruct, REAL *restrict xx[3],
                                        const REAL custom_wavespeed[NUM_EVOL_GFS],
                                        const REAL custom_f_infinity[NUM_EVOL_GFS],
                                        REAL *restrict gfs, REAL *restrict rhs_gfs"""
    body = r"""
  // Unpack bc_info from bcstruct
  const bc_info_struct *bc_info = &bcstruct->bc_info;

  ////////////////////////////////////////////////////////
  // STEP 1 of 2: Apply BCs to pure outer boundary points.
  //              By "pure" we mean that these points are
  //              on the outer boundary and not also on
  //              an inner boundary.
  //              Here we fill in the innermost ghost zone
  //              layer first and move outward. At each
  //              layer, we fill in +/- x0 faces first,
  //              then +/- x1 faces, finally +/- x2 faces,
  //              filling in the edges as we go.
  // Spawn N OpenMP threads, either across all cores, or according to e.g., taskset.
#pragma omp parallel
  {
    for(int which_gz=0;which_gz<NGHOSTS;which_gz++) for(int dirn=0;dirn<3;dirn++) {
        // This option results in about 1.6% slower runtime for SW curvilinear at 64x24x24 on 8-core Ryzen 9 4900HS
        //#pragma omp for collapse(2)
        //for(int which_gf=0;which_gf<NUM_EVOL_GFS;which_gf++) for(int idx2d=0;idx2d<bc_info->num_pure_outer_boundary_points[which_gz][dirn];idx2d++) {
        //  {
        // Don't spawn a thread if there are no boundary points to fill; results in a nice little speedup.
        if(bc_info->num_pure_outer_boundary_points[which_gz][dirn] > 0) {
#pragma omp for  // threads have been spawned; here we distribute across them
          for(int idx2d=0;idx2d<bc_info->num_pure_outer_boundary_points[which_gz][dirn];idx2d++) {
            const short i0 = bcstruct->pure_outer_bc_array[dirn + (3*which_gz)][idx2d].i0;
            const short i1 = bcstruct->pure_outer_bc_array[dirn + (3*which_gz)][idx2d].i1;
            const short i2 = bcstruct->pure_outer_bc_array[dirn + (3*which_gz)][idx2d].i2;
            const short FACEX0 = bcstruct->pure_outer_bc_array[dirn + (3*which_gz)][idx2d].FACEX0;
            const short FACEX1 = bcstruct->pure_outer_bc_array[dirn + (3*which_gz)][idx2d].FACEX1;
            const short FACEX2 = bcstruct->pure_outer_bc_array[dirn + (3*which_gz)][idx2d].FACEX2;
            const int idx3 = IDX3S(i0,i1,i2);
            for(int which_gf=0;which_gf<NUM_EVOL_GFS;which_gf++) {
              // *** Apply radiation BCs to all outer boundary points. ***
              rhs_gfs[IDX4ptS(which_gf, idx3)] = radiation_bcs(params, bcstruct, xx, gfs, rhs_gfs, which_gf,
                                                               custom_wavespeed[which_gf], custom_f_infinity[which_gf],
                                                               i0,i1,i2, FACEX0,FACEX1,FACEX2);
            }
          }
        }
      }
  }

  ///////////////////////////////////////////////////////
  // STEP 2 of 2: Apply BCs to inner boundary points.
  //              These map to either the grid interior
  //              ("pure inner") or to pure outer boundary
  //              points ("inner maps to outer"). Those
  //              that map to outer require that outer be
  //              populated first; hence this being
  //              STEP 2 OF 2.
  apply_bcs_inner_only(params, bcstruct, rhs_gfs); // <- apply inner BCs to RHS gfs only
"""
    add_to_Cfunction_dict(
        includes=includes,
        prefunc=prefunc,
        desc=desc,
        c_type=c_type, name=name, params=params,
        body=body,
        rel_path_to_Cparams=rel_path_to_Cparams)


# Only call this after ALL gridfunctions have been registered!
def CurviBoundaryConditions_register_NRPy_basic_defines(verbose=True):
    # Then set up the dictionary entry for CurviBC in NRPy_basic_defines
    Nbd_str  = NRPy_basic_defines_CurviBC_data_structures()
    Nbd_str += NRPy_basic_defines_set_gridfunction_defines_with_parity_types(verbose=verbose)
    outC_NRPy_basic_defines_h_dict["CurviBoundaryConditions"] = Nbd_str


def CurviBoundaryConditions_register_C_functions(rel_path_to_Cparams=os.path.join("./"),
                                                 radiation_BC_FD_order=4):
    add_to_Cfunction_dict_bcstruct_set_up(rel_path_to_Cparams=rel_path_to_Cparams)
    add_to_Cfunction_dict_apply_bcs_outerradiation_and_inner(rel_path_to_Cparams=rel_path_to_Cparams,
                                                             radiation_BC_FD_order=radiation_BC_FD_order)
    add_to_Cfunction_dict_apply_bcs_inner_only(rel_path_to_Cparams=rel_path_to_Cparams)
    add_to_Cfunction_dict_apply_bcs_outerextrap_and_inner(rel_path_to_Cparams=rel_path_to_Cparams)
