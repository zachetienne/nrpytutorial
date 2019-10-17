# The A-to-B driver
from outputC import *
import os
import grid as gri
import indexedexp as ixp
import NRPy_param_funcs as par
import finite_difference as fin

# Set spatial dimension (must be 3 for BSSN)
DIM = 3
par.set_parval_from_str("grid::DIM",DIM)

def GiRaFFE_HO_A2B(outdir): 
    # Register the gridfunction gammadet. This determinant will be calculated separately
    gammadet = gri.register_gridfunctions("AUXEVOL","gammadet")
    # Import the Levi-Civita symbol and build the corresponding tensor.
    # We already have a handy function to define the Levi-Civita symbol in WeylScalars
    import WeylScal4NRPy.WeylScalars_Cartesian as weyl
    LeviCivitaDDD = weyl.define_LeviCivitaSymbol_rank3()
    LeviCivitaUUU = ixp.zerorank3()
    for i in range(DIM):
        for j in range(DIM):
            for k in range(DIM):
                LCijk = LeviCivitaDDD[i][j][k]
                #LeviCivitaDDD[i][j][k] = LCijk * sp.sqrt(gho.gammadet)
                LeviCivitaUUU[i][j][k] = LCijk / sp.sqrt(gammadet)

    # We can use this function to compactly reset to expressions to print at each FD order.
    def set_BU_to_print():
        return [lhrh(lhs=gri.gfaccess("out_gfs","BU0"),rhs=BU[0]),\
                lhrh(lhs=gri.gfaccess("out_gfs","BU1"),rhs=BU[1]),\
                lhrh(lhs=gri.gfaccess("out_gfs","BU2"),rhs=BU[2])]    

    AD = ixp.register_gridfunctions_for_single_rank1("EVOL","AD")
    BU = ixp.register_gridfunctions_for_single_rank1("AUXEVOL","BU")
    AD_dD = ixp.declarerank2("AD_dD","nosym")
    BU = ixp.zerorank1() # BU is already registered as a gridfunction, but we need to zero its values and declare it in this scope.

    for i in range(DIM):
        for j in range(DIM):
            for k in range(DIM):
                BU[i] += LeviCivitaUUU[i][j][k] * AD_dD[k][j]

    par.set_parval_from_str("finite_difference::FD_CENTDERIVS_ORDER", 10)
    fin.FD_outputC(os.path.join(outdir,"B_from_A_order10.h"),set_BU_to_print(),params="outCverbose=False")

    par.set_parval_from_str("finite_difference::FD_CENTDERIVS_ORDER", 8)
    fin.FD_outputC(os.path.join(outdir,"B_from_A_order8.h"),set_BU_to_print(),params="outCverbose=False")

    par.set_parval_from_str("finite_difference::FD_CENTDERIVS_ORDER", 6)
    fin.FD_outputC(os.path.join(outdir,"B_from_A_order6.h"),set_BU_to_print(),params="outCverbose=False")

    par.set_parval_from_str("finite_difference::FD_CENTDERIVS_ORDER", 4)
    fin.FD_outputC(os.path.join(outdir,"B_from_A_order4.h"),set_BU_to_print(),params="outCverbose=False")

    par.set_parval_from_str("finite_difference::FD_CENTDERIVS_ORDER", 2)
    fin.FD_outputC(os.path.join(outdir,"B_from_A_order2.h"),set_BU_to_print(),params="outCverbose=False")

    # For the outermost points, we'll need a separate file for each face. 
    # These will correspond to an upwinded and a downwinded file for each direction.
    AD_ddnD = ixp.declarerank2("AD_ddnD","nosym")
    for i in range(DIM):
        BU[i] = 0
        for j in range(DIM):
            for k in range(DIM):
                if j==0:
                    BU[i] += LeviCivitaUUU[i][j][k] * AD_ddnD[k][j]
                else:
                    BU[i] += LeviCivitaUUU[i][j][k] * AD_dD[k][j]

    fin.FD_outputC(os.path.join(outdir,"B_from_A_order2_dirx0_dnwind.h"),set_BU_to_print(),params="outCverbose=False")

    AD_dupD = ixp.declarerank2("AD_dupD","nosym")
    for i in range(DIM):
        BU[i] = 0
        for j in range(DIM):
            for k in range(DIM):
                if j==0:
                    BU[i] += LeviCivitaUUU[i][j][k] * AD_dupD[k][j]
                else:
                    BU[i] += LeviCivitaUUU[i][j][k] * AD_dD[k][j]

    fin.FD_outputC(os.path.join(outdir,"B_from_A_order2_dirx0_upwind.h"),set_BU_to_print(),params="outCverbose=False")

    for i in range(DIM):
        BU[i] = 0
        for j in range(DIM):
            for k in range(DIM):
                if j==1:
                    BU[i] += LeviCivitaUUU[i][j][k] * AD_ddnD[k][j]
                else:
                    BU[i] += LeviCivitaUUU[i][j][k] * AD_dD[k][j]

    fin.FD_outputC(os.path.join(outdir,"B_from_A_order2_dirx1_dnwind.h"),set_BU_to_print(),params="outCverbose=False")
    for i in range(DIM):
        BU[i] = 0
        for j in range(DIM):
            for k in range(DIM):
                if j==1:
                    BU[i] += LeviCivitaUUU[i][j][k] * AD_dupD[k][j]
                else:
                    BU[i] += LeviCivitaUUU[i][j][k] * AD_dD[k][j]

    fin.FD_outputC(os.path.join(outdir,"B_from_A_order2_dirx1_upwind.h"),set_BU_to_print(),params="outCverbose=False")

    for i in range(DIM):
        BU[i] = 0
        for j in range(DIM):
            for k in range(DIM):
                if j==2:
                    BU[i] += LeviCivitaUUU[i][j][k] * AD_ddnD[k][j]
                else:
                    BU[i] += LeviCivitaUUU[i][j][k] * AD_dD[k][j]

    fin.FD_outputC(os.path.join(outdir,"B_from_A_order2_dirx2_dnwind.h"),set_BU_to_print(),params="outCverbose=False")
    for i in range(DIM):
        BU[i] = 0
        for j in range(DIM):
            for k in range(DIM):
                if j==2:
                    BU[i] += LeviCivitaUUU[i][j][k] * AD_dupD[k][j]
                else:
                    BU[i] += LeviCivitaUUU[i][j][k] * AD_dD[k][j]

    fin.FD_outputC(os.path.join(outdir,"B_from_A_order2_dirx2_upwind.h"),set_BU_to_print(),params="outCverbose=False")
