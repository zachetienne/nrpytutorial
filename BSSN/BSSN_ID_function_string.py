# This module sets up an initial data function meant to
# be called in a pointwise manner at all gridpoints.

# Author: Zachariah B. Etienne
#         zachetie **at** gmail **dot* com

from outputC import *

def BSSN_ID_function_string(cf,hDD,lambdaU,aDD,trK,alpha,vetU,betU):
    returnstring = "void BSSN_ID(REAL xx0,REAL xx1,REAL xx2,\n"
    returnstring += "\tREAL *hDD00,REAL *hDD01,REAL *hDD02,REAL *hDD11,REAL *hDD12,REAL *hDD22,\n"
    returnstring += "\tREAL *aDD00,REAL *aDD01,REAL *aDD02,REAL *aDD11,REAL *aDD12,REAL *aDD22,\n"
    returnstring += "\tREAL *trK,\n"
    returnstring += "\tREAL *lambdaU0,REAL *lambdaU1,REAL *lambdaU2,\n"
    returnstring += "\tREAL *vetU0,REAL *vetU1,REAL *vetU2,\n"
    returnstring += "\tREAL *betU0,REAL *betU1,REAL *betU2,\n"
    returnstring += "\tREAL *alpha,REAL *cf) {\n"
    returnstring += outputC([hDD[0][0], hDD[0][1], hDD[0][2], hDD[1][1], hDD[1][2], hDD[2][2],
                             aDD[0][0], aDD[0][1], aDD[0][2], aDD[1][1], aDD[1][2], aDD[2][2],
                             trK,
                             lambdaU[0], lambdaU[1], lambdaU[2],
                             vetU[0], vetU[1], vetU[2],
                             betU[0], betU[1], betU[2],
                             alpha, cf],
                            ["*hDD00", "*hDD01", "*hDD02", "*hDD11", "*hDD12", "*hDD22",
                             "*aDD00", "*aDD01", "*aDD02", "*aDD11", "*aDD12", "*aDD22",
                             "*trK",
                             "*lambdaU0", "*lambdaU1", "*lambdaU2",
                             "*vetU0", "*vetU1", "*vetU2",
                             "*betU0", "*betU1", "*betU2",
                             "*alpha", "*cf"], filename="returnstring",
                            params="preindent=1,CSE_enable=True,outCverbose=False",  # outCverbose=False to prevent
                            # enormous output files.
                            prestring="", poststring="")
    returnstring += "}\n"
    return returnstring