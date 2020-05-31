# This module sets up a standard initial data function used for
#   setting up SENR initial data at all gridpoints.

# Author: Zachariah B. Etienne
#         zachetie **at** gmail **dot* com

from outputC import outputC, add_to_Cfunction_dict # NRPy+: Core C code output module

def BSSN_ID_function_string(cf,hDD,lambdaU,aDD,trK,alpha,vetU,betU):
    rhss = [trK,alpha,cf]
    lhss = ["in_gfs[IDX4S(TRKGF,i0,i1,i2)]","in_gfs[IDX4S(ALPHAGF,i0,i1,i2)]","in_gfs[IDX4S(CFGF,i0,i1,i2)]"]
    for i in range(3):
        rhss.append(lambdaU[i])
        lhss.append("in_gfs[IDX4S(LAMBDAU"+str(i)+"GF,i0,i1,i2)]")
        rhss.append(vetU[i])
        lhss.append("in_gfs[IDX4S(VETU"+str(i)+"GF,i0,i1,i2)]")
        rhss.append(betU[i])
        lhss.append("in_gfs[IDX4S(BETU"+str(i)+"GF,i0,i1,i2)]")
        for j in range(i,3):
            rhss.append(hDD[i][j])
            lhss.append("in_gfs[IDX4S(HDD" + str(i) + str(j) + "GF,i0,i1,i2)]")
            rhss.append(aDD[i][j])
            lhss.append("in_gfs[IDX4S(ADD" + str(i) + str(j) + "GF,i0,i1,i2)]")

    # Sort the lhss list alphabetically, and rhss to match:
    lhss,rhss = [list(x) for x in zip(*sorted(zip(lhss, rhss), key=lambda pair: pair[0]))]

    body = outputC(rhss, lhss, filename="returnstring",
                   params="preindent=1,CSE_enable=True,outCverbose=False",  # outCverbose=False to prevent
                                                                            # enormous output files.
                   prestring="", poststring="")

    desc = "Set up the initial data at all points on the numerical grid."
    add_to_Cfunction_dict(
        desc    =desc,
        name    ="initial_data",
        params  ="const paramstruct *restrict params,REAL *restrict xx[3], REAL *restrict in_gfs",
        body    =body,
        loopopts="AllPoints,Read_xxs")
