# This module contains infrastructure for generating
#   C-code loops of arbitrary dimension

# Primary Author: Zachariah B. Etienne
#         zachetie **at** gmail **dot* com
# Secondary Author: Ken Sible
#         ksible **at** outlook **dot* com

import sys  # Standard Python modules for multiplatform OS-level functions

# loop1D() is just a special case of loop(), and in fact is called by loop().
#   For documentation on the inputs, see loop()'s documentation below.
def loop1D(idxvar="i0", lower="0", upper="Nx0", incr="1", OpenMPpragma="#pragma omp parallel for", tabprefix=""):
    if not (all(isinstance(i, str) for i in (idxvar, lower, upper, incr, OpenMPpragma))):
        raise Exception('all list parameters must be strings.')
    OMPheader = OpenMPpragma + '\n' if OpenMPpragma else ''
    incrstr = '+=' + incr if incr != '1' else '++'
    loopbody = tabprefix + 'for(int ' + idxvar + '=' + lower + '; ' + idxvar + '<' + upper + '; ' + idxvar + incrstr + ')'
    return OMPheader + loopbody + ' {\n', tabprefix + '} // END LOOP: ' + loopbody.replace(tabprefix, '') + '\n'

# loop() creates a C code loop, taking as input:
#    idxvar (string or list of strings): the index variable name or list of names.
#        In the case that idxvar is a list of N strings, we adopt the formulation:
#            idxvar[0]=outermost loop variable
#            idxvar[N-1]=innermost loop variable
#    lower (string or list of strings): lower loop index of idxvar or idxvar[i].
#         See definition of "idxvar" if lower is a list of strings. 
#         idxvar[] must have the same length as idxvar.
#    upper (string or list of strings): Defined similarly to "lower", except
#         this refers to the *upper* loop index of idxvar or idxvar[i].
#    incr (string or list of strings):  Defined similarly to "lower", except
#         this refers to the loop increment of idxvar or idxvar[i]
#         (incr[i] = 1 -> idxvar++)
#    OpenMPpragma (string or list of strings): Defined similarly to "lower", except
#         this refers to the OpenMP pragma corresponding to idxvar or idxvar[i].
#    tabprefix (optional; string): Sets the tab stop for the C code.
#    loopbody (optional; string): The loop interior.
#    tilesize (optional; string): The cache block size for loop tiling.
#
# Example: loop(["i0","i1"],["0","5"],["N","N-5"],["1","5"],["#pragma omp parallel for",""],
#               "double a=-2;\n gf[IDX3(GF,i0,i1)]=-2*a;\n"])
# Output:
# #pragma omp parallel for
# for(int i0=0;i0<N;i0++) {
#     for(int i1=5;i1<N-5;i1+=5) {
#         a=-2;
#         gf[IDX3(GF,i0,i1)]=-2*a;
#     } // END LOOP: for(int i1=5;i1<N-5;i1+=5)
# } // END LOOP: for(int i0=0;i0<N;i0++)
def loop(idxvar, lower, upper, incr, OpenMPpragma, tabprefix="", loopbody="", tilesize=""):
    if (all(isinstance(i, str) for i in (idxvar, lower, upper, incr, OpenMPpragma))):
        idxvar, lower, upper, incr, OpenMPpragma = [idxvar], [lower], [upper], [incr], [OpenMPpragma]
    list_len = len(idxvar)
    if any(len(i) != list_len for i in (lower, upper, incr, OpenMPpragma)):
        raise Exception('all list parameters must be the same length.')
    if tilesize:
        if isinstance(tilesize, str): tilesize = [tilesize]
        if len(tilesize) != list_len:
            raise Exception('tile size must be the same length as all other list parameters.')
    header_list, footer_list = [], []
    for i in range(list_len):
        if len(tilesize) > 0:
            ext_header, ext_footer = loop1D(idxvar[i] + 'B', lower[i], upper[i], tilesize[i], '', tabprefix + i*'    ')
            header, footer = loop1D(idxvar[i], idxvar[i] + 'B', 'MIN(%s,%s+%s)' % (upper[i], idxvar[i] + 'B', \
                tilesize[i]), incr[i], OpenMPpragma[i], tabprefix + (list_len + i)*'    ')
            header_list.insert(i, ext_header)
            footer_list.insert(i, ext_footer)
        else:
            header, footer = loop1D(idxvar[i], lower[i], upper[i], incr[i], OpenMPpragma[i], tabprefix + i*'    ')
        header_list.append(header)
        footer_list.append(footer)
    if loopbody:
        loopbody = [tabprefix + (list_len + len(tilesize))*'    ' + line + '\n' for line in loopbody.split('\n')]
    header = ''.join(header_list)
    footer = ''.join(footer_list[::-1])
    if not loopbody: return header, footer
    return header + ''.join(loopbody) + footer

# Automatic generation of C-code loops around an arbitrarily
#     defined loop body.
def simple_loop(loopopts, loopbody):
    if not loopopts: return loopbody

    if "AllPoints" in loopopts:
        i2i1i0_mins = ["0", "0", "0"]
        i2i1i0_maxs = ["Nxx_plus_2NGHOSTS2", "Nxx_plus_2NGHOSTS1", "Nxx_plus_2NGHOSTS0"]
        if "oldloops" in loopopts:
            i2i1i0_maxs = ["Nxx_plus_2NGHOSTS[2]", "Nxx_plus_2NGHOSTS[1]", "Nxx_plus_2NGHOSTS[0]"]
    elif "InteriorPoints" in loopopts:
        i2i1i0_mins = ["NGHOSTS","NGHOSTS","NGHOSTS"]
        i2i1i0_maxs = ["NGHOSTS+Nxx2","NGHOSTS+Nxx1","NGHOSTS+Nxx0"]
        if "oldloops" in loopopts:
            i2i1i0_maxs = ["NGHOSTS+Nxx[2]", "NGHOSTS+Nxx[1]", "NGHOSTS+Nxx[0]"]
    else: raise Exception('no points over which to loop were specified.')

    Read_1Darrays = ["", "", ""]
    if "Read_xxs" in loopopts:
        if not "EnableSIMD" in loopopts:
            Read_1Darrays = ["const REAL xx0 = xx[0][i0];",
                             "            const REAL xx1 = xx[1][i1];",
                             "        const REAL xx2 = xx[2][i2];", ]
        else: raise Exception('no SIMD support on Read_xxs (currently).')

    if "Enable_rfm_precompute" in loopopts:
        if "Read_xxs" in loopopts:
            raise Exception('Enable_rfm_precompute and Read_xxs cannot both be enabled.')
        if "EnableSIMD" in loopopts:
            Read_1Darrays = ["#include \"rfm_files/rfm_struct__SIMD_inner_read0.h\"",
                             "#include \"rfm_files/rfm_struct__SIMD_outer_read1.h\"",
                             "#include \"rfm_files/rfm_struct__SIMD_outer_read2.h\""]
        else:
            Read_1Darrays = ["#include \"rfm_files/rfm_struct__read0.h\"",
                             "#include \"rfm_files/rfm_struct__read1.h\"",
                             "#include \"rfm_files/rfm_struct__read2.h\""]

    OpenMPpragma   = "" if "DisableOpenMP" in loopopts else "#pragma omp parallel for"
    loopincrements = ["1", "1", "SIMD_width"] if "EnableSIMD" in loopopts else ["1","1","1"]
    tilesize       = ["16", "16", "16"] if "EnableLoopTiling" in loopopts else ""

    return loop(["i2","i1","i0"], i2i1i0_mins, i2i1i0_maxs, loopincrements, [OpenMPpragma, Read_1Darrays[2], Read_1Darrays[1]], \
                tabprefix="    ", loopbody=Read_1Darrays[0] + "\n" + loopbody, tilesize=tilesize)
