# loop1D() is just a special case of loop(), and in fact is called by loop().
#   For documentation on the inputs, see loop()'s documentation below.
def loop1D(idxvar="i0",lower="0",upper="Nx0",incr="1",OpenMPpragma="#pragma omp parallel for",tabprefix=""):
    if not (isinstance(idxvar, str) and isinstance(lower, str) and
            isinstance(upper, str) and isinstance(incr, str) and isinstance(OpenMPpragma, str)):
        print("Error: all inputs to loop1D() must be STRINGS, and to loop() must be LISTS OF STRINGS")
        exit(1)
    OMPheader = ""
    if OpenMPpragma != "":
        OMPheader = OpenMPpragma + "\n"
    incrstr = "++"
    if incr != "1":
        incrstr = "+="+incr
    loopbody = tabprefix+"for(int "+idxvar+"="+lower+"; "+idxvar+"<"+upper+"; "+idxvar+incrstr+")"
    return OMPheader+loopbody+" {\n", tabprefix+"} // END LOOP: "+loopbody.replace(tabprefix,"")+"\n"

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
#    loopguts (optional; string): The loop interior.
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
def loop(idxvar,lower,upper,incr,OpenMPpragma,tabprefix="",loopguts=""):
    # Step 1: Check and/or clean input. 
    # Step 1a: If only strings are passed, then create lists out of them:
    if (isinstance(idxvar,str) and isinstance(lower,str) and
        isinstance(upper, str) and isinstance(incr, str) and isinstance(OpenMPpragma, str)):
        idxvar = [idxvar]
        lower = [lower]
        upper = [upper]
        incr = [incr]
        OpenMPpragma = [OpenMPpragma]
    # Step 1b: At this point all inputs should be lists. If they are not, then exit.
    if not (isinstance(idxvar,list) and isinstance(lower,list) and
            isinstance(upper, list) and isinstance(incr, list) and isinstance(OpenMPpragma, list)):
        print("Error: loop(idxvar,lower,upper,incr,OpenMPpragma) requires all inputs be lists")
        sys.exit(1)
    # Step 1c: At this point all inputs should be lists. If they are not the same length, then exit.
    if len(idxvar) != len(lower) or len(lower) != len(upper) or len(upper) != len(incr) or len(incr) != len(
            OpenMPpragma):
        print("Error: loop(idxvar,lower,upper,incr,OpenMPpragma) requires all inputs be lists OF THE SAME LENGTH")
        sys.exit(1)
    # Step 2: tabprefix will be set according to the loop nesting, so the loop has proper tabination;
    #         one tab for each nesting of the loop.
    # Step 3: header will be the top of the loop
    header = ""
    # Step 4: footerarray 
    footerarray = []
    for i in range(len(idxvar)):
        headerstr,footerstr = loop1D(idxvar[i],lower[i],upper[i],incr[i],OpenMPpragma[i],tabprefix)
        header += headerstr
        footerarray.append(footerstr)
        tabprefix += "    "
    loopgutsout = ""
    if loopguts != "":
        loopgutsarray = loopguts.split("\n")
        for line in loopgutsarray:
            loopgutsout += tabprefix+line+"\n"
    footer = ""
    for i in range(len(idxvar)-1,-1,-1):
        footer += footerarray[i]
    if loopguts == "":
        return header,footer
    return header+loopgutsout+footer

#loopheader,loopfooter = loop(idxvar=["i0","i1"],lower=["0","0"],upper=["Nx0","Nx1"],incr=["1","1"],
#                             OpenMPpragma=["", "#pragma omp parallel for"])
#loopheader,loopfooter = loop(idxvar=["i0","i1"],lower=["0","0"],upper=["Nx0","Nx1"],incr=["1","1"],
#                             OpenMPpragma=["", "#pragma omp parallel for"])
#print(loopheader+loopfooter)
