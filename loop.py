def loop1D(idxvar="i0",lower="0",upper="Nx0",incr="1",OpenMPpragma="#pragma omp parallel for",prefix=""):
    OMPheader = ""
    if OpenMPpragma != "":
        OMPheader = OpenMPpragma + "\n"
    incrstr = "++"
    if incr != "1":
        incrstr = "+="+incr
    loopbody = prefix+"for(int "+idxvar+"="+lower+"; "+idxvar+"<"+upper+"; "+idxvar+incrstr+")"
    return OMPheader+loopbody+" {\n", prefix+"} // END LOOP: "+loopbody.replace(prefix,"")+"\n"

def loop(idxvar,lower,upper,incr,OpenMPpragma):
    prefix = ""
    header = ""
    footerarray = []
    for i in range(len(idxvar)):
        headerstr,footerstr = loop1D(idxvar[i],lower[i],upper[i],incr[i],OpenMPpragma[i],prefix)
        header += headerstr
        footerarray.append(footerstr)
        prefix += "    "
    footer = ""
    for i in range(len(idxvar)-1,-1,-1):
        footer += footerarray[i]
    return header,footer

loopheader,loopfooter = loop(idxvar=["i0","i1"],lower=["0","0"],upper=["Nx0","Nx1"],incr=["1","1"],
                             OpenMPpragma=["", "#pragma omp parallel for"])
print(loopheader+loopfooter)