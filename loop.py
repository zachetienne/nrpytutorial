def loop1D(idxvar="i0",loweridx="0",upperidx="Nx0",increment="1",OpenMPpragma="#pragma omp parallel for",prefix=""):
    str = OpenMPpragma + "\n"
    str += prefix+"for(int "+idxvar+"="+loweridx+"; "+idxvar+"<"+upperidx+"; "+idxvar+"+="+increment+") {"
    return str

def loop(idxvar=["i0","i1"],
         loweridx=["0","0"],
         upperidx=["Nx0","Nx1"],
         increment=["1","1"],
         OpenMPpragma=["","#pragma omp parallel for"]):
    prefix = ""
    str = ""
    for i in range(len(idxvar)):
        str += loop1D(idxvar[i],loweridx[i],upperidx[i],increment[i],OpenMPpragma[i],prefix) + "\n"
        prefix += "    "
    return str

loopstr = loop()
print(loopstr)