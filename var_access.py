import re

from_access = {}

def set_access(ret, varname):
    from_access[ret] = varname

def var_from_access(access):
    access = access.strip() # This really shouldn't be necessary
    v = from_access.get(access, None)
    if v is not None:
        return v
    g = re.match(r'^(\w+)(\[\w+\])*$', access)
    if g:
        return g.group(1)
    g = re.match(r'in_gfs\w*\[IDX4S\((\w+),i0,i1,i2\)\]', access)
    if g:
        return g.group(1)
    g = re.match(r'^const\s+(\w+)\s+(\w+)', access)
    if g:
        return g.group(2)
    g = re.match(r'^\*?([\w.]+?)[DU]*\d*$', access)
    if g:
        return g.group(1)
    g = re.match(r'^(\w+)\[CCTK_GFINDEX3D\(cctkGH,i0,i1,i2\)\]', access)
    if g:
        return g.group(1)
    return "?"
    #raise Exception("Could not identify a variable name from the access string '"+access+"'")
