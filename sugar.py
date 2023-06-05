import sys
import re
import sympy as sp
import indexedexp as ixp
from outputC import lhrh
from colored import colored
from inspect import currentframe

from nrpylatex import parse_latex as parse_latex_
from here import here
from fstr import f

import warnings
warnings.filterwarnings('ignore',r'some variable\(s\) in the namespace were overridden')

coords = None

class excepthook:
    def __init__(self):
        self.excepthook = sys.excepthook
    def __enter__(self):
        pass
    def __exit__(self, ty, val, tb):
        sys.excepthook = self.excepthook

def parse_latex(*args,**kwargs):
    with excepthook():
        args = (re.sub(r'dt2alp',r'dtdtalp',args[0]),)
        return parse_latex_(*args,**kwargs)

def set_coords(*c):
    global coords
    latex = "% coord ["+", ".join(c)+"]"
    parse_latex(latex)
    coords = c

def flatten(lists):
    """
    Convert lists of lists of list to one list.
    >>> flatten([[1,2,[3,4],5,[6,[7,[8]]]]])
    [1, 2, 3, 4, 5, 6, 7, 8]
    """
    new_list = []
    for item in lists:
        if isinstance(item, list):
            new_list += flatten(item)
        else:
            new_list += [item]
    return new_list

###
class Seq:
    """
    Convert a sequence of D's and U's signifying "down" and "up" indexes
    into a latex equivalent which, when combined with specific index values
    becomes a format sequence for latex. Thus, "DD" becomes "_{%s %s}",
    "U" becomes "^{%s}" and "UUD" becomes "^{%s %s}_{%s}".
    """
    def __init__(self,seq):
        self.fmt = ""
        self.prev = ""
        for letter in seq:
            self.add(letter)
        if self.fmt != "":
            self.fmt += "}"
    def add(self, letter):
        if letter in ["l","D"]:
            if self.prev == letter:
                self.fmt += " %s"
            elif self.fmt == "":
                self.fmt += "_{%s"
            else:
                self.fmt += "}_{%s"
        elif letter in ["u","U"]:
            if self.prev == letter:
                self.fmt += " %s"
            elif self.fmt == "":
                self.fmt += "^{%s"
            else:
                self.fmt += "}^{%s"
        else:
            assert False
        self.prev = letter


def _latex_def(basename, fmt, gf, args):
    assert coords is not None, "Please call set_coords()"
    if isinstance(gf, list):
        for i in range(len(gf)):
            _latex_def(basename, fmt, gf[i], args + [coords[i]])
    else:
        indexes = fmt % tuple(args)
        latex = f("% \\mathrm{{{basename}}}{indexes} = \\mathrm{{{gf}}}")
        parse_latex(latex)


def latex_def(basename, sig, gf):
    s = Seq(sig)
    _latex_def(basename, s.fmt, gf, [])

def match_expr(expr):
    """
    Parse a sequence of tensor expressions, e.g. T_{a b} beta^a, into a tuple sequence.
    """
    result = []
    while True:
        g = re.match(r'''(?x)
        \s+| # spaces
        %.*| # comments
        \\mathrm{([A-Za-z_]+)}| # var name, group 1
        \\(hat|tilde|bar){\\?([A-Za-z_]+)}| # var name, group 2 and 3
        \\?([A-Za-z]+)| # var name, group 4
        ([\^_]){\ *((?:[a-z]\ )*[a-z])\ *}| # multi-index, group 5 and 6
        ([\^_])([a-z])| # single index, group 7 and 8
        =
        ''',expr)
        if g:
            s = g.group(0)
            assert s != '', f("{g.group(0)}, {expr}")
            if g.group(1) is not None:
                result += [("symbol",g.group(1))]
            elif g.group(2) is not None:
                assert g.group(3) is not None
                result += [("symbol",g.group(3)+g.group(2))]
            elif g.group(4) is not None:
                result += [("symbol",g.group(4))]
            elif g.group(5) is not None:
                if g.group(5) == "^":
                    typestr = "upindex"
                else:
                    typestr = "downindex"
                answer = re.split(r'\s+',g.group(6))
                result += [(typestr,answer)]
            elif g.group(7) is not None:
                if g.group(7) == "^":
                    typestr = "upindex"
                else:
                    typestr = "downindex"
                answer = [g.group(8)]
                result += [(typestr,answer)]
            elif g.group(0) == "=":
                result += [("equals","=")]
                result += [("end",expr[g.end():])]
                break
            expr = expr[g.end():]
        else:
            result += [("end",expr)]
            break
    return result
###

def get_add_type():
    a, b = sp.symbols("a b")
    return type(a+b)

# Need to get type type of two symbols added together
# comparing with sp.core.add.Add seems to sometimes
# not work.
Add = get_add_type()

properties = {}
variants = {}
definitions = {}

verbose = False

sim_params = {}

def gfparams(**kw):
    global sim_params
    sim_params = {}
    for k in kw:
        if k == "symmetries":
            sim_params["symmetry_option"] = kw[k]
        else:
            sim_params[k] = kw[k]
    if "DIM" not in sim_params:
        sim_params["DIM"] = 3

def n(s):
    if s is None:
        return ""
    assert isinstance(s, str)
    return s

def numstr(n):
    assert n >= 0
    if n == 0:
        return ""
    return str(n)

def _deriv_decl(globs, tens, base, suffix, rank, div, indexes, rstart, sym_start, DIM):
    """
    >>> import NRPy_param_funcs as par
    >>> decl_indexes()
    >>> set_coords("x","y","z")
    >>> gfparams(gf_type="EXTERNAL",symmetries="sym01",centering="VVV",external_module="ADMBase",namefun=name_xyz)
    >>> gfdecl("chi",[])
    >>> deriv_decl(chi, ("_d",[li]), ("_dup",[li]))
    >>> chi_dD
    [chi_dD0, chi_dD1, chi_dD2]
    >>> chi_dupD
    [chi_dupD0, chi_dupD1, chi_dupD2]
    """
    # base += div
    suffix2 = ''
    for k in indexes:
        rank += 1
        if str(k)[0] == 'u':
            suffix2 += 'U'
        else:
            suffix2 += 'D'
    sym_list = list(sym_start)
    for i in range(rstart, rank-1):
        sym_list += [f("sym{i}{i+1}")]
    fullname = base + suffix + div + suffix2
    syms = "_".join(sym_list)
    result = ixp.declare_indexedexp(rank=rank, \
        symbol=fullname, dimension=DIM, \
        symmetry=syms)
    globs[fullname] = result

    nm = base + numstr(len(suffix)) + div + numstr(len(suffix2))
    globs[nm] = sp.IndexedBase(nm)
    definitions[nm + suffix + suffix2] = result

def deriv_decl(tens, *args):
    DIM = sim_params.get("DIM",3)
    globs = currentframe().f_back.f_globals
    rank = 0
    suffix = ''
    if isinstance(tens, sp.Symbol):
        base = str(tens)
    else:
        base = str(tens.base)
        for k in tens.args[1:]:
            rank += 1
            if str(k)[0] == 'u':
                suffix += 'U'
            else:
                suffix += 'D'
    rstart = rank
    props = properties.get(base+suffix, {})
    symmetries = getsyms(props.get("symmetry_option",None))
    sym_start = []
    for sym in symmetries:
        if sym[2] == 1:
            sym_start += [f("sym{sym[0]}{sym[1]}")]
        else:
            sym_start += [f("asym{sym[0]}{sym[1]}")]
    div = None
    indexes = None
    for arg in args:
        if isinstance(arg, str):
            div = arg
            if indexes is None:
                continue
        elif isinstance(arg, list):
            indexes = arg
            if div is None:
                continue
        elif isinstance(arg, tuple):
            assert isinstance(arg[0], str)
            assert isinstance(arg[1], list)
            div, indexes = arg
        _deriv_decl(globs, tens, base, suffix, rank, div, indexes, rstart, sym_start, DIM)

def gfdecl(*args):
    """
    >>> import NRPy_param_funcs as par
    >>> decl_indexes()
    >>> set_coords("x","y","z")
    >>> gfparams(gf_type="EXTERNAL",symmetries="sym01",centering="VVV",external_module="ADMBase",namefun=name_xyz)
    >>> gfdecl("eTtt",[],"eTt",[la],"eT",[la,lb])
    >>> (eTtt,eTtD,eTDD)
    (eTtt, [eTtx, eTty, eTtz], [[eTxx, eTxy, eTxz], [eTxy, eTyy, eTyz], [eTxz, eTyz, eTzz]])
    """
    if len(args)>0 and isinstance(args[-1], dict):
        globs = args[-1]
        args = args[:-1]
    else:
        globs = currentframe().f_back.f_globals
    namelist = []
    for arg in args:
        if isinstance(arg, str):
            namelist += [arg]
            #sim_params["gf_basename"] = name
            #globs[name] = ixp.register_gridfunctions_for_single_rankN(**sim_params)
        elif isinstance(arg, list):
            #gfparams(rank=len(arg))
            rank = len(arg)
            suffix = ""
            for k in arg:
                assert isinstance(k, sp.tensor.indexed.Idx)
                if str(k)[0] == "u":
                    suffix += "U"
                elif str(k)[0] == "l":
                    suffix += "D"
                else:
                    assert False,f("{k} {type(k)}")
            for basename in namelist:
                assert not re.match(r'^.*[DU]$',basename), f("Bad declaration for '{basename}'. Basenames should not end in D or U.")
                fullname = basename + suffix
                if basename not in globs:
                    if rank>0:
                        globs[basename] = sp.IndexedBase(basename,shape=tuple([sim_params["DIM"]]*len((arg))))
                name = basename + suffix
                copy = {}
                for k in sim_params:
                    copy[k] = sim_params[k]
                copy["gf_basename"] = name
                copy["rank"] = rank
                if rank < 2:
                    copy["symmetry_option"] = None
                if copy["gf_type"] != "EXTERNAL":
                    copy["external_module"] = None

                assert fullname not in properties, f("Redefinition of {fullname}")

                base_variants = variants.get(basename,set())
                sym1 = copy.get("symmetry_option", "")

                if verbose:
                    print(colored("Adding Definition for:","cyan"),basename)
                    for k in copy:
                        print("  ",colored(k+":","yellow"),copy[k])
                    if len(base_variants) > 0:
                        print("  ",colored("Previous Definitions:","yellow"),base_variants)
                    print()
                latex = f("% define {fullname} --dim {properties.get('DIM',3)}")
                if sym1 not in ["", None]:
                    latex += f(" --sym {sym1}")
                parse_latex(latex)

                base_variants.add(fullname)
                variants[basename] = base_variants

                properties[fullname] = copy

                if copy["gf_type"] == "DERIV":
                    s = []
                    for i in range(0,copy["rank"]-1):
                        assert i+1 <= 9
                        s += ["sym"+str(i)+str(i+1)]
                    gf = ixp.declare_indexedexp(rank=copy["rank"], \
                        symbol=copy["gf_basename"], dimension=copy["DIM"], \
                        symmetry="_".join(s), \
                        namefun=copy.get("namefun",None))
                else:
                    gf = ixp.register_gridfunctions_for_single_rankN(**copy)

                namefun = copy.get("namefun",None)
                if namefun is not None:
                    # - basename is something like "g" for the metric
                    # - suffix is something like DD
                    # - gf will be an array of values
                    #
                    # These will be assembled into a latex expression
                    # which will be passed off to parse_latex().
                    latex_def(basename, suffix, gf)

                globs[name] = gf
                definitions[name] = gf
            namelist = []
    assert len(namelist)==0, "Missing final index args"

indexdefs = {}

def latex_tensor(inp, globs = None):
    """
    Parse a sequence of tensor expressions, e.g. T_{a b} beta^a, into a tuple sequence.
    """
    if globs is None:
        globs = currentframe().f_back.f_globals
    assert isinstance(inp, str), "input shoud be str"
    assert "," not in inp, f("Commas are not valid in input: '{inp}'")
    args = match_expr(inp)
    assert len(args) > 0, f("Failed to parse {inp}")
    assert args[-1][0] == "end", args[-1]
    symbol = None
    indexes = []
    results = []
    for a in args:
        if a[0] == "symbol":
            if symbol is not None:
                results += [(symbol, indexes)]
                symbol = None
                indexes = []
            symbol = a[1]
        elif a[0] == "upindex":
            for letter in a[1]:
                pair = indexdefs[letter]
                up = pair[1]
                indexes += [ up ]
        elif a[0] == "downindex":
            for letter in a[1]:
                pair = indexdefs[letter]
                down = pair[0]
                indexes += [ down ]
    results += [(symbol, indexes)]
    return results

def symlatex(inp, globs=None):
    if globs is None:
        globs = currentframe().f_back.f_globals
    symbol, indexes = latex_tensor(inp,globs)[0]
    if len(indexes) == 0:
        estr = symbol
    else:
        estr = f("{symbol}{indexes}")
    try:
        return eval(estr,globs)
    except Exception as e:
        here("Bad expression:",estr)
        raise e

def gflatex(inp, globs = None):
    """
    >>> import NRPy_param_funcs as par
    >>> decl_indexes()
    >>> set_coords("x","y","z")
    >>> gfparams(gf_type="EXTERNAL",symmetries="sym01",centering="VVV",external_module="ADMBase",namefun=name_xyz)
    >>> gflatex(r"k_{i j} alpt beta^i")
    >>> (kDD,alpt,betaU)
    ([[kxx, kxy, kxz], [kxy, kyy, kyz], [kxz, kyz, kzz]], alpt, [betax, betay, betaz])
    """
    if globs is None:
        globs = currentframe().f_back.f_globals
    for symbol, indexes in latex_tensor(inp, globs):
        if symbol is not None:
            gfdecl(symbol, indexes, globs)

def decl_indexes():
    g = currentframe().f_back.f_globals
    for c in range(ord('a'),ord('z')+1):
        letter = chr(c)
        dn, up = sp.symbols(f("l{letter} u{letter}"), cls=sp.Idx)
        g[f("l{letter}")] = dn
        g[f("u{letter}")] = up
        indexdefs[letter] = (dn, up)

def ixnam(i):
    return ["x", "y", "z"][i]

def namefun(symbol, index, shape, prefix):
    symbol = prefix
    result = [sp.Symbol(symbol + ''.join(ixnam(n) for n in index + [i]))
            if symbol else sp.sympify(0) for i in range(shape[0])]
    return result

def name_xyz(sym,ind,shape):
    symbase = re.sub("[UD]+$","",sym)
    return namefun(sym,ind,shape,symbase)

def matchindex(s):
    return re.match(r'^([ul])([a-z])$', str(s))

UP_INDEX         = 1
DOWN_INDEX       = 2
CONTRACTED_INDEX = UP_INDEX | DOWN_INDEX

def getindexes(expr):
    """
    getindexes(expr) finds all the symbols in expr that
    represent up or down indexes.
    """
    indexes = {}
    for sym in expr.free_symbols:
        ssym = str(sym)
        g = matchindex(ssym)
        if g:
            updn = g.group(1)
            let = g.group(2)
            if updn == "u":
                mask = UP_INDEX
            else:
                mask = DOWN_INDEX
            indexes[let] = indexes.get(let,0) | mask
    return indexes

def incrindexes(indexes_input,dim,symmetries=None):
    """
    This function is designed to generate all permuations
    of values for a set of indexes. The `indexes_input`
    should be an array of zeros. Thus incrindexes(2,2) should yield
    the sequence [0,0], [0,1], [1,0], [1,1]. Symmetries
    is passed in as a triple of values which represent a pair indexes
    plus a sign. If the indices are switched, a symmetric (or antisymmetric)
    part of the matrix is identified.
    Thus, incrindexes(2,2,[0,1,1]) should yield [0,0], [1,0], and [1,1].
    By symmetry, the index [0,1] is not needed.
    """
    if symmetries is None:
        symmetries = []
    # Make a copy of the input
    indexes = [0]*indexes_input
    yield indexes
    while True:
        for i in range(len(indexes)):
            indexes[i] += 1
            max_val = dim
            for symmetry in symmetries:

                # Check the symmetry
                assert len(symmetry) == 3, f("The symmetry is {symmetry}, it should be (index1,index2,sign)")
                assert symmetry[2] in [1, -1] # Third index is the sign
                assert isinstance(symmetry, (list, tuple))

                # We want the symmetry ordered
                assert symmetry[0] < symmetry[1]

                if symmetry[0] == i:
                    j = symmetry[1]
                    max_val = min(max_val,indexes[j]+1)
            if indexes[i] >= max_val:
                if i+1 == len(indexes):
                    return
                indexes[i] = 0
            else:
                result = list(indexes)
                yield result
                break

def lookup(array, indexes, i=0):
    if i >= len(indexes):
        return array
    else:
        return lookup(array[indexes[i]], indexes, i+1)

def getsyms(syms):
    """
    Parse a symmetry string and return a triple of the
    form index1, index2, and sign.
    """
    li = []
    if syms in [None,""]:
        return li
    for sym in syms.split('_'):
        g = re.match(r'^(a?sym)(\d)(\d)', sym)
        assert g, f("Bad symmetry: '{sym}'")
        if g.group(1) == "sym":
            sign = 1
        else:
            sign = -1
        sym = (int(g.group(2)),int(g.group(3)),sign)
        assert sym[0] < sym[1], f("Symmetry indexes should be in order: {sym}")
        li += [sym]
    return li

def getsuffix(expr):
    suffix = ""
    # Can't use free_symbols because the order varies.
    # We need the order that the user supplied.
    for sym in expr.args[1:]: #expr.free_symbols:
        g = matchindex(sym)
        if g:
            if g.group(1) == "u":
                suffix += "U"
            else:
                suffix += "D"
    return suffix

def make_sum(expr, dim=3):
    expr = sp.expand(expr)
    if isinstance(expr, Add):
        sume = sp.sympify(0)
        for a in expr.args:
            sume += make_sum(a,dim)
        return sume
    elif isinstance(expr, sp.Piecewise):
        nargs = []
        for k in expr.args:
            narg = make_sum(k[0]), k[1]
            nargs += [narg]
        return sp.Piecewise(*nargs)
    elif expr.is_Function:
        if len(expr.args)==1:
            expr = expr.func(make_sum(expr.args[0]))
        elif len(expr.args)==2:
            expr = expr.func(make_sum(expr.args[0]),make_sum(expr.args[1]))
        else:
            raise Exception(f("len(expr.args)={len(expr.args)} not handled."))
        return expr

    indexes = {}
    for fsym in expr.free_symbols:
        fs = str(fsym)
        g = matchindex(fs)
        if g:
            # This is an index
            updown = g.group(1)
            letter = g.group(2)
            if letter not in indexes:
                indexes[letter] = 0
            if updown == "u":
                indexes[letter] |= 1
            else:
                indexes[letter] |= 2
    for index in indexes:
        if indexes[index] == 3:
            new_expr = sp.sympify(0)
            un, dn = sp.symbols(f("u{index} l{index}"), cls=sp.Idx)
            for d in range(dim):
                u1, d1 = sp.symbols(f("u{d} l{d}"), cls=sp.Idx)
                new_expr += expr.subs(un,u1).subs(dn,d1)
            expr = new_expr
    return expr

def eval_expression(expr):
    subs = {}
    for sym in expr.free_symbols:
        if isinstance(sym, sp.tensor.indexed.Indexed):
            nm = str(sym.base)
            indexes = []
            for k in sym.args[1:]:
                ks = str(k)
                assert not re.match(r'([ul])([a-z])$', ks), f("Unevaluated index '{ks}' in expression '{expr}'")
                g = re.match(r'([ul])(\d)$', ks)
                if not g:
                    here(f("`{ks}'"))
                    continue
                if g.group(1) == "u":
                    nm += "U"
                else:
                    nm += "D"
                indexes += [int(g.group(2))]
            assert nm in definitions, f("Missing defenition for '{nm}'.")
            subs[sym] = lookup(definitions[nm],indexes)
    return expr.subs(subs)

def eval_sum(expr,dim=3):
    expr = make_sum(expr, dim)
    return eval_expression(expr)

def geneqns3(eqn, DIM=3, globs=None, loop=False):
    if globs is None:
        globs = currentframe().f_back.f_globals
    m = match_expr(eqn)
    assert len(m) >= 2, f("match_expr failed for '{eqn}' -> {m}")
    sy = symlatex(eqn, globs)
    return geneqns2(lhs=sy, rhs=m[-1][1], loop=loop, globs=globs, DIM=DIM)

def geneqns2(lhs, rhs, DIM=3, globs=None, loop=False):
    if globs is None:
        globs = currentframe().f_back.f_globals
    lhs_str = r"\mathrm{result}"
    last = ""
    suffix = ""

    # We expect a tensor expression of the form Foo[la,lb,ua,ub]
    # We will create a latex expression for this to pass to nrpylatex.
    assert isinstance(lhs, sp.Indexed), f("Type of lhs was '{type(lhs)}', not Indexed")
    for a in lhs.args[1:]:
        sa = str(a)
        # Ensure that we have an index, ua, ub,... or la, lb, ...
        assert len(sa) == 2
        if sa[0] == 'u':
            if last != "u":
                if last != "":
                    lhs_str += "}"
                lhs_str += "^{"
            suffix += "U"
        else:
            assert sa[0] == 'l'
            if last != "l":
                if last != "":
                    lhs_str += "}"
                lhs_str += "_{"
            suffix += "D"
        lhs_str += sa[1] +  " "
        last = sa[0]
    latex = lhs_str + "}=" + rhs

    # The parse expression is of the form "result_{a b ...}^{c d ...} = foo_{a b ...}^{c d ...}"
    # When evaluated, it will assign to the "result" global variable.
    parse_latex(latex) #,verbose=True)

    return geneqns(lhs=lhs, values=globals()["result"+suffix], globs=globs, loop=loop, DIM=DIM)

def geneqns(lhs, rhs=None, values=None, DIM=3, globs=None, loop=False):
    """
    >>> #import NRPy_param_funcs as par
    >>> decl_indexes()
    >>> #set_coords("x","y","z")
    >>> gfparams(gf_type="EXTERNAL",symmetries="sym01",centering="VVV",external_module="ADMBase",namefun=name_xyz)
    >>> gfdecl("alp",[])
    >>> deriv_decl(alp, ("_d",[li]), ("_dup",[li]))
    >>> gfdecl("dalp",[la])
    >>> geneqns(lhs=dalp[li], values=alp_dD)
    [lhrh(lhs=dalpx, rhs=alp_dD0), lhrh(lhs=dalpy, rhs=alp_dD1), lhrh(lhs=dalpz, rhs=alp_dD2)]
    """
    if globs is None:
        globs = currentframe().f_back.f_globals
    if rhs is None and values is None and isinstance(lhs, str):
        return geneqns3(lhs, DIM=DIM, globs=globs, loop=loop)
    if values is None and isinstance(rhs, str):
        return geneqns2(lhs, rhs, DIM=DIM, globs=globs, loop=loop)
    if hasattr(lhs,"base"):
        nm = str(lhs.base)+getsuffix(lhs)
        props = properties.get(nm, {})
        symmetries = getsyms(props.get("symmetry_option",""))
    else:
        assert lhs, f("{lhs}, {type(lhs)}")
    if values is None:
        assert rhs is not None, "Must supply either values or rhs to geneqns"
        lhs_indexes = getindexes(lhs)
        rhs_indexes = getindexes(rhs)
        for k in lhs_indexes:
            assert lhs_indexes[k] != CONTRACTED_INDEX, f("Contracted indexes are not allowed on the left hand side: '{lhs}'")
            assert lhs_indexes.get(k,-1) == rhs_indexes.get(k,-1), f("Free index '{k}' does not match on the lhs and rhs.")
        for k in rhs_indexes:
            if rhs_indexes[k] == CONTRACTED_INDEX:
                continue
            assert lhs_indexes.get(k,-1) == rhs_indexes.get(k,-1), f("Free index '{k}' does not match on the rhs and lhs.")
        if len(lhs_indexes)==0:
            return [lhrh(lhs=lhs, rhs=eval_sum(rhs))]
        else:
            result = []
            indexes = lhs.args[1:] # These will be the indexes
            for index in incrindexes(len(indexes), DIM, symmetries):
                expr_lhs = lhs
                expr_rhs = rhs
                for i in range(len(indexes)):
                    ix = indexes[i]
                    letter = str(ix)[0]
                    assert letter in "ul"
                    rix = sp.symbols(letter+str(index[i]))
                    expr_lhs = expr_lhs.subs(ix, rix)
                    expr_rhs = expr_rhs.subs(ix, rix)
                expr_lhs = eval_sum(expr_lhs)
                expr_rhs = eval_sum(expr_rhs)
                result += [lhrh(lhs=expr_lhs, rhs=expr_rhs)]
                if loop:
                    result += [lhrh(lhs=None, rhs=None)]
            return result
    elif rhs is None:
        result = []
        assert values is not None, "Must supply either values or rhs to geneqns"
        indexes = getindexes(lhs)
        for index in indexes:
            assert indexes[index] != CONTRACTED_INDEX, f("Error, contracted index in lhs: '{index}'")
        for index in incrindexes(len(indexes), DIM, symmetries):
            result += [lhrh(lhs=lookup(definitions[nm],index),rhs=lookup(values, index))]
            if loop:
                loop += [lhrh(lhs=None, rhs=None)]
        return result #generatevalues(lhs,rhs,[0]*len(indexes),symmetries)
    else:
        assert False, "Must supply either values or rhs to geneqns"

if __name__ == "__main__":
    import doctest
    doctest.testmod()
