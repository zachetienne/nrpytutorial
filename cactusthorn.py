# An attempt at a WaveToy using the NRPy+ infrastructure
# TODO: Parity on grid functions
from __future__ import unicode_literals
from __future__ import print_function
from doc_init import doc_init
import os, re
from datetime import date
from sympy import symbols
import finite_difference as fin
from outputC import lhrh
import grid as gri
import NRPy_param_funcs as par
from safewrite import SafeWrite
from sympy.core.symbol import Symbol
import indexedexp as ixp
from fstr import f

def makedirs(dname, exist_ok):
    if exist_ok:
        try:
            # Exist OK flag is not in Python2
            os.makedirs(dname)
        except:
            pass
    else:
        os.makedirs(dname)

def get_user():
    try:
        import pwd
        return pwd.getpwuid(os.getuid()).pw_name
    except:
        return os.environ.get("USER","jovyan")

def flatten(lists):
    new_list = []
    for item in lists:
        if isinstance(item, list):
            new_list += flatten(item)
        else:
            new_list += [item]
    return new_list

msgs = {}

def check_msg(vtype, v, scalar_writes, tile_writes, forgotten):
    if v in forgotten:
        msg = " Was forgotten because of loop break."
    else:
        msg = ""
    assert v in scalar_writes or v in tile_writes, f("{vtype} variable '{v}' was read without being written.{msg}")

class sortedset(set):
    def __init__(self):
        set.__init__(self)
    def __iter__(self):
        li = list(set.__iter__(self))
        li = sorted(li)
        return li.__iter__()

def check_eqns(eqns):
    forgotten = sortedset()
    scalar_reads = sortedset()
    scalar_writes = sortedset()
    tile_writes = sortedset()
    tmp_tile_writes = sortedset()
    for lr in eqns:
        if lr.lhs is None:
            for v in scalar_reads:
                check_msg("SCALAR_TMP", v, scalar_writes, tile_writes, forgotten)
            for v in tmp_tile_writes:
                tile_writes.add(v)
            tmp_tile_writes.clear()
            scalar_reads.clear()
            for v in scalar_writes:
                forgotten.add(v)
            scalar_writes.clear()
            continue
        for vr in lr.rhs.free_symbols:
            v = str(vr)
            gftype = gri.find_gftype(v, fail_on_missing=False)
            if gftype == "SCALAR_TMP":
                scalar_reads.add(v)
            elif gftype == "TILE_TMP":
                check_msg("TILE_TMP", v, scalar_writes, tile_writes, forgotten)
        v = str(lr.lhs)
        gftype = gri.find_gftype(v, fail_on_missing=False)
        if gftype == "SCALAR_TMP":
            scalar_writes.add(v)
        elif gftype == "TILE_TMP":
            tmp_tile_writes.add(v)
    for v in scalar_reads:
        check_msg("SCALAR_TMP", v, scalar_writes, tile_writes, forgotten)

today = date.today().strftime("%B %d, %Y")

makefile_init="""\
# Main make.code.defn file for thorn {thornname}

# Subdirectories containing source files
SUBDIRS =

# Source files in this directory
SRCS ="""


def check_centering(centering):
    if centering is None:
        return
    assert len(centering)==3, f("The centering string '{centering}' should be 3 characters long")
    for c in centering:
        assert c in "cCvV", f("Invalid Centering Character: {c}")

# Divider used to break loops
loop = lhrh(lhs=None, rhs=None)

def typeof(*args):
    t = None
    for a in args:
        if a is None:
            continue
        ta = type(a)
        assert ta in [int, float, bool]
        if t is None:
            t = ta
        if t == int and ta == float:
            t = float
        elif t == float and ta == int:
            pass
        else:
            assert t == ta, f("{t} and {ta}")
    assert t is not None
    return t

class CactusFunc:
    def __init__(self, name, body, schedule_bin, doc, centering):
        self.name = name
        self.body = body
        self.writegfs = sortedset()
        self.readgfs = sortedset()
        self.schedule_bin = schedule_bin
        self.doc = doc
        self.centering = centering

class CactusSrc:
    def __init__(self, name):
        self.name = name
        self.funcs = []
        self.where = None
    def add_func(self, func, where):
        self.funcs += [func]
        self.where = where

class CactusThorn:

    def _add_src(self, src):
        if isinstance(src, str) or isinstance(src, unicode):  # unicode is a valid type only in Python2
            csrc = CactusSrc(src)
        elif isinstance(src, CactusSrc):
            csrc = src
        else:
            assert False, "Bad type for src"
        self.src_files[csrc.name] = csrc
        self.last_src_file = csrc

    def add_func(self, name, body, schedule_bin, doc, where='interior', centering=None, sync=None):
        body = flatten(body)
        check_eqns(body)
        self.sync[name] = sync
        check_centering(centering)
        if gri.ET_driver == "Carpet":
            schedule_bin = re.sub(r'\bRHS\b','MoL_CalcRHS',schedule_bin)
        elif gri.ET_driver == "CarpetX":
            schedule_bin = re.sub(r'\bRHS\b','ODESolvers_RHS',schedule_bin)
        else:
            assert False
        self._add_src(name + ".cc")
        writegfs = sortedset()
        readgfs = sortedset()
        tmps = sortedset()
        if isinstance(body, list):
            new_body = []
            for item in body:
                assert type(item) == lhrh, "Equations should be stored in outputC.lhrh objects."

                if hasattr(item.rhs, "free_symbols"):
                    for sym in item.rhs.free_symbols:
                        rdsym = str(sym)
                        gftype = gri.find_gftype(rdsym,fail_on_missing=False)
                        if gftype == "TILE_TMP":
                            tmps.add(rdsym)
                            continue
                        elif gftype == "SCALAR_TMP":
                            continue
                        # Check if the symbol name is a derivative
                        g = re.match(r'(.*)_dD+\d+$', rdsym)
                        if g:
                            rdsym = g.group(1)
                        if rdsym not in writegfs:
                            readgfs.add(rdsym)

                writem = str(item.lhs)
                gftype = gri.find_gftype(writem,fail_on_missing=False)
                if gftype == "TILE_TMP":
                    tmps.add(writem)
                elif gftype == "SCALAR_TMP":
                    pass
                elif item.lhs is None:
                    pass
                else:
                    writegfs.add(writem)
                    if writem == "regrid_error":
                        cent_gf = "CCC"
                    else:
                        cent_gf = gri.find_centering(writem)
                    assert cent_gf is not None, f("Centering for grid function '{writem}' is unknown")
                    if centering is None:
                        centering = cent_gf
                    assert cent_gf == centering, f("Centering of '{writem}' is '{cent_gf},' but it must match the loop centering of function '{name}': '{centering}'")
                if isinstance(item.lhs, Symbol):
                    new_lhs=gri.gfaccess(varname=str(item.lhs),context="USE")
                    new_body += [lhrh(lhs=new_lhs, rhs=item.rhs)]
                else:
                    new_body += [item]

            assert centering is not None, "The centering for loop '{name}' is none"
            body = new_body
            new_body = []
            body_str = ""
            for i in range(len(body)+1):
                if i == len(body) or (body[i].lhs is None and body[i].rhs is None):
                    if len(new_body) > 0:
                        body_str += self.do_body(new_body,where,centering) + "\n"
                    new_body = []
                else:
                    new_body += [body[i]]
        elif isinstance(body, str):
            # Pass the body through literally
            pass
        else:
            assert False, "Body must be list or str"

        centerings_used = sortedset()
        for gf in readgfs:
            gtype = gri.find_gftype(gf,fail_on_missing=False)
            c = gri.find_centering(gf)
            if gtype in ["AUX","EVOL","AUXEVOL"] and centering is not None:
                centerings_used.add(c)
        for gf in writegfs:
            gtype = gri.find_gftype(gf,fail_on_missing=False)
            c = gri.find_centering(gf)
            if gtype in ["AUX","EVOL","AUXEVOL"] and centering is not None:
                centerings_used.add(c)
        centerings_needed = sortedset()
        centerings_needed.add(centering)
        layout_decls = ""
        tmp_centerings = {}
        for gf in tmps:
            c = gri.find_centering(gf)
            if c not in tmp_centerings:
                tmp_centerings[c] = sortedset()
            tmp_centerings[c].add(gf)
        if gri.ET_driver == "CarpetX":
            if where == 'interior':
                layout_decls += f("  // Allocate temporary grid functions without ghost zones\n")
                layout_decls += f("  const GF3D5layout CCTK_ATTRIBUTE_UNUSED VVV_tmp_layout(cctkGH, {{0,0,0}}, {{0,0,0}});\n")
                layout_decls += f("  const GF3D5layout CCTK_ATTRIBUTE_UNUSED VVC_tmp_layout(cctkGH, {{0,0,1}}, {{0,0,0}});\n")
                layout_decls += f("  const GF3D5layout CCTK_ATTRIBUTE_UNUSED VCV_tmp_layout(cctkGH, {{0,1,0}}, {{0,0,0}});\n")
                layout_decls += f("  const GF3D5layout CCTK_ATTRIBUTE_UNUSED VCC_tmp_layout(cctkGH, {{0,1,1}}, {{0,0,0}});\n")
                layout_decls += f("  const GF3D5layout CCTK_ATTRIBUTE_UNUSED CVV_tmp_layout(cctkGH, {{1,0,0}}, {{0,0,0}});\n")
                layout_decls += f("  const GF3D5layout CCTK_ATTRIBUTE_UNUSED CVC_tmp_layout(cctkGH, {{1,0,1}}, {{0,0,0}});\n")
                layout_decls += f("  const GF3D5layout CCTK_ATTRIBUTE_UNUSED CCV_tmp_layout(cctkGH, {{1,1,0}}, {{0,0,0}});\n")
                layout_decls += f("  const GF3D5layout CCTK_ATTRIBUTE_UNUSED CCC_tmp_layout(cctkGH, {{1,1,1}}, {{0,0,0}});\n")
            else:
                layout_decls += f("  const GF3D5layout CCTK_ATTRIBUTE_UNUSED VVV_tmp_layout(cctkGH, {{0,0,0}});\n")
                layout_decls += f("  const GF3D5layout CCTK_ATTRIBUTE_UNUSED VVC_tmp_layout(cctkGH, {{0,0,1}});\n")
                layout_decls += f("  const GF3D5layout CCTK_ATTRIBUTE_UNUSED VCV_tmp_layout(cctkGH, {{0,1,0}});\n")
                layout_decls += f("  const GF3D5layout CCTK_ATTRIBUTE_UNUSED VCC_tmp_layout(cctkGH, {{0,1,1}});\n")
                layout_decls += f("  const GF3D5layout CCTK_ATTRIBUTE_UNUSED CVV_tmp_layout(cctkGH, {{1,0,0}});\n")
                layout_decls += f("  const GF3D5layout CCTK_ATTRIBUTE_UNUSED CVC_tmp_layout(cctkGH, {{1,0,1}});\n")
                layout_decls += f("  const GF3D5layout CCTK_ATTRIBUTE_UNUSED CCV_tmp_layout(cctkGH, {{1,1,0}});\n")
                layout_decls += f("  const GF3D5layout CCTK_ATTRIBUTE_UNUSED CCC_tmp_layout(cctkGH, {{1,1,1}});\n")

        tmp_decls = ""
        for c in tmp_centerings:
            ctmps = tmp_centerings[c]
            tmp_decls=f("  const GF3D5vector<CCTK_REAL> tiles_{c}({c}_tmp_layout, {len(ctmps)});\n")
            tmp_list = sorted(list(ctmps))
            for ix in range(len(tmp_list)):
                tname = tmp_list[ix]
                c = gri.find_centering(tname)
                tmp_decls += f("  const GF3D5 {tname}(tiles_{c}({ix}));\n")

        body_str = layout_decls + tmp_decls + body_str

        func = CactusFunc(name, body_str, schedule_bin, doc, centering)
        self.last_src_file.add_func(func, where)
        func.readgfs = readgfs
        func.writegfs = writegfs

    def do_body(self, body, where, centering):
        # idxs is a list of integer offsets for stencils used by the FD routine
        idxs = sortedset()
        kernel = fin.FD_outputC("returnstring",body,idxs=idxs)
        kernel = f("""
        // Begin NRPy+ Kernel
        {kernel}
        // End NRPy+ Kernel\n""")
        # Compute the max offset from 0 for the stencils
        maxx = 0
        maxy = 0
        maxz = 0
        for idx4s in idxs:
            idx4 = [int(ii) for ii in idx4s.split(",")]
            maxx = max(abs(idx4[0]),maxx)
            maxy = max(abs(idx4[1]),maxy)
            maxz = max(abs(idx4[2]),maxz)

        # Generate code to check the stencils
        # have enough ghost zones
        checkbounds = f("""
  if(cctk_nghostzones[0] < {maxx}) CCTK_ERROR("cctk_nghostzones[0] must be at least {maxx}");
  if(cctk_nghostzones[1] < {maxy}) CCTK_ERROR("cctk_nghostzones[1] must be at least {maxy}");
  if(cctk_nghostzones[2] < {maxz}) CCTK_ERROR("cctk_nghostzones[2] must be at least {maxz}");
  if(2*cctk_nghostzones[0] >= cctk_lsh[0]) CCTK_ERROR("cctk_nghostzones[0] is too large");
  if(2*cctk_nghostzones[1] >= cctk_lsh[1]) CCTK_ERROR("cctk_nghostzones[1] is too large");
  if(2*cctk_nghostzones[2] >= cctk_lsh[2]) CCTK_ERROR("cctk_nghostzones[2] is too large");""")
        decl = ""
        if gri.ET_driver == "Carpet":
            if where == "interior":
                body = f("""
  {decl}
  {checkbounds}
  for(int i2=cctk_nghostzones[2];i2<cctk_lsh[2]-cctk_nghostzones[2];i2++) {{
  for(int i1=cctk_nghostzones[1];i1<cctk_lsh[1]-cctk_nghostzones[1];i1++) {{
  for(int i0=cctk_nghostzones[0];i0<cctk_lsh[0]-cctk_nghostzones[0];i0++) {{
  {kernel}
  }} }} }}
                    """).strip()
            elif where == "everywhere":
                body = f("""
  {decl}
  for(int i2=0;i2<cctk_lsh[2];i2++) {{
  for(int i1=0;i1<cctk_lsh[1];i1++) {{
  for(int i0=0;i0<cctk_lsh[0];i0++) {{
  {kernel}
  }} }} }}
                    """).strip()
            elif where == "boundary":
                body = f("""
  {decl}
  for(int i2=0;i2<cctk_nghostzones[2];i2++) {{
  for(int i1=0;i1<cctk_lsh[1];i1++) {{
  for(int i0=0;i0<cctk_lsh[0];i0++) {{
  {kernel}
  }} }} }}
  for(int i2=0;i2<cctk_lsh[2];i2++) {{
  for(int i1=0;i1<cctk_nghostzones[1];i1++) {{
  for(int i0=0;i0<cctk_lsh[0];i0++) {{
  {kernel}
  }} }} }}
  for(int i2=0;i2<cctk_lsh[2];i2++) {{
  for(int i1=0;i1<cctk_lsh[1];i1++) {{
  for(int i0=0;i0<cctk_nghostzones[0];i0++) {{
  {kernel}
  }} }} }}

  for(int i2=cctk_lsh[2]-cctk_nghostzones[2];i2<cctk_lsh[2];i2++) {{
  for(int i1=0;i1<cctk_lsh[1];i1++) {{
  for(int i0=0;i0<cctk_lsh[0];i0++) {{
  {kernel}
  }} }} }}
  for(int i2=0;i2<cctk_lsh[2];i2++) {{
  for(int i1=cctk_lsh[1]-cctk_nghostzones[1];i1<cctk_lsh[1];i1++) {{
  for(int i0=0;i0<cctk_lsh[0];i0++) {{
  {kernel}
  }} }} }}
  for(int i2=0;i2<cctk_lsh[2];i2++) {{
  for(int i1=0;i1<cctk_lsh[1];i1++) {{
  for(int i0=cctk_lsh[0]-cctk_nghostzones[0];i0<cctk_lsh[0];i0++) {{
  {kernel}
  }} }} }}
                    """).strip()
            else:
                assert False, f("where={where} is not supported")
        elif gri.ET_driver == "CarpetX":
            if where == 'everywhere':
                wtag = 'all'
            elif where == 'interior':
                wtag = 'int'
            elif where == 'boundary':
                wtag = 'bnd'
            else:
                assert False, "where should be in ['interior', 'everywhere', 'boundary']"
            body = f("""
  {decl}
  {checkbounds}
  grid.loop_{wtag}_device<{centering}_centered[0], {centering}_centered[1], {centering}_centered[2], CCTK_VECSIZE>(
  grid.nghostzones, [=] CCTK_DEVICE (
  const PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {{
  const GF3D5index CCTK_ATTRIBUTE_UNUSED {centering}_index({centering}_layout, p.I);
  const GF3D5index CCTK_ATTRIBUTE_UNUSED {centering}_tmp_index({centering}_tmp_layout, p.I);
  const CCTK_BOOLVEC mask CCTK_ATTRIBUTE_UNUSED = mask_for_loop_tail<CCTK_BOOLVEC>(p.i, p.imax);
  {kernel}
  }});
                """).strip()
        return body

    def declare_param(self, name, default, doc, vmin=None, vmax=None, options=None):
        self.params += [(name, default, doc, vmin, vmax, options)]
        ty = type(default)
        if ty == bool:
            c_type = "bool"
        elif ty == int:
            c_type = "int"
        elif ty == float:
            c_type = "REAL"
        else:
            raise Exception("Unkown parameter type: "+str(ty))
        par.Cparameters(c_type,self.thornname,[name],default)
        return symbols(name)

    def use_coords(self):
        return self.coords

    def register_gridfunctions(self, gtype, gf_names, centering=None, external_module=None):
        check_centering(centering)
        if external_module is not None:
            assert gtype == "EXTERNAL"
            thorn = external_module
        else:
            thorn = self.thornname
        if external_module is not None:
            self.external_modules[external_module]=1
        return gri.register_gridfunctions(gtype, gf_names, external_module=external_module,centering=centering)

    def __init__(self, arrangement, thornname, author=None, email=None, license='BSD'):
        self.sync = {}
        self.ET_driver = gri.ET_driver
        self.arrangement = arrangement
        self.thornname = thornname
        self.thorn_dir = os.path.join("arrangements", self.arrangement, self.thornname)

        self.src_files = {}
        self.last_src_file = None
        self.params = []
        self.external_modules = {}
        if self.ET_driver == "CarpetX":
            self.external_modules["CarpetX"]=1

        self.param_ccl = os.path.join(self.thorn_dir, "param.ccl")
        self.interface_ccl = os.path.join(self.thorn_dir, "interface.ccl")
        self.schedule_ccl = os.path.join(self.thorn_dir, "schedule.ccl")
        self.configuration_ccl = os.path.join(self.thorn_dir, "configuration.ccl")

        self.register_cc = os.path.join(self.thorn_dir, "src", "register.cc")

        self.src_dir = os.path.join(self.thorn_dir, "src")
        self.makefile = os.path.join(self.src_dir, "make.code.defn")

        self.test_dir = os.path.join(self.thorn_dir, "test")
        self.test_ccl = os.path.join(self.test_dir, "test.ccl")

        self.doc_dir = os.path.join(self.thorn_dir, "doc")
        self.doc_tex = os.path.join(self.doc_dir, "documentation.tex")
        self.author = author
        self.email = email
        self.license = license

        par.set_parval_from_str("grid::DIM",3)
        #par.set_parval_from_str("reference_metric::CoordSystem","Cartesian")
        par.set_parval_from_str("grid::GridFuncMemAccess","ETK")

        if self.author is None or self.email is None:
            # Does gitconfig exist? Look up author and email there...
            home = os.environ["HOME"]
            git_config = os.path.join(home, ".gitconfig")
            if os.path.exists(git_config):
                import configparser
                config = configparser.ConfigParser(strict=False)
                config.read(git_config)
                if 'user' in config:
                    user = config['user']
                    if 'name' in user and self.author is None:
                        self.author = user['name']
                    if 'email' in user and self.email is None:
                        self.email = user['email']
        if self.author is None:
            self.author = get_user()
        if self.email is None:
            self.email = get_user()

    def get_full_name(self,gf_name):
        if gf_name == "regrid_error":
            return "CarpetX::regrid_error"
        gfthorn = gri.find_gfmodule(gf_name,fail_on_missing=False)
        if gfthorn is None:
            gfthorn = self.thornname
        if gfthorn == self.thornname:
            return f("{gfthorn}::{gf_name}GF")
        else:
            return f("{gfthorn}::{gf_name}")

    def get_full_group_name(self,gf_name):
        gf_group = ixp.get_group_name(gf_name)
        gfthorn = gri.find_gfmodule(gf_name,fail_on_missing=False)
        if gfthorn is None:
            gfthorn = self.thornname
        if gfthorn == self.thornname:
            return f("{gfthorn}::{gf_group}GF")
        else:
            return f("{gfthorn}::{gf_name}")

    def generate(self,dirname=None,cactus_config="sim",cactus_thornlist=None,schedule_raw=""):
        assert self.ET_driver == gri.ET_driver
        if self.ET_driver == "CarpetX":
            assert os.path.exists(os.path.join(dirname,"arrangements","CarpetX","CarpetX")), \
                "Generating for CarpetX, but the CarpetX driver is not present."
        else:
            assert os.path.exists(os.path.join(dirname,"arrangements","Carpet","Carpet")), \
                "Generating for Carpet, but the Carpet driver is not present."
        rhs_pairs = sortedset()
        cwd = None
        try:
            if dirname is not None:
                cwd = os.getcwd()
                os.chdir(dirname)
            assert os.path.isdir("arrangements"), "Please run this script from a Cactus root directory"
            makedirs(self.src_dir, exist_ok=True)
            makedirs(self.test_dir, exist_ok=True)
            makedirs(self.doc_dir, exist_ok=True)
            with SafeWrite(self.configuration_ccl) as fd:
                print(f(u"# Configuration definitions for thorn {self.thornname}"),file=fd)
                if gri.ET_driver == "CarpetX":
                    print(f(u"REQUIRES Arith Loop"),file=fd)
            with SafeWrite(self.param_ccl) as fd:
                print(f(u"# Parameter definitions for thorn {self.thornname}"),file=fd)
                for name, default, doc, vmin, vmax, options in self.params:
                    t = typeof(default, vmin, vmax)
                    vmin = vmin if vmin is not None else '*'
                    vmax = vmax if vmax is not None else '*'
                    if t == bool:
                        print(f(u'BOOLEAN {name} "{doc}" {{}} {default}'), file=fd)
                    elif t == int:
                        print(f(u'CCTK_INT {name} "{doc}" {{ {vmin}:{vmax} :: "" }} {default}'), file=fd)
                    elif t == float:
                        print(f(u'CCTK_REAL {name} "{doc}" {{ {vmin}:{vmax} :: "" }} {default}'), file=fd)
            with SafeWrite(self.interface_ccl) as fd:
                print(f(u"# Interface definitions for thorn {self.thornname}"),file=fd)
                print(f(u"IMPLEMENTS: {self.thornname}"),file=fd)
                for v in gri.glb_gridfcs_list:
                    if v.external_module is not None:
                        self.external_modules[v.external_module]=1
                print(f(u"INHERITS: {', '.join(sorted(list(self.external_modules.keys())))}"),file=fd)
                if gri.ET_driver == "CarpetX":
                    print(f(u"USES INCLUDE HEADER: loop_device.hxx"),file=fd)
                    print(f(u"USES INCLUDE HEADER: simd.hxx"),file=fd)
                elif gri.ET_driver == "Carpet":
                    print("CCTK_INT FUNCTION MoLRegisterEvolved(CCTK_INT IN EvolvedIndex, CCTK_INT IN RHSIndex)",file=fd)
                    print("USES FUNCTION MoLRegisterEvolved",file=fd)

                all_group_names = ixp.get_all_group_names()

                for gf_group in all_group_names:
                    #TODO: change lines with "if gf_name in ["x","y","z"]:" to check if gftype=="CORE"
                    if gf_group in ["x","y","z"]:
                        continue
                    if ixp.find_gftype_for_group(gf_group,fail_on_missing=False) in ["TILE_TMP","SCALAR_TMP"]:
                        continue

                    rhs="rhs_" + gf_group
                    if re.match(r'rhs_.*',gf_group):
                        tag = "TAGS='checkpoint=\"no\"'"
                        rhs_pairs.add(gf_group)
                    elif rhs in all_group_names:
                        tag = f("TAGS='rhs=\"{self.thornname}::{rhs}GF\"'")
                    else:
                        # We assume that variables without RHS are not part of the state vector
                        tag = f("TAGS='checkpoint=\"no\"'")

                    module = ixp.find_gfmodule_for_group(gf_group,fail_on_missing=False)
                    if module is None or module == self.thornname:
                        gfs = ixp.get_gfnames_for_group(gf_group)
                        gfs_str = "GF, ".join(list(gfs))
                        if gri.ET_driver == "Carpet":
                            print(f(u'REAL {gf_group}GF TYPE=gf TIMELEVELS=3 {{ {gfs_str}GF }} "{gf_group}"'), file=fd)
                        elif gri.ET_driver == "CarpetX":
                            # CarpetX does not use sub-cycling in time, so it only needs one time level
                            _centering = ixp.find_centering_for_group(gf_group)
                            if _centering is None:
                                _centering = ''
                            else:
                                _centering=f("CENTERING={{ {_centering} }}")
                            print(f(u'REAL {gf_group}GF TYPE=gf TIMELEVELS=1 {tag} {_centering} {{ {gfs_str}GF }} "{gf_group}"'), file=fd)

            if gri.ET_driver == "Carpet":
                with SafeWrite(self.register_cc,do_format=True) as fd:
                    print(f(u"""
#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include "cctk_Parameters.h"

void {self.thornname}_RegisterVars(CCTK_ARGUMENTS)
{{
  DECLARE_CCTK_ARGUMENTS_{self.thornname}_RegisterVars;
  DECLARE_CCTK_PARAMETERS;
  int ierr, var, rhs;""").strip(),file=fd)
                    for rhs_var in sorted(rhs_pairs):
                        var = rhs_var[4:]
                        assert "rhs_" + var == rhs_var, f("rhs_{var} != {rhs_var}")
                        print("   ",f("""
  var   = CCTK_VarIndex("{self.thornname}::{var}GF");
  rhs   = CCTK_VarIndex("{self.thornname}::{rhs_var}GF");
  ierr += MoLRegisterEvolved(var, rhs);
                        """).strip(),file=fd)
                    print("}",file=fd)
            with SafeWrite(self.schedule_ccl) as fd:
                print(f(u"# Schedule definitions for thorn {self.thornname}"),file=fd)
                for gf_name in sorted(gri.find_gfnames()):
                    if gri.ET_driver == "CarpetX" and gf_name in ["x","y","z"]:
                        continue
                    storage_written = sortedset()
                    module = gri.find_gfmodule(gf_name)
                    if module is None or module == self.thornname:
                        gf_group = ixp.rev_index_group.get(gf_name, gf_name)
                        if gf_group in storage_written:
                            continue
                        storage_written.add(gf_group)
                        if gri.find_gftype(gf_name,fail_on_missing=False) in ["TILE_TMP","SCALAR_TMP"]:
                            continue
                        if gri.ET_driver == "Carpet":
                            print(f(u"STORAGE: {gf_group}GF[3]"), file=fd)
                        elif gf_group == "regrid_error":
                            pass
                        else:
                            print(f(u"STORAGE: {gf_group}GF[1]"), file=fd)
                if gri.ET_driver == "Carpet":
                    print(f(u"""
schedule {self.thornname}_RegisterVars in MoL_Register
{{
  LANG: C
  OPTIONS: META
}} "Register variables for MoL"
                    """).strip(), file=fd)
                for src in self.src_files.values():
                    fsrc = os.path.join(self.src_dir, src.name)

                    readgfs = sortedset()
                    writegfs = sortedset()

                    for func in src.funcs:
                        # ignore empty line:
                        # print(file=fd)
                        # We split the `schedule_bin` to allow conditions such as "initial AFTER ADMBase_Initial"
                        if func.schedule_bin.split()[0].lower() in [
                            "basegrid", "initial", "postinitial", "prestep", "evol", "poststep", "analysis"]:
                            atin = "AT"
                        else:
                            atin = "IN"
                        print(f(u"SCHEDULE {func.name} {atin} {func.schedule_bin} {{"),file=fd)
                        print(f(u"    LANG: C"),file=fd)
                        for readgf in sorted(func.readgfs):
                            if gri.ET_driver == "CarpetX" and readgf in ["x","y","z"]:
                                continue
                            # The symbols in readgfs might not actually be grid
                            # functions. Make sure that they are registered as
                            # such before generating read/write decls.
                            if readgf in gri.find_gfnames():
                                full_name = self.get_full_name(readgf)
                                print(f(u"    READS: {full_name}(everywhere)"),file=fd)
                        for writegf in sorted(func.writegfs):
                            if writegf in gri.find_gfnames():
                                full_name = self.get_full_name(writegf)
                                print(f(u"    WRITES: {full_name}({src.where})"),file=fd)
                        sync = self.sync.get(func.name, None)
                        if sync is not None:
                            print(f(u"    SYNC: {sync}"), file=fd)
                        print(f(u'}} "{func.doc}"'), file=fd)
                print(schedule_raw, end='', file=fd)

            if not os.path.exists(self.doc_tex):
                doc_src = doc_init
                doc_src = re.sub(r'{thorn}',self.thornname,doc_src)
                doc_src = re.sub(r'{author}',self.author,doc_src)
                doc_src = re.sub(r'{email}',self.email,doc_src)
                doc_src = re.sub(r'{date}',today,doc_src)
                with open(self.doc_tex, "w") as fd:
                    print(doc_src, file=fd)
            with SafeWrite(self.makefile) as fd:
                # for item in outC_function_master_list
                print(re.sub(r'{thornname}',self.thornname,makefile_init), end='', file=fd)
                for src in self.src_files.values():
                    print(" ",src.name,sep="",end="",file=fd)
                if gri.ET_driver == "Carpet":
                    print(" ","register.cc",file=fd)
                print(u'',file=fd)

            for src in self.src_files.values():
                fsrc = os.path.join(self.src_dir, src.name)
                with SafeWrite(fsrc,do_format=True) as fd:
                    if gri.ET_driver == "Carpet":
                        print("#include <cctk.h>", file=fd)
                        print("#include <cctk_Arguments.h>", file=fd)
                        print("#include <cctk_Parameters.h>", file=fd)
                        print("#include <iostream>", file=fd)
                    elif gri.ET_driver == "CarpetX":
                        print("#include <fixmath.hxx>", file=fd)
                        print("#include <cctk.h>", file=fd)
                        print("#include <cctk_Arguments.h>", file=fd)
                        print("#include <cctk_Parameters.h>", file=fd)
                        print("", file=fd)
                        print("#define CARPETX_GF3D5", file=fd)
                        print("#include <loop_device.hxx>", file=fd)
                        # Activate this line to disable SIMD parallelization
                        # print("#define SIMD_CPU", file=fd)
                        print("#include <simd.hxx>", file=fd)
                        print("", file=fd)
                        print("#include <cmath>", file=fd)
                        print("#include <tuple>", file=fd)
                        print("", file=fd)
                        print("using namespace Arith;", file=fd)
                        print("using namespace Loop;", file=fd)
                        print("using std::cbrt, std::fmax, std::fmin, std::sqrt;", file=fd)
                    else:
                        assert "Bad value for gri.ET_driver={gri.ET_driver}"
                    for func in src.funcs:
                        print('',file=fd)
                        print(f(u"void {func.name}(CCTK_ARGUMENTS) {{"),file=fd)
                        if gri.ET_driver == "Carpet":
                            print(f(u"  DECLARE_CCTK_ARGUMENTS_{func.name};"),file=fd)
                        elif gri.ET_driver == "CarpetX":
                            print(f(u"  DECLARE_CCTK_ARGUMENTSX_{func.name};"),file=fd)
                        print("  DECLARE_CCTK_PARAMETERS;",file=fd)
                        if gri.ET_driver == "CarpetX":
                            print("  using CCTK_BOOLVEC = simdl<CCTK_REAL>;",file=fd)
                            print("  using CCTK_REALVEC = simd<CCTK_REAL>;",file=fd)
                            print("  constexpr std::size_t CCTK_VECSIZE CCTK_ATTRIBUTE_UNUSED = std::tuple_size_v<CCTK_REALVEC>;",file=fd)
                            print("  const Loop::GF3D5layout CCTK_ATTRIBUTE_UNUSED VVV_layout(cctkGH, {0,0,0});",file=fd)
                            print("  const Loop::GF3D5layout CCTK_ATTRIBUTE_UNUSED VVC_layout(cctkGH, {0,0,1});",file=fd)
                            print("  const Loop::GF3D5layout CCTK_ATTRIBUTE_UNUSED VCV_layout(cctkGH, {0,1,0});",file=fd)
                            print("  const Loop::GF3D5layout CCTK_ATTRIBUTE_UNUSED VCC_layout(cctkGH, {0,1,1});",file=fd)
                            print("  const Loop::GF3D5layout CCTK_ATTRIBUTE_UNUSED CVV_layout(cctkGH, {1,0,0});",file=fd)
                            print("  const Loop::GF3D5layout CCTK_ATTRIBUTE_UNUSED CVC_layout(cctkGH, {1,0,1});",file=fd)
                            print("  const Loop::GF3D5layout CCTK_ATTRIBUTE_UNUSED CCV_layout(cctkGH, {1,1,0});",file=fd)
                            print("  const Loop::GF3D5layout CCTK_ATTRIBUTE_UNUSED CCC_layout(cctkGH, {1,1,1});",file=fd)
                        for ii in range(3):
                            print(f(u"  const CCTK_REAL invdx{ii} CCTK_ATTRIBUTE_UNUSED = 1/CCTK_DELTA_SPACE({ii});"),file=fd)
                        print(f("  {func.body}"),file=fd)
                        print(f('}}'),file=fd)

            # Create the thorn_list entry if needed
            entry = f("{self.arrangement}/{self.thornname}")
            thorn_list = os.path.join("configs", cactus_config, "ThornList")
            self.update_thornlist(entry, thorn_list)
            if cactus_thornlist is not None:
                self.update_thornlist(entry, cactus_thornlist)
        finally:
            if cwd is not None:
                os.chdir(cwd)

    def get_xyz(self):
        assert self.ET_driver == gri.ET_driver
        if gri.ET_driver == "Carpet":
            x, y, z = self.register_gridfunctions("EXTERNAL", ["x", "y", "z"], external_module="grid")
        elif gri.ET_driver == "CarpetX":
            x, y, z = self.register_gridfunctions("CORE", ["x", "y", "z"])
        else:
            assert "Bad value for gri.ET_driver={gri.ET_driver}"
        return x, y, z

    def get_regrid_error(self):
        assert self.ET_driver == gri.ET_driver and self.ET_driver == "CarpetX"
        return self.register_gridfunctions("CORE", ["regrid_error"])

    def update_thornlist(self, entry, thorn_list):
        if not thorn_list.startswith("/"):
            thorn_list = os.path.join(os.getcwd(), thorn_list)
        if os.path.exists(thorn_list):
            with open(thorn_list, "r") as fd:
                contents = fd.read()
            if re.search(f('^{entry}\\b'), contents, re.MULTILINE):
                print(f(u"Thorn {entry} is already in {thorn_list}"))
            else:
                if not thorn_list.startswith("/"):
                    thorn_list = os.path.join(os.getcwd(), thorn_list)
                print(f(u"Appending {entry} to {thorn_list}"))
                with open(thorn_list, "a") as fd:
                    print(entry, file=fd)
        else:
            print(f(u"Thornlist {thorn_list} does not exist. Entry {entry} not added."))
