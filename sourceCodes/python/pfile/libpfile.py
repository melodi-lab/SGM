# This file was automatically generated by SWIG (http://www.swig.org).
# Version 2.0.1
#
# Do not make changes to this file unless you know what you are doing--modify
# the SWIG interface file instead.
# This file is compatible with both classic and new-style classes.

from sys import version_info
if version_info >= (2,6,0):
    def swig_import_helper():
        from os.path import dirname
        import imp
        fp = None
        try:
            fp, pathname, description = imp.find_module('_libpfile', [dirname(__file__)])
        except ImportError:
            import _libpfile
            return _libpfile
        if fp is not None:
            try:
                _mod = imp.load_module('_libpfile', fp, pathname, description)
            finally:
                fp.close()
            return _mod
    _libpfile = swig_import_helper()
    del swig_import_helper
else:
    import _libpfile
del version_info
try:
    _swig_property = property
except NameError:
    pass # Python < 2.2 doesn't have 'property'.
def _swig_setattr_nondynamic(self,class_type,name,value,static=1):
    if (name == "thisown"): return self.this.own(value)
    if (name == "this"):
        if type(value).__name__ == 'SwigPyObject':
            self.__dict__[name] = value
            return
    method = class_type.__swig_setmethods__.get(name,None)
    if method: return method(self,value)
    if (not static) or hasattr(self,name):
        self.__dict__[name] = value
    else:
        raise AttributeError("You cannot add attributes to %s" % self)

def _swig_setattr(self,class_type,name,value):
    return _swig_setattr_nondynamic(self,class_type,name,value,0)

def _swig_getattr(self,class_type,name):
    if (name == "thisown"): return self.this.own()
    method = class_type.__swig_getmethods__.get(name,None)
    if method: return method(self)
    raise AttributeError(name)

def _swig_repr(self):
    try: strthis = "proxy of " + self.this.__repr__()
    except: strthis = ""
    return "<%s.%s; %s >" % (self.__class__.__module__, self.__class__.__name__, strthis,)

try:
    _object = object
    _newclass = 1
except AttributeError:
    class _object : pass
    _newclass = 0


PFILE_HEADER_SIZE = _libpfile.PFILE_HEADER_SIZE
PRIst = _libpfile.PRIst
PRIsst = _libpfile.PRIsst
SEGID_UNKNOWN = _libpfile.SEGID_UNKNOWN
SEGID_BAD = _libpfile.SEGID_BAD
SIZET_BAD = _libpfile.SIZET_BAD
SIZET_MAX = _libpfile.SIZET_MAX
ALL = _libpfile.ALL
OK = _libpfile.OK
BAD = _libpfile.BAD
ENDOF = _libpfile.ENDOF
class PFile_Val(_object):
    __swig_setmethods__ = {}
    __setattr__ = lambda self, name, value: _swig_setattr(self, PFile_Val, name, value)
    __swig_getmethods__ = {}
    __getattr__ = lambda self, name: _swig_getattr(self, PFile_Val, name)
    __repr__ = _swig_repr
    __swig_setmethods__["l"] = _libpfile.PFile_Val_l_set
    __swig_getmethods__["l"] = _libpfile.PFile_Val_l_get
    if _newclass:l = _swig_property(_libpfile.PFile_Val_l_get, _libpfile.PFile_Val_l_set)
    __swig_setmethods__["f"] = _libpfile.PFile_Val_f_set
    __swig_getmethods__["f"] = _libpfile.PFile_Val_f_get
    if _newclass:f = _swig_property(_libpfile.PFile_Val_f_get, _libpfile.PFile_Val_f_set)
    def __init__(self): 
        this = _libpfile.new_PFile_Val()
        try: self.this.append(this)
        except: self.this = this
    __swig_destroy__ = _libpfile.delete_PFile_Val
    __del__ = lambda self : None;
PFile_Val_swigregister = _libpfile.PFile_Val_swigregister
PFile_Val_swigregister(PFile_Val)

PF_LLU = _libpfile.PF_LLU
PF_LLD = _libpfile.PF_LLD
PFILE_LARGEFILES = _libpfile.PFILE_LARGEFILES
class InFtrLabStream_PFile(_object):
    __swig_setmethods__ = {}
    __setattr__ = lambda self, name, value: _swig_setattr(self, InFtrLabStream_PFile, name, value)
    __swig_getmethods__ = {}
    __getattr__ = lambda self, name: _swig_getattr(self, InFtrLabStream_PFile, name)
    __repr__ = _swig_repr
    def __init__(self, *args): 
        this = _libpfile.new_InFtrLabStream_PFile(*args)
        try: self.this.append(this)
        except: self.this = this
    __swig_destroy__ = _libpfile.delete_InFtrLabStream_PFile
    __del__ = lambda self : None;
    def num_labs(self): return _libpfile.InFtrLabStream_PFile_num_labs(self)
    def num_ftrs(self): return _libpfile.InFtrLabStream_PFile_num_ftrs(self)
    def num_segs(self): return _libpfile.InFtrLabStream_PFile_num_segs(self)
    def read_ftrslabs(self, *args): return _libpfile.InFtrLabStream_PFile_read_ftrslabs(self, *args)
    def read_ftrs(self, *args): return _libpfile.InFtrLabStream_PFile_read_ftrs(self, *args)
    def read_labs(self, *args): return _libpfile.InFtrLabStream_PFile_read_labs(self, *args)
    def get_pos(self, *args): return _libpfile.InFtrLabStream_PFile_get_pos(self, *args)
    def nextseg(self): return _libpfile.InFtrLabStream_PFile_nextseg(self)
    def rewind(self): return _libpfile.InFtrLabStream_PFile_rewind(self)
    def num_frames(self, *args): return _libpfile.InFtrLabStream_PFile_num_frames(self, *args)
    def set_pos(self, *args): return _libpfile.InFtrLabStream_PFile_set_pos(self, *args)
InFtrLabStream_PFile_swigregister = _libpfile.InFtrLabStream_PFile_swigregister
InFtrLabStream_PFile_swigregister(InFtrLabStream_PFile)

class OutFtrLabStream_PFile(_object):
    __swig_setmethods__ = {}
    __setattr__ = lambda self, name, value: _swig_setattr(self, OutFtrLabStream_PFile, name, value)
    __swig_getmethods__ = {}
    __getattr__ = lambda self, name: _swig_getattr(self, OutFtrLabStream_PFile, name)
    __repr__ = _swig_repr
    def __init__(self, *args): 
        this = _libpfile.new_OutFtrLabStream_PFile(*args)
        try: self.this.append(this)
        except: self.this = this
    __swig_destroy__ = _libpfile.delete_OutFtrLabStream_PFile
    __del__ = lambda self : None;
    def num_labs(self): return _libpfile.OutFtrLabStream_PFile_num_labs(self)
    def num_ftrs(self): return _libpfile.OutFtrLabStream_PFile_num_ftrs(self)
    def write_ftrslabs(self, *args): return _libpfile.OutFtrLabStream_PFile_write_ftrslabs(self, *args)
    def write_ftrs(self, *args): return _libpfile.OutFtrLabStream_PFile_write_ftrs(self, *args)
    def write_labs(self, *args): return _libpfile.OutFtrLabStream_PFile_write_labs(self, *args)
    def doneseg(self, *args): return _libpfile.OutFtrLabStream_PFile_doneseg(self, *args)
OutFtrLabStream_PFile_swigregister = _libpfile.OutFtrLabStream_PFile_swigregister
OutFtrLabStream_PFile_swigregister(OutFtrLabStream_PFile)


def new_doubleArray(*args):
  return _libpfile.new_doubleArray(*args)
new_doubleArray = _libpfile.new_doubleArray

def delete_doubleArray(*args):
  return _libpfile.delete_doubleArray(*args)
delete_doubleArray = _libpfile.delete_doubleArray

def doubleArray_getitem(*args):
  return _libpfile.doubleArray_getitem(*args)
doubleArray_getitem = _libpfile.doubleArray_getitem

def doubleArray_setitem(*args):
  return _libpfile.doubleArray_setitem(*args)
doubleArray_setitem = _libpfile.doubleArray_setitem

def new_uintArray(*args):
  return _libpfile.new_uintArray(*args)
new_uintArray = _libpfile.new_uintArray

def delete_uintArray(*args):
  return _libpfile.delete_uintArray(*args)
delete_uintArray = _libpfile.delete_uintArray

def uintArray_getitem(*args):
  return _libpfile.uintArray_getitem(*args)
uintArray_getitem = _libpfile.uintArray_getitem

def uintArray_setitem(*args):
  return _libpfile.uintArray_setitem(*args)
uintArray_setitem = _libpfile.uintArray_setitem

def fopen(*args):
  return _libpfile.fopen(*args)
fopen = _libpfile.fopen

def fclose(*args):
  return _libpfile.fclose(*args)
fclose = _libpfile.fclose


