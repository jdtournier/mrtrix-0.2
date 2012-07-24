from sysconf.common import *

obj_suffix = '.o'
exe_suffix = '.exe'
lib_prefix = ''
lib_suffix = '.dll'

cpp = [ 'g++', '-c', '$flags$', '$gtk$', '$path$', '$src$', '-o', '$obj$' ]
cpp_flags = [ '-Wall', '-march=i686', '-fno-strict-aliasing', '-DGL_GLEXT_PROTOTYPES', '-DUSE_TR1' ]

windres = [ 'windres' ]

ld = [ 'g++', '$flags$', '$path$', '$obj$', '$mrtrix$', '$gsl$', '$gtk$', '$lz$', '-o', '$bin$' ]
ld_flags = []
ld_flags_lib_prefix = '-l'

ld_lib = [ 'g++', '-shared', '$flags$', '$obj$', '-o', '$lib$' ]
ld_lib_flags = []

cpp_flags_debug = cpp_flags + [ '-g' ]
ld_flags_debug = ld_flags + [ '-g' ]
ld_lib_flags_debug = ld_lib_flags + [ '-g' ]

cpp_flags_profile = [ '-pg' ] + cpp_flags_debug
ld_flags_profile = ld_flags_debug + [ '-pg' ]
ld_lib_flags_profile = ld_lib_flags_debug + [ '-pg' ]

cpp_flags += [ '-O2' ]

cpp_flags_release = cpp_flags + [ '-DNDEBUG' ]

cpp_flags_gsl = [ '-IC:/MinGW/msys/1.0/local/include', '-DGSL_DLL' ] 
ld_flags_gsl = [ '-lgsl', '-lgslcblas' ]
ld_flags_gl = [ '-lglu32', '-mwindows' ]
pkgconfig = [ 'pkg-config' ]
pkgconfig_env = None

ld_flags_zlib = [ '-lz' ]

default_installto = ''
default_linkto = ''

