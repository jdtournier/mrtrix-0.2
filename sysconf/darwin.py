from sysconf.common import *

obj_suffix = '.o'
exe_suffix = ''
lib_prefix = 'lib'
lib_suffix = '.dylib'

cpp = [ 'g++', '-c', '$flags$', '$gtk$', '$path$', '$src$', '-o', '$obj$' ]
cpp_flags = [ '-Wall', '-mtune=native', '-fPIC', '-fno-strict-aliasing', '-DGL_GLEXT_PROTOTYPES', '-DUSE_TR1' ]

ld = [ 'g++', '$flags$', '$path$', '$obj$', '$mrtrix$', '$gsl$', '$gtk$', '-o', '$bin$' ]
ld_flags = []
ld_flags_lib_prefix = '-l'

ld_lib = [ 'g++', '-shared', '$flags$', '$obj$', '-o', '$lib$' ]
ld_lib_flags = []

# look for FINK dependencies, and act accordingly if found:
fink_dir = [ '/sw64', '/sw' ]
for entry in fink_dir:
  if os.path.isdir (entry):
    cpp_flags += [ '-I' + entry + '/include' ]
    ld_flags += [ '-L' + entry + '/lib' ]
    ld_lib_flags += [ '-L' + entry + '/lib' ]
    break

cpp_flags_debug = cpp_flags + [ '-g' ]
ld_flags_debug = ld_flags + [ '-g' ]
ld_lib_flags_debug = ld_lib_flags + [ '-g' ]

cpp_flags_profile = [ '-pg' ] + cpp_flags_debug
ld_flags_profile = ld_flags_debug + [ '-pg' ]
ld_lib_flags_profile = ld_lib_flags_debug + [ '-pg' ]

cpp_flags += [ '-O2' ]

cpp_flags_release = [ '-DNDEBUG' ]

cpp_flags_gsl = [] 
# uncomment this line for default GSL BLAS implementation, or the next line for the optimised ATLAS libraries (recommended for performance):
ld_flags_gsl = [ '-lgsl', '-lgslcblas' ]
#ld_flags_gsl = [ '-lgsl', '-lcblas', '-latlas', '-llapack', '-lf77blas' ]
ld_flags_gl = [] 

pkgconfig = [ 'pkg-config' ]
pkgconfig_env = None


