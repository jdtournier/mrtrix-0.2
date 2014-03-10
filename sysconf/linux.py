from sysconf.common import *

obj_suffix = '.o'
exe_suffix = ''
lib_prefix = 'lib'
lib_suffix = '.so'

cpp = [ 'g++', '-c', '$flags$', '$gtk$', '$path$', '$src$', '-o', '$obj$' ]
cpp_flags = [ '-Wall', '-pedantic', '-fPIC', '-fno-strict-aliasing', '-DGL_GLEXT_PROTOTYPES', '-DUSE_TR1' ]
# set CPP_MARCH environment variable to override -march=native option
# set to empty string to remove option altogether
if 'CPP_MARCH' in os.environ.keys():
  march=os.environ['CPP_MARCH']
  if len(march):
    cpp_flags += [ '-march='+march ]
else:
  cpp_flags += [ '-march=native' ]



ld = [ 'g++', '$flags$', '$path$', '$obj$', '$mrtrix$', '$gsl$', '$gtk$', '$lz$', '-o', '$bin$' ]
ld_flags = [ '-Wl,-rpath,$ORIGIN/../lib' ]
ld_flags_lib_prefix = '-l'

ld_lib = [ 'g++', '-shared', '$flags$', '$obj$', '-o', '$lib$' ]
ld_lib_flags = []

if 'LDFLAGS' in os.environ.keys(): 
  env_flags = shlex.split (os.environ['LDFLAGS'])
  ld_flags.extend (env_flags)
  ld_lib_flags.extend (env_flags)

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

ld_flags_zlib = [ '-lz' ]

default_installto = '/opt/mrtrix'
default_linkto = '/usr'

