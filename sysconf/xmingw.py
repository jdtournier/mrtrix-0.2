from sysconf.common import *

obj_suffix = '.o'
exe_suffix = '.exe'
lib_prefix = ''
lib_suffix = '.dll'

cpp = [ 'i586-mingw32msvc-g++', '-c', '$flags$', '$gtk$', '$path$', '$src$', '-o', '$obj$' ]
cpp_flags = [ '-std=c++11', '-Wall', '-march=i686', '-fno-strict-aliasing', '-DGL_GLEXT_PROTOTYPES', '-mno-cygwin', '-mms-bitfields' ]
windres = [ 'i586-mingw32msvc-windres' ]

ld_use_shell = True
ld = [ 'i586-mingw32msvc-g++', '--no-undefined', '--enable-runtime-pseudo-reloc', '$flags$', '$obj$', '$path$', '$gsl$', '$gtk$', '$mrtrix$', '-o', '$bin$' ]
ld_flags = [ '-mno-cygwin', '-Wl,-subsystem,windows' ]
ld_flags_lib_prefix = '-l'

ld_lib = [ 'i586-mingw32msvc-g++', '--no-undefined', '--enable-runtime-pseudo-reloc', '-shared', '$flags$', '$obj$', '-o', '$lib$' ]
ld_lib_flags = [ '-mno-cygwin' ]

cpp_flags_debug = cpp_flags + [ '-g' ]
ld_flags_debug = ld_flags + [ '-g' ]
ld_lib_flags_debug = ld_lib_flags + [ '-g' ]

cpp_flags_profile = [ '-pg' ] + cpp_flags_debug
ld_flags_profile = ld_flags_debug + [ '-pg' ]
ld_lib_flags_profile = ld_lib_flags_debug + [ '-pg' ]

cpp_flags += [ '-O2' ]

cpp_flags_release = cpp_flags + [ '-DNDEBUG' ]

cpp_flags_gsl = [ '-I/target/include', '-DGSL_DLL' ]
ld_flags_gsl = [ '-L/target/lib -lgsl', '-lgslcblas' ]
ld_flags_gl = [ '-lopengl32', '-lglu32' ]
pkgconfig = [ 'pkg-config' ]
pkgconfig_env = { 'PKG_CONFIG_PATH': '/target/lib/pkgconfig' }



