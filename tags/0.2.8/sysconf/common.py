import os

bin_dir = 'bin'
cmd_dir = 'cmd'
lib_dir = 'lib'
misc_dir = 'src'
doc_dir = 'doc'
dev_dir = 'dev'
gl = os.path.join ('src', 'use_gl.h')

cpp_suffix = '.cpp'
h_suffix = '.h'

libname = 'mrtrix'
pkgconfig_gtk = 'gtkmm-2.4'
pkgconfig_glib = 'glibmm-2.4,gthread-2.0'
pkgconfig_gl = 'gtkglext-1.0'

icon = 'icons/icon.o'
icon_dep = 'icons/icon.rc'
