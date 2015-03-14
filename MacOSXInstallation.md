Thanks to Michael Zeineh for providing these instructions.

# Installation #

Installation on MacOS X relies on the packages made available by the [Fink project](http://www.finkproject.org). This project provides Mac versions of the dependencies required to compile MRtrix. Once the dependencies are installed, MRtrix itself can be compiled and installed.


Here are the complete instructions:

  1. Install the Fink applications, by following [these instructions](http://www.finkproject.org/download/index.php).
  1. Use Fink to install dependencies: <br> <code>fink install python glib gtk+ glibmm2.4 gtkmm2.4 gtkglext1 gsl</code> <br> You may need to repeat the command above multiple times for it complete successfully.<br>
<ol><li>build MRtrix: <br> <code>$ cd mrtrix-0.2.X</code> <br> <code>$ ./build</code>
</li><li>once the compilation stage is complete, copy the contents to a more suitable folder, and change permissions to make it accessible for all users: <br> <code>$ cp -r * /usr/local/mrtrix</code> <br> <code>$ chmod -R a+rX /usr/local/mrtrix</code>
</li><li>change the path to include /usr/local/mrtrix<br>
</li><li>set the library search path appropriately: <br> <code>export DYLD_LIBRARY_PATH="/usr/local/mrtrix/lib"</code></li></ol>

At this point, you should be able to start using MRtrix. You may also want to set up <a href='UnixInstallation#Enable_multi-threading.md'>multi-threading</a>, as per the UnixInstallation instructions, to make full use of a multi-core processor.<br>
<br>
<br>
You may encounter problems related to the freetype library. In this case, you may need to copy the /usr/x11/bin/libfreetype.6.dylib to /sw/lib/freetype219/lib.