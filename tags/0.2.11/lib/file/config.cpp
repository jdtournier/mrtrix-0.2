/*
    Copyright 2008 Brain Research Institute, Melbourne, Australia

    Written by J-Donald Tournier, 27/06/08.

    This file is part of MRtrix.

    MRtrix is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    MRtrix is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with MRtrix.  If not, see <http://www.gnu.org/licenses/>.


    08-07-2008 J-Donald Tournier <d.tournier@brain.org.au>
    * fixed get_int() & get_float()
      They were previously defined as returning bool

*/

#include <glib/gstdio.h>
#include <glibmm/miscutils.h>
#include <glibmm/fileutils.h>

#include "file/config.h"

#define CONFIG_FILE "mrtrix.conf"

#ifdef G_OS_WIN32
#define SYS_CONFIG_FILE "C:" G_DIR_SEPARATOR_S CONFIG_FILE
#define USER_CONFIG_FILE CONFIG_FILE
#else 
#define SYS_CONFIG_FILE "/etc/" CONFIG_FILE
#define USER_CONFIG_FILE "." CONFIG_FILE
#endif



namespace MR {
  namespace File {

    std::map<String, String> Config::config;

    void Config::init ()
    {
      if (Glib::file_test (SYS_CONFIG_FILE, Glib::FILE_TEST_IS_REGULAR)) {
        try { 
          KeyValue kv (SYS_CONFIG_FILE);
          while (kv.next()) { config[kv.key()] = kv.value(); }
        }
        catch (...) { }
      }

      String path = Glib::build_filename (Glib::get_home_dir(), USER_CONFIG_FILE);
      if (Glib::file_test (path, Glib::FILE_TEST_IS_REGULAR)) {
        try {
          KeyValue kv (path);
          while (kv.next()) { config[kv.key()] = kv.value(); }
        }
        catch (...) { }
      }
    }



    bool Config::get_bool (const String& key, bool default_value) 
    {
      String value = get (key); 
      if (value.empty()) return (default_value);
      value = lowercase (value);
      if (value == "true") return (true);
      if (value == "false") return (false);
      error ("malformed boolean entry \"" + value + "\" for key \"" + key + "\" in configuration file - ignored");
      return (default_value);
    }


    int Config::get_int (const String& key, int default_value) 
    {
      String value = get (key); 
      if (value.empty()) return (default_value);
      try { return (to<int> (value)); }
      catch (...) { 
        error ("malformed integer entry \"" + value + "\" for key \"" + key + "\" in configuration file - ignored");
        return (default_value); 
      }
    }


    float Config::get_float (const String& key, float default_value) 
    {
      String value = get (key); 
      if (value.empty()) return (default_value);
      try { return (to<float> (value)); }
      catch (...) { 
        error ("malformed floating-point entry \"" + value + "\" for key \"" + key + "\" in configuration file - ignored");
        return (default_value); 
      }
    }


  }
}

