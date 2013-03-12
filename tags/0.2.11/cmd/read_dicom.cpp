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

*/

#include <glibmm/fileutils.h>
#include <glibmm/stringutils.h>

#include "app.h"
#include "file/dicom/image.h"

using namespace std; 
using namespace MR; 

SET_VERSION_DEFAULT;

DESCRIPTION = {
  "output DICOM fields in human-readable format.",
  NULL
};


ARGUMENTS = {
  Argument ("file", "DICOM file", "the DICOM file to be scanned.", true, true).type_file (),
  Argument::End
};


OPTIONS = {
  Option ("all", "all DICOM fields", "print all DICOM fields."),
  Option ("csa", "Siemens CSA fields", "print all Siemens CSA fields"),
  Option::End
};


EXECUTE {
  bool print_DICOM_fields = false;
  bool print_CSA_fields = false;

  if (get_options (0).size()) print_DICOM_fields = true;
  if (get_options (1).size()) print_CSA_fields = true;

  for (guint n = 0; n < argument.size();  n++) {

    if (Glib::file_test (argument[n].get_string(), Glib::FILE_TEST_IS_DIR)) {
      Glib::Dir* dir;
      try { dir = new Glib::Dir (argument[n].get_string()); }
      catch (...) { throw Exception (String ("error opening folder \"") + argument[n].get_string() 
          + "\": " + Glib::strerror (errno)); }
      
      String entry;
      while ((entry = dir->read_name()).size()) {
        String filename = Glib::build_filename (argument[n].get_string(), entry);
        try {
          File::Dicom::Image image;
          image.filename = filename;
          image.read (print_DICOM_fields, print_CSA_fields);
          if (!print_DICOM_fields && !print_CSA_fields) 
            cout << image << "\n";
        }
        catch (Exception) {
          throw Exception ("error reading file \"" + filename + "\"");
        }
      }

    }
    else {

      try {
        File::Dicom::Image image;
        image.filename = argument[n].get_string();
        image.read (print_DICOM_fields, print_CSA_fields);
        // if (!print_DICOM_fields && !print_CSA_fields) 
          // cout << image << "\n";
        std::sort (image.frames.begin(), image.frames.end());
        cout << image << "\n";
      }
      catch (Exception) {
        throw Exception (String("error reading file \"") + argument[n].get_string() + "\"");
      }
      
    }

  }

}
  
