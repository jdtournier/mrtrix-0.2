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
#include "file/dicom/quick_scan.h"

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
  File::Dicom::QuickScan reader;

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
        if (reader.read (Glib::build_filename (argument[n].get_string(), entry), print_DICOM_fields, print_CSA_fields))
          error ("error reading file \"" + reader.filename + "\"");
        else cout << reader << "\n";
      }
    }
    else if (reader.read (argument[n].get_string(), print_DICOM_fields, print_CSA_fields))
      error ("error reading file \"" + reader.filename + "\"");

    else if (!print_DICOM_fields) cout << reader << "\n";
  }

}
  
