/*
    Copyright 2008 Brain Research Institute, Melbourne, Australia

    Written by J-Donald Tournier, 11/05/09.

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

#include <fstream>
#include <glibmm/stringutils.h>

#include "app.h"
#include "get_set.h"
#include "dwi/tractography/file.h"
#include "dwi/tractography/properties.h"

using namespace MR; 
using namespace MR::DWI; 
using namespace std; 

SET_VERSION_DEFAULT;

DESCRIPTION = {
  "truncate a tracks file by selecting only the first N tracks.",
  NULL
};

ARGUMENTS = {
  Argument ("tracks", "track file", "the input track file.").type_file (),
  Argument ("N", "number of tracks", "the number of tracks to include.").type_integer (1,INT_MAX, 1000),
  Argument ("output", "output file", "the output track file").type_file(),
  Argument::End
};



OPTIONS = { Option::End };




EXECUTE {
  Tractography::Properties properties;
  Tractography::Reader file;
  file.open (argument[0].get_string(), properties);

  int N = argument[1].get_int();

  Tractography::Writer writer;
  writer.create (argument[2].get_string(), properties);

  std::vector<Point> tck;

  ProgressBar::init (N, "truncating tracks...");

  int num = 0;
  while (file.next (tck) && num < N) {
    writer.append (tck);
    writer.total_count++;
    ProgressBar::inc();
    num++;
  }

  file.close();
  writer.close();
  ProgressBar::done();
}

