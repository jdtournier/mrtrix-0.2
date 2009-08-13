/*
    Copyright 2009 Brain Research Institute, Melbourne, Australia

    Written by J-Donald Tournier, 13/08/09.

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
  "sample image intensity values along the tracks, producing one intensity value per point along each track.",
  "the track file should generally have been produced by resample_tracks to ensure even sampling.",
  NULL
};

ARGUMENTS = {
  Argument ("tracks", "track file", "the input track file.").type_file (),
  Argument ("image", "sampled image", "the image to be sampled.").type_image_in(),
  Argument ("output", "output file", "the output text file containing the intensity values").type_file(),
  Argument::End
};



OPTIONS = { Option::End };


EXECUTE {
  Tractography::Properties properties;
  Tractography::Reader file;
  file.open (argument[0].get_string(), properties);

  Image::Interp interp (*argument[1].get_image());

  std::ofstream out (argument[2].get_string());

  std::vector<Point> tck;

  ProgressBar::init (0, "sampling tracks...");

  while (file.next (tck)) {
    for (std::vector<Point>::iterator i = tck.begin(); i != tck.end(); ++i) {
      interp.R (*i);
      out << interp.value() << " ";
    }
    out << "\n";
    ProgressBar::inc();
  }
  ProgressBar::done();
  out.close();
}



