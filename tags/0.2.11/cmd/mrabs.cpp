/*
    Copyright 2009 Brain Research Institute, Melbourne, Australia

    Written by J-Donald Tournier, 29/10/09.

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

#include "app.h"
#include "image/position.h"

using namespace std; 
using namespace MR; 

SET_VERSION_DEFAULT;

DESCRIPTION = {
  "take absolute value of image intensity",
  NULL
};

ARGUMENTS = {
  Argument ("input", " input image", "the  input image.").type_image_in (),
  Argument ("output", "output image", "the output image.").type_image_out (),
  Argument::End
};

OPTIONS = { Option::End };




EXECUTE {
  Image::Position in (*argument[0].get_image());
  Image::Header header (in.image.header());

  Image::Position out (*argument[1].get_image (header));

  ProgressBar::init (out.voxel_count(), "taking absolute value...");
  do {
    out.value (fabs(in.value()));
    ProgressBar::inc();
    in++;
  } while (out++);
  ProgressBar::done();
}

