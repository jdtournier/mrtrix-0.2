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

#include "app.h"
#include "image/position.h"

using namespace std; 
using namespace MR; 

SET_VERSION_DEFAULT;

DESCRIPTION = {
  "compute images statistics.",
  NULL
};

ARGUMENTS = {
  Argument ("image", "image", "the input image from which statistics will be computed.").type_image_in (),
  Argument::End
};


OPTIONS = { 
  Option ("mask", "brain mask", "only perform computation within the specified binary brain mask image.")
    .append (Argument ("image", "image", "the mask image to use.").type_image_in ()),

  Option::End 
};





EXECUTE {
  Image::Position ima (*argument[0].get_image());

  RefPtr<Image::Position> mask;
  std::vector<OptBase> opt = get_options (0); // mask
  if (opt.size()) {
    mask = new Image::Position (*opt[0][0].get_image ());
    if (mask->dim(0) != ima.dim(0) || mask->dim(1) != ima.dim(1) || mask->dim(2) != ima.dim(2)) 
      throw Exception ("dimensions of mask image do not match that of data image - aborting");
  }

  bool header_shown = false;
  do {
    double mean = 0.0, std = 0.0;
    float min = INFINITY, max = -INFINITY;
    size_t count = 0;

    if (mask) mask->set(2,0); 
    for (ima.set(2,0); ima[2] < ima.dim(2); ima.inc(2)) {
      if (mask) mask->set(1,0); 
      for (ima.set(1,0); ima[1] < ima.dim(1); ima.inc(1)) {
        if (mask) mask->set(0,0); 
        for (ima.set(0,0); ima[0] < ima.dim(0); ima.inc(0)) {

          bool skip = false;
          if (mask) if (mask->value() < 0.5) skip = true;
          if (!skip) {
            float val = ima.value();
            if (isfinite (val)) {
              mean += val;
              std += val*val;
              if (min > val) min = val;
              if (max < val) max = val;
              count++;
            }
          }
          if (mask) mask->inc (0);
        }
        if (mask) mask->inc (1);
      }
      if (mask) mask->inc (2);
    }

    if (count == 0) throw Exception ("no voxels in mask - aborting");

    mean /= double(count);
    std = sqrt(std/double(count) - mean*mean);

    String s = "[ ";
    for (int n = 3; n < ima.ndim(); n++) s += str(ima[n]) + " ";
    s += "] ";

    if (!header_shown) print ("channel         mean        std. dev.   min         max         count\n");
    header_shown = true;
    print (MR::printf ("%-15s %-11g %-11g %-11g %-11g %-11d\n", s.c_str(), mean, std, min, max, count));

  } while (ima++);

}

