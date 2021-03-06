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


    04-11-2009 J-Donald Tournier <d.tournier@brain.org.au>
    * fix -axis option

    09-12-2010 J-Donald Tournier <d.tournier@brain.org.au>
    * fix -axis option again

*/

#include "app.h"
#include "image/position.h"

using namespace std; 
using namespace MR; 

SET_VERSION_DEFAULT;

DESCRIPTION = {
  "clean up warp field generated by ANTS (or more specifically, WarpMultiImageTransform, "
  "replacing background voxels with NaN when they map to regions outside the original image.",
  NULL
};

ARGUMENTS = {
  Argument ("unit_warp", "unit warp", "the original unit warp image").type_image_in (),
  Argument ("ants_warp", "ANTS warp", "the warp image as produced by WarpMultiImageTransform.").type_image_in (),
  Argument ("clean_warp", "clean warp", "the output cleaned-up warp image.").type_image_out (),
  Argument::End
};

OPTIONS = { 
  Option::End 
};




EXECUTE {

  float background[3];
  {
    Image::Object &unit_warp_obj (*argument[0].get_image());
    if (unit_warp_obj.ndim() < 4) 
      throw Exception ("expected 4D image for unit_warp - aborting");
    if (unit_warp_obj.dim(3) < 3) 
      throw Exception ("expected unit_warp image to contain 3 volumes - aborting");

    Image::Position bg (unit_warp_obj);
    bg.set (0,0);
    bg.set (1,0);
    bg.set (2,0);
    bg.set (3,0);
    background[0] = bg.value();
    bg.inc (3);
    background[1] = bg.value();
    bg.inc (3);
    background[2] = bg.value();
  }


  Image::Object &warp_obj (*argument[1].get_image());
  Image::Header header (warp_obj);

  Image::Position in (warp_obj);
  Image::Position out (*argument[2].get_image (header));
  
  ProgressBar::init (out.dim(0)*out.dim(1)*out.dim(2), "cleaning up ANTS transform...");

  for (out.set(2,0), in.set(2,0); out[2] < out.dim(2); out.inc(2), in.inc(2)) {
    for (out.set(1,0), in.set(1,0); out[1] < out.dim(1); out.inc(1), in.inc(1)) {
      for (out.set(0,0), in.set(0,0); out[0] < out.dim(0); out.inc(0), in.inc(0)) {
        float val[3];

        in.set(3,0);
        val[0] = in.value();
        in.inc(3);
        val[1] = in.value();
        in.inc(3);
        val[2] = in.value();

        if (val[0] == background[0] && val[1] == background[1] && val[2] == background[2]) {
          out.set (3,0);
          out.value (GSL_NAN);
          out.inc(3);
          out.value (GSL_NAN);
          out.inc(3);
          out.value (GSL_NAN);
        }
        else {
          out.set (3,0);
          out.value (val[0]);
          out.inc(3);
          out.value (val[1]);
          out.inc(3);
          out.value (val[2]);
        }
      }
    }
  }

  ProgressBar::done();
}
