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
#include "math/linalg.h"

using namespace std; 
using namespace MR; 

SET_VERSION_DEFAULT;

DESCRIPTION = {
  "generate map of the major eigenvector.",
  NULL
};

ARGUMENTS = {
  Argument ("tensor", "input tensor image", "the input diffusion tensor image.").type_image_in (),
  Argument ("vector", "output vector image", "the output image of the major eigenvector of the diffusion tensor.").type_image_out (),
  Argument::End
};


OPTIONS = { Option::End };


EXECUTE {
  Image::Object &dt_obj (*argument[0].get_image());
  Image::Header header (dt_obj);

  if (header.ndim() != 4) 
    throw Exception ("base image should contain 4 dimensions");

  if (header.dim(3) != 6) 
    throw Exception ("expecting dimension 3 of image \"" + header.name + "\" to be 6");

  header.axes.dim[3] = 3;
  header.data_type = DataType::Float32;

  Math::Matrix V(3,3), M(3,3);
  double ev[3];

  Math::eig_init (M, true);

  Image::Position dt (dt_obj);
  Image::Position vec (*argument[1].get_image (header));

  ProgressBar::init (vec.dim(0)*vec.dim(1)*vec.dim(2), "generating major eigenvector map...");

  for (dt.set(2,0), vec.set(2,0); dt[2] < dt.dim(2); dt.inc(2), vec.inc(2)) {
    for (dt.set(1,0), vec.set(1,0); dt[1] < dt.dim(1); dt.inc(1), vec.inc(1)) {
      for (dt.set(0,0), vec.set(0,0); dt[0] < dt.dim(0); dt.inc(0), vec.inc(0)) {

        dt.set(3,0);

        M(0,0) = dt.value();          dt.inc(3);
        M(1,1) = dt.value();          dt.inc(3);
        M(2,2) = dt.value();          dt.inc(3);
        M(0,1) = M(1,0) = dt.value(); dt.inc(3);
        M(0,2) = M(2,0) = dt.value(); dt.inc(3);
        M(1,2) = M(2,1) = dt.value(); 

        Math::eig (M, ev, V);

        vec.set(3,0);
        vec.value (V(0,2)); vec.inc(3);
        vec.value (V(1,2)); vec.inc(3);
        vec.value (V(2,2));

        ProgressBar::inc();
      }
    }
  }
  ProgressBar::done();
  Math::eig_end();
}
