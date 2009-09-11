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
#include "dwi/tensor.h"

using namespace std; 
using namespace MR; 

SET_VERSION_DEFAULT;

DESCRIPTION = {
 "calculate map of mean apparent diffusion coefficient (ADC) from diffusion tensor image.",
  NULL
};

ARGUMENTS = {
  Argument ("tensor", "input tensor image", "the input diffusion tensor image.").type_image_in (),
  Argument ("ADC", "output ADC image", "the output mean apparent diffusion coefficient image.").type_image_out (),
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

  header.axes.set_ndim (3);
  header.data_type = DataType::Float32;

  Image::Position dt (dt_obj);
  Image::Position adc (*argument[1].get_image (header));

  float buf[6];

  ProgressBar::init(adc.voxel_count(), "generating ADC map...");

  for (adc.set(2,0), dt.set(2,0); adc[2] < adc.dim(2); adc.inc(2), dt.inc(2)) {
    for (adc.set(1,0), dt.set(1,0); adc[1] < adc.dim(1); adc.inc(1), dt.inc(1)) {
      for (adc.set(0,0), dt.set(0,0); adc[0] < adc.dim(0); adc.inc(0), dt.inc(0)) {

        for (dt.set(3,0); dt[3] < 6; dt.inc(3)) 
          buf[dt[3]] = dt.value();

        adc.value (DWI::tensor2ADC (buf));
      }
    }
    ProgressBar::inc();
  }
  ProgressBar::done();
}

