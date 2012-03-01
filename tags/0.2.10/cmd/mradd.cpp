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


    05-11-2010 J-Donald Tournier <d.tournier@brain.org.au>
    * fix incorrect progressbar count, and handling of large numbers of input
    * files.

*/

#include "app.h"
#include "image/position.h"

using namespace std; 
using namespace MR; 

SET_VERSION_DEFAULT;

DESCRIPTION = {
  "add or subtract images",
  NULL
};

ARGUMENTS = {
  Argument ("image1", "first input image", "the first input image.").type_image_in (),
  Argument ("image2", "second input image", "the second input image.", true, true).type_image_in (),
  Argument ("output", "output image", "the output image.").type_image_out (),
  Argument::End
};

OPTIONS = { Option::End };




EXECUTE {
  guint num_images = argument.size()-1;

  RefPtr<Image::Object> in[num_images];
  in[0] = argument[0].get_image();
  Image::Header header (in[0]->header());

  header.data_type = DataType::Float32;

  for (guint i = 1; i < num_images; i++) {
    in[i] = argument[i].get_image();

    if (in[i]->is_complex()) header.data_type = DataType::CFloat32;

    if (in[i]->ndim() > header.axes.ndim()) 
      header.axes.set_ndim (in[i]->ndim());

    for (int n = 0; n < header.axes.ndim(); n++) { 
      if (header.axes.dim[n] != in[i]->dim(n)) {
        if (header.axes.dim[n] < 2) header.axes.copy (n, in[i]->header().axes, n);
        else if (in[i]->dim(n) > 1) throw Exception ("dimension mismatch between input files");
      }
    }
  }



  Image::Position out (*argument[num_images].get_image (header));


  ProgressBar::init (num_images * out.voxel_count(), "adding...");

  for (guint i = 0; i < num_images; i++) {
    Image::Position y (*in[i]);

    do {

      for (int n = 0; n < y.ndim(); n++)
        y.set (n, y.dim(n) > 1 ? out[n] : 0);

      float val = i ? out.re() : 0.0;
      out.re (val + y.re());
      
      if (out.is_complex()) {
        val = i ? out.im() : 0.0;
        if (y.is_complex()) val += y.im();
        out.im (val);
      }

      ProgressBar::inc();

    } while (out++);

    in[i] = NULL;
  }

  ProgressBar::done();
}

