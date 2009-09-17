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


    28-07-2008 J-Donald Tournier <d.tournier@brain.org.au>
    * fix option parsing to allow multiple ignoreslices and ignorestudies instances
    
    17-09-2009 J-Donald Tournier <d.tournier@brain.org.au>
    * improved support for images that might use a different dimension for the DWI
    
*/

#include "app.h"
#include "image/position.h"
#include "math/matrix.h"
#include "math/linalg.h"
#include "dwi/gradient.h"

using namespace std; 
using namespace MR; 

SET_VERSION_DEFAULT;

DESCRIPTION = {
 "convert diffusion-weighted images to tensor images.",
  NULL
};

ARGUMENTS = {
  Argument ("dwi", "input DW image", "the input diffusion-weighted image.").type_image_in (),
  Argument ("tensor", "output tensor image", "the output diffusion tensor image.").type_image_out (),
  Argument::End
};


OPTIONS = {
  Option ("grad", "supply gradient encoding", "specify the diffusion-weighted gradient scheme used in the acquisition. The program will normally attempt to use the encoding stored in image header.")
    .append (Argument ("encoding", "gradient encoding", "the gradient encoding, supplied as a 4xN text file with each line is in the format [ X Y Z b ], where [ X Y Z ] describe the direction of the applied gradient, and b gives the b-value in units (1000 s/mm^2).").type_file ()),

  Option ("ignoreslices", "ignore image slices", "ignore the image slices specified when computing the tensor.", false, true)
    .append (Argument ("slice", "slice z coordinate", "the z coordinate of the slice to be ignored").type_integer (0, G_MAXINT, 0))
    .append (Argument ("volumes", "study numbers", "the volume numbers of the slice to be ignored").type_sequence_int ()),

  Option ("ignorevolumes", "ignore image volumes", "ignore the image volumes specified when computing the tensor.", false, true)
    .append (Argument ("volumes", "volume numbers", "the volumes to be ignored").type_sequence_int ()),

  Option::End
};


EXECUTE {
  Image::Object &dwi_obj (*argument[0].get_image());
  Image::Header header (dwi_obj);

  if (header.ndim() < 4) 
    throw Exception ("dwi image should contain at least 4 dimensions");

  int axis = 3;
  while (header.dim(axis) <= 1 && axis < header.ndim()) axis++;

  Math::Matrix grad, bmat, binv;

  std::vector<OptBase> opt = get_options (0);
  if (opt.size()) grad.load (opt[0][0].get_string());
  else {
    if (!header.DW_scheme.is_valid()) 
      throw Exception ("no diffusion encoding found in image \"" + header.name + "\"");
    grad.copy (header.DW_scheme);
  }

  if (grad.rows() < 7 || grad.columns() != 4) 
    throw Exception ("unexpected diffusion encoding matrix dimensions");

  info ("found " + str(grad.rows()) + "x" + str(grad.columns()) + " diffusion-weighted encoding");

  if (header.dim(axis) != (int) grad.rows()) 
    throw Exception ("number of volumes in base image does not match that in encoding file");

  DWI::normalise_grad (grad);
  DWI::grad2bmatrix (bmat, grad);

  grad.allocate (bmat.rows(), bmat.columns());
  binv.allocate (bmat.columns(), bmat.rows());

  Math::PseudoInverter pinverter (binv, grad);


  std::vector<std::vector<int> > islc (header.dim(2));
  std::vector<int> ivol;

  opt = get_options (1);
  for (guint n = 0; n < opt.size(); n++) {
    int z = opt[n][0].get_int();
    if (z >= (int) islc.size()) throw Exception ("slice number out of bounds");
    islc[z] = parse_ints (opt[n][1].get_string());
  }

  opt = get_options (2);
  for (guint n = 0; n < opt.size(); n++) {
    std::vector<int> v = parse_ints (opt[n][0].get_string());
    ivol.insert (ivol.end(), v.begin(), v.end());
  }


  header.axes.set_ndim (4);
  header.axes.dim[3] = 6;
  header.data_type = DataType::Float32;
  header.DW_scheme.reset();

  Image::Position dwi (dwi_obj);
  Image::Position dt (*argument[1].get_image (header));

  info ("converting base image \"" + dwi.name() + " to tensor image \"" + dt.name() + "\"");

  ProgressBar::init (dwi.dim(0)*dwi.dim(1)*dwi.dim(2), "converting DW images to tensor image..."); 

  float d[dwi.dim(axis)];
  for (dwi.set(2,0), dt.set(2,0); dwi[2] < dwi.dim(2); dwi.inc(2), dt.inc(2)) {

    grad.copy (bmat);
    for (guint i = 0; i < ivol.size(); i++)
      for (int j = 0; j < 7; j++)
        grad (ivol[i],j) = 0.0;

    for (guint i = 0; i < islc[dt[2]].size(); i++)
      for (int j = 0; j < 7; j++)
        grad (islc[dt[2]][i],j) = 0.0;

    pinverter.invert (binv, grad);

    for (dwi.set(1,0), dt.set(1,0); dwi[1] < dwi.dim(1); dwi.inc(1), dt.inc(1)) {
      for (dwi.set(0,0), dt.set(0,0); dwi[0] < dwi.dim(0); dwi.inc(0), dt.inc(0)) {
        for (dwi.set(axis,0); dwi[axis] < dwi.dim(axis); dwi.inc(axis)) {
          d[dwi[axis]] = dwi.value();
          d[dwi[axis]] = d[dwi[axis]] > 0.0 ? -log (d[dwi[axis]]) : 1e-12;
        }
        for (dt.set(3,0); dt[3] < dt.dim(3); dt.inc(3)) {
          float val = 0.0;
          for (int i = 0; i < dwi.dim(axis); i++)
            val += (float) (binv(dt[3], i)*d[i]);
          dt.value (val);
        }
        ProgressBar::inc();
      }
    }
  }

  ProgressBar::done();
}

