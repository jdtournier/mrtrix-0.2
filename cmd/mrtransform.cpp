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
#include "image/interp.h"
#include "math/linalg.h"

using namespace std; 
using namespace MR; 

SET_VERSION_DEFAULT;

DESCRIPTION = {
 "apply spatial transformations or reslice images.",
  NULL
};

ARGUMENTS = {
  Argument ("input", "input image", "input image to be transformed.", true, true).type_image_in (),
  Argument ("output", "output image", "the output image.").type_image_out (),
  Argument::End
};


OPTIONS = {

  Option ("transform", "the transform to use", "specified the 4x4 transform to apply.")
    .append (Argument ("transform", "transform", "the transform to apply, in the form of a 4x4 ascii file.").type_file ()),

  Option ("replace", "replace transform", "replace the current transform by that specified, rather than applying it to the current transform."),

  Option ("inverse", "use inverse transform", "invert the specified transform before using it."),

  Option ("template", "template image", "reslice the input image to match the specified template image.")
    .append (Argument ("image", "image", "the template image.").type_image_in ()),

  Option ("reference", "reference image for transform", "in case the transform supplied maps from the input image onto a reference image, use this option to specify the reference. Noe that this implicitly sets the -replace option.")
    .append (Argument ("image", "image", "the reference image.").type_image_in ()),

  Option ("flipx", "assume x-flipped transform", "assume the transform is supplied assuming a coordinate system with the x-axis reversed relative to the MRtrix convention (i.e. x increases from right to left). This is required to handle transform matrices produced by FSL's FLIRT command. This is only used in conjunction with the -reference option."),

  Option::End 
};



EXECUTE {
  Math::Matrix T(4,4);
  T.identity();
  bool transform_supplied = false;

  std::vector<OptBase> opt = get_options (0); // transform
  if (opt.size()) {
    transform_supplied = true;
    T.load (opt[0][0].get_string());
    if (T.rows() != 4 || T.columns() != 4) 
      throw Exception (String("transform matrix supplied in file \"") + opt[0][0].get_string() + "\" is not 4x4");
  }

  bool replace = get_options(1).size(); // replace
  bool inverse = get_options(2).size(); // inverse

  Image::Object& in_obj (*argument[0].get_image());
  Image::Header header (in_obj);

  if (transform_supplied) {

    if (inverse) {
      Math::Matrix I;
      Math::invert (I, T);
      T = I;
    }

    opt = get_options(4); // reference 
    if (opt.size() && transform_supplied) {
      replace = true;
      const Image::Header& ref_header (opt[0][0].get_image()->header());

      if (get_options(5).size()) { // flipx 
        Math::Matrix R(4,4);
        R.identity();
        R(0,0) = -1.0;
        R(0,3) = (ref_header.axes.dim[0]-1) * ref_header.axes.vox[0];

        Math::Matrix M;
        M.multiply (R, T);

        R(0,3) = (header.axes.dim[0]-1) * header.axes.vox[0];

        T.multiply (M, R);
      }

      Math::Matrix M;
      M.multiply (ref_header.transform(), T);
      T = M;
    }



    if (replace) header.set_transform (T);
    else {
      Math::Matrix M;
      M.multiply (T, header.transform());
      header.set_transform (M);
    }

    header.comments.push_back ("transform modified");
  }


  opt = get_options(3); // template : need to reslice
  if (opt.size()) {
    Math::Matrix Mi (header.R2P());

    Image::Header template_header (opt[0][0].get_image()->header());
    header.axes.dim[0] = template_header.axes.dim[0];
    header.axes.dim[1] = template_header.axes.dim[1];
    header.axes.dim[2] = template_header.axes.dim[2];
    header.axes.vox[0] = template_header.axes.vox[0];
    header.axes.vox[1] = template_header.axes.vox[1];
    header.axes.vox[2] = template_header.axes.vox[2];
    header.set_transform (template_header.transform());
    header.comments.push_back ("resliced to reference image \"" + template_header.name + "\"");

    Math::Matrix M;
    M.multiply (Mi, header.P2R());
    VAR (M);

    float R[] = { 
      M(0,0), M(0,1), M(0,2), M(0,3), 
      M(1,0), M(1,1), M(1,2), M(1,3), 
      M(2,0), M(2,1), M(2,2), M(2,3)
    };

    in_obj.optimise();
    Image::Interp in (in_obj);
    Image::Position out (*argument[1].get_image (header));
    Point pos;

    ProgressBar::init (out.voxel_count(), "reslicing image...");
    do { 
      for (out.set(2,0); out[2] < out.dim(2); out.inc(2)) {
        for (out.set(1,0); out[1] < out.dim(1); out.inc(1)) {
          pos[0] = R[1]*out[1] + R[2]*out[2] + R[3];
          pos[1] = R[5]*out[1] + R[6]*out[2] + R[7];
          pos[2] = R[9]*out[1] + R[10]*out[2] + R[11];
          for (out.set(0,0); out[0] < out.dim(0); out.inc(0)) {
            in.P (pos);
            if (!in) out.value (0.0);
            else out.value (in.value());
            pos[0] += R[0];
            pos[1] += R[4];
            pos[2] += R[8];
            ProgressBar::inc();
          }
        }
      }
    } while (out++);
    ProgressBar::done();
  }
  else {
    Image::Position in (in_obj);
    Image::Position out (*argument[1].get_image (header));
    ProgressBar::init (out.voxel_count(), "copying image data...");
    do { 
      out.value (in.value());
      in++;
      ProgressBar::inc();
    } while (out++);
    ProgressBar::done();
  }
}

