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


    15-10-2008 J-Donald Tournier <d.tournier@brain.org.au>
    * fix -prs option handling
    * remove MR::DICOM_DW_gradients_PRS flag

    15-10-2008 J-Donald Tournier <d.tournier@brain.org.au>
    * add -layout option to manipulate data ordering within the image file

    14-02-2010 J-Donald Tournier <d.tournier@brain.org.au>
    * fix -coord option so that the "end" keyword can be used


*/

#include "app.h"
#include "image/position.h"
#include "image/axis.h"
#include "math/linalg.h"

using namespace std; 
using namespace MR; 

SET_VERSION_DEFAULT;

DESCRIPTION = {
  "perform conversion between different file types and optionally extract a subset of the input image.",
  "If used correctly, this program can be a very useful workhorse. In addition to converting images between different formats, it can be used to extract specific studies from a data set, extract a specific region of interest, flip the images, or to scale the intensity of the images.",
  NULL
};

ARGUMENTS = {
  Argument ("input", "input image", "the input image.").type_image_in (),
  Argument ("ouput", "output image", "the output image.").type_image_out (),
  Argument::End
};


const gchar* type_choices[] = { "REAL", "IMAG", "MAG", "PHASE", "COMPLEX", NULL };
const gchar* data_type_choices[] = { "FLOAT32", "FLOAT32LE", "FLOAT32BE", "FLOAT64", "FLOAT64LE", "FLOAT64BE", 
    "INT32", "UINT32", "INT32LE", "UINT32LE", "INT32BE", "UINT32BE", 
    "INT16", "UINT16", "INT16LE", "UINT16LE", "INT16BE", "UINT16BE", 
    "CFLOAT32", "CFLOAT32LE", "CFLOAT32BE", "CFLOAT64", "CFLOAT64LE", "CFLOAT64BE", 
    "INT8", "UINT8", "BIT", NULL };

OPTIONS = {
  Option ("coord", "select coordinates", "extract data only at the coordinates specified.", false, true)
    .append (Argument ("axis", "axis", "the axis of interest").type_integer (0, INT_MAX, 0))
    .append (Argument ("coord", "coordinates", "the coordinates of interest").type_sequence_int()),

  Option ("vox", "voxel size", "change the voxel dimensions.")
    .append (Argument ("sizes", "new dimensions", "A comma-separated list of values. Only those values specified will be changed. For example: 1,,3.5 will change the voxel size along the x & z axes, and leave the y-axis voxel size unchanged.")
        .type_sequence_float ()),

  Option ("datatype", "data type", "specify output image data type.")
    .append (Argument ("spec", "specifier", "the data type specifier.").type_choice (data_type_choices)),

  Option ("scale", "scaling factor", "apply scaling to the intensity values.")
    .append (Argument ("factor", "factor", "the factor by which to multiply the intensities.").type_float (NAN, NAN, 1.0)),

  Option ("offset", "offset", "apply offset to the intensity values.")
    .append (Argument ("bias", "bias", "the value of the offset.").type_float (NAN, NAN, 0.0)),

  Option ("zero", "replace NaN by zero", "replace all NaN values with zero."),

  Option ("output", "output type", "specify type of output")
    .append (Argument ("type", "type", "type of output.")
        .type_choice (type_choices)),

  Option ("layout", "data layout", "specify the layout of the data in memory. The actual layout produced will depend on whether the output image format can support it.")
    .append (Argument ("spec", "specifier", "the data layout specifier.").type_string ()),

  Option ("prs", "DW gradient specified as PRS", "assume that the DW gradients are specified in the PRS frame (Siemens DICOM only)."),

  Option ("nocomment", "remove comments", "remove all comments from the header. Useful to anonymise data."),

  Option ("addcomment", "add comments", "add a new comment to the header.", false, true)
    .append (Argument ("comment", "comment", "the text to add as a comment.").type_string()),

  Option::End
};



inline bool next (Image::Position& ref, Image::Position& other, const std::vector<int>* pos)
{
  int axis = 0;
  do {
    ref.inc (axis);
    if (ref[axis] < ref.dim(axis)) {
      other.set (axis, pos[axis][ref[axis]]);
      return (true);
    }
    ref.set (axis, 0);
    other.set (axis, pos[axis][0]);
    axis++;
  } while (axis < ref.ndim());
  return (false);
}





EXECUTE {
  std::vector<OptBase> opt = get_options (1); // vox
  std::vector<float> vox;
  if (opt.size()) 
    vox = parse_floats (opt[0][0].get_string());


  opt = get_options (3); // scale
  float scale = 1.0;
  if (opt.size()) scale = opt[0][0].get_float();

  opt = get_options (4); // offset
  float offset = 0.0;
  if (opt.size()) offset = opt[0][0].get_float();

  opt = get_options (5); // zero
  bool replace_NaN = opt.size();

  opt = get_options (6); // output
  Image::OutputType output_type = Image::Default;
  if (opt.size()) {
    switch (opt[0][0].get_int()) {
      case 0: output_type = Image::Real; break;
      case 1: output_type = Image::Imaginary; break;
      case 2: output_type = Image::Magnitude; break;
      case 3: output_type = Image::Phase; break;
      case 4: output_type = Image::RealImag; break;
    }
  }




  Image::Object &in_obj (*argument[0].get_image());

  Image::Header header (in_obj);

  if (output_type == 0) {
    if (in_obj.is_complex()) output_type = Image::RealImag;
    else output_type = Image::Default;
  }

  if (output_type == Image::RealImag) header.data_type = DataType::CFloat32;
  else if (output_type == Image::Phase) header.data_type = DataType::Float32;
  else header.data_type.unset_flag (DataType::ComplexNumber);

  
  opt = get_options (2); // datatype
  if (opt.size()) header.data_type.parse (data_type_choices[opt[0][0].get_int()]);

  for (guint n = 0; n < vox.size(); n++) 
    if (isfinite (vox[n])) header.axes.vox[n] = vox[n];

  opt = get_options (7); // layout
  if (opt.size()) {
    std::vector<Image::Axis> ax = parse_axes_specifier (header.axes, opt[0][0].get_string());
    if (ax.size() != (guint) header.axes.ndim()) 
      throw Exception (String("specified layout \"") + opt[0][0].get_string() + "\" does not match image dimensions");

    for (guint i = 0; i < ax.size(); i++) {
      header.axes.axis[i] = ax[i].axis;
      header.axes.forward[i] = ax[i].forward;
    }
  }


  opt = get_options (8); // prs
  if (opt.size() && header.DW_scheme.rows() && header.DW_scheme.columns()) {
    for (guint row = 0; row < header.DW_scheme.rows(); row++) {
      double tmp = header.DW_scheme(row, 0);
      header.DW_scheme(row, 0) = header.DW_scheme(row, 1);
      header.DW_scheme(row, 1) = tmp;
      header.DW_scheme(row, 2) = -header.DW_scheme(row, 2);
    }
  }

  opt = get_options (9); // nocomment
  if (opt.size()) 
    header.comments.clear();

  opt = get_options (10); // addcomment
  for (guint n = 0; n < opt.size(); ++n) 
    header.comments.push_back (opt[n][0].get_string());


  std::vector<int> pos[in_obj.ndim()];

  opt = get_options (0); // coord
  for (guint n = 0; n < opt.size(); n++) {
    int axis = opt[n][0].get_int();
    if (pos[axis].size()) throw Exception ("\"coord\" option specified twice for axis " + str (axis));
    pos[axis] = parse_ints (opt[n][1].get_string(), header.dim(axis)-1);
    header.axes.dim[axis] = pos[axis].size();
  }

  for (int n = 0; n < in_obj.ndim(); n++) {
    if (pos[n].empty()) { 
      pos[n].resize (in_obj.dim(n));
      for (guint i = 0; i < pos[n].size(); i++) pos[n][i] = i;
    }
  }


  in_obj.apply_scaling (scale, offset);






  Image::Position in (in_obj);
  Image::Position out (*argument[1].get_image (header));

  for (int n = 0; n < in.ndim(); n++) in.set (n, pos[n][0]);

  ProgressBar::init (out.voxel_count(), "copying data...");

  do { 

    float re, im = 0.0;
    in.get (output_type, re, im);
    if (replace_NaN) if (gsl_isnan (re)) re = 0.0;
    out.re (re);

    if (output_type == Image::RealImag) {
      if (replace_NaN) if (gsl_isnan (im)) im = 0.0;
      out.im (im);
    }

    ProgressBar::inc();
  } while (next (out, in, pos));

  ProgressBar::done();
}





