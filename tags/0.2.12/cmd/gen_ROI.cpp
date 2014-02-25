/*
    Copyright 2009 Brain Research Institute, Melbourne, Australia

    Written by J-Donald Tournier, 12/08/09.

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
#include "image/axis.h"

using namespace std; 
using namespace MR; 

SET_VERSION_DEFAULT;

DESCRIPTION = {
  "generate a mask ROI with the same dimensions as the input dataset, and all voxel values set to one in the regions specified",
  NULL
};

ARGUMENTS = {
  Argument ("input", "input image", "the input image.").type_image_in (),
  Argument ("ouput", "output image", "the output image.").type_image_out (),
  Argument::End
};


const gchar* data_type_choices[] = { "FLOAT32", "FLOAT32LE", "FLOAT32BE", "FLOAT64", "FLOAT64LE", "FLOAT64BE", 
    "INT32", "UINT32", "INT32LE", "UINT32LE", "INT32BE", "UINT32BE", 
    "INT16", "UINT16", "INT16LE", "UINT16LE", "INT16BE", "UINT16BE", 
    "CFLOAT32", "CFLOAT32LE", "CFLOAT32BE", "CFLOAT64", "CFLOAT64LE", "CFLOAT64BE", 
    "INT8", "UINT8", "BIT", NULL };

OPTIONS = {
  Option ("coord", "select coordinates", "extract data only at the coordinates specified.", false, true)
    .append (Argument ("axis", "axis", "the axis of interest").type_integer (0, INT_MAX, 0))
    .append (Argument ("coord", "coordinates", "the coordinates of interest").type_sequence_int()),

  Option ("datatype", "data type", "specify output image data type.")
    .append (Argument ("spec", "specifier", "the data type specifier.").type_choice (data_type_choices)),

  Option::End
};



inline bool next (Image::Position& ref, Image::Position& other, const std::vector<int>* pos)
{
  int axis = 0;
  do {
    ref.inc (axis);
    other.set (axis, pos[axis][ref[axis]]);
    if (ref[axis] < ref.dim(axis)) return (true);
    ref.set (axis, 0);
    other.set (axis, pos[axis][0]);
    axis++;
  } while (axis < ref.ndim());
  return (false);
}





EXECUTE {
  Image::Object &in_obj (*argument[0].get_image());
  Image::Header header (in_obj);
  header.axes.set_ndim (3);
  header.data_type = DataType::Bit;
  header.offset = 0.0;
  header.scale = 1.0;

  std::vector< std::vector<int> > pos (in_obj.ndim());

  std::vector<OptBase> opt = get_options (0); // coord
  for (guint n = 0; n < opt.size(); n++) {
    int axis = opt[n][0].get_int();
    if (axis >= 3) throw Exception ("axis must be between 0 & 2 for option argument 1 of option \"coord\"");
    if (pos[axis].size()) throw Exception ("\"coord\" option specified twice for axis " + str (axis));
    pos[axis] = parse_ints (opt[n][1].get_string());
  }

  for (int n = 0; n < in_obj.ndim(); n++) {
    if (pos[n].empty()) { 
      pos[n].resize (in_obj.dim(n));
      for (guint i = 0; i < pos[n].size(); i++) pos[n][i] = i;
    }
  }

  opt = get_options (1); // datatype
  if (opt.size()) header.data_type.parse (data_type_choices[opt[0][0].get_int()]);


  Image::Position in (in_obj);
  Image::Position out (*argument[1].get_image (header));

  ProgressBar::init (out.voxel_count(), "generating ROI...");

  for (out.set(2,0); out[2] < out.dim(2); out.inc(2)) 
    for (out.set(1,0); out[1] < out.dim(1); out.inc(1)) 
      for (out.set(0,0); out[0] < out.dim(0); out.inc(0)) 
        out.value (0.0);

  for (guint z = 0; z < pos[2].size(); z++) {
    out.set(2, pos[2][z]);
    for (guint y = 0; y < pos[1].size(); y++) {
      out.set(1, pos[1][y]);
      for (guint x = 0; x < pos[0].size(); x++) {
        out.set(0, pos[0][x]);
        out.value (1.0);
        ProgressBar::inc();
      }
    }
  }

  ProgressBar::done();
}






