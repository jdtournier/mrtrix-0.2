/*

    This source has been copied from track_info and then been modified so that
    it converts the mrtrix Tracks to a vtk readble format.
    the source has been rewritten by
    Philip Broser 1/06/2011 Univeristy of Tuebingen, philip.broser@me.com


    The original copyright reads:

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

#include <fstream>
#include <glibmm/stringutils.h>

#include "app.h"
#include "dwi/tractography/file.h"
#include "dwi/tractography/properties.h"


//To get Transformation Matrix

#include "image/object.h"
#include "image/axis.h"

using namespace MR;
using namespace MR::DWI;
using namespace std;



Point scanner_to_image_space (const Image::Object& H, const Point& pos)
{
  const Math::Matrix& M (H.R2I());
  return Point (
      M(0,0)*pos[0] + M(0,1)*pos[1] + M(0,2)*pos[2] + M(0,3),
      M(1,0)*pos[0] + M(1,1)*pos[1] + M(1,2)*pos[2] + M(1,3),
      M(2,0)*pos[0] + M(2,1)*pos[1] + M(2,2)*pos[2] + M(2,3) );
}

Point scanner_to_voxel_space (const Image::Object& H, const Point& pos)
{
  const Math::Matrix& M (H.R2P());
  return Point (
      M(0,0)*pos[0] + M(0,1)*pos[1] + M(0,2)*pos[2] + M(0,3),
      M(1,0)*pos[0] + M(1,1)*pos[1] + M(1,2)*pos[2] + M(1,3),
      M(2,0)*pos[0] + M(2,1)*pos[1] + M(2,2)*pos[2] + M(2,3) );
}


SET_VERSION_DEFAULT;

DESCRIPTION = {
  "convert a track file to a vtk format, cave: coordinates are in XYZ coordinates not reference",
  NULL
};

ARGUMENTS = {
  Argument ("tracks.tck", "track file", "the input track file.").type_file (),
  Argument ("vtkoutputfile.vtk", "track file", "the output vtk file name (use .vtk as suffix)").type_file (),
  Argument::End
};



OPTIONS = {

  Option ("voxel", "voxel space", 
      "if specified, the properties of this image will be used to convert "
      "track point positions from real (scanner) coordinates into voxel coordinates.")
    .append (Argument ("image", " image", "the reference image.").type_image_in ()),

  Option ("image", "image space", 
      "if specified, the properties of this image will be used to convert "
      "track point positions from real (scanner) coordinates into image coordinates (in mm).")
    .append (Argument ("image", " image", "the reference image.").type_image_in ()),

  Option::End
};





EXECUTE {
  RefPtr<Image::Object> reference;
  bool to_voxel = false;

  std::vector<OptBase> opt = get_options (0);
  if (opt.size()) {
    to_voxel = true;
    reference = opt[0][0].get_image();
  }

  opt = get_options (1);
  if (opt.size()) {
    if (reference)
      throw Exception ("options \"-voxel\" and \"-image\" are mutually exclusive");
    to_voxel = false;
    reference = opt[0][0].get_image();
  }

  Tractography::Properties properties;
  Tractography::Reader file;
  file.open (argument[0].get_string(), properties);


  // create and write header of VTK output file:
  String VTKFileName (argument[1].get_string());
  std::ofstream VTKout (VTKFileName.c_str(), std::ios_base::out | std::ios_base::binary);
  if (!VTKout) 
    throw Exception ("error opening file \"" + VTKFileName + "\": " + Glib::strerror (errno));
  
  VTKout << 
    "# vtk DataFile Version 1.0\n"
    "Data values for Tracks\n"
    "ASCII\n"
    "\n" 
    "DATASET POLYDATA\n"
    "POINTS ";
  // keep track of offset to write proper value later:
  guint offset_num_points = VTKout.tellp();
  VTKout << "XXXXXXXXXX float\n";


  ProgressBar::init (0, "writing track data to VTK file");
  std::vector<Point> tck;
  std::vector< std::pair<guint,guint> > track_list;
  guint current_index = 0;

  // write out points, and build index of tracks:
  while (file.next (tck)) {

    guint start_index = current_index;
    current_index += tck.size();
    track_list.push_back (std::pair<guint,guint> (start_index, current_index));

    for (std::vector<Point>::iterator i = tck.begin(); i != tck.end(); ++i) {
      Point pos (*i);
      if (reference) {
        if (to_voxel) pos = scanner_to_voxel_space (*reference, pos);
        else pos = scanner_to_image_space (*reference, pos);
      }
      VTKout << pos[0] << " " << pos[1] << " " << pos[2] << "\n";
    }

    ProgressBar::inc();
  }
  ProgressBar::done();

  // write out list of tracks:
  VTKout << "\nLINES " << track_list.size() << " " << track_list.size() + current_index << "\n";
  for (guint n = 0; n < track_list.size(); ++n) {
    const std::pair<guint,guint>& track (track_list[n]);
    VTKout << track.second - track.first << "\n";
    for (guint i = track.first; i < track.second; ++i)
      VTKout << i << "\n";
    VTKout << "\n";

  };

  // write back total number of points:
  VTKout.seekp (offset_num_points);
  String num_points (str (current_index));
  num_points.resize (10, ' ');
  VTKout.write (num_points.c_str(), 10);

  VTKout.close();
}

