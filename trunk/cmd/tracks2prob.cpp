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


    31-10-2008 J-Donald Tournier <d.tournier@brain.org.au>
    * various optimisations to improve performance

*/

#include <fstream>
#include <glibmm/stringutils.h>

#include "app.h"
#include "get_set.h"
#include "dwi/tractography/file.h"
#include "dwi/tractography/properties.h"

using namespace MR; 
using namespace MR::DWI; 
using namespace std; 

SET_VERSION_DEFAULT;

DESCRIPTION = {
  "convert a tracks file into a map of the fraction of tracks to enter each voxel.",
  NULL
};

ARGUMENTS = {
  Argument ("tracks", "track file", "the input track file.").type_file (),
  Argument ("template", "template image", "an image file to be used as a template for the output (the output image wil have the same voxel size and dimensions).").type_image_in(),
  Argument ("output", "output image", "the output fraction image").type_image_out(),
  Argument::End
};



OPTIONS = {
  Option ("count", "output fibre count", "produce an image of the fibre count through each voxel, rather than the fraction."),
  Option::End
};





EXECUTE {
  bool fibre_count = get_options(0).size();

  Image::Object& template_image (*argument[1].get_image());
  Image::Header header (template_image.header());



  Tractography::Properties properties;
  Tractography::Reader file;
  file.open (argument[0].get_string(), properties);

  header.data_type = fibre_count ? DataType::UInt32 : DataType::Float32;
  header.axes.set_ndim (3);
  header.comments.push_back (String ("track ") + (fibre_count ? "count" : "fraction") + " map");
  for (Tractography::Properties::iterator i = properties.begin(); i != properties.end(); ++i) 
    header.comments.push_back (i->first + ": " + i->second);
  for (std::vector<RefPtr<Tractography::ROI> >::iterator i = properties.roi.begin(); i != properties.roi.end(); ++i)
    header.comments.push_back ("ROI: " + (*i)->specification());
  for (std::vector<String>::iterator i = properties.comments.begin(); i != properties.comments.end(); ++i)
    header.comments.push_back ("comment: " + *i);

  guint total_count = properties["total_count"].empty() ? 0 : to<guint> (properties["total_count"]);
  guint num_tracks = properties["count"].empty() ? 0 : to<guint> (properties["count"]);
  guint count = 0;

  guint voxel_count = header.voxel_count(3);
  int xmax = header.dim(0);
  int ymax = header.dim(1);
  int zmax = header.dim(2);

  guint yskip = xmax;
  guint zskip = xmax*ymax;

  guint* countbuf = new guint [voxel_count];
  memset (countbuf, 0, voxel_count*sizeof(guint));

  voxel_count = (voxel_count+7)/8;
  guint8* visited = new guint8 [voxel_count];
  std::vector<Point> tck;

  Image::Interp interp (template_image);
  ProgressBar::init (num_tracks, "generating track count image...");

  while (file.next (tck)) {
    memset (visited, 0, voxel_count*sizeof(guint8));
    for (std::vector<Point>::iterator i = tck.begin(); i != tck.end(); ++i) {
      Point p (interp.R2P (*i));
      int x = round (p[0]);
      int y = round (p[1]);
      int z = round (p[2]);
      if (x >= 0 && y >= 0 && z >= 0 && x < xmax && y < ymax && z < zmax) {
        guint offset = x + yskip*y + zskip*z;
        if (!get<bool>(visited, static_cast<gsize>(offset))) {
          put<bool> (true, visited, static_cast<gsize>(offset));
          countbuf[offset]++;
        }
      }
    }
    count++;
    ProgressBar::inc();
  }
  ProgressBar::done();

  delete [] visited;


  Image::Position pos (*argument[2].get_image (header));
  if (total_count > 0) count = total_count;
  float multiplier = 1.0;
  if (!fibre_count) multiplier /= float(count);

  ProgressBar::init (xmax*ymax*zmax, "writing track count image...");
  for (pos.set (2,0); pos[2] < zmax; pos.inc(2)) {
    for (pos.set(1,0); pos[1] < ymax; pos.inc(1)) {
      for (pos.set(0,0); pos[0] < xmax; pos.inc(0)) {
        pos.value (multiplier * countbuf[pos[0] + yskip*pos[1] + zskip*pos[2]]);
        ProgressBar::inc();
      }
    }
  }
  ProgressBar::done();

  delete [] countbuf;
}

