/*
    Copyright 2009 Brain Research Institute, Melbourne, Australia

    Written by J-Donald Tournier, 13/08/09.

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
#include "get_set.h"
#include "dwi/tractography/file.h"
#include "dwi/tractography/properties.h"

using namespace MR; 
using namespace MR::DWI; 
using namespace std; 

SET_VERSION_DEFAULT;

DESCRIPTION = {
  "Resample tracks so they all have the same number of samples, spaced at constant intervals within the range specified.",
  "Use the -extent option to specify the range within which the resampling should take place. For highly curved tracks, you may need to also use the -direction option to ensure the resampling is done correctly. In this case, the start direction specifies the plane through the start point from which the resampling will take place, and the end direction specifies the plane through the end point where the resampling will stop. In both cases, the direction vector should point in the direction of increasing sample point index.",
  NULL
};

ARGUMENTS = {
  Argument ("tracks", "track file", "the input track file.").type_file (),
  Argument ("output", "output file", "the output track file").type_file(),
  Argument::End
};



OPTIONS = {

  Option ("extent", "extent", "specify the start and end of the resampling range (default is 0,0,-50 to 0,0,50).")
    .append (Argument ("start", "start point", "a point lying within the start of the resampling range.").type_sequence_float ())
    .append (Argument ("end", "end point", "a point lying within the end of the resampling range.").type_sequence_float ()),

  Option ("direction", "direction", "specify the direction of the resampling range at the start and end points (default is the vector pointing from the start point to the end point for both start and end directions).")
    .append (Argument ("start", "start direction", "a vector pointing along the direction of the resampling at the start of the range.").type_sequence_float ())
    .append (Argument ("end", "end direction", "a vector pointing along the direction of the resampling at the end of the range.").type_sequence_float ()),

  Option ("warp", "warp image", "specify an image containing the warp field to the space in which the resampling is to take place.")
    .append (Argument ("image", "warp image", "the warp image.").type_image_in()),
 
  Option ("num", "number of samples", "specify the number of samples to take along each track.")
    .append (Argument ("nsamples", "number of samples", "the number of samples.").type_integer(1, 100000, 100)),
 
  Option::End 
};




inline Point get_pos (const Point& pos, RefPtr<Image::Interp>& warp)
{
  if (!warp) return (pos);
  warp->R (pos);
  Point ret;
  if (!(*warp)) return (ret);
  warp->set(3,0);
  ret[0] = warp->value();
  warp->inc(3);
  ret[1] = warp->value();
  warp->inc(3);
  ret[2] = warp->value();
  return (ret);
}




inline int get_limits (
    const std::vector<Point>& tck, 
    RefPtr<Image::Interp>& warp,
    const Point& start, 
    const Point& mid,
    const Point& end,
    const Point& start_dir, 
    const Point& mid_dir,
    const Point& end_dir, 
    guint& idx_start, 
    guint& idx_end) 
{
  idx_start = idx_end = 0;
  guint i = 0;
  
  while (( start_dir.dot (get_pos (tck[i], warp) - start) > 0.0 || 
        mid_dir.dot (get_pos (tck[i], warp) - mid) > 0.0 ) && i < tck.size()) i++;
  if (i >= tck.size()) return (1);
  idx_start = i;

  while (( end_dir.dot (get_pos (tck[i], warp) - end) < 0.0 || 
        mid_dir.dot (get_pos (tck[i], warp) - mid) < 0.0 ) && i < tck.size()) {
    if (start_dir.dot (get_pos (tck[i], warp) - start) <= 0.0 && mid_dir.dot (get_pos (tck[i], warp) - mid) < 0.0) idx_start = i;
    i++;
  }
  if (i >= tck.size()) return (1);
  idx_end = i;

  return (0);
}




EXECUTE {
  Tractography::Properties properties;
  Tractography::Reader file;
  file.open (argument[0].get_string(), properties);
  
  Point start (0.0, 0.0, -50.0);
  Point mid (0.0, 0.0, 0.0);
  Point end (0.0, 0.0, 50.0);
  Point start_dir (0.0, 0.0, 1.0);
  Point mid_dir (0.0, 0.0, 1.0);
  Point end_dir (0.0, 0.0, 1.0);
  guint nsamples = 100;
  RefPtr<Image::Interp> warp;

  std::vector<OptBase> opt = get_options (0); // extent
  if (opt.size()) {
    std::vector<float> v = parse_floats (opt[0][0].get_string());
    if (v.size() != 3) throw Exception ("start point should be specified as a comma-separated vector of 3 coordinates");
    start.set (v[0], v[1], v[2]);
    v = parse_floats (opt[0][1].get_string());
    if (v.size() != 3) throw Exception ("end point should be specified as a comma-separated vector of 3 coordinates");
    end.set (v[0], v[1], v[2]);

    start_dir = end - start;
    start_dir.normalise();
    mid_dir = end_dir = start_dir;
    mid = 0.5 * (start + end);
  }
 
  opt = get_options (1); // direction
  if (opt.size()) {
    std::vector<float> v = parse_floats (opt[0][0].get_string());
    if (v.size() != 3) throw Exception ("start direction should be specified as a comma-separated vector of 3 coordinates");
    start_dir.set (v[0], v[1], v[2]);
    v = parse_floats (opt[0][1].get_string());
    if (v.size() != 3) throw Exception ("end direction should be specified as a comma-separated vector of 3 coordinates");
    end_dir.set (v[0], v[1], v[2]);
  }

  opt = get_options (2); // warp
  if (opt.size()) {
    warp = new Image::Interp (*opt[0][0].get_image());
    if (warp->ndim() < 4) throw Exception ("warp image should contain at least 4 dimensions");
    if (warp->dim(3) < 3) throw Exception ("4th dimension of warp image should have length 3");
  }

  opt = get_options (3); // num
  if (opt.size()) nsamples = opt[0][0].get_int();

  Tractography::Writer writer;
  writer.create (argument[1].get_string(), properties);

  std::vector<Point> tck;
  guint skipped = 0;

  ProgressBar::init (0, "resampling tracks...");

  while (file.next (tck)) {

    guint idx_start, idx_end;
    if (get_limits (tck, warp, start, mid, end, start_dir, mid_dir, end_dir, idx_start, idx_end)) {
      std::reverse (tck.begin(), tck.end());
      if (get_limits (tck, warp, start, mid, end, start_dir, mid_dir, end_dir, idx_start, idx_end)) {
        skipped++;
        continue;
      }
    }

    std::vector<Point> tck_out;
    for (guint n = 0; n < nsamples; n++) {
      float f = float(n) / float (nsamples-1);
      Point pos = (1.0-f) * start + f * end;
      Point dir = (1.0-f) * start_dir + f * end_dir;

      for (guint i = idx_start; i <= idx_end; i++) {
        float d = dir.dot (get_pos (tck[i], warp) - pos);
        if (d > 0.0) {
          float f = d / (d - dir.dot (get_pos (tck[i-1], warp) - pos));
          tck_out.push_back (f*tck[i-1] + (1.0-f)*tck[i]);
          break;
        }
      }
    }

    writer.append (tck_out);
    writer.total_count++;
    
    ProgressBar::inc();
  }

  ProgressBar::done();

  writer.close();
  print ("skipped " + str(skipped) + " tracks that did not fit within range\n");
}


