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
  "resample tracks so they all have the same number of samples, spaced at constant intervals along the axis specified.",
  NULL
};

ARGUMENTS = {
  Argument ("tracks", "track file", "the input track file.").type_file (),
  Argument ("output", "output file", "the output track file").type_file(),
  Argument::End
};



OPTIONS = {

  Option ("start", "start point", "specify the start of the resampling range (default is 0,0,-50).")
    .append (Argument ("point", "point", "a point lying within the start of the resampling range.").type_sequence_float ()),

  Option ("end", "end point", "specify the end of the resampling range (default is 0,0,50).")
    .append (Argument ("point", "point", "a point lying within the end of the resampling range.").type_sequence_float ()),

  Option ("transform", "transform image", "specify an image containing the transform to the space in which the resampling is to take place.")
    .append (Argument ("transform", "transform image", "the transform image.").type_image_in()),
 
  Option ("num", "number of samples", "specify the number of samples to take along each track.")
    .append (Argument ("nsamples", "number of samples", "the number of samples.").type_integer(1, 100000, 100)),
 
  Option::End 
};


inline Point get_pos (const Point& pos, RefPtr<Image::Interp>& trans)
{
  if (!trans) return (pos);
  trans->R (pos);
  Point ret;
  if (!(*trans)) return (ret);
  trans->set(3,0);
  ret[0] = trans->value();
  trans->inc(3);
  ret[1] = trans->value();
  trans->inc(3);
  ret[2] = trans->value();
  return (ret);
}


EXECUTE {
  Tractography::Properties properties;
  Tractography::Reader file;
  file.open (argument[0].get_string(), properties);
  
  Point start (0.0, 0.0, -50.0);
  Point end (0.0, 0.0, 50.0);
  guint nsamples = 100;
  RefPtr<Image::Interp> trans;

  std::vector<OptBase> opt = get_options (0); // start
  if (opt.size()) {
    std::vector<float> v = parse_floats (opt[0][0].get_string());
    if (v.size() != 3) throw Exception ("midpoint should be specified as a comma-separated vector of 3 coordinates");
    start.set (v[0], v[1], v[2]);
  }
 
  opt = get_options (1); // end
  if (opt.size()) {
    std::vector<float> v = parse_floats (opt[0][0].get_string());
    if (v.size() != 3) throw Exception ("midpoint should be specified as a comma-separated vector of 3 coordinates");
    end.set (v[0], v[1], v[2]);
  }

  opt = get_options (2); // transform
  if (opt.size()) {
    trans = new Image::Interp (*opt[0][0].get_image());
    if (trans->ndim() < 4) throw Exception ("transform image should contain at least 4 dimensions");
    if (trans->dim(3) < 3) throw Exception ("4th dimension of transform image should have length 3");
  }

  opt = get_options (3); // num
  if (opt.size()) nsamples = opt[0][0].get_int();

  Tractography::Writer writer;
  writer.create (argument[1].get_string(), properties);

  std::vector<Point> tck;
  guint skipped = 0;

  Point dir = end - start;
  float extent = dir.norm();
  dir.normalise();

  float inc = extent / (nsamples-1);

  ProgressBar::init (0, "resampling tracks...");

  while (file.next (tck)) {
    float d1 = dir.dot (get_pos (tck.front(), trans) - start);
    float d2 = dir.dot (get_pos (tck.back(), trans) - start);
    if ( !( d1 > 0.0 && d1 < extent ) && 
         !( d2 > 0.0 && d2 < extent ) &&
          (d1 - extent/2.0) * (d2 - extent/2.0) < 0.0) { 
      if (d2 < 0.0) std::reverse (tck.begin(), tck.end());

      guint i = 0;
      std::vector<Point> tck_out;
      for (guint n = 0; n < nsamples; n++) {
        float loc = n*inc;
        while (i < tck.size()) {
          d2 = dir.dot (get_pos (tck[i], trans) - start);
          if (d2 > loc) break;
          i++;
        }
        d1 = dir.dot (get_pos (tck[i-1], trans) - start) - loc;
        d2 -= loc;

        float f = d2 / (d2-d1);
        tck_out.push_back (f*tck[i-1] + (1.0-f)*tck[i]);
      }
        
      writer.append (tck_out);
      writer.total_count++;
    }
    else skipped++;
    ProgressBar::inc();
  }

  ProgressBar::done();

  writer.close();
  print ("skipped " + str(skipped) + " tracks that did not fit within range\n");
}


