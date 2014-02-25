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

using namespace std; 
using namespace MR; 

SET_VERSION_DEFAULT;

DESCRIPTION = {
 "smooth images using a 3x3x3 median filter.",
  NULL
};

ARGUMENTS = {
  Argument ("input", "input image", "input image to be median-filtered.", true, true).type_image_in (),
  Argument ("output", "output image", "the output image.").type_image_out (),
  Argument::End
};


OPTIONS = { Option::End };

EXECUTE {
  Image::Object& in_obj (*argument[0].get_image());
  in_obj.optimise();

  Image::Position in (in_obj);
  Image::Header header (in.image.header());

  Image::Position out (*argument[1].get_image (header));

  int from[3], to[3], n, n1, n2, nc, i;
  float val, v[14], cm, t;
  bool avg;

  ProgressBar::init (out.voxel_count(), "median filtering...");

  do { 
    for (out.set(2,0); out[2] < out.dim(2); out.inc(2)) {
      from[2] = out[2] > 0 ? out[2]-1 : 0;
      to[2] = out[2] < out.dim(2)-1 ? out[2]+2 : out.dim(2);
      n2 = to[2]-from[2];
      for (out.set(1,0); out[1] < out.dim(1); out.inc(1)) {
        from[1] = out[1] > 0 ? out[1]-1 : 0;
        to[1] = out[1] < out.dim(1)-1 ? out[1]+2 : out.dim(1);
        n1 = n2*(to[1]-from[1]);
        for (out.set(0,0); out[0] < out.dim(0); out.inc(0)) {
          from[0] = out[0] > 0 ? out[0]-1 : 0;
          to[0] = out[0] < out.dim(0)-1 ? out[0]+2 : out.dim(0);
          n = n1*(to[0]-from[0]);
          avg = (n+1)%2;
          n = (n/2)+1;
          nc = 0;
          cm = GSL_NEGINF;

          for (in.set(2,from[2]); in[2] < to[2]; in.inc(2)) {
            for (in.set(1,from[1]); in[1] < to[1]; in.inc(1)) {
              for (in.set(0,from[0]); in[0] < to[0]; in.inc(0)) {
                val = in.value();
                if (nc < n) {
                  v[nc] = val;
                  if (v[nc] > cm) cm = v[nc];
                  nc++;
                }
                else if (val < cm) {
                  for (i = 0; v[i] != cm; i++);
                  v[i] = val;
                  cm = GSL_NEGINF;
                  for (i = 0; i < n; i++)
                    if (v[i] > cm) cm = v[i];
                }
              }
            }
          }

          if (avg) {
            t = cm = GSL_NEGINF;
            for (i = 0; i < n; i++) {
              if (v[i] > cm) {
                t = cm;
                cm = v[i];
              }
              else if (v[i] > t) t = v[i];
            }
            cm = (cm+t)/2.0;
          }

          out.value (cm);
          ProgressBar::inc();
        }
      }
    }
  } while (out++);

  ProgressBar::done();
}

