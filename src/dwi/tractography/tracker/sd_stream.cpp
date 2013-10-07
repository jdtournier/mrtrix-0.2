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

    03-03-2010 J-Donald Tournier <d.tournier@brain.org.au>
    * tracking now stops immediately before the track leaves the mask, rather
    * than immediately after.

*/

#include "dwi/tractography/tracker/sd_stream.h"
#include "dwi/SH.h"

namespace MR {
  namespace DWI {
    namespace Tractography {
      namespace Tracker {


        SDStream::SDStream (Image::Object& source_image, Properties& properties) : 
          Base (source_image, properties),
          lmax (SH::LforN (source.dim(3))),
          precomputed (true)
        {
          float min_curv = step_size / ( 2.0 * sin (0.5 * M_PI_2));

          properties["method"] = "SD_STREAM";
          if (props["min_curv"].empty()) props["min_curv"] = str (min_curv); else min_curv = to<float> (props["min_curv"]);
          if (props["max_num_tracks"].empty()) props["max_num_tracks"] = "100";

          if (props["lmax"].empty()) props["lmax"] = str (lmax); else lmax = to<int> (props["lmax"]);
          if (props["sh_precomputed"].empty()) props["sh_precomputed"] = ( precomputed ? "1" : "0" ); else precomputed = to<int> (props["sh_precomputed"]);

          min_dp = cos (curv2angle (step_size, min_curv));
          if (precomputed) SH::precompute (lmax);
        }





        bool SDStream::init_direction (const Point& seed_dir)
        {
          if (get_source_data (pos)) return (true);

          if (!seed_dir) dir.set (rng.normal(), rng.normal(), rng.normal());
          else dir = seed_dir;
          dir.normalise();
          float val = SH::get_peak (&values[0], lmax, dir, precomputed);
          if (gsl_finite (val)) if (val > init_threshold) return (false);
          return (true);
        }




        bool SDStream::next_point ()
        {
          if (get_source_data (pos)) return (true);

          Point prev_dir (dir);
          dir.normalise ();
          float val = SH::get_peak (&values[0], lmax, dir, precomputed);

          if (!gsl_finite (val)) return (true);
          if (val < threshold) return (true);
          if (dir.dot (prev_dir) < min_dp) return (true);

          return (false);
        }



      }
    }
  }
}


