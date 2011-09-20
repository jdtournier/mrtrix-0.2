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

#ifndef __dwi_tractography_tracker_sd_stream_h__
#define __dwi_tractography_tracker_sd_stream_h__

#include "dwi/tractography/tracker/base.h"
#include "math/linalg.h"
#include "dwi/SH.h"

namespace MR {
  namespace DWI {
    namespace Tractography {
      namespace Tracker {

        class SDStream : public Base {
          public:
            SDStream (Image::Object& source_image, Properties& properties);

          protected:
            int   lmax;
            bool  precomputed;

            virtual bool  init_direction (const Point& seed_dir);
            virtual bool  next_point ();

            bool init_direction (const Point& seed_dir, const float* values)
            {
              if (!seed_dir) dir.set (rng.normal(), rng.normal(), rng.normal());
              else dir = seed_dir;
              dir.normalise();
              float val = SH::get_peak (values, lmax, dir, precomputed);
              if (gsl_finite (val)) if (val > init_threshold) return (false);
              return (true);
            }

            bool next_point (const float* values)
            {
              Point prev_dir (dir);
              dir.normalise ();
              float val = SH::get_peak (values, lmax, dir, precomputed);

              if (!gsl_finite (val)) return (true);
              if (val < threshold) return (true);
              if (dir.dot (prev_dir) < min_dp) return (true);

              return (false);
            }

            float         min_dp;
        };




      }
    }
  }
}

#endif

