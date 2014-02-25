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

#include "math/fft.h"
#include "image/fft.h"
#include "image/position.h"

namespace MR {
  namespace Image {


    namespace {
      inline bool next (Position& pos, const int limits[MRTRIX_MAX_NDIMS])
      {
        int axis = 0;
        do {
          pos.inc (axis);
          if (pos[axis] < limits[axis]) return (true);
          pos.set (axis, 0);
          axis++;
        } while (axis < pos.ndim());
        return (false);
      }
    }




    void FFT::fft (Position& dest, Position& source, int axis, bool inverse, bool shift)
    {
      int shift_dist = (source.dim(axis)+1)/2;
      int shift_up   = source.dim(axis)/2;

      std::vector<Math::ComplexNumber<double> > array (source.dim (axis));

      guint count = 1;
      int limits[MRTRIX_MAX_NDIMS];
      for (int n = 0; n < source.ndim(); n++) {
        if (n == axis) limits[n] = 1;
        else { limits[n] = source.dim(n); count *= limits[n]; }
      }

      ProgressBar::init (count, String ("performing ") + ( shift ? "shifted " : "" ) +  ( inverse ? "inverse " : "" ) 
          + "FFT along axis " + str (axis) +"...");

      do {
        for (int n = 0; n < source.dim(axis); n++) {
          if (shift && inverse) source.set (axis, ( n >= shift_dist ? n - shift_dist : n + shift_up ));
          else source.set (axis, n);
          array[n].re() = source.re();
          array[n].im() = source.im();
        }

        ft.fft (array, inverse);

        for (int n = 0; n < source.dim(axis); n++) {
          if (shift && !inverse) dest.set (axis, ( n >= shift_dist ? n - shift_dist : n + shift_up ));
          else dest.set (axis, n);
          if (dest.is_complex()) {
            dest.re(array[n].re());
            dest.im(array[n].im());
          }
          else dest.value (sqrt(array[n].re()*array[n].re() + array[n].im()*array[n].im()));
        }

        ProgressBar::inc();  
      } while (next (source, limits));

      ProgressBar::done();
    }

  }
}


