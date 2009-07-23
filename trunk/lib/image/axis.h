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


    01-10-2008 J-Donald Tournier <d.tournier@brain.org.au>
    * added sanitise() method to correct inconsistent axes ordering

*/

#ifndef __image_axis_h__
#define __image_axis_h__

#include "mrtrix.h"

namespace MR {
  namespace Image {

    class Axis {
      public:
        int   axis;
        bool  forward;

        static const gchar*  left_to_right;
        static const gchar*  posterior_to_anterior;
        static const gchar*  inferior_to_superior;
        static const gchar*  time;
        static const gchar*  real_imag;
        static const gchar*  millimeters;
        static const gchar*  milliseconds;

        static const int     undefined = G_MAXINT;
    };

    class Axes {
      public:
        Axes () : size_p (0) { }

        int             dim[MRTRIX_MAX_NDIMS];
        float           vox[MRTRIX_MAX_NDIMS];
        String          desc[MRTRIX_MAX_NDIMS];
        String          units[MRTRIX_MAX_NDIMS];
        int             axis[MRTRIX_MAX_NDIMS];
        bool            forward[MRTRIX_MAX_NDIMS];

        int             ndim () const { return (size_p); }
        void            set_ndim (int num);
        int             direction (int index) const { return (forward[index] ? 1 : -1); }
        void            sanitise ();

        void            copy (guint dest_axis, const Axes& src, guint src_axis);

      protected:
        int             size_p;

        int             find_free_axis () const;
    };

    std::ostream& operator<< (std::ostream& stream, const Axes& axes);
    std::vector<Axis> parse_axes_specifier (const Axes& original, const String& specifier);
    void check_axes_specifier (const std::vector<Axis>& parsed, int ndims);











    inline void Axes::set_ndim (int new_size)
    {
      for (int n = GSL_MIN (size_p, new_size); n < MRTRIX_MAX_NDIMS; n++) {
        dim[n] = 0;
        vox[n] = GSL_NAN;
        axis[n] = Axis::undefined;
        forward[n] = true;
        desc[n].clear();
        units[n].clear();
      }
      size_p = new_size;
    }



    inline int Axes::find_free_axis () const
    {
      for (int a = 0; a < size_p; a++) {
        int m = 0;
        for (; m < size_p; m++) if (axis[m] == a) break; 
        if (m >= size_p) return (a);
      }
      return (Axis::undefined);
    }





    inline void Axes::sanitise ()
    {
      // remove unset/invalid axes:
      for (int a = 0; a < size_p; a++) 
        if (axis[a] >= size_p) 
          axis[a] = find_free_axis();

      // remove duplicates:
      for (int a = 1; a < size_p; a++) {
        for (int n = 0; n < a; n++) {
          if (axis[a] == axis[n]) { 
            axis[a] = find_free_axis();
              break; 
          }
        }
      }   
            
    }





    inline void Axes::copy (guint dest_axis, const Axes& src, guint src_axis)
    {
      dim[dest_axis] = src.dim[src_axis];
      vox[dest_axis] = src.vox[src_axis];
      desc[dest_axis] = src.desc[src_axis];
      units[dest_axis] = src.units[src_axis];
      axis[dest_axis] = src.axis[src_axis];
      forward[dest_axis] = src.forward[src_axis];
    }

  }
}

#endif


