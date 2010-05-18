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

#ifndef __image_header_h__
#define __image_header_h__

#include "data_type.h"
#include "image/header.h"
#include "image/axis.h"
#include "math/matrix.h"

namespace MR {
  namespace Image {

    class Object;

    class Header {
      public:
        Header ();
        Header (const Object& H);

        void                  reset ();

        Axes                  axes;

        std::vector<String>   comments;
        DataType              data_type;
        Math::Matrix          DW_scheme;
        float                 offset, scale;
        String                name;
        
        bool                  read_only;
        const gchar*          format;

        int                   ndim () const     { return (axes.ndim()); }
        int                   dim (int axis) const     { return (axes.dim[axis]); }
        float                 vox (int axis) const     { return (axes.vox[axis]); }

        const Math::Matrix&   transform () const { return (trans_I2R); }
        void                  set_transform (const Math::Matrix& M);
        const Math::Matrix&   I2R () const       { return (trans_I2R); }
        const Math::Matrix&   R2I () const       { return (trans_R2I); }
        const Math::Matrix&   P2R () const       { return (trans_P2R); }
        const Math::Matrix&   R2P () const       { return (trans_R2P); }

        gsize                memory_footprint (guint up_to_dim = MRTRIX_MAX_NDIMS) const;
        gsize                memory_footprint (const gchar* axes_spec) const;

        gsize                voxel_count (int up_to_dim = MRTRIX_MAX_NDIMS) const;
        gsize                voxel_count (const gchar* axes_spec) const;

        String               description () const;
        void                 sanitise_transform ();

        void                 merge (const Header& H);

      protected:
        Math::Matrix         trans_I2R, trans_R2I, trans_P2R, trans_R2P;
    };

    std::ostream& operator<< (std::ostream& stream, const Header& H);













    inline Header::Header() :
      offset (0.0),
      scale (1.0),
      read_only (true),
      format (NULL)
    {
    }


    inline gsize Header::voxel_count (int up_to_dim) const
    { 
      if (up_to_dim > axes.ndim()) up_to_dim = axes.ndim();
      gsize fp = 1;
      for (int n = 0; n < up_to_dim; n++) fp *= axes.dim[n];
      return (fp); 
    }


    inline gsize Header::voxel_count (const gchar* axes_spec) const
    { 
      gsize fp = 1;
      for (int n = 0; n < axes.ndim() && axes_spec[n]; n++) 
        if (axes_spec[n] != '0') fp *= axes.dim[n];
      return (fp); 
    }


    inline gsize Header::memory_footprint (guint up_to_dim) const
    { 
      if (data_type.bits() < 8) return ((voxel_count (up_to_dim)+7)/8);
      else return (data_type.bytes() * voxel_count (up_to_dim));
    }


    inline gsize Header::memory_footprint (const gchar* axes_spec) const
    { 
      if (data_type.bits() < 8) return ((voxel_count (axes_spec)+7)/8);
      else return (data_type.bytes() * voxel_count (axes_spec));
    }



  }
}

#endif

