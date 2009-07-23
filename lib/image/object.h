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

#ifndef __image_object_h__
#define __image_object_h__

#include "ptr.h"
#include "image/header.h"
#include "image/mapper.h"
#include "image/format/base.h"
#include "math/complex_number.h"

namespace MR {
  namespace Dialog { class File; }
  namespace Image {

    class Interp;
    class Position;

    /*! \defgroup Image Image access
     * \brief Classes and functions providing access to image data. */

    //! \addtogroup Image 
    // @{
    
    class Object {
      public:
        Object () : start (0) { memset (stride, 0, MRTRIX_MAX_NDIMS*sizeof(gssize)); }
        ~Object ()
        { 
          info ("closing image \"" + H.name + "\"...");
          M.unmap (H); 
        }

        const Header&        header () const         { return (H); }

        void                 open (const String& imagename, bool is_read_only = true);
        void                 create (const String& imagename, Header &template_header);
        void                 concatenate (std::vector<RefPtr<Object> >& images);

        void                 map ()                  { if (!is_mapped()) M.map (H); }
        void                 unmap ()                { if (is_mapped()) M.unmap (H); }
        bool                 is_mapped () const      { return (M.is_mapped()); }

        int                  dim (guint index) const { return (H.axes.dim[index]); } 
        int                  ndim () const           { return (H.axes.ndim()); }
        float                vox (guint index) const { return (H.axes.vox[index]); }
        const String&        name () const           { return (H.name); }
        bool                 is_complex () const     { return (H.data_type.is_complex()); }

        const std::vector<String>&  comments () const { return (H.comments); }
        DataType             data_type () const       { return (H.data_type); }
        const Math::Matrix&  DW_scheme () const       { return (H.DW_scheme); }
        gfloat               offset () const          { return (H.offset); }
        gfloat               scale () const           { return (H.scale); }
        
        bool                 output_name () const     { return (M.output_name.size()); }
        void                 no_output_name ()        { M.output_name.clear(); }

        void                 set_temporary (bool yesno = true) {
          M.temporary = yesno; 
          if (M.temporary) { for (guint n = 0; n < M.list.size(); n++) M.list[n].fmap.mark_for_deletion(); }
        }

        bool                 read_only () const   { return (H.read_only); }
        void                 set_read_only (bool read_only) { M.set_read_only (read_only); H.read_only = read_only; }

        const gchar*         format () const      { return (H.format); }


        const Axes&          axes () const;

        gsize                memory_footprint (guint up_to_dim = MRTRIX_MAX_NDIMS) const { return (H.memory_footprint (up_to_dim)); }
        gsize                memory_footprint (const gchar* axes_spec) const { return (H.memory_footprint (axes_spec)); }

        gsize                voxel_count (guint up_to_dim = MRTRIX_MAX_NDIMS) const  { return (H.voxel_count (up_to_dim)); }
        gsize                voxel_count (const gchar* axes_spec) const { return (H.voxel_count (axes_spec)); }

        String               description () const { return (H.description()); }



        const Math::Matrix&  I2R () const            { return (H.I2R()); }
        const Math::Matrix&  R2I () const            { return (H.R2I()); }
        const Math::Matrix&  P2R () const            { return (H.P2R()); }
        const Math::Matrix&  R2P () const            { return (H.R2P()); }

        void                 apply_scaling (float scale, float bias = 0.0) { H.scale *= scale; H.offset = scale * H.offset + bias; }
        void                 set_transform (const Math::Matrix& T) { H.set_transform (T); }

        void                 optimise () { M.optimised = true; }

        friend std::ostream& operator<< (std::ostream& stream, const Object& obj);

      protected:
        static const Format::Base* handlers[];

        Header               H;
        Mapper               M;
        gsize                start;
        gssize               stride[MRTRIX_MAX_NDIMS];

        void                 setup ();

        float32              scale_from_storage (float32 val) const { return (H.offset + H.scale * val); }
        float32              scale_to_storage (float32 val) const   { return ((val - H.offset) / H.scale); }

        float                re (gsize offset) const                { return (scale_from_storage (M.re (offset))); }
        void                 re (gsize offset, float val)           { M.re (scale_to_storage (val), offset); }

        float                im (gsize offset) const                { return (scale_from_storage (M.im (offset))); }
        void                 im (gsize offset, float val)           { M.im (scale_to_storage (val), offset); }

        friend class Interp;
        friend class Dialog::File;
        friend class Position;
        friend class Header;
    };

    //! @}

    inline Header::Header (const Object& H) { *this = H.H; }

  }
}

#endif

