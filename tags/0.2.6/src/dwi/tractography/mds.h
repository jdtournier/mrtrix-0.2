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

#ifndef __dwi_tractography_mds_h__
#define __dwi_tractography_mds_h__

#include "point.h"
#include "get_set.h"
#include "file/mmap.h"
#include "dwi/tractography/mds_tags.h"
#include "dwi/tractography/mds.h"
#include "dwi/tractography/properties.h"

#define MDS_SIZE_INC 4096

/*
   MRI format:
   Magic number:              MDS#         (4 bytes)
   Byte order specifier:    guint16 = 1     (2 bytes)
   ...
Elements:
ID:                      guint32       (4 bytes)
size:                    guint32       (4 bytes)
contents:              unspecified  ('size' bytes)
...

Datatype:      ( ID & 0x000000FFU ) 
*/

namespace MR {
  namespace DWI {
    namespace Tractography {

      class MDS {
        public:
          class Track {
            public:
              MDS*       parent;
              gsize      offset;
              guint      count;
              bool       is_BE;
              std::vector<Point> next () const
              {
                std::vector<Point> V;
                float32* data = (float32*) ((guint8*) parent->get_mmap().address() + offset);
                for (guint n = 0; n < count; n++) 
                  V.push_back (Point (get<float32> (data+3*n, is_BE), get<float32> (data+3*n+1, is_BE), get<float32> (data+3*n+2, is_BE) ));
                return (V);
              }
          };



          MDS () { current_offset = next = 0; }

          void         read (const String& filename, Properties& properties);
          std::vector<Track>  tracks;
          bool         changed () const     { return (mmap.changed()); }

          void         create (const String& filename, const Properties& properties);
          void         append (const std::vector<Point>& points);
          void         finalize ();

          String       name () const    { return (mmap.name()); }

        private:
          MR::File::MMap  mmap;
          gsize        current_offset, next;
          bool         is_BE;
          std::vector<Tag>   stack;

          guint32       read_uint32 (int idx) const;

        protected:

          ROI::Type type;
          guint8 shape;
          Point sphere_pos;
          float sphere_rad;
          String mask_name;
          Properties* prop;


          bool         BE () const { return (is_BE); }
          void         interpret ();

          const MR::File::MMap&  get_mmap () const { return (mmap); }

          guint32      size () const { return (read_uint32 (1)); }
          Tag          tag () const { Tag T (read_uint32(0)); return (T); }
          guint32      count () const;
          gsize        offset (guint index = 0) const;
          const guint8* data (guint index = 0) const;
          guint8*       data (guint index = 0);

          String       get_string () const;
          gint8        get_int8 (guint index = 0) const;
          guint8       get_uint8 (guint index = 0) const;
          gint16       get_int16 (guint index = 0) const;
          guint16      get_uint16 (guint index = 0) const;
          gint32       get_int32 (guint index = 0) const;
          guint32      get_uint32 (guint index = 0) const;
          float32      get_float32 (guint index = 0) const;
          float64      get_float64 (guint index = 0) const;

          void         append (Tag tag_id, guint32 nbytes = 0);

          void         append_string (const Tag& T, const String& str);
          void         append_int8 (const Tag& T, gint8 val);
          void         append_uint8 (const Tag& T, guint8 val);
          void         append_int32 (const Tag& T, gint32 val);
          void         append_uint32 (const Tag& T, guint32 val);
          void         append_float32 (const Tag& T, float32 val);
          void         append_float32 (const Tag& T, const float32* val, int num);
          void         append_float64 (const Tag& T, const float64* val, int num);

          const std::vector<Tag>& containers () const;
      };












      /************************************************************************
       *            MDS inline function definitions             *
       ************************************************************************/



      inline guint32 MDS::read_uint32 (int idx) const 
      { return (get<guint32> ((guint8*) mmap.address() + current_offset + idx*sizeof (guint32), is_BE)); }





      inline guint32 MDS::count () const
      {
        if (size() == 0) return (0);
        if (tag().type() == DataType::Text) return (1);
        assert (tag().type() != DataType::Bit);
        return (size() / tag().type().bytes());
      }

      inline gsize MDS::offset (guint index) const
      {
        if (index == 0) return (current_offset + 2*sizeof (guint32));
        assert (tag().type() != DataType::Bit && tag().type() != DataType::Text);
        return (current_offset + 2*sizeof (guint32) + index*tag().type().bytes());
      }

      inline const guint8* MDS::data (guint index) const
      { 
        if (index == 0) return ((guint8*) mmap.address() + current_offset + 2*sizeof (guint32));
        assert (tag().type() != DataType::Bit && tag().type() != DataType::Text);
        return ((guint8*) mmap.address() + current_offset + 2*sizeof (guint32) + index*tag().type().bytes());
      }

      inline guint8* MDS::data (guint index)
      {
        if (index == 0) return ((guint8*) mmap.address() + current_offset + 2*sizeof (guint32));
        assert (tag().type() != DataType::Bit && tag().type() != DataType::Text);
        return ((guint8*) mmap.address() + current_offset + 2*sizeof (guint32) + index*tag().type().bytes());
      }

      inline const std::vector<Tag>& MDS::containers () const   { return (stack); }





      inline String MDS::get_string () const { return (String ((const gchar*) data(), size())); }
      inline gint8 MDS::get_int8 (guint index) const { return (*((gint8*) data (index))); }
      inline guint8 MDS::get_uint8 (guint index) const { return (*((guint8*) data (index))); }
      inline gint16 MDS::get_int16 (guint index) const { return (get<gint16> (data(index), is_BE)); }
      inline guint16 MDS::get_uint16 (guint index) const { return (get<guint16> (data(index), is_BE)); }
      inline gint32 MDS::get_int32 (guint index) const { return (get<gint32> (data(index), is_BE)); }
      inline guint32 MDS::get_uint32 (guint index) const { return (get<guint32> (data(index), is_BE)); }
      inline float32 MDS::get_float32 (guint index) const { return (get<float32> (data(index), is_BE)); }
      inline float64 MDS::get_float64 (guint index) const { return (get<float64> (data(index), is_BE)); }




      inline void MDS::append_string (const Tag& T, const String& str)
      {
        guint s = strlen (str.c_str ()); 
        append (T, s);
        memcpy (data(), str.c_str (), s); 
      }



      inline void MDS::append_int8 (const Tag& T, gint8 val)
      {
        append (T, sizeof (gint8));
        *((gint8*) data()) = val;
      }




      inline void MDS::append_uint8 (const Tag& T, guint8 val)
      {
        append (T, sizeof (guint8));
        *((guint8*) data()) = val;
      }




      inline void MDS::append_int32 (const Tag& T, gint32 val)
      {
        append (T, sizeof (gint32));
        put<gint32> (val, data(), is_BE);
      }




      inline void MDS::append_uint32 (const Tag& T, guint32 val)
      {
        append (T, sizeof (guint32));
        put<guint32> (val, data(), is_BE);
      }




      inline void MDS::append_float32 (const Tag& T, float32 val)
      {
        append (T, sizeof (float32));
        put<float32> (val, data(), is_BE);
      }




      inline void MDS::append_float32 (const Tag& T, const float32* val, int num)
      {
        append (T, num * sizeof (float32));
        for (int n = 0; n < num; n++) put<float32> (val[n], data(n), is_BE);
      }





      inline void MDS::append_float64 (const Tag& T, const float64* val, int num)
      {
        append (T, num * sizeof (float64));
        for (int n = 0; n < num; n++) put<float64> (val[n], data(n), is_BE);
      }



    }
  }
}

#endif



