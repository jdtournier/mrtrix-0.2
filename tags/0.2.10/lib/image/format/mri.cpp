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


    02-09-2008 J-Donald Tournier <d.tournier@brain.org.au>
    * remove temporary file creation capacity

*/

#include <glibmm/stringutils.h>

#include "file/config.h"
#include "image/mapper.h"
#include "file/mmap.h"
#include "image/format/list.h"
#include "image/header.h"
#include "get_set.h"

/*
 MRI format:
     Magic number:              MRI#         (4 bytes)
     Byte order specifier:    guint16 = 1     (2 bytes)
     ...
     Elements:
       ID specifier:            guint32       (4 bytes)
       size:                    guint32       (4 bytes)
       contents:              unspecified  ('size' bytes)
     ...

*/

#define MRI_DATA       0x01
#define MRI_DIMENSIONS 0x02
#define MRI_ORDER      0x03
#define MRI_VOXELSIZE  0x04
#define MRI_COMMENT    0x05
#define MRI_TRANSFORM  0x06
#define MRI_DWSCHEME   0x07

namespace MR {
  namespace Image {
    namespace Format {

      namespace {

        const gchar* FormatMRI = "MRTools (legacy format)";

        inline guint char2order (gchar item, bool& forward)
        {
          switch (item) {
            case 'L': forward = true;  return (0);
            case 'R': forward = false; return (0);
            case 'P': forward = true;  return (1);
            case 'A': forward = false; return (1);
            case 'I': forward = true;  return (2);
            case 'S': forward = false; return (2);
            case 'B': forward = true;  return (3);
            case 'E': forward = false; return (3);
          }
          return (Axis::undefined);
        }


        inline gchar order2char (guint axis, bool forward)
        {
          switch (axis) {
            case 0: if (forward) return ('L'); else return ('R');
            case 1: if (forward) return ('P'); else return ('A');
            case 2: if (forward) return ('I'); else return ('S');
            case 3: if (forward) return ('B'); else return ('E');
          }
          return ('\0');
        }


        inline guint type (const guint8* pos, bool is_BE) { return (get<guint32> (pos, is_BE)); }
        inline guint size (const guint8* pos, bool is_BE) { return (get<guint32> (pos + sizeof (guint32), is_BE)); }
        inline guint8* data (guint8* pos)                 { return (pos + 2*sizeof (guint32)); }
        inline const guint8* data (const guint8* pos)     { return (pos + 2*sizeof (guint32)); }

        inline guint8* next (guint8* current_pos, bool is_BE) 
        {
          return (current_pos + 2*sizeof (guint32) + size (current_pos, is_BE));
        }
        inline const guint8* next (const guint8* current_pos, bool is_BE) 
        {
          return (current_pos + 2*sizeof (guint32) + size (current_pos, is_BE));
        }
        inline void write (guint8* pos, guint32 Type, guint32 Size, bool is_BE)
        {
          put<guint32> (Type, pos, is_BE);
          put<guint32> (Size, pos + sizeof (guint32), is_BE);
        }

      }








      bool MRI::read (Mapper& dmap, Header& H) const
      {
        if (!Glib::str_has_suffix (H.name, ".mri")) return (false);

        H.format = FormatMRI;

        File::MMap fmap (H.name);
        fmap.map();

        if (memcmp ((guint8*) fmap.address(), "MRI#", 4)) 
          throw Exception ("file \"" + H.name + "\" is not in MRI format (unrecognised magic number)");

        bool is_BE = false;
        if (get<guint16> ((guint8*) fmap.address() + sizeof (gint32), is_BE) == 0x0100U) is_BE = true;
        else if (get<guint16> ((guint8*) fmap.address() + sizeof (guint32), is_BE) != 0x0001U) 
          throw Exception ("MRI file \"" + H.name + "\" is badly formed (invalid byte order specifier)");

        H.axes.set_ndim (4);

        guint data_offset = 0;
        gchar* c;
        Math::Matrix M (4,4);
        const guint8* current = (guint8*) fmap.address() + sizeof (gint32) + sizeof (guint16);
        const guint8* last = (guint8*) fmap.address() + fmap.size() - 2*sizeof (guint32);

        while (current <= last) {
          switch (type (current, is_BE)) {
            case MRI_DATA:
              H.data_type = (((const gchar*) data (current)) - 4)[0];
              data_offset = current + 5 - (guint8*) fmap.address();
              fmap.unmap();
              break;
            case MRI_DIMENSIONS:
              H.axes.dim[0] = get<guint32> (data (current), is_BE);
              H.axes.dim[1] = get<guint32> (data (current) + sizeof (guint32), is_BE);
              H.axes.dim[2] = get<guint32> (data (current) + 2*sizeof (guint32), is_BE);
              H.axes.dim[3] = get<guint32> (data (current) + 3*sizeof (guint32), is_BE);
              break;
            case MRI_ORDER:
              c = (gchar*) data (current);
              for (guint n = 0; n < 4; n++) {
                bool forward = true;
                guint ax = char2order (c[n], forward);
                H.axes.axis[ax] = n;
                H.axes.forward[ax] = forward;
              }
              break;
            case MRI_VOXELSIZE:
              H.axes.vox[0] = get<float32> (data (current), is_BE);
              H.axes.vox[1] = get<float32> (data (current) + sizeof (float32), is_BE);
              H.axes.vox[2] = get<float32> (data (current) + 2*sizeof (float32), is_BE);
              break;
            case MRI_COMMENT:
              H.comments.push_back (String ((const gchar*) data (current), size (current, is_BE)));
              break;
            case MRI_TRANSFORM:
              for (guint i = 0; i < 4; i++)
                for (guint j = 0; j < 4; j++)
                  M(i,j) = get<float32> (data (current) + ( i*4 + j )*sizeof (float32), is_BE);
              H.set_transform (M);
              break;
            case MRI_DWSCHEME:
              H.DW_scheme.allocate (size (current, is_BE)/(4*sizeof (float32)), 4);
              for (guint i = 0; i < H.DW_scheme.rows(); i++)
                for (guint j = 0; j < 4; j++)
                  H.DW_scheme(i,j) = get<float32> (data (current) + ( i*4 + j )*sizeof (float32), is_BE);
              break;
            default:
              error ("unknown header entity (" + str (type (current, is_BE)) 
                  + ", offset " + str (current - (guint8*) fmap.address()) 
                  + ") in image \"" + H.name + "\" - ignored");
              break;
          }
          if (data_offset) break;
          current = next (current, is_BE);
        }


        if (!data_offset) throw Exception ("no data field found in MRI image \"" + H.name + "\"");


        if (!H.axes.desc[0].size()) H.axes.desc[0] = Axis::left_to_right;
        if (!H.axes.units[0].size()) H.axes.units[0] = Axis::millimeters;
        if (H.axes.ndim() > 1) {
          if (!H.axes.desc[1].size()) H.axes.desc[1] = Axis::posterior_to_anterior;
          if (!H.axes.units[1].size()) H.axes.units[1] = Axis::millimeters;
          if (H.axes.ndim() > 2) {
            if (!H.axes.desc[2].size()) H.axes.desc[2] = Axis::inferior_to_superior;
            if (!H.axes.units[2].size()) H.axes.units[2] = Axis::millimeters;
          }
        }

        dmap.add (fmap, data_offset);

        return (true);
      }





      bool MRI::check (Header& H, int num_axes) const
      {
        if (!Glib::str_has_suffix (H.name, ".mri")) return (false);
        if ((int) H.axes.ndim() > num_axes && num_axes != 4) throw Exception ("MRTools format can only support 4 dimensions");

        H.format = FormatMRI;

        H.axes.set_ndim (num_axes);

        if (H.axes.desc[0].empty()) H.axes.desc[0]= Axis::left_to_right;
        if (H.axes.units[0].empty()) H.axes.units[0] = Axis::millimeters;
        if (H.axes.ndim() > 1) {
          if (H.axes.desc[1].empty()) H.axes.desc[1] = Axis::posterior_to_anterior;
          if (H.axes.units[1].empty()) H.axes.units[1] = Axis::millimeters;
          if (H.axes.ndim() > 2) {
            if (H.axes.desc[2].empty()) H.axes.desc[2] = Axis::inferior_to_superior;
            if (H.axes.units[2].empty()) H.axes.units[2] = Axis::millimeters;
          }
        }

        return (true);
      }







      void MRI::create (Mapper& dmap, const Header& H) const
      {
        File::MMap fmap (H.name, 65536, "mri");
        fmap.map();

#if G_BYTE_ORDER == G_BIG_ENDIAN
        bool is_BE = true;
#else
        bool is_BE = false;
#endif

        memcpy ((guint8*) fmap.address(), "MRI#", 4);
        put<guint16> (0x01U, (guint8*) fmap.address() + sizeof (guint32), is_BE);

        guint8* current = (guint8*) fmap.address() + sizeof (guint32) + sizeof (guint16);

        write (current, MRI_DIMENSIONS, 4*sizeof (guint32), is_BE);
        put<guint32> (H.axes.dim[0], data (current), is_BE);
        put<guint32> (( H.axes.ndim() > 1 ? H.axes.dim[1] : 1 ), data (current) + sizeof (guint32), is_BE);
        put<guint32> (( H.axes.ndim() > 2 ? H.axes.dim[2] : 1 ), data (current) + 2*sizeof (guint32), is_BE);
        put<guint32> (( H.axes.ndim() > 3 ? H.axes.dim[3] : 1 ), data (current) + 3*sizeof (guint32), is_BE);

        current = next (current, is_BE);
        write (current, MRI_ORDER, 4*sizeof (guchar), is_BE);
        int n;
        for (n = 0; n < H.axes.ndim(); n++) 
          ((gchar*) data (current))[H.axes.axis[n]] = order2char (n, H.axes.forward[n]);
        for (; n < 4; n++) ((gchar*) data (current))[n] = order2char (n, true);

        current = next (current, is_BE);
        write (current, MRI_VOXELSIZE, 3*sizeof (float32), is_BE);
        put<float32> (H.axes.vox[0], data (current), is_BE);
        put<float32> (( H.axes.ndim() > 1 ? H.axes.vox[1] : 2.0 ), data (current)+sizeof (float32), is_BE);
        put<float32> (( H.axes.ndim() > 2 ? H.axes.vox[2] : 2.0 ), data (current)+2*sizeof (float32), is_BE);

        for (guint n = 0; n < H.comments.size(); n++) {
          guint l = H.comments[n].size();
          if (l) {
            current = next (current, is_BE);
            write (current, MRI_COMMENT, l, is_BE);
            memcpy (data (current), H.comments[n].c_str(), l);
          }
        }

        if (H.transform().is_valid()) {
          current = next (current, is_BE);
          write (current, MRI_TRANSFORM, 16*sizeof (float32), is_BE);
          for (guint i = 0; i < 4; i++)
            for (guint j = 0; j < 4; j++)
              put<float32> (H.transform()(i,j), data (current) + ( i*4 + j )*sizeof (float32), is_BE);
        }

        if (H.DW_scheme.is_valid()) {
          current = next (current, is_BE);
          write (current, MRI_DWSCHEME, 4*H.DW_scheme.rows()*sizeof (float32), is_BE);
          for (guint i = 0; i < H.DW_scheme.rows(); i++)
            for (guint j = 0; j < 4; j++)
              put<float32> (H.DW_scheme(i,j), data (current) + ( i*4 + j )*sizeof (float32), is_BE);
        }

        current = next (current, is_BE);
        write (current, MRI_DATA, 0, is_BE);
        ((gchar*) current)[4] = H.data_type();

        guint data_offset = current + 5 - (guint8*) fmap.address();
        fmap.resize (data_offset + H.memory_footprint());
        dmap.add (fmap, data_offset);
      }

    }
  }
}



