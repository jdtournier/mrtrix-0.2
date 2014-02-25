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

#include <fstream>
#include <glibmm/fileutils.h>
#include <glibmm/stringutils.h>

#include "image/mapper.h"
#include "image/format/list.h"
#include "image/name_parser.h"
#include "get_set.h"

namespace MR {
  namespace Image {
    namespace Format {

      namespace {

        const gchar* FormatBFloat = "XDS (floating point)";
        const gchar* FormatBShort = "XDS (integer)";

      }






      bool XDS::read (Mapper& dmap, Header& H) const
      {
        if (!Glib::str_has_suffix (H.name, ".bfloat") && !Glib::str_has_suffix (H.name, ".bshort")) return (false);

        H.axes.set_ndim (4);
        int BE;

        String name (H.name);
        name.replace (name.size()-6, 6, "hdr");

        std::ifstream in (name.c_str());
        if (!in) throw Exception ("error reading header file \"" + name + "\": " + Glib::strerror (errno));
        in >> H.axes.dim[1] >> H.axes.dim[0] >> H.axes.dim[3] >> BE;
        in.close();

        if (Glib::str_has_suffix (H.name, ".bfloat")) {
          H.data_type = DataType::Float32;
          H.format = FormatBFloat;
        }
        else {
          H.data_type = DataType::UInt16;
          H.format = FormatBShort;
        }

        if (BE) H.data_type.set_flag (DataType::LittleEndian);
        else H.data_type.set_flag (DataType::BigEndian);

        H.axes.dim[2] = 1;

        H.axes.vox[0] = H.axes.vox[1] = 3.0;
        H.axes.vox[2] = 10.0;
        H.axes.vox[3] = 1.0;

        H.axes.axis[0] = 0;                H.axes.forward[0] = false;
        H.axes.axis[1] = 1;                H.axes.forward[1] = false;
        H.axes.axis[2] = Axis::undefined;  H.axes.forward[2] = true;
        H.axes.axis[3] = 2;                H.axes.forward[3] = true;

        H.axes.desc[0] = Axis::left_to_right;
        H.axes.desc[1] = Axis::posterior_to_anterior;
        H.axes.desc[2] = Axis::inferior_to_superior;
        H.axes.desc[3] = Axis::time;

        H.axes.units[0] = Axis::millimeters; 
        H.axes.units[1] = Axis::millimeters; 
        H.axes.units[2] = Axis::millimeters;
        H.axes.units[3] = Axis::milliseconds;

        dmap.add (H.name, 0);

        return (true);
      }







      bool XDS::check (Header& H, int num_axes) const
      {
        if (!Glib::str_has_suffix (H.name, ".bfloat") && !Glib::str_has_suffix (H.name, ".bshort")) return (false);
        if (num_axes > 4) throw Exception ("cannot create XDS image with more than 4 dimensions");
        if (num_axes == 4 && H.axes.dim[2] > 1) throw Exception ("cannot create multi-slice XDS image with a single file");
        if (num_axes < 2) throw Exception ("cannot create XDS image with less than 2 dimensions");

        H.axes.set_ndim (4);

        H.axes.dim[2] = 1;
        for (guint n = 0; n < 4; n++)
          if (H.axes.dim[n] < 1) H.axes.dim[n] = 1;

        H.axes.vox[0] = H.axes.vox[1] = 3.0;
        H.axes.vox[2] = 10.0;
        H.axes.vox[3] = 1.0;

        H.axes.axis[0] = 0;                H.axes.forward[0] = false;
        H.axes.axis[1] = 1;                H.axes.forward[1] = false;
        H.axes.axis[2] = Axis::undefined;  H.axes.forward[2] = true;
        H.axes.axis[3] = 2;                H.axes.forward[3] = true;

        H.axes.desc[0] = Axis::left_to_right;
        H.axes.desc[1] = Axis::posterior_to_anterior;
        H.axes.desc[2] = Axis::inferior_to_superior;
        H.axes.desc[3] = Axis::time;

        H.axes.units[0] = Axis::millimeters; 
        H.axes.units[1] = Axis::millimeters; 
        H.axes.units[2] = Axis::millimeters;
        H.axes.units[3] = Axis::milliseconds;

        bool is_BE = H.data_type.is_big_endian();

        if (Glib::str_has_suffix (H.name, ".bfloat")) {
          H.data_type = DataType::Float32;
          H.format = FormatBFloat;
        }
        else {
          H.data_type = DataType::UInt16;
          H.format = FormatBShort;
        }

        if (is_BE) H.data_type.set_flag (DataType::BigEndian);
        else H.data_type.set_flag (DataType::LittleEndian);

        return (true);
      }




      void XDS::create (Mapper& dmap, const Header& H) const
      {
        guint msize = H.memory_footprint ("1101");

        String header_name (H.name);
        header_name.replace (header_name.size()-6, 6, "hdr");

        std::ofstream out(header_name.c_str());
        if (!out) throw Exception ("error writing header file \"" + header_name + "\": " + Glib::strerror (errno));

        out << H.axes.dim[1] << " " << H.axes.dim[0] << " " << H.axes.dim[3] 
            << " " << ( H.data_type.is_little_endian() ? 1 : 0 ) << "\n";
        out.close();

        dmap.add (H.name, 0, msize);
      }

    } 
  }  
}

