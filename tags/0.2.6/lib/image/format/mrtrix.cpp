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
    * add capacity to create temporary files for use through pipes

*/

#include <unistd.h>
#include <fcntl.h>
#include <glib/gstdio.h>
#include <glibmm/stringutils.h>

#include <fstream>
#include "image/header.h"
#include "image/mapper.h"
#include "file/key_value.h"
#include "image/format/list.h"
#include "image/name_parser.h"

namespace MR {
  namespace Image {
    namespace Format {

      namespace {

        const gchar* FormatMRtrix = "MRtrix";

      }



      // extensions are: 
      // mih: MRtrix Image Header
      // mif: MRtrix Image File

      bool MRtrix::read (Mapper& dmap, Header& H) const
      { 
        if (!Glib::str_has_suffix (H.name, ".mih") && !Glib::str_has_suffix (H.name, ".mif")) return (false);

        File::KeyValue kv (H.name, "mrtrix image");

        H.format = FormatMRtrix;

        String dtype, layout, file;
        std::vector<int> dim;
        std::vector<float> transform, dw_scheme, vox, scaling;
        std::vector<String> units, labels;

        while (kv.next()) {
          String key = lowercase (kv.key());
          if (key == "dim") dim = parse_ints (kv.value());
          else if (key == "vox") vox = parse_floats (kv.value());
          else if (key == "layout") layout = kv.value();
          else if (key == "datatype") dtype = kv.value();
          else if (key == "file") file = kv.value();
          else if (key == "scaling") scaling = parse_floats (kv.value());
          else if (key == "comments") H.comments.push_back (kv.value());
          else if (key == "units") units = split (kv.value(), "\\");
          else if (key == "labels") labels = split (kv.value(), "\\");
          else if (key == "transform") { std::vector<float> V (parse_floats (kv.value())); transform.insert (transform.end(), V.begin(), V.end()); }
          else if (key == "dw_scheme") { std::vector<float> V (parse_floats (kv.value())); dw_scheme.insert (dw_scheme.end(), V.begin(), V.end()); }
          else error ("WARNING: invalid key \"" + kv.key() + " in generic image header \"" + H.name + "\" - ignored");
        }

        if (dim.empty()) throw Exception ("missing \"dim\" specification for generic image \"" + H.name + "\"");
        H.axes.set_ndim (dim.size());
        for (guint n = 0; n < dim.size(); n++) {
          if (dim[n] < 1) throw Exception ("invalid dimensions for generic image \"" + H.name + "\"");
          H.axes.dim[n] = dim[n];
        }

        if (vox.empty()) throw Exception ("missing \"vox\" specification for generic image \"" + H.name + "\"");
        for (int n = 0; n < H.axes.ndim(); n++) {
          if (vox[n] < 0.0) throw Exception ("invalid voxel size for generic image \"" + H.name + "\"");
          H.axes.vox[n] = vox[n];
        }


        if (dtype.empty()) throw Exception ("missing \"datatype\" specification for generic image \"" + H.name + "\"");
        H.data_type.parse (dtype);


        if (layout.empty()) throw Exception ("missing \"layout\" specification for generic image \"" + H.name + "\"");
        std::vector<Axis> ax = parse_axes_specifier (H.axes, layout);
        if (ax.size() != (guint) H.axes.ndim()) 
          throw Exception ("specified layout does not match image dimensions for generic image \"" + H.name + "\"");

        for (guint i = 0; i < ax.size(); i++) {
          H.axes.axis[i] = ax[i].axis;
          H.axes.forward[i] = ax[i].forward;
        }

        for (guint n = 0; n < MIN ((guint) H.axes.ndim(), labels.size()); n++) H.axes.desc[n] = labels[n];
        for (guint n = 0; n < MIN ((guint) H.axes.ndim(), units.size()); n++) H.axes.units[n] = units[n];

        if (transform.size()) {
          if (transform.size() < 9) throw Exception ("invalid \"transform\" specification for generic image \"" + H.name + "\"");
          Math::Matrix T (4,4);
          int count = 0;
          for (int row = 0; row < 3; row++) 
            for (int col = 0; col < 4; col++) 
              T(row,col) = transform[count++];
          T(3,0) = T(3,1) = T(3,2) = 0.0; T(3,3) = 1.0;
          H.set_transform (T);
        }


        if (dw_scheme.size()) {
          if (dw_scheme.size() % 4) info ("WARNING: invalid \"dw_scheme\" specification for generic image \"" + H.name + "\" - ignored");
          else {
            Math::Matrix M (dw_scheme.size()/4, 4);
            int count = 0;
            for (guint row = 0; row < M.rows(); row++) 
              for (guint col = 0; col < 4; col++) 
                M(row,col) = dw_scheme[count++];
            H.DW_scheme = M;
          }
        }


        if (scaling.size()) {
          if (scaling.size() != 2) throw Exception ("invalid \"scaling\" specification for generic image \"" + H.name + "\"");
          H.offset = scaling[0];
          H.scale = scaling[1];
        }

        if (file.empty()) throw Exception ("missing \"file\" specification for generic image \"" + H.name + "\"");
        std::istringstream files_stream (file);
        String fname;
        files_stream >> fname;
        guint offset = 0;
        if (files_stream.good()) {
          try { files_stream >> offset; }
          catch (...) { throw Exception ("invalid offset specified for file \"" + fname + "\" in generic image header \"" + H.name + "\""); }
        }

        if (fname == ".") {
          if (offset == 0) throw Exception ("invalid offset specified for embedded generic image \"" + H.name + "\""); 
          fname = H.name;
        }
        else fname = Glib::build_filename (Glib::path_get_dirname (H.name), fname);

        ParsedNameList list;
        std::vector<int> num = list.parse_scan_check (fname);

        for (guint n = 0; n < list.size(); n++) 
          dmap.add (list[n]->name(), offset);

        return (true);
      }





      bool MRtrix::check (Header& H, int num_axes) const
      {
        if (H.name.size() && !Glib::str_has_suffix (H.name, ".mih") && !Glib::str_has_suffix (H.name, ".mif")) return (false);

        H.format = FormatMRtrix;

        H.axes.set_ndim (num_axes);
        for (int i = 0; i < H.axes.ndim(); i++) 
          if (H.axes.dim[i] < 1) H.axes.dim[i] = 1;

        return (true);
      }




      void MRtrix::create (Mapper& dmap, const Header& H) const
      {
        if (!is_temporary (H.name) && Glib::file_test (H.name, Glib::FILE_TEST_IS_REGULAR)) 
          throw Exception ("cannot create generic image file \"" + H.name + "\": file exists");

        std::ofstream out (H.name.c_str(), std::ios::out | std::ios::binary);
        if (!out) throw Exception ("error creating file \"" + H.name + "\":" + Glib::strerror(errno));

        out << "mrtrix image\n";
        out << "dim: " << H.axes.dim[0];
        for (int n = 1; n < H.axes.ndim(); n++) out << "," << H.axes.dim[n];

        out << "\nvox: " << H.axes.vox[0];
        for (int n = 1; n < H.axes.ndim(); n++) out << "," << H.axes.vox[n];

        out << "\nlayout: " << ( H.axes.forward[0] ? "+" : "-" ) << H.axes.axis[0];
        for (int n = 1; n < H.axes.ndim(); n++) out << "," << ( H.axes.forward[n] ? "+" : "-" ) << H.axes.axis[n];

        out << "\ndatatype: " << H.data_type.specifier();

        out << "\nlabels: " << H.axes.desc[0];
        for (int n = 1; n < H.axes.ndim(); n++) out << "\\" << H.axes.desc[n];

        out << "\nunits: " <<  H.axes.units[0];
        for (int n = 1; n < H.axes.ndim(); n++) out << "\\" << H.axes.units[n];

        for (std::vector<String>::const_iterator i = H.comments.begin(); i != H.comments.end(); i++) 
          out << "\ncomments: " << *i;


        if (H.transform().is_valid()) {
          out << "\ntransform: " << H.transform() (0,0) << "," <<  H.transform() (0,1) << "," << H.transform() (0,2) << "," << H.transform() (0,3);
          out << "\ntransform: " << H.transform() (1,0) << "," <<  H.transform() (1,1) << "," << H.transform() (1,2) << "," << H.transform() (1,3);
          out << "\ntransform: " << H.transform() (2,0) << "," <<  H.transform() (2,1) << "," << H.transform() (2,2) << "," << H.transform() (2,3);
        }

        if (H.offset != 0.0 || H.scale != 1.0) 
          out << "\nscaling: " << H.offset << "," << H.scale;

        if (H.DW_scheme.is_valid()) {
          for (guint i = 0; i < H.DW_scheme.rows(); i++)
            out << "\ndw_scheme: " << H.DW_scheme (i,0) << "," << H.DW_scheme (i,1) << "," << H.DW_scheme (i,2) << "," << H.DW_scheme (i,3);
        }

        bool single_file = Glib::str_has_suffix (H.name, ".mif");

        size_t offset = 0;
        out << "\nfile: ";
        if (single_file) {
          offset = out.tellp();
          offset += 14;
          out << ". " << offset << "\nEND\n";
        }
        else out << Glib::path_get_basename (H.name.substr (0, H.name.size()-4) + ".dat") << "\n";

        out.close();

        if (single_file) {
          int fd = g_open (H.name.c_str(), O_RDWR, 0755);
          if (fd < 0) throw Exception ("error opening file \"" + H.name + "\" for resizing: " + Glib::strerror(errno));
          int status = ftruncate (fd, offset + H.memory_footprint());
          close (fd);
          if (status) throw Exception ("cannot resize file \"" + H.name + "\": " + Glib::strerror(errno));
          dmap.add (H.name, offset);
        }
        else dmap.add (H.name.substr (0, H.name.size()-4) + ".dat", 0, H.memory_footprint());
      }


    }
  }
}
