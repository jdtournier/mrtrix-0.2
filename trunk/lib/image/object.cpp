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


    29-08-2008 J-Donald Tournier <d.tournier@brain.org.au>
    * in create(), finalise the byte order after the handler's check() method 
      to allow different file formats to override the data type more easily.

    02-09-2008 J-Donald Tournier <d.tournier@brain.org.au>
    * update temporary file handling (i.e. those sent via pipes)
    - switch to MRtrix format as the standard format for temporary files
    - handle any type of image supplied through the standard input

    01-10-2008 J-Donald Tournier <d.tournier@brain.org.au>
    * sanitise axes prior to creating an image 
*/

#include "app.h"
#include "image/object.h"
#include "image/format/list.h"
#include "image/name_parser.h"
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>

namespace MR {
  namespace Image {


    void Object::open (const String& imagename, bool is_read_only)
    {
      M.reset();
      H.read_only = is_read_only;

      if (imagename == "-") getline (std::cin, H.name);
      else H.name = imagename;

      if (H.name.empty()) throw Exception ("no name supplied to open image!");

      info ("opening image \"" + H.name + "\"...");

      ParsedNameList list;
      std::vector<int> num = list.parse_scan_check (H.name);

      try { 
        const Format::Base** handler = handlers;
        std::vector< RefPtr<ParsedName> >::iterator item = list.begin();
        Header header;
        header.name = (*item)->name();

        for (; *handler; handler++) if ((*handler)->read (M, header)) break;
        if (!*handler) throw Exception ("unknown format for image \"" + header.name + "\"");

        String old_name (H.name);

        H = header;
        if (header.name == (*item)->name()) H.name = old_name;

        while (++item != list.end()) {
          header.name = (*item)->name();
          if (!(*handler)->read (M, header)) throw Exception ("image specifier contains mixed format files");
          H.merge (header);
        }
      }
      catch (...) { throw Exception ("error opening image \"" + H.name + "\""); }

      if (num.size()) {
        int a = 0, n = 0;
        for (int i = 0; i < H.axes.ndim(); i++) if (H.axes.axis[i] != Axis::undefined) n++;
        H.axes.set_ndim (n + num.size());

        for (std::vector<int>::const_iterator item = num.begin(); item != num.end(); item++) {
          while (H.axes.axis[a] != Axis::undefined) a++;
          H.axes.dim[a] = *item;
          H.axes.axis[a] = n++;
        }
      }

      if (is_temporary (H.name)) M.set_temporary (true);
      setup();
    }





    void Object::create (const String& imagename, Image::Header&template_header)
    {
      M.reset();
      H = template_header;
      H.read_only = false;
      H.axes.sanitise();

      if (imagename.empty()) {

        H.name = "scratch image";
        try { M.add (new guint8 [H.memory_footprint()]); }
        catch (...) { throw Exception ("failed to allocate memory for scratch image data!"); }

      }
      else {

        if (imagename == "-") {
          File::MMap fmap ("", 1024, "mif");
          H.name = fmap.name();
        }
        else H.name = imagename;

        info ("creating image \"" + name() + "\"...");

        NameParser parser;
        parser.parse (H.name);
        std::vector<int> dim (parser.ndim());
        const Format::Base** handler = handlers;

        Axes axes = H.axes;

        try {
          for (; *handler; handler++) 
            if ((*handler)->check (H, H.axes.ndim() - dim.size())) break;
          if (!*handler) throw Exception ("unknown format for image \"" + H.name + "\"");
        }
        catch (...) { throw Exception ("error creating image \"" + H.name + "\""); }

        H.data_type.set_byte_order_native();
        int a = 0;
        for (int n = 0; n < (int) dim.size(); n++) {
          while (H.axes.axis[a] != Axis::undefined) a++;
          dim[n] = axes.dim[a];
        }
        parser.calculate_padding (dim);

        try { 
          std::vector<int> num (dim.size());
          do {
            H.name = parser.name (num);
            (*handler)->create (M, H);
          } while (get_next (num, dim));
        }
        catch (...) { throw Exception ("error creating image \"" + H.name + "\""); }

        if (dim.size()) {
          int a = 0, n = 0;
          for (int i = 0; i < H.axes.ndim(); i++) if (H.axes.axis[i] != Axis::undefined) n++;
          H.axes.set_ndim (n + dim.size());

          for (std::vector<int>::const_iterator item = dim.begin(); item != dim.end(); item++) {
            while (H.axes.axis[a] != Axis::undefined) a++;
            H.axes.dim[a] = *item;
            H.axes.axis[a] = n++;
          }
        }

        if (is_temporary (H.name)) M.output_name = H.name;
      }

      setup();
    }







    void Object::concatenate (std::vector<RefPtr<Object> >& images)
    {
      M.reset();
      if (!images.front() || ! images.back()) throw Exception ("cannot concatenate images: some images are NULL");
      debug ("concatenating images \"" + images.front()->name() + " -> " + images.back()->name() + "\"...");
      M.optimised = false;

      Image::Object& ref (*images[0]);
      for (std::vector<RefPtr<Object> >::iterator it = images.begin()+1; it != images.end(); ++it) {
        if (!*it) throw Exception ("cannot concatenate images: some images are NULL");
        Image::Object& ima (**it);
        if (ima.ndim() != ref.ndim()) throw Exception ("cannot concatenate images: number of dimensions do not match");
        if (ima.header().data_type != ref.header().data_type) throw Exception ("cannot concatenate images: data types do not match");
        for (int n = 0; n < ref.ndim(); n++) {
          if (ima.header().axes.dim[n] != ref.header().axes.dim[n]) throw Exception ("cannot concatenate images: dimensions do not match");
          if (ima.header().axes.axis[n] != ref.header().axes.axis[n] || 
              ima.header().axes.forward[n] != ref.header().axes.forward[n]) throw Exception ("cannot concatenate images: data layouts do not match");
        }
        if (ima.M.list.size() != ref.M.list.size()) throw Exception ("cannot concatenate images: number of files do not match");
      }

      H = ref.header();
      H.name = "{ " + images.front()->name() + " -> " + images.back()->name() + " }";
      H.axes.set_ndim (ref.ndim()+1);
      H.axes.dim[ref.ndim()] = images.size();
      M.optimised = false;
      M.temporary = false;
      M.set_data_type (H.data_type);

      for (std::vector<RefPtr<Object> >::iterator it = images.begin(); it != images.end(); ++it) 
        for (std::vector<Mapper::Entry>::iterator f = (*it)->M.list.begin(); f != (*it)->M.list.end(); ++f) 
          M.add (f->fmap, f->offset);


      start = ref.start;
      memcpy (stride, ref.stride, MRTRIX_MAX_NDIMS*sizeof(gssize));
      stride[ref.ndim()] = ref.voxel_count();
      if (H.data_type.is_complex()) stride[ref.ndim()] *= 2; 

      if (App::log_level > 2) {
        String string ("data increments initialised with start = " + str (start) + ", stride = [ ");
        for (int i = 0; i < ndim(); i++) string += str (stride[i]) + " "; 
        debug (string + "]");
      }
    }







    void Object::setup ()
    {
      if (H.name == "-") H.name = M.list[0].fmap.name();

      debug ("setting up image \"" + H.name + "\"...");
      M.optimised = false;

      set_temporary (M.temporary);

      M.set_read_only (H.read_only);
      M.set_data_type (H.data_type);
      H.sanitise_transform();

      if (M.list.size() == 1 && H.data_type == DataType::Native) M.optimised = true;

      debug ("setting up data increments for \"" + H.name + "\"...");

      start = 0;
      memset (stride, 0, MRTRIX_MAX_NDIMS*sizeof(gssize));

      guint order[ndim()];
      guint last = ndim()-1;
      for (int i = 0; i < ndim(); i++) {
        if (H.axes.axis[i] != Axis::undefined) order[H.axes.axis[i]] = i; 
        else order[last--] = i;
      }

      guint axis;
      gssize mult = 1;

      for (int i = 0; i < ndim(); i++) {
        axis = order[i];
        assert (axis < guint (ndim()));
        if (stride[axis]) throw Exception ("invalid data order specifier for image \"" + H.name + "\": same dimension specified twice");

        stride[axis] = mult * gssize(H.axes.direction(axis));
        if (stride[axis] < 0) start += gsize(-stride[axis]) * gsize(H.axes.dim[axis]-1);
        mult *= gssize(H.axes.dim[axis]);
      }

      if (H.data_type.is_complex()) {
        start *= 2;
        for (int i = 0; i < ndim(); i++) stride[i] *= 2;
      }


      if (App::log_level > 2) {
        String string ("data increments initialised with start = " + str (start) + ", stride = [ ");
        for (int i = 0; i < ndim(); i++) string += str (stride[i]) + " "; 
        debug (string + "]");
      }
    }




    std::ostream& operator<< (std::ostream& stream, const Object& obj)
    {
      stream << "Image object: \"" << obj.name() << "\" [ ";
      for (int n = 0; n < obj.ndim(); n++) stream << obj.dim(n) << " ";
      stream << "]\n Offset: start = " << obj.start << ", stride = [ ";
      for (int n = 0; n < obj.ndim(); n++) stream << obj.stride[n] << " ";
      stream << "]\nHeader:\n" << obj.H << obj.M;
      return (stream);
    }


  }
}
