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

#ifndef __image_mapper_h__
#define __image_mapper_h__

#include "data_type.h"
#include "file/mmap.h"
#include "image/header.h"
#include "image/format/base.h"
#include "math/complex_number.h"

namespace MR {
  namespace Dialog { class File; }
  namespace Image {

    class Object;
    class Header;

    class Mapper {
      public:
        void                   reset ();
        void                   add (const String& filename, gsize offset = 0, gsize desired_size_if_inexistant = 0);
        void                   add (const File::MMap& fmap, gsize offset = 0);
        void                   add_gz (const String& filename, const String& gz_filename, gsize offset = 0, gsize desired_size_if_inexistant = 0);
        void                   add_gz (const File::MMap& fmap, const String& gz_filename, gsize offset = 0);
        void                   add (guint8* memory_buffer);


        float32                re (gsize offset) const;
        void                   re (float32 val, gsize offset);
        float32                im (gsize offset) const;
        void                   im (float32 val, gsize offset); 

        void                   set_temporary (bool temp);
        String                 output_name;



        static void gzip (const String& original, const String& gzfile);
        static String gunzip (const String& gzfile, const char* suffix);

      protected:
        Mapper ();
        ~Mapper ();

        class Entry {
          public:
            File::MMap fmap;
            gsize    offset;
            guint8*    start () const;
            String gzfilename;
            friend std::ostream& operator<< (std::ostream& stream, const Entry& m)
            {
              stream << "Mapper::Entry: offset = " << m.offset << ", " << m.fmap;
              return (stream);
            }
        };

        std::vector<Entry>    list;
        guint8*               mem;
        guint8**              segment;
        gsize                 segsize;

        bool                  optimised, temporary, files_new;

        void                  set_data_type (DataType dt);
        void                  set_read_only (bool read_only);

        Entry&                operator[] (guint index);
        
        float32                (*get_func) (const void* data, gsize i);
        void                   (*put_func) (float32 val, void* data, gsize i);

        static float32         getBit       (const void* data, gsize i);
        static float32         getInt8      (const void* data, gsize i);
        static float32         getUInt8     (const void* data, gsize i);
        static float32         getInt16LE   (const void* data, gsize i);
        static float32         getUInt16LE  (const void* data, gsize i);
        static float32         getInt16BE   (const void* data, gsize i);
        static float32         getUInt16BE  (const void* data, gsize i);
        static float32         getInt32LE   (const void* data, gsize i);
        static float32         getUInt32LE  (const void* data, gsize i);
        static float32         getInt32BE   (const void* data, gsize i);
        static float32         getUInt32BE  (const void* data, gsize i);
        static float32         getFloat32LE (const void* data, gsize i);
        static float32         getFloat32BE (const void* data, gsize i);
        static float32         getFloat64LE (const void* data, gsize i);
        static float32         getFloat64BE (const void* data, gsize i);

        static void            putBit       (float32 val, void* data, gsize i);
        static void            putInt8      (float32 val, void* data, gsize i);
        static void            putUInt8     (float32 val, void* data, gsize i);
        static void            putInt16LE   (float32 val, void* data, gsize i);
        static void            putUInt16LE  (float32 val, void* data, gsize i);
        static void            putInt16BE   (float32 val, void* data, gsize i);
        static void            putUInt16BE  (float32 val, void* data, gsize i);
        static void            putInt32LE   (float32 val, void* data, gsize i);
        static void            putUInt32LE  (float32 val, void* data, gsize i);
        static void            putInt32BE   (float32 val, void* data, gsize i);
        static void            putUInt32BE  (float32 val, void* data, gsize i);
        static void            putFloat32LE (float32 val, void* data, gsize i);
        static void            putFloat32BE (float32 val, void* data, gsize i);
        static void            putFloat64LE (float32 val, void* data, gsize i);
        static void            putFloat64BE (float32 val, void* data, gsize i);


        void                   map (const Header& H);
        void                   unmap (const Header& H);
        bool                   is_mapped () const { return (segment); }


        friend class Object;
        friend class Dialog::File;
        friend class Value;
        friend std::ostream& operator<< (std::ostream& stream, const Mapper& dmap);
    };













    inline guint8* Mapper::Entry::start () const { return ((guint8*) fmap.address() + offset); }


    inline Mapper::Mapper () :
      mem (NULL),
      segment (NULL),
      segsize (0),
      optimised (false),
      temporary (false),
      files_new (true),
      get_func (NULL),
      put_func (NULL)
    { 
    }





    inline void Mapper::reset ()
    {
      list.clear();
      segsize = 0; 
      get_func = NULL;
      put_func = NULL;
      optimised = temporary = false;
      files_new = true;
      output_name.clear();
      delete [] mem;
      delete [] segment;
      mem = NULL;
      segment = NULL;
    }





    inline void Mapper::add (const String& filename, gsize offset, gsize desired_size_if_inexistant)
    {
      Entry entry;
      entry.fmap.init (filename, desired_size_if_inexistant, "tmp"); 
      if (entry.fmap.is_read_only()) 
        files_new = false;
      entry.offset = offset;
      list.push_back (entry);
    }









    inline void Mapper::add (const File::MMap& fmap, gsize offset)
    {
      assert (!fmap.is_mapped());
      Entry entry;
      entry.fmap = fmap;
      if (entry.fmap.is_read_only()) 
        files_new = false;
      entry.offset = offset;
      list.push_back (entry);
    }




    inline void Mapper::add_gz (const String& filename, const String& gz_filename, gsize offset, gsize desired_size_if_inexistant)
    {
      add (filename, offset, desired_size_if_inexistant);
      list.back().gzfilename = gz_filename;
    }



    inline void Mapper::add_gz (const File::MMap& fmap, const String& gz_filename, gsize offset)
    {
      add (fmap, offset);
      list.back().gzfilename = gz_filename;
    }



    inline void Mapper::add (guint8* memory_buffer)
    {
      assert (mem == NULL);
      assert (list.size() == 0);
      mem = memory_buffer;
    }








    /** \brief make all files in the list read-only or read/write.
     *
     * Makes all files currently in the list read-only or read/write.
     */
    inline void Mapper::set_read_only (bool read_only)
    {
      for (guint s = 0; s < list.size(); s++) {
        list[s].fmap.set_read_only (read_only); 
        if (segment) segment[s] = list[s].start();
      }
    }






    inline Mapper::Entry& Mapper::operator[] (guint index) { return (list[index]); }




    inline float32 Mapper::re (gsize offset) const 
    {
      if (optimised) return (((float32*) segment[0])[offset]);
      gssize nseg (offset/segsize);
      return (get_func (segment[nseg], offset - nseg*segsize)); 
    }





    inline void Mapper::re (float32 val, gsize offset)
    { 
      if (optimised) ((float32*) segment[0])[offset] = val;
      gssize nseg (offset/segsize);
      put_func (val, segment[nseg], offset - nseg*segsize); 
    }





    inline float32 Mapper::im (gsize offset) const
    { 
      if (optimised) return (((float32*) segment[0])[offset+1]);
      gssize nseg (offset/segsize);
      return (get_func (segment[nseg], offset - nseg*segsize + 1)); 
    }




    inline void Mapper::im (float32 val, gsize offset)
    { 
      if (optimised) ((float32*) segment[0])[offset+1] = val;
      gssize nseg (offset/segsize);
      put_func (val, segment[nseg], offset - nseg*segsize + 1); 
    }



    inline void Mapper::set_temporary (bool temp)                { temporary = temp; }

  }
}

#endif


