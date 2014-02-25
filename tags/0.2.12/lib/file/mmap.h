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
    * tidied up class structure (MMap::Base is now private to MMap)

*/

#ifndef __file_mmap_h__
#define __file_mmap_h__

#include "ptr.h"

namespace MR {
  namespace File {

    class MMap {
      public:

        MMap () : base (NULL) { }
        MMap (const String& fname, gsize desired_size_if_inexistant = 0, const gchar* suffix = NULL) : base (NULL) { init (fname, desired_size_if_inexistant, suffix); }

        void              init (const String& fname, gsize desired_size_if_inexistant = 0, const gchar* suffix = NULL);
        String            name () const                     { return (base ? base->filename : ""); }
        gsize             size () const                     { return (base ? base->msize : 0); }
        void              resize (gsize  new_size);

        void              map ();
        void*             address () const                  { return (base ? base->addr : NULL); }
        void              unmap ();

        void              set_read_only (bool is_read_only);
        void              mark_for_deletion ()           { if (base) base->delete_after = true; }

        bool              is_ready () const                  { return (base ? base->msize : false); }
        bool              is_mapped () const                 { return (base ? ( base->addr != NULL ) : false); }
        bool              is_read_only () const              { return (base ? base->read_only : true); }
        bool              is_marked_for_deletion () const    { return (base ? base->delete_after : false); }

        bool              changed () const;

        friend std::ostream& operator<< (std::ostream& stream, const MMap& m)
        {
          stream << "MMap: ";
          if (!m.base) stream << "(null)";
          else stream << *m.base;
          return (stream);
        }

      private:
        class Base {
          private:
            Base () : fd (-1), addr (NULL), msize (0), read_only (true), delete_after (false), mtime (0) { }
            ~Base ();

            int               fd;
            String            filename;     /**< The name of the file. */
            void*             addr;         /**< The address in memory where the file has been mapped. */
            gsize             msize;        /**< The size of the memory-mapping, should correspond to the size of the file. */
            bool              read_only;    /**< A flag to indicate whether the file is mapped as read-only. */
            bool              delete_after;
            time_t            mtime;

            void              map ();
            void              unmap ();
            void              resize (gsize  new_size);

            friend class RefPtr<Base>;
            friend class MMap;
            friend std::ostream& operator<< (std::ostream& stream, const Base& m)
            {
              stream << "{ " << m.filename << " [" << m.fd << "] mapped at " << m.addr << ", size " << m.msize << " }";
              return (stream);
            }
        };


        RefPtr<Base>  base;
    };












    inline void MMap::set_read_only (bool is_read_only)   
    { 
      if (!base) return;
      if (base->read_only == is_read_only) return; 
      bool was_mapped = (base->addr != NULL);
      base->unmap(); 
      base->read_only = is_read_only; 
      if (was_mapped) base->map();
    }




    inline void MMap::map ()
    {
      if (!base) throw Exception ("MMap not initialised!");
      if (!base->addr) base->map();
    }



    inline void MMap::unmap ()
    {
      if (!base) throw Exception ("MMap not initialised!");
      if (base->addr) base->unmap();
    }


    inline void MMap::resize (gsize new_size)
    {
      if (!base) throw Exception ("MMap not initialised!");
      base->resize (new_size);
    }





  }
}

#endif

