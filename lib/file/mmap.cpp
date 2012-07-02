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
 

    11-07-2008 J-Donald Tournier <d.tournier@brain.org.au>
    * fixed TMPFILE_ROOT_LEN - now set to 7
    
*/

#include <glib/gstdio.h>
#include <glibmm/stringutils.h>
#include <fcntl.h>
#include <unistd.h>

#ifdef G_OS_WIN32
#include <windows.h>
#else 
#include <sys/mman.h>
#endif

#include "file/mmap.h"
#include "file/config.h"




namespace MR {
  namespace File {

    namespace {

      gchar random_char ()
      {
        gchar c = rand () % 62;
        if (c < 10) return (c+48);
        if (c < 36) return (c+55);
        return (c+61);
      }

    }

    MMap::Base::~Base ()
    {
      unmap(); 
      if (delete_after) {
        debug ("deleting file \"" + filename + "\"...");
        if (g_unlink (filename.c_str())) 
          error ("WARNING: error deleting file \"" + filename + "\": " + Glib::strerror(errno));
      }
    }


    void MMap::Base::map()
    {
      if (msize == 0) throw Exception ("attempt to map file \"" + filename + "\" using invalid mmap!");
      if (addr) return;

      if ((fd = g_open (filename.c_str(), (read_only ? O_RDONLY : O_RDWR), 0755)) < 0) 
        throw Exception ("error opening file \"" + filename + "\": " + Glib::strerror(errno));

      try {
#ifdef G_OS_WIN32
        HANDLE handle = CreateFileMapping ((HANDLE) _get_osfhandle(fd), NULL, 
            (read_only ? PAGE_READONLY : PAGE_READWRITE), 0, msize, NULL);
        if (!handle) throw 0;
        addr = (void*) MapViewOfFile (handle, (read_only ? FILE_MAP_READ : FILE_MAP_ALL_ACCESS), 0, 0, msize);
        if (!addr) throw 0;
        CloseHandle (handle);
#else 
        addr = (void *) mmap((char*)0, msize, (read_only ? PROT_READ : PROT_READ | PROT_WRITE), MAP_SHARED, fd, 0);
        if (addr == MAP_FAILED) throw 0;
#endif
        debug ("file \"" + filename + "\" mapped at " + str (addr) 
            + ", size " + str (msize) 
            + " (read-" + ( read_only ? "only" : "write" ) + ")"); 
      }
      catch (...) {
        close (fd);
        addr = NULL;
        throw Exception ("memory-mapping failed for file \"" + filename + "\": " + Glib::strerror(
#ifdef G_OS_WIN32
              GetLastError()
#else
              errno
#endif
              ));
      }
    }





    void MMap::Base::unmap()
    {
      if (!addr) return;

      debug ("unmapping file \"" + filename + "\"");

#ifdef G_OS_WIN32
      if (!UnmapViewOfFile ((LPVOID) addr))
#else 
        if (munmap (addr, msize))
#endif
          error ("error unmapping file \"" + filename + "\": " + Glib::strerror(
#ifdef G_OS_WIN32
              GetLastError()
#else
              errno
#endif
              ));

      close (fd);
      fd = -1;
      addr = NULL;
    }









    void MMap::Base::resize (gsize new_size)
    {
      debug ("resizing file \"" + filename + "\" to " + str (new_size) + "...");

      if (read_only) throw Exception ("attempting to resize read-only file \"" + filename + "\"");
      unmap();

      if ((fd = g_open (filename.c_str(), O_RDWR, 0755)) < 0) 
        throw Exception ("error opening file \"" + filename + "\" for resizing: " + Glib::strerror(errno));

      int status = ftruncate (fd, new_size);

      close (fd);
      fd = -1;
      if (status) throw Exception ("cannot resize file \"" + filename + "\": " + Glib::strerror(errno));

      msize = new_size;
    }









    void MMap::init (const String& fname, gsize desired_size_if_inexistant, const gchar* suffix)
    {
      base = new Base;

      try {
        if (fname.size()) {
          debug ("preparing file \"" + fname + "\"");

          base->filename = fname;
          struct_stat64 sbuf;
          if (STAT64 (base->filename.c_str(), &sbuf)) {

            if (errno != ENOENT) 
              throw Exception ("cannot stat file \"" + base->filename + "\": " + Glib::strerror(errno));

            if (desired_size_if_inexistant == 0) 
              throw Exception ("cannot access file \"" + base->filename + "\": " + Glib::strerror(errno));

            int fid = g_open (base->filename.c_str(), O_CREAT | O_RDWR | O_EXCL, 0755);
            if (fid < 0) throw Exception ("error creating file \"" + base->filename + "\": " + Glib::strerror(errno));

            int status = ftruncate (fid, desired_size_if_inexistant);
            close (fid);
            if (status) throw Exception ("WARNING: cannot resize file \"" + base->filename + "\": " + Glib::strerror(errno));

            base->read_only = false;
            base->msize = desired_size_if_inexistant;

            return;
          }
          else {
            if (desired_size_if_inexistant) 
              throw Exception ("cannot create file \"" + base->filename + "\": it already exists");

            base->msize = sbuf.st_size;
            base->mtime = sbuf.st_mtime;

            return;
          }
        }
      }
      catch (Exception) {
        base = NULL;
        throw;
      }




      if (!desired_size_if_inexistant) throw Exception ("cannot create empty scratch file");

      debug ("creating and mapping scratch file");

      assert (suffix);
      base->filename = String (TMPFILE_ROOT) + "XXXXXX." + suffix; 

      int fid;
      do {
        for (int n = 0; n < 6; n++) 
          base->filename[TMPFILE_ROOT_LEN+n] = random_char();
      } while ((fid = g_open (base->filename.c_str(), O_CREAT | O_RDWR | O_EXCL, 0755)) < 0);


      int status = ftruncate (fid, desired_size_if_inexistant);
      close (fid);
      if (status) throw Exception ("cannot resize file \"" + base->filename + "\": " + Glib::strerror(errno));

      base->msize = desired_size_if_inexistant;
      base->read_only = false;
    }




    bool MMap::changed () const
    { 
      if (!base) return (false);
      struct_stat64 sbuf;
      if (STAT64 (base->filename.c_str(), &sbuf)) return (false);
      if (off_t (base->msize) != sbuf.st_size) return (true);
      if (base->mtime != sbuf.st_mtime) return (true);
      return (false);
    }




  }
}


