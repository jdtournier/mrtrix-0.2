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

#ifndef __file_ximg_h__
#define __file_ximg_h__

#include "get_set.h"
#include "file/mmap.h"


namespace MR {
  namespace File {

    class XImg {

      private:
        MMap           mmap;
        const guint8*  bof () const;

      public:
        void           read (const String& filename);
        String         name () const;

        const guint8*  pixel_data () const;
        gint           width () const;
        gint           height () const;
        gint           depth () const;

    };

    std::ostream& operator<< (std::ostream& stream, const XImg& X);











    inline const guint8*  XImg::bof () const { return ((guint8*) mmap.address()); }

    inline void XImg::read (const String& filename) 
    {
      mmap.init (filename);
      mmap.map();
    }

    inline String XImg::name () const { return (mmap.name()); }

    inline const guint8* XImg::pixel_data () const { return (bof() + getBE<gint32> (bof() + 0x4)); }
    inline gint XImg::width () const      { return (getBE<gint32> (bof() + 0x8)); }
    inline gint XImg::height () const     { return (getBE<gint32> (bof() + 0xc)); }
    inline gint XImg::depth () const      { return (getBE<gint32> (bof() + 0x10)); }

    inline std::ostream& operator<< (std::ostream& stream, const XImg& X) 
    {
      stream << "name: \"" << X.name() << ", pixel_data at " << X.pixel_data() << ", dim: [ " << X.width() << " " << X.height() << " ]\n";
      return (stream);
    }

  }
}

#endif



