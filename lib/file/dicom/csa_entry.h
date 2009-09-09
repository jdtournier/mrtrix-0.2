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


    23-09-2008 J-Donald Tournier <d.tournier@brain.org.au>
    * fix handling of ascii floating-point numbers on locales that use a comma as the decimal separator.

    09-09-2009 J-Donald Tournier <d.tournier@brain.org.au>
    * include <cstdio> to allow compilation on Fedora 11

*/

#ifndef __file_dicom_csa_entry_h__
#define __file_dicom_csa_entry_h__

#include <cstdio>

#include "data_type.h"
#include "get_set.h"
#include "file/dicom/element.h"

namespace MR {
  namespace File {
    namespace Dicom {

        class CSAEntry {
          protected:
            const guint8* start;
            const guint8* next;
            const guint8* end;
            bool          print;
            gchar         name[65], vr[5];
            gint32        nitems, num, cnum;

          public:
            CSAEntry (const guint8* start_p, const guint8* end_p, bool output_fields = false);
            bool         parse ();
            const gchar* key () const;
            gint         get_int () const;
            gfloat       get_float () const;  
            void         get_float (float* v) const;

            friend std::ostream& operator<< (std::ostream& stream, const CSAEntry& item);
        };










      inline CSAEntry::CSAEntry (const guint8* start_p, const guint8* end_p, bool output_fields) :
        start (start_p),
        end (end_p),
        print (output_fields)
      {
        if (strncmp ("SV10", (const gchar*) start, 4)) 
          throw Exception ("CSA data is not in SV10 format");

        cnum = 0;
        num = getLE<guint32> (start+8);
        next = start + 16;
      }





      inline bool CSAEntry::parse ()
      {
        if (cnum >= num) return (false);
        start = next;
        if (start >= end + 84) return (false);
        strncpy (name, (const gchar*) start, 64);
        getLE<gint32> (start+64); // vm
        strncpy (vr, (const gchar*) start+68, 4);
        getLE<gint32> (start+72); // syngodt
        nitems = getLE<gint32> (start+76);
        if (print) fprintf (stdout, "    [CSA] %s: ", name);
        next = start + 84;
        if (next + 4 >= end) return (false);

        for (gint32 m = 0; m < nitems; m++) {
          gint32 length = getLE<gint32> (next);
          guint size = 16 + 4*((length+3)/4);
          if (next + size > end) return (false);
          if (print) fprintf(stdout, "%.*s ", length, (const gchar*) next+16);
          next += size;
        }
        if (print) fprintf(stdout, "\n");

        cnum++;
        return (true);
      }




      
      inline const gchar* CSAEntry::key () const        { return (name); }
      
      
      
      
      inline int CSAEntry::get_int () const
      { 
        const guint8* p = start + 84;
        for (gint32 m = 0; m < nitems; m++) {
          gint32 length = getLE<gint32> (p);
          if (length) return (to<int> (String ((const char*) p+16, 4*((length+3)/4))));
          p += 16 + 4*((length+3)/4);
        }
        return (0);
      }




      inline float CSAEntry::get_float () const  
      { 
        const guint8* p = start + 84;
        for (gint32 m = 0; m < nitems; m++) {
          gint32 length = getLE<gint32> (p);
          if (length) return (to<float> (String ((const char*) p+16, 4*((length+3)/4))));
          p += 16 + 4*((length+3)/4);
        }
        return (GSL_NAN);
      }





      inline void CSAEntry::get_float (float* v) const
      { 
        const guint8* p = start + 84;
        for (gint32 m = 0; m < nitems; m++) {
          gint32 length = getLE<gint32> (p);
          if (length) v[m] = to<float> (String ((const char*) p+16, 4*((length+3)/4)));
          p += 16 + 4*((length+3)/4);
        }
      }





      inline std::ostream& operator<< (std::ostream& stream, const CSAEntry& item)
      {
        stream << "[CSA] " << item.name << ":";
        const guint8* next = item.start + 84;

        for (gint32 m = 0; m < item.nitems; m++) {
          gint32 length = getLE<gint32> (next);
          guint size = 16 + 4*((length+3)/4);
          while (length > 0 && !next[16+length-1]) length--;
          stream << " ";
          stream.write ((const gchar*) next+16, length);
          next += size;
        }
         
        return (stream);
      }



    }
  }
}

#endif


