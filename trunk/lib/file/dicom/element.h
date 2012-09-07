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


    17-03-2009 J-Donald Tournier <d.tournier@brain.org.au>
    * modify to allow use of either TR1 unordered map or SGI hash_map for the DICOM dictionary
    
*/

#ifndef __file_dicom_element_h__
#define __file_dicom_element_h__

#include "hash_map.h"
#include "file/mmap.h"
#include "file/dicom/definitions.h"

namespace MR {
  namespace File {
    namespace Dicom {

      typedef enum _ElementType {
        INVALID,
        INT,
        UINT,
        FLOAT,
        STRING,
        SEQ,
        OTHER
      } ElementType;


      class Sequence {
        public:
          Sequence (guint16 group, guint16 element, guint8* end) : group (group), element (element), end (end) { }
          guint16 group, element;
          guint8* end;
      };



      class Element {
        protected:
          static UnorderedMap<guint32, const gchar*>::Type dict;
          static void           init_dict();

          File::MMap            fmap;
          void                  set_explicit_encoding();
          bool                  read_GR_EL();

          guint8*               next;
          guint8*               start;
          bool                  is_explicit;
          bool                  is_BE;
          bool                  previous_BO_was_BE;


        public:

          guint16               group, element, VR;
          guint32               size;
          guint8*               data;
          std::vector<Sequence>  sequence;

          void                  set (const String& filename);
          bool                  read ();

          bool                  is (guint16 group, guint16 element) const;

          String                tag_name () const;
          guint32               tag () const;
          ElementType           type () const;
          bool                  is_big_endian () const;
          std::vector<gint32>   get_int () const;
          std::vector<guint32>  get_uint () const;
          std::vector<double>   get_float () const;
          std::vector<String>   get_string () const;
          guint                 offset (guint8* address) const;

          guint                 level () const { return sequence.size(); }

          void                  print () const;

          friend std::ostream& operator<< (std::ostream& stream, const Element& item);
      };



















      inline bool Element::is (guint16 Group, guint16 Element) const
      {
        if (group != Group) return (false);
        return (element == Element);
      }





      inline bool Element::is_big_endian () const { return (is_BE); }




      inline guint32 Element::tag () const
      {
        union __DICOM_group_element_pair__ { guint16 s[2]; guint32 i; } val = { {
#if G_BYTE_ORDER == G_BIG_ENDIAN
          group, element
#else
            element, group
#endif 
        } };
        return (val.i);
      }





      inline String Element::tag_name() const
      {
        if (dict.empty()) init_dict();
        const gchar* s = dict[tag()];
        return (s ? s : "");
      }





      inline guint Element::offset (guint8* address) const { return (address - (guint8*) fmap.address()); }






    }
  }
}

#endif

