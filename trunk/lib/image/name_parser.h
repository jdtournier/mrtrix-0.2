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

#ifndef __image_name_parser_h__
#define __image_name_parser_h__

#include <glibmm/fileutils.h>
#include <glibmm/refptr.h>

#include "ptr.h"

namespace MR {
  namespace Image {

    class NameParserItem {
      protected:
        guint                  seq_length;
        String                 str;
        std::vector<int>       seq;

      public:
        NameParserItem ();

        void                   set_str (const String& s);
        void                   set_seq (const String& s);
        void                   clear ();

        String                 string () const;
        const std::vector<int>&  sequence () const;
        std::vector<int>&      sequence ();
        bool                   is_string () const;
        bool                   is_sequence () const;
        guint                  size () const;

        void                   calc_padding (guint maxval = 0);

        friend std::ostream& operator<< (std::ostream& stream, const NameParserItem& item);
    };




    class NameParser {
      private:
        std::vector<NameParserItem> array;
        std::vector<guint>          seq_index;
        String                      folder_name;
        String                      specification;
        String                      current_name;
        Glib::Dir*                  folder;

        void                        insert_str (const String& str);
        void                        insert_seq (const String& str);

      public:
        NameParser ();
        void                       parse (const String& imagename, guint max_num_sequences = UINT_MAX);
        guint                      num () const;
        String                     spec () const;
        const NameParserItem&      operator[] (guint i) const;
        const std::vector<int>&    sequence (guint index) const;
        guint                      ndim () const;
        guint                      index_of_sequence (guint number = 0) const;

        bool                       match (const String& file_name, std::vector<int>& indices) const;
        void                       calculate_padding (const std::vector<int>& maxvals);
        String                     name (const std::vector<int>& indices);

        String                     get_next_match (std::vector<int>& indices, bool return_seq_index = false);

        friend std::ostream& operator<< (std::ostream& stream, const NameParser& parser);
    };





    class ParsedName {
      protected:
        std::vector<int>    indices;
        String              filename;

      public:
        ParsedName (const String& name, const std::vector<int>& index);

        String               name () const;
        guint                ndim () const;
        int                  index (guint num) const;

        bool                 operator< (const ParsedName& pn) const;
        friend std::ostream& operator<< (std::ostream& stream, const ParsedName& pin);
    };





    class ParsedNameList : public std::vector< RefPtr<ParsedName> > {
      protected:
        void                count_dim (std::vector<int>& dim, guint& current_entry, guint current_dim) const;
        guint               max_name_size;

      public:
        std::vector<int>    parse_scan_check (const String& specifier, guint max_num_sequences = UINT_MAX);
        void                scan (NameParser& parser);

        std::vector<int>    count () const;
        guint               biggest_filename_size () const;
    };
















    /****************************************************************
     *               NameParserItem inline functions                *
     ***************************************************************/

    inline NameParserItem::NameParserItem () : seq_length (0) { }



    inline void NameParserItem::clear () 
    {
      str.clear();
      seq.clear();
      seq_length = 0;
    }



    inline void NameParserItem::set_str (const String& s)
    {
      clear ();
      str = s;
    }



    inline void NameParserItem::set_seq (const String& s)
    { 
      clear ();
      if (s.size()) seq = parse_ints (s);
      seq_length = 1; 
    }



    inline guint NameParserItem::size () const                         { return (seq_length ? seq_length : str.size()); }
    inline String NameParserItem::string () const                      { return (str); }
    inline const std::vector<int>& NameParserItem::sequence () const   { return (seq); }
    inline std::vector<int>& NameParserItem::sequence ()               { return (seq); }
    inline bool NameParserItem::is_string () const                     { return (seq_length == 0); }
    inline bool NameParserItem::is_sequence () const                   { return (seq_length != 0); }




    /****************************************************************
     *                 NameParser inline functions                  *
     ***************************************************************/

    inline NameParser::NameParser () : folder (NULL) { } 
    inline String NameParser::spec () const                                  { return (specification); }
    inline guint NameParser::num () const                                    { return (array.size()); }
    inline guint NameParser::ndim () const                                   { return (seq_index.size()); }
    inline const NameParserItem& NameParser::operator[] (guint i) const      { return (array[i]); }
    inline const std::vector<int>& NameParser::sequence (guint index) const { return (array[seq_index[index]].sequence()); }
    inline guint NameParser::index_of_sequence (guint number) const          { return (seq_index[number]); }

    inline void NameParser::insert_str (const String& str)
    {
      NameParserItem item;
      item.set_str (str);
      array.insert (array.begin(), item);
    }
    inline void NameParser::insert_seq (const String& str)
    {
      NameParserItem item;
      item.set_seq (str);
      array.insert (array.begin(), item);
      seq_index.push_back (array.size()-1);
    }





    /****************************************************************
     *                 ParsedName inline functions                  *
     ***************************************************************/

    inline ParsedName::ParsedName (const String& name, const std::vector<int>& index) : indices (index), filename (name) { }
    inline String ParsedName::name () const                    { return (filename); }
    inline guint ParsedName::ndim () const                            { return (indices.size()); }
    inline int ParsedName::index (guint num) const                   { return (indices[num]); }
    inline guint ParsedNameList::biggest_filename_size () const       { return (max_name_size); }




  }
}

#endif

