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

#include <algorithm>
#include <glib/gprintf.h>
#include <glibmm/miscutils.h>
#include <glibmm/stringutils.h>

#include "image/name_parser.h"

namespace MR {
  namespace Image {

    namespace {

      inline bool is_seq_char (gchar c) { return (g_ascii_isdigit (c) || c == ',' || c == ':' || c == '%'); } 

      inline int last_bracket (const String& str, int from, gchar c)
      {
        while (str[from] != c && from >= 0) from--;
        return (from);
      }

      inline bool in_seq (const std::vector<int>& seq, int val) 
      {
        if (seq.size() == 0) return (true);
        for (guint i = 0; i < seq.size(); i++) 
          if (seq[i] == val) return (true);
        return (false);
      }


    }






    void NameParser::parse (const String& imagename, guint max_num_sequences)
    {
      specification = imagename;
      if (Glib::file_test (imagename, Glib::FILE_TEST_IS_DIR)) {
        array.resize (1);
        array[0].set_str (imagename);
        return;
      }

      folder_name = Glib::path_get_dirname (specification);


      try {
        String::size_type pos;
        String basename = Glib::path_get_basename (specification);
        guint num = 0;

        while ((pos = basename.find_last_of (']')) < String::npos && num < max_num_sequences) {
          insert_str (basename.substr (pos+1));
          basename = basename.substr (0, pos);
          if ((pos = basename.find_last_of ('[')) == String::npos) 
          throw Exception ("malformed image sequence specifier for image \"" + specification + "\"");

          insert_seq (basename.substr (pos+1));
          num++;
          basename = basename.substr (0, pos);
        }

        insert_str (basename);


        for (guint i = 0; i < array.size(); i++) 
          if (array[i].is_sequence()) 
            if (array[i].sequence().size()) 
              for (guint n = 0; n < array[i].sequence().size()-1; n++) 
                for (guint m = n+1; m < array[i].sequence().size(); m++) 
                  if (array[i].sequence()[n] == array[i].sequence()[m]) 
                    throw Exception ("malformed image sequence specifier for image \"" + specification + "\" (duplicate indices)");
      }
      catch (...) {
        array.resize (1);
        array[0].set_str (imagename);
        throw;
      }
    }





    std::ostream& operator<< (std::ostream& stream, const NameParserItem& item)
    {
      if (item.is_string()) stream << "\"" << item.string() << "\"";
      else {
        if (item.sequence().size()) stream << item.sequence();
        else stream << "[ any ]";
      }
      return (stream);
    }





    std::ostream& operator<< (std::ostream& stream, const NameParser& parser)
    {
      stream << "Image::NameParser: " << parser.specification << "\n";
      for (guint i = 0; i < parser.array.size(); i++) 
        stream << "  " << i << ": " << parser.array[i] << "\n";
      return (stream);
    }










    bool NameParser::match (const String& file_name, std::vector<int>& indices) const
    {
      int current = 0;
      guint num = 0;
      indices.resize (seq_index.size());

      for (guint i = 0; i < array.size(); i++) {
        if (array[i].is_string()) {
          if (file_name.substr(current, array[i].string().size()) != array[i].string()) return (false);
          current += array[i].string().size();
        }
        else {
          int x = current;
          while (isdigit (file_name[current])) current++;
          x = to<int> (file_name.substr (x, current-x));
          if (!in_seq (array[i].sequence(), x)) return (false);
          indices[num] = x;
          num++;
        }
      }
      return (true);
    }




    void NameParser::calculate_padding (const std::vector<int>& maxvals)
    {
      assert (maxvals.size() == seq_index.size());
      for (guint n = 0; n < seq_index.size(); n++) 
        assert (maxvals[n] > 0);

      for (guint n = 0; n < seq_index.size(); n++) {
        guint m = seq_index.size() - 1 - n;
        NameParserItem& item (array[seq_index[n]]);
        if (item.sequence().size()) {
          if (maxvals[m]) 
            if (item.sequence().size() != (guint) maxvals[m]) 
              throw Exception ("dimensions requested in image specifier \"" + specification 
                  + "\" do not match supplied header information");
        }
        else {
          item.sequence().resize (maxvals[m]);
          for (guint i = 0; i < item.sequence().size(); i++)
            item.sequence()[i] = i;
        }

        item.calc_padding (maxvals[m]);
      }
    }





    void NameParserItem::calc_padding (guint maxval)
    {
      for (guint i = 0; i < sequence().size(); i++) {
        assert (sequence()[i] >= 0);
        if (maxval < (guint) sequence()[i]) maxval = sequence()[i];
      }

      seq_length = 1;
      for (guint num = 10; maxval >= num; num *= 10) 
        seq_length += 1; 
    }





    String NameParser::name (const std::vector<int>& indices)
    {
      if (!seq_index.size()) 
        return (Glib::build_filename (folder_name, array[0].string()));

      assert (indices.size() == seq_index.size());

      String str;
      guint n = seq_index.size()-1;
      for (guint i = 0; i < array.size(); i++) {
        if (array[i].is_string()) str += array[i].string();
        else { 
          guint last = str.size();
          str.resize (last + array[i].size() + 1);
          g_sprintf (&str[last], "%*.*d", array[i].size(), array[i].size(), array[i].sequence()[indices[n]]);
          n--; 
        }
      }

      return (Glib::build_filename (folder_name, str));
    }







    String NameParser::get_next_match (std::vector<int>& indices, bool return_seq_index)
    {
      if (!folder) {
        try { folder = new Glib::Dir (folder_name); }
        catch (...) {
          folder = NULL;
          throw Exception ("failed to open folder \"" + folder_name + "\": " + Glib::strerror(errno));
        }
      }

      String fname;
      while ((fname = folder->read_name()).size()) {
        if (match (fname, indices)) {
          if (return_seq_index) {
            for (guint i = 0; i < ndim(); i++) {
              if (sequence(i).size()) {
                guint n = 0;
                while (indices[i] != sequence(i)[n]) n++; 
                indices[i] = n;
              }
            }
          }
          return (Glib::build_filename (folder_name, fname));
        }
      }

      return ("");
    }


    bool ParsedName::operator< (const ParsedName& pn) const
    {
      for (guint i = 0; i < ndim(); i++) 
        if (index(i) != pn.index(i)) 
          return (index(i) < pn.index(i));
      return (false);
    }




    std::vector<int> ParsedNameList::parse_scan_check (const String& specifier, guint max_num_sequences)
    {
      NameParser parser;
      parser.parse (specifier);

      scan (parser);
      std::sort (begin(), end());
      std::vector<int> dim = count();

      for (guint n = 0; n < dim.size(); n++) 
        if (parser.sequence(n).size()) 
          if (dim[n] != (int) parser.sequence(n).size()) 
            throw Exception ("number of files found does not match specification \"" + specifier + "\"");

      return (dim);
    }






    void ParsedNameList::scan (NameParser& parser)
    {
      std::vector<int> index;
      if (parser.ndim() == 0) {
        push_back (RefPtr<ParsedName> (new ParsedName (parser.name (index), index)));
        return;
      }

      String entry;

      while ((entry = parser.get_next_match (index, true)).size()) 
        push_back (RefPtr<ParsedName> (new ParsedName (entry, index)));

      if (!size()) 
        throw Exception ("no matching files found for image specifier \"" + parser.spec() + "\"");
    }






    std::vector<int> ParsedNameList::count () const   
    { 
      if (!(*this)[0]->ndim()) {
        if (size() == 1) return (std::vector<int>());
        else throw Exception ("image number mismatch");
      }

      std::vector<int> dim ((*this)[0]->ndim(), 0);
      guint current_entry = 0;
      
      count_dim (dim, current_entry, 0);

      return (dim);
    }




    void ParsedNameList::count_dim (std::vector<int>& dim, guint& current_entry, guint current_dim) const
    {
      int n;
      bool stop = false;
      RefPtr<const ParsedName> first_entry ((*this)[current_entry]);

      for (n = 0; current_entry < size(); n++) {
        for (guint d = 0; d < current_dim; d++)
          if ((*this)[current_entry]->index(d) != first_entry->index(d)) stop = true;
        if (stop) break;

        if (current_dim < (*this)[0]->ndim()-1) count_dim (dim, current_entry, current_dim+1);
        else current_entry++;
      }

      if (dim[current_dim] && dim[current_dim] != n) 
        throw Exception ("number mismatch between number of images along different dimensions");

      dim[current_dim] = n;
    }













    std::ostream& operator<< (std::ostream& stream, const ParsedName& pin)
    {
      stream << "[ ";
      for (guint n = 0; n < pin.ndim(); n++) 
        stream << pin.index(n) << " ";
      stream << "] " << pin.name();
      return (stream);
    }




  }
}
