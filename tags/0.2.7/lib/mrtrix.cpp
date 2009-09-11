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

#include <gsl/gsl_version.h>


#include "mrtrix.h"
#include "image/object.h"



namespace MR {


  /************************************************************************
   *                     MRtrix version information                      *
   ************************************************************************/

  const guint mrtrix_major_version = MRTRIX_MAJOR_VERSION;
  const guint mrtrix_minor_version = MRTRIX_MINOR_VERSION;
  const guint mrtrix_micro_version = MRTRIX_MICRO_VERSION;


  void (*print) (const String& msg) = NULL;
  void (*error) (const String& msg) = NULL;
  void (*info) (const String& msg) = NULL;
  void (*debug) (const String& msg) = NULL;

  int Exception::level_offset = 0;






  /************************************************************************
   *                         ProgressBar functions                        *
   ************************************************************************/

  namespace ProgressBar {

    String message; 
    guint current_val = 0, percent = 0;
    float multiplier = 0.0;
    bool display = true;
    bool stop = false;
    Glib::Timer stop_watch;

    void (*init_func) () = NULL;
    void (*display_func) () = NULL;
    void (*done_func) () = NULL;





    void init (guint target, const String& msg)
    {
      stop = false;
      message = msg;
      if (target) multiplier = 100.0/((float) target);
      else multiplier = NAN;
      current_val = percent = 0;
      if (gsl_isnan (multiplier)) stop_watch.start();
      init_func ();
      if (display) display_func ();
    }



  }








  /************************************************************************
   *                       Miscellaneous functions                        *
   ************************************************************************/

  std::vector<gfloat> parse_floats (const String& spec)
  {
    std::vector<float> V;
    try {
      if (!spec.size()) throw (0);
      String::size_type start = 0, end;
      do {
        end = spec.find_first_of (',', start);
        V.push_back (to<gfloat> (spec.substr (start, end-start)));
        start = end+1;
      } while (end < String::npos);
    }
    catch (...) { throw Exception ("can't parse floating-point sequence specifier \"" + spec + "\""); }
    return (V);
  }




  std::vector<gint> parse_ints (const String& spec, gint last)
  {
    std::vector<gint> V;
    try {
      if (!spec.size()) throw (0);
      String::size_type start = 0, end;
      int num[3];
      int i = 0;
      do {
        end = spec.find_first_of (",:", start);
        String str (strip (spec.substr (start, end-start)));
        lowercase (str);
        if (str == "end") {
          if (last == INT_MAX) throw (0);
          num[i] = last;
        }
        else num[i] = to<gint> (spec.substr (start, end-start));

        if (spec[end] == ':') { i++; if (i > 2) throw (0); }
        else {
          if (i) {
            gint inc, last;
            if (i == 2) { inc = num[1]; last = num[2]; }
            else { inc = 1; last = num[1]; }
            if (inc * (last - num[0]) < 0) inc = -inc;
            for (; ( inc > 0 ? num[0] <= last : num[0] >= last ) ; num[0] += inc) V.push_back (num[0]);
          }
          else V.push_back (num[0]);
          i = 0;
        }

        start = end+1;
      } while (end < String::npos);
    }
    catch (...) { throw Exception ("can't parse integer sequence specifier \"" + spec + "\""); }
    return (V);
  }









  std::vector<String> split (const String& string, const gchar* delimiters, bool ignore_empty_fields)
  {
    std::vector<String> V;
    String::size_type start = 0, end;
    try {
      do {
        end = string.find_first_of (delimiters, start);
        V.push_back (string.substr (start, end-start));
        start = ignore_empty_fields ? string.find_first_not_of (delimiters, end+1) : end+1;
      } while (end < String::npos);
    }
    catch (...)  { throw Exception ("can't split string \"" + string + "\""); }
    return (V);
  }



}
