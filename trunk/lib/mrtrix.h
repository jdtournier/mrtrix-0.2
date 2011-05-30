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


    22-07-2008 J-Donald Tournier <d.tournier@brain.org.au>
    * fix va_list handing in printf()

    09-07-2008 J-Donald Tournier <d.tournier@brain.org.au>
    * added getline() method to handle Unix/DOS end-of-line

    02-09-2008 J-Donald Tournier <d.tournier@brain.org.au>
    * added is_temporary() method to identify temporary files

    02-10-2008 J-Donald Tournier <d.tournier@brain.org.au>
    * change Exception silencing to a priority level changing approach

    15-10-2008 J-Donald Tournier <d.tournier@brain.org.au>
    * remove MR::DICOM_DW_gradients_PRS flag

    31-10-2008 J-Donald Tournier <d.tournier@brain.org.au>
    * include <cstring> & <cstdlib> to allow compilation on Fedora 9

    15-03-2010 J-Donald Tournier <d.tournier@brain.org.au>
    * add shorten() function to reduce long filenames 

    15-07-2010 J-Donald Tournier <d.tournier@brain.org.au>
    * version upped to 9

*/

#ifndef __mrtrix_h__
#define __mrtrix_h__

#include <limits.h>
#include <stdarg.h>
#include <errno.h>
#include <assert.h>
#include <sstream>
#include <algorithm>
#include <cctype>
#include <iostream>
#include <vector>
#include <string>
#include <cstring>
#include <cstdlib>
#include <cmath>

#include <glib/gtypes.h>
#include <glib/gutils.h>
#include <glibmm/timer.h>
#include <glibmm/miscutils.h>

#include <gsl/gsl_version.h>
#include <gsl/gsl_math.h>



#define MRTRIX_MAJOR_VERSION 0
#define MRTRIX_MINOR_VERSION 2
#define MRTRIX_MICRO_VERSION 10


/** Prints the file and line number. Useful for debugging purposes. */
#define TEST std::cerr << Glib::get_application_name() << ": line " << __LINE__ << " in " << __func__ << "() from file " << __FILE__ << "\n"


/** Prints a variable name and its value, followed by the file and line number. Useful for debugging purposes. */
#define VAR(variable) std::cerr << Glib::get_application_name() << ": " << #variable << " = " << (variable) << " (in " << __func__ << "() from " << __FILE__  << ": " << __LINE__ << ")\n"

#define TRACE std::cerr << "in " << __PRETTY_FUNCTION__ << " (" << __FILE__ << ": " << __LINE__ << ")\n"

#define MRTRIX_MAX_NDIMS 16
#define BUSY_INTERVAL 0.1
#define GUI_SPACING 5

#define MODIFIERS ( \
        GDK_SHIFT_MASK | \
        GDK_CONTROL_MASK | \
        GDK_BUTTON1_MASK | \
        GDK_BUTTON2_MASK | \
        GDK_BUTTON3_MASK | \
        GDK_BUTTON4_MASK | \
        GDK_BUTTON5_MASK | \
        GDK_SUPER_MASK | \
        GDK_HYPER_MASK | \
        GDK_META_MASK )

#ifdef __APPLE__
# define CTRL_CMD_MASK GDK_META_MASK
#else
# define CTRL_CMD_MASK GDK_CONTROL_MASK
#endif

#define TMPFILE_ROOT "mrtrix-"
#define TMPFILE_ROOT_LEN 7


#ifdef G_OS_WIN32
typedef struct __stat64 struct_stat64;
# define STAT64 _stat64
#else
# ifdef __APPLE__
typedef struct stat struct_stat64;
#  define STAT64 stat
# else
typedef struct stat64 struct_stat64;
#  define STAT64 stat64
# endif
#endif

namespace std {

  template <class T> inline ostream& operator<< (ostream& stream, const vector<T>& V)
  {
    stream << "[ ";
    for (guint n = 0; n < V.size(); n++) stream << V[n] << " ";
    stream << "]";
    return (stream);
  }

}


namespace MR {

  extern const guint mrtrix_major_version;
  extern const guint mrtrix_minor_version;
  extern const guint mrtrix_micro_version;


  typedef std::string String;
  typedef gfloat  float32;
  typedef gdouble float64;




  extern void (*print) (const String& msg);
  extern void (*error) (const String& msg);
  extern void (*info)  (const String& msg);
  extern void (*debug) (const String& msg);








  class Exception {
    public:
      Exception (const String& msg, int log_level = 1) : 
        description (msg),
        level (log_level) { display(); }

      const String description;
      const int level;

      void  display () const {
        if (level + level_offset < 2) error (description);
        else if (level + level_offset == 2) info (description);
        else debug (description);
      }

      class Lower {
        public:
          Lower (int amount = 1)  { level_offset = amount; }
          ~Lower () { level_offset = 0; }
          friend class Exception;
      };

    private:
      static int level_offset;
  };




  namespace ProgressBar {

    extern String message; 
    extern guint current_val, percent;
    extern float multiplier;
    extern bool display, stop;
    extern Glib::Timer stop_watch;

    extern void (*init_func)();
    extern void (*display_func)();
    extern void (*done_func)();


    void init (guint target, const String& msg);
    int inc ();
    void done ();
  }


  namespace Image {
    typedef enum {
      Default,
      Magnitude,
      Real,
      Imaginary,
      Phase,
      RealImag
    } OutputType;
  }


  namespace ByteOrder {

    gint16 swap (gint16 v);
    guint16 swap (guint16 v);
    gint32 swap (gint32 v);
    guint32 swap (guint32 v);
    float32 swap (float32 v);
    float64 swap (float64 v);

    gint16 LE (gint16 v);
    gint16 BE (gint16 v);
    guint16 LE (guint16 v);
    guint16 BE (guint16 v);
    gint32 LE (gint32 v);
    gint32 BE (gint32 v);
    guint32 LE (guint32 v);
    guint32 BE (guint32 v);
    float32 LE (float32 v);
    float32 BE (float32 v);
    float64 LE (float64 v);
    float64 BE (float64 v);

  }






  inline bool is_temporary (const String& file) { return (Glib::path_get_basename(file).compare (0, TMPFILE_ROOT_LEN, TMPFILE_ROOT) == 0); }






  inline std::istream& getline (std::istream& stream, String& string) 
  {
    std::getline (stream, string);
    if (string[string.size()-1] == 015) string.resize (string.size()-1);
    return (stream);
  }



  template <class T> inline String str (const T& value)
  { 
    std::ostringstream stream; 
    try { stream << value; }
    catch (...) { throw Exception ("error converting value to string"); }
    return (stream.str());
  }



  inline String shorten (const String& string, guint longest = 40, guint prefix = 10) 
  {
    if (string.size() > longest) 
      return (string.substr (0,prefix) + "..." + string.substr (string.size()-longest+prefix-3));
    else return (string);
  }



  template <class T> inline T to (const String& string)
  {
    std::istringstream stream (string);
    T value;
    try { stream >> value; }
    catch (...) { throw Exception ("error converting string \"" + string + "\""); }
    return (value);
  }


  inline String printf (const gchar* format, ...)
  {
    va_list list;
    va_start (list, format);
    guint len = g_vsnprintf (NULL, 0, format, list) + 1;
    va_end (list);
    gchar buf[len];
    va_start (list, format);
    g_vsnprintf (buf, len, format, list);
    va_end (list);
    return (buf);
  }


  inline String strip (const String& string, const gchar* ws = " \t\n", bool left = true, bool right = true)
  {
    std::string::size_type start = ( left ? string.find_first_not_of (ws) : 0 );
    if (start == String::npos) return ("");
    std::string::size_type end   = ( right ? string.find_last_not_of (ws) + 1 : String::npos );
    return (string.substr (start, end - start));
  }



  inline void  replace (String& string, gchar orig, gchar final)
  {
    for (String::iterator i = string.begin(); i != string.end(); ++i)
      if (*i == orig) *i = final;
  }


  inline String& lowercase (String& string)
  {
    for_each (string.begin(), string.end(), tolower);
    return (string);
  }

  inline String& uppercase (String& string)
  {
    for_each (string.begin(), string.end(), toupper);
    return (string);
  }

  inline String lowercase (const String& string)
  {
    String ret; ret.resize (string.size());
    transform (string.begin(), string.end(), ret.begin(), tolower);
    return (ret);
  }

  inline String uppercase (const String& string)
  {
    String ret; ret.resize (string.size());
    transform (string.begin(), string.end(), ret.begin(), toupper);
    return (ret);
  }

  std::vector<String> split (const String& string, const gchar* delimiters = " \t\n", bool ignore_empty_fields = false); 

  inline String join (std::vector<String>& V, const String& delimiter)
  {
    String ret;
    if (V.empty()) return (ret);
    ret = V[0];
    for (std::vector<String>::iterator i = V.begin()+1; i != V.end(); ++i)
      ret += delimiter + *i;
    return (ret);
  }

  std::vector<gfloat> parse_floats (const String& spec);
  std::vector<int>   parse_ints (const String& spec, int last = G_MAXINT);

  inline int round (float x) { return (int (x + (x > 0.0 ? 0.5 : -0.5))); }




  template <class T> inline T maxvalue (const T& v0, const T& v1, const T& v2)
  {
    return (v0 > v1 ? ( v0 > v2 ? v0 : ( v1 > v2 ? v1 : v2 ) ) : ( v1 > v2 ? v1 : ( v0 > v2 ? v0 : v2 ))); 
  }

  template <class T> inline int maxindex (const T& v0, const T& v1, const T& v2)
  {
    return (v0 > v1 ? ( v0 > v2 ? 0 : ( v1 > v2 ? 1 : 2 ) ) : ( v1 > v2 ? 1 : ( v0 > v2 ? 0 : 2 ))); 
  }

  template <class T> inline T maxvalue (const T v[3]) { return (maxvalue (v[0], v[1], v[2])); }
  template <class T> inline int maxindex (const T v[3]) { return (maxindex (v[0], v[1], v[2])); }




  template <class T> inline T minvalue (const T& v0, const T& v1, const T& v2)
  {
    return (v0 < v1 ? ( v0 < v2 ? v0 : ( v1 < v2 ? v1 : v2 ) ) : ( v1 < v2 ? v1 : ( v0 < v2 ? v0 : v2 ))); 
  }

  template <class T> inline int minindex (const T& v0, const T& v1, const T& v2)
  {
    return (v0 < v1 ? ( v0 < v2 ? 0 : ( v1 < v2 ? 1 : 2 ) ) : ( v1 < v2 ? 1 : ( v0 < v2 ? 0 : 2 ))); 
  }

  template <class T> inline T minvalue (const T v[3]) { return (minvalue (v[0], v[1], v[2])); }
  template <class T> inline int minindex (const T v[3]) { return (minindex (v[0], v[1], v[2])); }




  template <class T> inline void set_all (std::vector<T>& V, const T& value)
  {
    for (guint n = 0; n < V.size(); n++) V[n] = value;
  }




  template <class T> inline bool get_next (std::vector<T>& pos, const std::vector<T>& limits)
  {
    guint axis = 0;
    while (axis < limits.size()) {
      pos[axis]++;
      if (pos[axis] < limits[axis]) return (true);
      pos[axis] = 0;
      axis++;
    }
    return (false);
  }








  namespace ProgressBar {

    inline int inc ()
    {
      current_val++; 
      if (!display) return (stop);
      guint t = gsl_isnan (multiplier) ? (guint) (stop_watch.elapsed() / BUSY_INTERVAL) : (guint) (multiplier*current_val);
      if (percent != t) {
        percent = t;
        display_func();
      }
      return (stop);
    }

    inline void done () { if (display) done_func(); }
  }



}

#endif

