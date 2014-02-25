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


    17-10-2008 J-Donald Tournier <d.tournier@brain.org.au>
    * add SET_VERSION_DEFAULT, SET_VERSION, SET_AUTHOR and SET_COPYRIGHT macros

*/

#ifndef __app_h__
#define __app_h__

#include "args.h"


#define DESCRIPTION const gchar* __command_description[]
#define ARGUMENTS   const MR::Argument __command_arguments[]
#define OPTIONS const MR::Option __command_options[]
#define SET_VERSION(a, b, c) const guint __command_version[3] = { (a), (b), (c) }
#define SET_AUTHOR(a) const char* __command_author = (a)
#define SET_COPYRIGHT(a) const char* __command_copyright = (a)

#define SET_VERSION_DEFAULT \
  SET_VERSION (MRTRIX_MAJOR_VERSION, MRTRIX_MINOR_VERSION, MRTRIX_MICRO_VERSION); \
  SET_COPYRIGHT ("Copyright 2008 Brain Research Institute, Melbourne, Australia"); \
  SET_AUTHOR ("J-Donald Tournier (d.tournier@brain.org.au)")

#define EXECUTE class MyApp : public MR::App { \
  public: \
  MyApp (int argc, gchar** argv, const gchar** cmd_desc, const MR::Argument* cmd_args, const MR::Option* cmd_opts, \
      const guint* cmd_version, const char* cmd_author, const char* cmd_copyright) : \
  App (argc, argv, cmd_desc, cmd_args, cmd_opts, cmd_version, cmd_author, cmd_copyright) { } \
  void execute (); \
}; \
int main (int argc, gchar* argv[]) { \
  try { \
    MyApp app (argc, argv, __command_description, __command_arguments, __command_options, __command_version, __command_author, __command_copyright); \
    app.run (argc, argv); \
  } \
  catch (Exception) { return (1); } \
  catch (int ret) { return (ret); } \
  return (0); } \
void MyApp::execute ()

#define DEFAULT_OPTIONS_OFFSET 65536U


namespace MR {
  
  namespace Viewer { class Window; }

  class ParsedOption {
    public:
      guint index;
      std::vector<const gchar*> args;
  };


  class App
  {
    public:
      App (int argc, gchar** argv, const gchar** cmd_desc, const MR::Argument* cmd_args, const MR::Option* cmd_opts,
          const guint* cmd_version, const char* cmd_author, const char* cmd_copyright);
      virtual ~App ();

      void   run (int argc, gchar** argv);

      static int    log_level;

      static const guint*           version;
      static const char*            copyright;
      static const char*            author;

      friend std::ostream& operator<< (std::ostream& stream, const App& app);
      friend std::ostream& operator<< (std::ostream& stream, const OptBase& opt);

    protected:
      void                          parse_arguments ();

      void                          print_help () const;
      void                          print_full_usage () const;
      void                          print_full_argument_usage (const Argument& arg) const;
      void                          print_full_option_usage (const Option& opt) const;
      guint                         match_option (const gchar* stub) const;
      void                          sort_arguments (int argc, gchar** argv);

      static const char**           command_description;
      static const Argument*        command_arguments;
      static const Option*          command_options;

      static const Option           default_options[];
      std::vector<const gchar*>     parsed_arguments;
      std::vector<ParsedOption>     parsed_options;

      const gchar*                  option_name (guint num) const { 
        return (num < DEFAULT_OPTIONS_OFFSET ? command_options[num].sname : default_options[num - DEFAULT_OPTIONS_OFFSET].sname ); 
      }

      virtual void                  execute () = 0;

      std::vector<ArgBase>          argument;
      std::vector<OptBase>          option;

      std::vector<OptBase>          get_options (guint index);
  };


  void cmdline_print (const String& msg);
  void cmdline_error (const String& msg);
  void cmdline_info  (const String& msg);
  void cmdline_debug (const String& msg);




  inline std::vector<OptBase> App::get_options (guint index)
  {
    std::vector<OptBase> a;
    for (guint n = 0; n < option.size(); n++) 
      if (option[n].index == index)
        a.push_back (option[n]);
    return (a);
  }

}

#endif
