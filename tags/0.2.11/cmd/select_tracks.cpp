/*
    Copyright 2008 Brain Research Institute, Melbourne, Australia

    Written by J-Donald Tournier, 03/06/2011.

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

#include "app.h"
#include "math/matrix.h"
#include "point.h"
#include "dwi/tractography/file.h"

using namespace MR; 

SET_VERSION_DEFAULT;

DESCRIPTION = {
  "select tracks from one or several track files.",
  "Note that the track properties (as listed in the header of the file, "
    "and displayed using track_info) will be taken from the first input file only.",
  NULL
};

ARGUMENTS = {
  Argument ("input", "input track", "The input track file(s).", true, true).type_file (),
  Argument ("output", "output tracks file", "The output tracks file").type_file(),
  Argument::End
};


OPTIONS = { 
  Option ("number", "index numbers", 
      "specify a sequence of indices into the track file(s) corresponding "
      "to the indices of the tracks to include in the output.")
    .append (Argument ("sequence", "sequence", "the sequence of index numbers").type_sequence_int ()),
  
  Option::End 
};


class Range {
  public:
    Range (int num) : A (num), B (num), inc (0) { }
    Range (int start, int end, int increment) : 
      A (start), B (end), inc (increment) { 
      }
    Range (const String& specifier) : inc (1) {
      std::vector<String> tokens = split (specifier, ":");
      if (tokens.empty()) 
        throw Exception ("error: empty number range");
      A = B = to<int> (tokens[0]);
      if (tokens.size() > 1) {
        B = to<int> (tokens[1]);
        if (tokens.size() > 2) 
          inc = to<int> (tokens[2]);
      }
      check_inc ();
    }

    int start () const { return A; }
    int end () const { return B; }
    int increment () const { return inc; }

    static std::vector<Range> parse (const String& specifier) 
    {
      std::vector<Range> list;
      std::vector<String> bits = split (specifier, ",");
      for (guint n = 0; n < bits.size(); ++n)
        list.push_back (Range (bits[n]));
      return list;
    }

    static void check_increasing (const std::vector<Range>& list) 
    {
      int last = -1;
      for (guint n = 0; n < list.size(); ++n) 
        if (list[n].start() <= last || ( list[n].start() != list[n].end() && list[n].increment() < 0 ))
          throw Exception ("number sequence should increase monotonically");
    }

  private:
    int A, B, inc;

    void check_inc () { 
      if (!inc) inc = 1;
      if ((B-A) * inc < 0) 
        inc = -inc;
      if (inc * ((B-A)/inc) != (B-A))
        throw Exception ("error: end value not reachable with specified increment");
    }
};



EXECUTE {
  DWI::Tractography::Properties properties;
  DWI::Tractography::Reader file;
  file.open (argument[0].get_string(), properties);
  file.close();

  DWI::Tractography::Writer writer;
  writer.create (argument.back().get_string(), properties);

  std::vector<Point> tck;

  std::vector<Range> list;
  
  std::vector<OptBase> opt = get_options (0); // number
  if (opt.size()) {
    list = Range::parse (opt[0][0].get_string());
    Range::check_increasing (list);
  }

  for (guint n = 0; n < list.size(); ++n)
    std::cout << "  " << list[n].start() << ":" << list[n].increment() << ":" << list[n].end() << "\n";

  ProgressBar::init (0, "selecting tracks...");

  try {
    int num = 0, target = list.size() ? list[0].start() : 0;
    guint n = 0;
    for (guint nfile = 0; nfile < argument.size()-1; ++nfile) {
      file.open (argument[nfile].get_string(), properties);
      writer.total_count += to<guint> (properties["total_count"]);
      while (file.next (tck)) {
        if (list.empty() || num == target) {
          writer.append (tck);
          if (list.size()) {
            if (target == list[n].end()) {
              ++n;
              if (n >= list.size())
                throw 1;
              target = list[n].start();
            }
            else target += list[n].increment();
          }
        }
        ++num;
        ProgressBar::inc();
      }
      file.close();
    }
  }
  catch (int) { }

  ProgressBar::done();
  writer.close();
}

