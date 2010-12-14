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

#include "app.h"
#include "dwi/SH.h"

using namespace MR; 

SET_VERSION_DEFAULT;

DESCRIPTION = {
  "this is used to test stuff.",
  NULL
};


ARGUMENTS = { 
  Argument ("SH", "SH", "SH").type_file(),
  Argument ("start", "start", "start").type_sequence_float(),
  Argument::End 
}; 
OPTIONS = { Option::End };

EXECUTE {
  Math::Vector SHv;
  SHv.load (argument[0].get_string());
  float SH [SHv.size()];
  for (size_t n = 0; n < SHv.size(); n++)
    SH[n] = SHv[n];

  std::vector<float> start_v = parse_floats (argument[1].get_string());

  Point start (start_v[0], start_v[1], start_v[2]);
  start.normalise();
  int lmax = DWI::SH::LforN (SHv.size());
  VAR (lmax);

  float peak = DWI::SH::get_peak (SH, lmax, start, false);
  VAR (peak);
}

