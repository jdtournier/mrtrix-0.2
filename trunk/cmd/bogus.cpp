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
#include "ptr.h"

using namespace MR; 

SET_VERSION_DEFAULT;

DESCRIPTION = {
  "this is used to test stuff.",
  NULL
};


ARGUMENTS = { Argument::End }; 
OPTIONS = { Option::End };

EXECUTE {
  std::vector<RefPtr<int> > V;
  for (int n = 0; n < 100; n++) V.push_back (RefPtr<int> (new int (rand())));
  std::sort (V.begin(), V.end());

  for (std::vector<RefPtr<int> >::iterator i = V.begin(); i != V.end(); ++i)
    **i = rand();

  std::sort (V.begin(), V.end());

  std::vector<RefPtr<int> > V2 (V);

  std::sort (V2.begin(), V2.end());
}

