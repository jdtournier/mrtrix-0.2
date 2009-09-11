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

#ifndef __histogram_h__
#define __histogram_h__

#include "mrtrix.h"

namespace MR {
  namespace Image {
    class Position;
  }

  class Histogram {
    public:
      Histogram (Image::Position& ima, guint num_buckets=100);

      guint   frequency (guint index) const    { return (list[index].frequency); }
      float   value (guint index) const        { return (list[index].value); }
      guint   num () const                    { return (list.size()); }
      float   first_min () const;

    protected:
      class Entry {
        public:
          Entry () : frequency (0), value (0.0) { }
          guint  frequency;
          float  value;
      };

      std::vector<Entry> list;

  };

}

#endif
