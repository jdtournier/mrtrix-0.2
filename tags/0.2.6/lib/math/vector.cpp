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

#include <glibmm/stringutils.h>
#include <fstream>
#include "math/vector.h"

namespace MR {
  namespace Math {

    void Vector::load (const String& filename)
    {
      std::ifstream in (filename.c_str());
      if (!in) throw Exception ("cannot open file \"" + filename + "\": " + Glib::strerror (errno));

      std::vector<double>   data;
      double value;

      while (true) {
        in >> value;
        if (in.eof()) break;
        data.push_back (value);
      }
      in.close();

      allocate (data.size());
      for (guint n = 0; n < size(); n++)
        (*this)[n] = data[n];
    }





    void Vector::save (const String& filename) const
    {
      std::ofstream out (filename.c_str());
      if (!out) throw Exception ("cannot open file \"" + filename + "\": " + Glib::strerror (errno));

      for (guint i = 0; i < size(); i++)
        out << (*this)[i] << "\n";
    }





    void Vector::print() const
    {
      gchar buf[12];
      for (guint i = 0; i < size(); i++) {
        sprintf (buf, "%.4g", (*this)[i]);
        fprintf (stderr, "%11.10s\n", buf);
      }
    }



    std::ostream& operator<<(std::ostream& stream, const Vector& vec)
    {
      stream << "[ ";
      for (guint i = 0; i < vec.size(); i++)
        stream << vec[i] << " ";
      stream << "]";
      return(stream);
    }

  }
}

