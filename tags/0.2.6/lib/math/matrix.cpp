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

#include "ptr.h"
#include "math/matrix.h"

namespace MR {
  namespace Math {


    void Matrix::load (const String& filename)
    {
      std::ifstream in (filename.c_str());
      if (!in) throw Exception ("cannot open matrix file \"" + filename + "\": " + Glib::strerror (errno));

      std::vector< RefPtr< std::vector<double> > >   data;

      do {
        String sbuf;
        getline (in, sbuf);
        if (in.bad()) throw Exception ("error reading matrix file \"" + filename + "\": " + Glib::strerror (errno));
        if (in.eof()) break;

        sbuf = strip (sbuf.substr (0, sbuf.find_first_of ('#')));
        if (sbuf.size()) {
          data.push_back (RefPtr< std::vector<double> > (new std::vector<double>));

          std::istringstream stream (sbuf);
          do {
            double val;
            stream >> val;
            data.back()->push_back (val);
          } while (stream.good());

          if (data.size() > 1) 
            if (data.back()->size() != data[0]->size())
              throw Exception ("uneven rows in matrix file \"" + filename + "\"");
        }
      } while (in.good());




      allocate (data.size(), data[0]->size());

      for (guint i = 0; i < rows(); i++)
        for (guint j = 0; j < columns(); j++)
          (*this)(i, j) = (*data[i])[j];
    }





    void Matrix::save (const String& filename) const
    {

      std::ofstream out (filename.c_str());
      if (!out) throw Exception ("cannot open matrix file \"" + filename + "\": " + Glib::strerror (errno));

      for (guint i = 0; i < rows(); i++) {
        for (guint j = 0; j < columns(); j++) 
          out << (*this)(i,j) << "\t";
        out << "\n";
      }

    }





    std::ostream& operator<< (std::ostream& stream, const Matrix& M)
    {
      for (guint i = 0; i < M.rows(); i++) {
        for (guint j = 0; j < M.columns(); j++) 
          stream << MR::printf ("%11.4g ", M(i,j));
        stream << "\n";
      }
      return (stream);
    }



  }
}
