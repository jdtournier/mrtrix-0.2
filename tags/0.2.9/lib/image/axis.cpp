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

#include "image/axis.h"

namespace MR {
  namespace Image {

    const gchar* Axis::left_to_right = "left->right";
    const gchar* Axis::posterior_to_anterior = "posterior->anterior";
    const gchar* Axis::inferior_to_superior = "inferior->superior";
    const gchar* Axis::time = "time";
    const gchar* Axis::real_imag = "real-imaginary";
    const gchar* Axis::millimeters = "mm";
    const gchar* Axis::milliseconds = "ms";


    std::vector<Axis> parse_axes_specifier (const Axes& original, const String& specifier)
    {
      std::vector<Axis> parsed (original.ndim());

      int str = 0;
      int lim = 0;
      int end = specifier.size();
      int current = 0;

      try {
        while (str <= end) {
          parsed[current].forward = original.forward[current];
          if (specifier[str] == '+') { parsed[current].forward = true; str++; }
          else if (specifier[str] == '-') { parsed[current].forward = false; str++; }
          else if (!( specifier[str] == '\0' || specifier[str] == ',' ) && !isdigit (specifier[str])) throw 0;

          lim = str;

          if (specifier[str] == '\0' || specifier[str] == ',') parsed[current].axis = original.axis[current];
          else {
            while (isdigit (specifier[lim])) lim++;
            if (specifier[lim] != ',' && specifier[lim] != '\0') throw 0;
            parsed[current].axis = to<guint> (specifier.substr (str, lim-str));
          }

          str = lim+1;
          current++;
        }
      }
      catch (int) { throw Exception ("malformed axes specification \"" + specifier + "\""); }
      catch (...) { throw; }

      if (current != original.ndim()) 
        throw Exception ("incorrect number of axes in axes specification \"" + specifier + "\"");

      check_axes_specifier (parsed, original.ndim());

      return (parsed);
    }





    void check_axes_specifier (const std::vector<Axis>& parsed, int ndims)
    {
      for (guint n = 0; n < parsed.size(); n++) {
        if (parsed[n].axis >= ndims) 
          throw Exception ("axis " + str (parsed[n].axis) + " out of range");

        for (guint i = 0; i < n; i++) 
          if (parsed[i].axis == parsed[n].axis) 
            throw Exception ("duplicate axis (" + str (parsed[n].axis) + ")");
      }
    }






    std::ostream& operator<< (std::ostream& stream, const Axes& axes)
    {
      stream << "dim [ ";
      for (int n = 0; n < axes.ndim(); n++) stream << axes.dim[n] << " ";
      stream << "], vox [ ";
      for (int n = 0; n < axes.ndim(); n++) stream << axes.vox[n] << " ";
      stream << "], axes [ ";
      for (int n = 0; n < axes.ndim(); n++) stream << ( axes.forward[n] ? '+' : '-' ) << axes.axis[n] << " ";
      stream << "], desc [ ";
      for (int n = 0; n < axes.ndim(); n++) stream << "\n" << axes.desc[n] << "\" ";
      stream << "], units [ ";
      for (int n = 0; n < axes.ndim(); n++) stream << "\n" << axes.units[n] << "\" ";

      return (stream);
    }




  }
}

