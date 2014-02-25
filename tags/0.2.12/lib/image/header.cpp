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


    17-03-2009 J-Donald Tournier <d.tournier@brain.org.au>
    * additional sanity checks in sanitise_transform(): 
    * - make sure voxel sizes are finite numbers
    * - make sure all entries in the transform matrix are finite.
    * use sane defaults otherwise.

    14-12-2010 J-Donald Tournier <d.tournier@brain.org.au>
    * fix bug in transform re-jigging code to handle 45° oblique orientations

*/

#include "image/header.h"
#include "math/matrix.h"
#include "math/linalg.h"

namespace MR {
  namespace Image {


    void Header::reset ()
    {
      comments.clear();
      axes.set_ndim (0);
      name.clear();
      data_type = DataType();
      offset = 0.0;
      scale = 1.0;
      read_only = true;
      format = NULL;
      trans_I2R.reset();
      trans_R2I.reset();
      trans_P2R.reset();
      trans_R2P.reset();
      DW_scheme.reset();
    }





    void Header::set_transform (const Math::Matrix& M)
    {
      if (M.rows() != 4 || M.columns() != 4) 
        throw Exception ("invalid transform specified for image \"" + name + "\"");

      trans_I2R.copy (M);
      trans_I2R(3,0) = trans_I2R(3,1) = trans_I2R(3,2) = 0.0; 
      trans_I2R(3,3) = 1.0;

      sanitise_transform();
    }




    namespace {

      guint find_max_in_row (const Math::Matrix& M, guint row)
      {
        if (fabs (M (row,0)) > fabs (M (row,1))) {
          if (fabs (M (row,0)) > fabs (M (row,2))) return (0);
          else return (2);
        }
        else {
          if (fabs (M (row,1)) > fabs (M (row,2))) return (1);
          else return (2);
        }
      }

      inline guint not_any_of (guint a, guint b) 
      {
        for (guint i = 0; i < 3; ++i) {
          if (a == i || b == i) 
            continue;
          return (i);
        }
        assert (0);
        return (UINT_MAX);
      }

      inline void disambiguate_reorientation (guint* permutation) 
      {
        if (permutation[0] == permutation[1]) 
          permutation[0] = not_any_of (permutation[0], permutation[2]);

        if (permutation[0] == permutation[2]) 
          permutation[0] = not_any_of (permutation[0], permutation[1]);

        if (permutation[1] == permutation[2]) 
          permutation[1] = not_any_of (permutation[0], permutation[1]);
      }
    }




    void Header::sanitise_transform ()
    {
      float default_vox_size = 0.0;
      debug ("sanitising transformation matrix...");
      {
        int count = 0;

        for (int n = 0; n < std::min (ndim(), 3); ++n) {
          if (gsl_finite (axes.vox[n])) {
            default_vox_size += axes.vox[n];
            ++count;
          }
        }

        default_vox_size = count ? default_vox_size / count : 1.0;

        bool issue_warning = false;
        for (int n = 0; n < std::min (ndim(), 3); ++n) {
          if (!gsl_finite (axes.vox[n])) {
            axes.vox[n] = default_vox_size;
            issue_warning = true;
          }
        }
        if (issue_warning) 
          error ("invalid voxel sizes - resetting to sane values");
      }

      if (trans_I2R.is_valid()) {
        if (trans_I2R.rows() != 4 || trans_I2R.columns() != 4) {
          trans_I2R.reset();
          error ("transform matrix is not 4x4 - resetting to sane defaults");
        }
        else {
          for (guint i = 0; i < 3; i++) {
            for (guint j = 0; j < 4; j++) {
              if (!gsl_finite (trans_I2R(i,j))) {
                trans_I2R.reset();
                error ("transform matrix contains invalid entries - resetting to sane defaults");
                break;
              }
            }
            if (!trans_I2R.is_valid()) break;
          }
        }
      }

      float vox[3];
      int dim[3];
      for (int n = 0; n < 3; ++n) {
        if (n < ndim()) { vox[n] = axes.vox[n]; dim[n] = axes.dim[n]; }
        else { vox[n] = default_vox_size; dim[n] = 1; }
      }

      if (!trans_I2R.is_valid()) {
        trans_I2R.allocate (4,4);
        trans_I2R.identity();
        trans_I2R(0,3) = -0.5 * dim[0] * vox[0];
        trans_I2R(1,3) = -0.5 * dim[1] * vox[1];
        trans_I2R(2,3) = -0.5 * dim[2] * vox[2];
      }

      trans_I2R(3,0) = trans_I2R(3,1) = trans_I2R(3,2) = 0.0; trans_I2R(3,3) = 1.0;

      guint permutation[3] = { 
        find_max_in_row (trans_I2R, 0),
        find_max_in_row (trans_I2R, 1),
        find_max_in_row (trans_I2R, 2) 
      };

      disambiguate_reorientation (permutation);

      bool flip[3] = {
        ( trans_I2R(0, permutation[0]) < 0.0 ),
        ( trans_I2R(1, permutation[1]) < 0.0 ),
        ( trans_I2R(2, permutation[2]) < 0.0 )
      };

      if (permutation[0] != 0 || permutation[1] != 1 || permutation[2] != 2 || flip[0] || flip[1] || flip[2]) {
        if (ndim() < 3) 
          axes.set_ndim (3);
        
        bool forward[] = { axes.forward [permutation[0]], axes.forward [permutation[1]], axes.forward [permutation[2]] };
        int newdim[] = { dim [permutation[0]], dim [permutation[1]], dim [permutation[2]] };
        int axis[] = { axes.axis [permutation[0]], axes.axis [permutation[1]], axes.axis [permutation[2]] };
        float newvox[] = { vox [permutation[0]], vox [permutation[1]], vox [permutation[2]] };
        String desc[] = { axes.desc [permutation[0]], axes.desc [permutation[1]], axes.desc [permutation[2]] };
        String units[] = { axes.units [permutation[0]], axes.units [permutation[1]], axes.units [permutation[2]] };

        Math::Matrix T(trans_I2R);

        for (guint i = 0; i < 3; i++) {
          for (guint n = 0; n < 3; n++) trans_I2R(n, i) = T (n, permutation[i]);
          if (flip[i]) {
            forward[i] = !forward[i];
            float length = (newdim[i]-1) * newvox[i];
            for (guint n = 0; n < 3; n++) {
              trans_I2R (n, i) = -trans_I2R (n, i);
              trans_I2R (n, 3) += length * T (n, permutation[i]);
            }
          }
          axes.dim[i] = newdim[i];
          axes.vox[i] = newvox[i];
          axes.forward[i] = forward[i];
          axes.axis[i] = axis[i];
          axes.desc[i] = desc[i];
          axes.units[i] = units[i];
        }

      }

      for (int n = 0; n < 3; ++n) {
        if (n < ndim()) vox[n] = axes.vox[n]; 
        else vox[n] = default_vox_size;
      }

      Math::PseudoInverter pinvert (trans_R2I, trans_I2R);
      pinvert.invert (trans_R2I, trans_I2R);

      Math::Matrix V(4,4);
      V.zero();
      V(0,0) = vox[0];
      V(1,1) = vox[1];
      V(2,2) = vox[2];
      V(3,3) = 1.0;

      trans_P2R.multiply (trans_I2R, V);
      V(0,0) = 1.0/V(0,0);
      V(1,1) = 1.0/V(1,1);
      V(2,2) = 1.0/V(2,2);
      trans_R2P.multiply (V, trans_R2I);
    }







    void Header::merge (const Header& H)
    {
      if (data_type != H.data_type) 
        throw Exception ("data types differ between image files for \"" + name + "\"");

      if (offset != H.offset || scale != H.scale) 
        throw Exception ("scaling coefficients differ between image files for \"" + name + "\"");

      if (axes.ndim() != H.axes.ndim()) 
        throw Exception ("dimension mismatch between image files for \"" + name + "\"");

      for (int n = 0; n < axes.ndim(); n++) {
        if (axes.dim[n] != H.axes.dim[n]) 
          throw Exception ("dimension mismatch between image files for \"" + name + "\"");

        if (axes.axis[n] != H.axes.axis[n] || axes.forward[n] != H.axes.forward[n])
          throw Exception ("data layout differs image files for \"" + name + "\"");

        if (axes.vox[n] != H.axes.vox[n])
          error ("WARNING: voxel dimensions differ between image files for \"" + name + "\"");
      }

      
      for (std::vector<String>::const_iterator item = H.comments.begin(); item != H.comments.end(); item++)
        if (find (comments.begin(), comments.end(), *item) == comments.end())
          comments.push_back (*item);

      if (!trans_I2R.is_valid() && H.trans_I2R.is_valid()) set_transform (H.trans_I2R);
      if (!DW_scheme.is_valid() && H.DW_scheme.is_valid()) DW_scheme = H.DW_scheme; 
    }





    String Header::description() const
    {
      String desc ( 
            "************************************************\n"
            "Image:               \"" + name + "\"\n"
            "************************************************\n"
            "  Format:            " + ( format ? format : "undefined" ) + "\n"
            "  Dimensions:        ");

      int i;
      for (i = 0; i < axes.ndim(); i++) {
        if (i) desc += " x ";
        desc += str (axes.dim[i]);
      }



      desc += "\n  Voxel size:        ";

      for (i = 0; i < axes.ndim(); i++) {
        if (i) desc += " x ";
        desc += gsl_isnan (axes.vox[i]) ? "?" : str (axes.vox[i]);
      }




      desc += "\n  Dimension labels:  ";

      for (i = 0; i < axes.ndim(); i++)  
        desc += ( i ? "                     " : "" ) + str (i) + ". " 
          + ( axes.desc[i].size() ? axes.desc[i] : "undefined" ) + " ("
          + ( axes.units[i].size() ? axes.units[i] : "?" ) + ")\n";



      desc += String ("  Data type:         ") + ( data_type.description() ? data_type.description() : "invalid" ) + "\n"
            "  Data layout:       [ ";


      for (i = 0; i < axes.ndim(); i++) 
        desc += axes.axis[i] == Axis::undefined ? "? " : ( axes.forward[i] ? '+' : '-' ) + str (axes.axis[i]) + " ";



      desc += "]\n"
            "  Data scaling:      offset = " + str (offset) + ", multiplier = " + str (scale) + "\n"
            "  Comments:          " + ( comments.size() ? comments[0] : "(none)" ) + "\n";



      for (i = 1; i < (int) comments.size(); i++)
        desc += "                     " + comments[i] + "\n";

      if (trans_I2R.is_valid()) {
        desc += "  Transform:         ";
        guint i, j;
        for (i = 0; i < trans_I2R.rows(); i++) {
          if (i) desc +=  "                     ";
          for (j = 0; j < trans_I2R.columns(); j++) {
            gchar buf[14], buf2[14];
            g_snprintf (buf, 14, "%.4g", trans_I2R(i,j));
            g_snprintf (buf2, 14, "%12.10s", buf);
            desc += buf2;
          }
          desc += "\n";

        }
      }

      if (DW_scheme.is_valid()) 
        desc += "  DW scheme:         " + str (DW_scheme.rows()) + " x " + str (DW_scheme.columns()) + "\n";

      return (desc);
    }



    std::ostream& operator<< (std::ostream& stream, const Header& H)
    {
      stream << H.description();
      return (stream);
    }

  }
}
