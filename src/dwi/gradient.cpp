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

#include "point.h"
#include "math/matrix.h"
#include "dwi/gradient.h"
#include "dwi/shells.h"

namespace MR {
  namespace DWI {

    void normalise_grad (Math::Matrix &grad)
    {
      if (grad.columns() != 4)
          throw Exception ("invalid gradient matrix dimensions");
      for (guint i = 0; i < grad.rows(); i++) {
        double norm = sqrt(grad(i,0)*grad(i,0)+grad(i,1)*grad(i,1)+grad(i,2)*grad(i,2));
        if (norm) {
          grad(i,0) /= norm;
          grad(i,1) /= norm;
          grad(i,2) /= norm;
        } else {
          grad(i,3) = 0; 
        }
      }
    }





    void grad2bmatrix (Math::Matrix &bmat, const Math::Matrix &grad)
    {
      bmat.allocate (grad.rows(),7);
      for (guint i = 0; i < grad.rows(); i++) {
        bmat(i,0) = grad(i,3) * grad(i,0)*grad(i,0);
        bmat(i,1) = grad(i,3) * grad(i,1)*grad(i,1);
        bmat(i,2) = grad(i,3) * grad(i,2)*grad(i,2);
        bmat(i,3) = grad(i,3) * 2*grad(i,0)*grad(i,1);
        bmat(i,4) = grad(i,3) * 2*grad(i,0)*grad(i,2);
        bmat(i,5) = grad(i,3) * 2*grad(i,1)*grad(i,2);
        bmat(i,6) = -1.0;
      }
    }







    void guess_DW_directions (std::vector<int>& dwi, std::vector<int>& bzero, const Math::Matrix& grad)
    {
      if (grad.columns() != 4)
        throw Exception ("invalid gradient encoding matrix: expecting 4 columns.");
      Shells shells(grad);
      int shell_count = shells.count();
      if (shell_count < 1 || shell_count > sqrt(grad.rows()))
        throw Exception ("Gradient encoding matrix does not represent a HARDI sequence!");
      info ("found " + str (shell_count) + " shells");
      Shell bzeroShell;
      Shell dwiShell;
      if (shell_count>1)
        bzeroShell = shells.first();
      dwiShell = shells.last();
      if (shell_count>1)
        info ("using " + str (bzeroShell.count()) + " volumes with b-value " + str (bzeroShell.avg_bval()) + " +/-" + str (bzeroShell.std_bval()) + " as b=0 volumes");
      info ("using " + str (dwiShell.count()) + " volumes with b-value " + str (dwiShell.avg_bval()) + " +/-" + str (dwiShell.std_bval()) + " as diffusion-weighted volumes");
      bzero = bzeroShell.idx();
      dwi = dwiShell.idx();
    }




    void gen_direction_matrix (Math::Matrix& dirs, const Math::Matrix& grad, const std::vector<int>& dwi)
    {
      dirs.allocate (dwi.size(), 2);
      for (guint i = 0; i < dwi.size(); i++) {
        double norm = Point (grad(dwi[i],0), grad(dwi[i],1), grad(dwi[i],2)).norm();

        double z = grad(dwi[i],2)/norm;
        if (z >= 1.0) z = 1.0;
        else if (z <= -1.0) z = -1.0;

        dirs(i,0) = atan2 (grad(dwi[i],1), grad(dwi[i],0));
        dirs(i,1) = acos (z);
      }
    }


  }
}
