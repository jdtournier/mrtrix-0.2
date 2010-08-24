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

namespace MR {
  namespace DWI {

    void normalise_grad (Math::Matrix &grad)
    {
      if (grad.columns() != 4) throw Exception ("invalid gradient matrix dimensions");

      double norm;
      for (guint i = 0; i < grad.rows(); i++) {
        norm = grad(i,3) ? 1.0/sqrt(grad(i,0)*grad(i,0)+grad(i,1)*grad(i,1)+grad(i,2)*grad(i,2)) : 0.0;
        grad(i,0) *= norm;
        grad(i,1) *= norm;
        grad(i,2) *= norm;
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
      if (grad.columns() != 4) throw Exception ("invalid gradient encoding matrix: expecting 4 columns.");

      dwi.clear();
      bzero.clear();

      for (int i = 0; i < (int) grad.rows(); i++) {
        if (grad(i,3)) dwi.push_back (i);
        else bzero.push_back (i);
      }
    }




    void gen_direction_matrix (Math::Matrix& dirs, const Math::Matrix& grad, const std::vector<int>& dwi)
    {
      dirs.allocate (dwi.size(), 2);
      for (guint i = 0; i < dwi.size(); i++) {
        double norm = Point (grad(dwi[i],0), grad(dwi[i],1), grad(dwi[i],2)).norm();
            /*sqrt (
            grad(dwi[i],0) * grad(dwi[i],0) + 
            grad(dwi[i],1) * grad(dwi[i],1) + 
            grad(dwi[i],2) * grad(dwi[i],2) );*/
        dirs(i,0) = atan2 (grad(dwi[i],1), grad(dwi[i],0));
        dirs(i,1) = acos (grad(dwi[i],2)/norm);
      }
    }


  }
}
