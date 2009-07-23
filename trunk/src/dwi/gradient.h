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

#ifndef __dwi_gradient_h__
#define __dwi_gradient_h__

namespace MR {
  namespace Math { class Matrix; }
  namespace DWI {

    void normalise_grad (Math::Matrix &grad);
    void grad2bmatrix (Math::Matrix &bmat, const Math::Matrix &grad);
    void guess_DW_directions (std::vector<int>& dwi, std::vector<int>& bzero, const Math::Matrix& grad);
    void gen_direction_matrix (Math::Matrix& dirs, const Math::Matrix& grad, const std::vector<int>& dwi);

  }
}

#endif

