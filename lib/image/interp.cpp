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

#include "image/interp.h"
#include "math/matrix.h"

namespace MR {
  namespace Image {

    Interp::Interp (Object& parent) : 
      Position (parent)
    { 
      bounds[0] = dim(0) - 0.5;
      bounds[1] = dim(1) - 0.5;
      bounds[2] = dim(2) - 0.5;

      out_of_bounds = true;

      
      PR[0][0] = image.P2R()(0,0); 
      PR[0][1] = image.P2R()(0,1); 
      PR[0][2] = image.P2R()(0,2); 

      PR[0][3] = image.P2R()(0,3);
      PR[1][0] = image.P2R()(1,0); 
      PR[1][1] = image.P2R()(1,1); 

      PR[1][2] = image.P2R()(1,2); 
      PR[1][3] = image.P2R()(1,3);
      PR[2][0] = image.P2R()(2,0);

      PR[2][1] = image.P2R()(2,1); 
      PR[2][2] = image.P2R()(2,2); 
      PR[2][3] = image.P2R()(2,3);


      RP[0][0] = image.R2P()(0,0); 
      RP[0][1] = image.R2P()(0,1); 
      RP[0][2] = image.R2P()(0,2); 

      RP[0][3] = image.R2P()(0,3);
      RP[1][0] = image.R2P()(1,0); 
      RP[1][1] = image.R2P()(1,1); 

      RP[1][2] = image.R2P()(1,2); 
      RP[1][3] = image.R2P()(1,3);
      RP[2][0] = image.R2P()(2,0);

      RP[2][1] = image.R2P()(2,1); 
      RP[2][2] = image.R2P()(2,2); 
      RP[2][3] = image.R2P()(2,3);
    }

  }
}
