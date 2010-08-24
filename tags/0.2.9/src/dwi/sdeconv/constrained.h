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

/** \file
  This file describes a number of functions related to performing
  Monte Carlo Markhov Chain simulations of regularised spherical deconvolution.
 */

#ifndef __dwi_SH_sdeconv_constrained_h__
#define __dwi_SH_sdeconv_constrained_h__

#include "dwi/SH.h"

namespace MR {
  namespace DWI {
    namespace SH {

      class CSDeconv
      {
        public:
          class Common 
          {
            public:
              Common (const Math::Vector& response, const Math::Vector& init_filter, const Math::Matrix& DW_dirs, const Math::Matrix& HR_dirs, int lmax = 8);
              Math::Matrix   fconv, rconv, HR_trans;
              double         lambda, threshold;
          };



          CSDeconv (const Common& common);
          ~CSDeconv() { }

          void      set (const Math::Vector& DW_signals);
          bool      iterate();

          const Math::Vector& FOD () const         { return (F); }
          const Math::Vector& signals () const     { return (S); }


        protected:
          const Common&      P;
          double             threshold;
          Math::Matrix       M2;
          Math::Vector       S, F, init_F, S_padded, HR_amps, buf, vec;
          std::vector<int>   neg;

      };

    }
  }
}

#endif


