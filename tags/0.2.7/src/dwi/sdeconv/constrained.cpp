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

#include "dwi/sdeconv/constrained.h"
#include "dwi/gradient.h"
#include "math/linalg.h"

#define TOLERANCE 1e-4
#define NUM_IT 100
#define STEP_SIZE 0.1
#define SEARCH_STEP 50

namespace MR {
  namespace DWI {
    namespace SH {

      CSDeconv::Common::Common (const Math::Vector& response, const Math::Vector& init_filter, const Math::Matrix& DW_dirs, const Math::Matrix& HR_dirs, int lmax) :
        lambda (1.0),
        threshold (0.1)
      {
        int lmax_data = (response.size()-1)*2;
        int n = LforN (DW_dirs.rows());
        if (lmax_data > n) lmax_data = n;
        if (lmax_data > lmax) lmax_data = lmax;
        info ("calculating even spherical harmonic components up to order " + str(lmax_data) + " for initialisation");

        if (init_filter.size() < (guint) (lmax_data/2)+1) 
          throw Exception ("not enough initial filter coefficients supplied for lmax = " + str (lmax_data));

        Math::Vector RH;
        SH2RH (RH, response);

        init_transform (fconv, DW_dirs, lmax_data);
        Math::invert (rconv, fconv);

        int l = 0;
        guint nl = 1;
        for (guint row = 0; row < rconv.rows(); row++) {
          if (row >= nl) { l++; nl = NforL (2*l); }
          for (guint col = 0; col < rconv.columns(); col++) {
            rconv(row, col) *= init_filter[l] / RH[l];
            fconv(col,row) *= RH[l];
          }
        }

        init_transform (HR_trans, HR_dirs, lmax);
        HR_trans.multiply(((float) fconv.rows())*response[0]/((float) HR_trans.rows()));

        info ("constrained spherical deconvolution initiated successfully");
      }




      CSDeconv::CSDeconv (const CSDeconv::Common& common) : P (common) 
      {
        M2.allocate (P.fconv.rows() + P.HR_trans.rows(), P.HR_trans.columns());
        for (guint row = 0; row < P.fconv.rows(); row++) {
          for (guint col = 0; col < P.fconv.columns(); col++)
            M2(row,col) = P.fconv(row,col);
          for (guint col = P.fconv.columns(); col < P.HR_trans.columns(); col++)
            M2(row,col) = 0.0;
        }

        S.allocate (P.fconv.rows());
        S_padded.zero (M2.rows());
        init_F.allocate (P.fconv.columns());
        F.allocate (P.HR_trans.columns());
      }




/*
      void CSDeconv::init (Math::Vector& response, Math::Vector& init_filter, Math::Matrix& DW_dirs, Math::Matrix& HR_dirs, int lmax)
      {
        int lmax_data = (response.size()-1)*2;
        int n = LforN (DW_dirs.rows());
        if (lmax_data > n) lmax_data = n;
        if (lmax_data > lmax) lmax_data = lmax;
        info ("calculating even spherical harmonic components up to order " + str(lmax_data) + " for initialisation");

        if (init_filter.size() < (guint) (lmax_data/2)+1) 
          throw Exception ("not enough initial filter coefficients supplied for lmax = " + str (lmax_data));

        Math::Vector RH;
        SH2RH (RH, response);

        init_transform (fconv, DW_dirs, lmax_data);
        Math::invert (rconv, fconv);

        int l = 0;
        guint nl = 1;
        for (guint row = 0; row < rconv.rows(); row++) {
          if (row >= nl) { l++; nl = NforL (2*l); }
          for (guint col = 0; col < rconv.columns(); col++) {
            rconv(row, col) *= init_filter[l] / RH[l];
            fconv(col,row) *= RH[l];
          }
        }

        init_transform (HR_trans, HR_dirs, lmax);
        HR_trans.multiply(((float) fconv.rows())*response[0]/((float) HR_trans.rows()));

        M2.allocate (fconv.rows() + HR_trans.rows(), HR_trans.columns());
        for (guint row = 0; row < fconv.rows(); row++) {
          for (guint col = 0; col < fconv.columns(); col++)
            M2(row,col) = fconv(row,col);
          for (guint col = fconv.columns(); col < HR_trans.columns(); col++)
            M2(row,col) = 0.0;
        }

        S.allocate (DW_dirs.rows());
        S_padded.zero (M2.rows());
        init_F.allocate (fconv.columns());
        F.allocate (HR_trans.columns());
        lambda = 1.0;
        threshold = 0.0;

        info ("constrained spherical deconvolution initiated successfully");
      }
*/




      void CSDeconv::set (const Math::Vector& DW_signals)
      {
        guint n;
        for (n = 0; n < S.size(); n++)
          S[n] = S_padded[n] = DW_signals[n];

        init_F.multiply (P.rconv, S);

        for (n = 0; n < init_F.size(); n++) F[n] = init_F[n];
        for (; n < F.size(); n++) F[n] = 0.0;

        HR_amps.multiply (P.HR_trans, F);
        threshold = P.threshold * HR_amps.mean();
      }






      bool CSDeconv::iterate()
      {
        neg.clear();
        HR_amps.multiply (P.HR_trans, F);
        for (guint n = 0; n < HR_amps.size(); n++)
          if (HR_amps[n] < threshold)
            neg.push_back (n);

        if (P.fconv.rows() + neg.size() < P.HR_trans.columns()) {
          error ("not enough negative directions! failed to converge.");
          F.set_all (GSL_NAN);
          return (true);
        }

        Math::Matrix M (P.fconv.rows()+neg.size(), P.HR_trans.columns());
        for (guint n = 0; n < P.fconv.columns(); n++) 
          for (guint m = 0; m < P.fconv.rows(); m++) 
            M(m,n) = P.fconv(m,n);

        for (guint n = P.fconv.columns(); n < P.HR_trans.columns(); n++) 
          for (guint m = 0; m < P.fconv.rows(); m++) 
            M(m,n) = 0.0;

        for (guint n = 0; n < P.HR_trans.columns(); n++)
          for (guint m = 0; m < neg.size(); m++)
            M(m+P.fconv.rows(),n) = P.lambda * P.HR_trans(neg[m],n);      

        vec.zero (M.rows());
        for (guint n = 0; n < S.size(); n++)
          vec[n] = S[n];

        HR_amps.copy (F);
        Math::QR_LS_solve (M, vec, F, buf);
        HR_amps.sub (F);

        return (HR_amps.norm2() == 0.0);
      }





    }
  }
}


