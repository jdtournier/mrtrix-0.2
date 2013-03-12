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



    21-07-2010 J-Donald Tournier <d.tournier@brain.org.au>
    * improved SH::delta() function

*/

#ifndef __spherical_deconvolution_h__
#define __spherical_deconvolution_h__

#include "point.h"
#include "math/linalg.h"

namespace MR {
  namespace DWI {
    namespace SH {

      inline int NforL (int lmax) { return ((lmax+1)*(lmax+2)/2); }
      inline int NforL_mpos (int lmax) { return ((lmax/2+1)*(lmax/2+1)); }

      inline int index (int l, int m) { return (l*(l+1)/2 + m); }
      inline int index_mpos (int l, int m) { return (l*l/4 + m); }

      inline int LforN (int N) { return (2*(((int) (sqrt((float) (1+8*N)))-3)/4)); }

      class Coefs {
        public:
          Coefs (int lmax = 8) : V (NforL (lmax)), _lmax (lmax) { }

          Math::Vector V;

          void        lmax (int l)                            { _lmax = l; V.allocate (NforL (l)); }
          int         lmax ()                                 { return (_lmax); }
          double&     operator() (int l, int m)               { return (V[index (l,m)]); }

        protected:
          int _lmax;
      };



      float value (float azimuth, float elevation, int l, int m);
      float value (Coefs &SH, float azimuth, float elevation, int lmax);
      float value (Coefs &SH, const Point& unit_dir);
      float value (const float *values, const Point& unit_dir, int lmax);

      void precompute (int lmax, int num = 256);
      float value_precomputed (const float *values, const Point& unit_dir);

      void delta (Coefs& SH, float azimuth, float elevation, int lmax);

      void init_transform (Math::Matrix& SHT, const Math::Matrix& dirs, int lmax);

      class Transform {
        public:
          Transform (const Math::Matrix& dirs, int lmax)                { init_transform (SHT, dirs, lmax); Math::invert (iSHT, SHT); }
          void set_filter (const Math::Vector& filter)
          { 
            int l = 0;
            guint nl = 1;
            for (guint n = 0; n < iSHT.rows(); n++) {
              if (n >= nl) { l++; nl = NforL(2*l); }
              for (guint i = 0; i < iSHT.columns(); i++) {
                iSHT(n,i) *= filter[l];
                SHT(i,n) = filter[l] == 0.0 ? 0.0 : SHT(i,n)/filter[l];
              }
            }
          }
          void A2SH (Math::Vector& SH, const Math::Vector& amplitudes)  { SH.multiply (iSHT, amplitudes); }
          void A2SH (Coefs& SH, const Math::Vector& amplitudes)         { A2SH (SH.V, amplitudes); }
          void SH2A (Math::Vector& amplitudes, const Math::Vector& SH)  { amplitudes.multiply (SHT, SH); }
          void SH2A (Math::Vector& amplitudes, const Coefs& SH)         { SH2A (amplitudes, SH.V); }

          int n_SH () const { return (SHT.columns()); }
          int n_amp () const  { return (SHT.rows()); }

          const Math::Matrix& mat_A2SH () const { return (iSHT); }
          const Math::Matrix& mat_SH2A () const { return (SHT); }

        protected:
          Math::Matrix SHT, iSHT;
      };


      void FA2SH (Coefs& SH, float FA, float ADC, float bvalue, int lmax, int precision = 100);
      void SH2RH (Math::Vector& RH, const Math::Vector& SH);

      float get_peak (const float* SH, int lmax, Point& unit_init_dir, bool precomputed = false);

      void derivatives (
          const float *SH,
          int   lmax,
          float elevation,
          float azimuth,
          float &amplitude,
          float &dSH_del,
          float &dSH_daz,
          float &d2SH_del2,
          float &d2SH_deldaz,
          float &d2SH_daz2,
          bool  precomputed = false);

    }
  }
}

#endif
