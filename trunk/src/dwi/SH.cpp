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


    29-10-2008 J-Donald Tournier <d.tournier@brain.org.au>
    * fix precomputed value calculation to handle rounding errors in the angle

    18-12-2008 J-Donald Tournier <d.tournier@brain.org.au>
    * modify precomputation to allow thread-safe operation
    
    11-09-2009 J-Donald Tournier <d.tournier@brain.org.au>
    * fix bug in computation of second derivative of SH series (pointed out by
    * Ben Jeurissen).
    
*/

#include <gsl/gsl_sf_legendre.h>
#include "image/position.h"
#include "image/interp.h"
#include "math/matrix.h"
#include "dwi/SH.h"

#define MAX_DIR_CHANGE 0.2
#define ANGLE_TOLERANCE 1e-4



namespace MR {
  namespace DWI {
    namespace SH {

      namespace {

        int       lmax_legendre = 0;
        int       num_legendre_coefs = 0;
        int       num_legendre_dirs = 0;
        float     *precomp_legendre = NULL;
        float     inc_legendre = 0.0;

        class PrecomputedFraction {
          public:
            PrecomputedFraction () : f1 (0.0), f2 (0.0), p1 (NULL), p2 (NULL) { }
            float     f1, f2;
            float     *p1, *p2;
        };


        inline void calc_index_fractions (PrecomputedFraction& f, float elevation)
        {
          f.f2 = elevation / inc_legendre;
          int index = (int) f.f2;
          if (index < 0) { index = 0; f.f1 = 1.0; f.f2 = 0.0; }
          else if (index >= num_legendre_dirs-1) { index = num_legendre_dirs-1; f.f1 = 1.0; f.f2 = 0.0; }
          else { f.f2 -= index; f.f1 = 1.0 - f.f2; }

          f.p1 = precomp_legendre + index*num_legendre_coefs;
          f.p2 = f.p1 + num_legendre_coefs;

          assert (f.p1 >= precomp_legendre);
          assert (f.f2 == 0.0 || f.p2 < precomp_legendre + num_legendre_coefs*num_legendre_dirs);
        }




        inline float legendre_precomputed (const PrecomputedFraction& f, int l, int m)
        {
          int i (index_mpos (l,m));
          float retval = f.f1*f.p1[i];
          if (f.f2) retval += f.f2*f.p2[i];
          return (retval);
        }


      }




      void delta (Coefs& SH, const Math::Matrix& dirs, float azimuth, float elevation, int lmax)
      {
        SH.lmax (lmax);
        SH.V.zero();
        for (int n = 0; n <= lmax; n+=2)
          SH (n,0) = value (0.0, 0.0, n, 0);

        Math::Matrix cart(3, dirs.rows());
        for (guint n = 0; n < dirs.rows(); n++) {
          double d = sin (dirs(n,1));
          cart(0,n) = d*cos (dirs(n,0));
          cart(1,n) = d*sin (dirs(n,0));
          cart(2,n) = cos (dirs(n,1));
        }

        Math::Matrix A(3,3), B(3,3);
        A.identity();
        A(0,0) = cos(azimuth);
        A(0,1) = -sin(azimuth);
        A(1,0) = -A(0,1);
        A(1,1) = A(0,0);

        B.identity();
        B(0,0) = cos(elevation);
        B(0,2) = sin(elevation);
        B(2,0) = -B(0,2);
        B(2,2) = B(0,0);

        Math::Matrix R;
        R.multiply (A, B);
        B.multiply (R, cart);

        A.allocate (dirs);
        for (guint n = 0; n < dirs.rows(); n++) {
          A(n,1) = acos (B(2,n));
          A(n,0) = atan2 (B(1,n), B(0,n));
        }

        Transform rotated (A, lmax);
        Transform non_rotated (dirs, lmax);

        Math::Vector tmp (NforL (lmax));
        rotated.SH2A (tmp, SH);
        non_rotated.A2SH (SH, tmp);
      }





      float value (Coefs& SH, float azimuth, float elevation, int lmax)
      {
        float cel = cos (elevation);
        float val = 0.0;
        for (int l = 0; l <= lmax; l+=2)
          val += SH(l,0) * gsl_sf_legendre_sphPlm (l, 0, cel);

        for (int m = 1; m <= lmax; m++) {
          float caz = cos (m*azimuth);
          float saz = sin (m*azimuth);
          for (int l = 2*((m+1)/2); l <= lmax; l+=2) {
            float buf = gsl_sf_legendre_sphPlm (l, m, cel);
            val += SH (l,m)  * buf * caz;
            val += SH (l,-m) * buf * saz;
          }
        }

        return (val);
      }




      float value (float azimuth, float elevation, int l, int m)
      {
        elevation = gsl_sf_legendre_sphPlm (l, abs(m), cos(elevation));
        if (!m) return (elevation);
        if (m > 0) return (elevation * cos (m*azimuth));
        return (elevation * sin (-m*azimuth));
      }





      float value (Coefs& SH, const Point& unit_dir)
      {
        int lmax = SH.lmax();
        float az = atan2 (unit_dir[1], unit_dir[0]);

        float val = 0.0;
        for (int l = 0; l <= lmax; l+=2)
          val += SH(l,0) * gsl_sf_legendre_sphPlm (l, 0, unit_dir[2]);

        for (int m = 1; m <= lmax; m++) {
          float caz = cos(m*az);
          float saz = sin(m*az);
          for (int l = 2*((m+1)/2); l <= lmax; l+=2) {
            float buf = gsl_sf_legendre_sphPlm (l, m, unit_dir[2]);
            val += SH(l,m)*buf*caz;
            val += SH(l,-m)*buf*saz;
          }
        }

        return (val);
      }






      float value (const float *values, const Point& unit_dir, int lmax)
      {
        float az = atan2(unit_dir[1], unit_dir[0]);

        float val = 0.0;
        for (int l = 0; l <= lmax; l+=2)
          val += values[index(l,0)] * gsl_sf_legendre_sphPlm (l, 0, unit_dir[2]);

        for (int m = 1; m <= lmax; m++) {
          float caz = cos(m*az);
          float saz = sin(m*az);
          for (int l = 2*((m+1)/2); l <= lmax; l+=2) {
            float buf = gsl_sf_legendre_sphPlm (l, m, unit_dir[2]);
            val += values[index(l,m)]*buf*caz;
            val += values[index(l,-m)]*buf*saz;
          }
        }

        return (val);
      }








      void init_transform (Math::Matrix& SHT, const Math::Matrix& dirs, int lmax)
      {
        if (dirs.columns() != 2) throw Exception ("direction matrix should have 2 columns: [ azimuth elevation ]");
        SHT.allocate (dirs.rows(), NforL (lmax));

        for (int l = 0; l <= lmax; l+=2) {
          for (int m = 0; m <= l; m++) {
            for (guint i = 0; i < dirs.rows(); i++) {
              double s = gsl_sf_legendre_sphPlm (l, m, cos(dirs(i, 1)));
              if (m) {
                SHT(i,index(l, m)) = s*cos(m*dirs(i, 0));
                SHT(i,index(l, -m)) = s*sin(m*dirs(i, 0));
              }
              else SHT(i,index(l, 0)) = s;
            }
          }
        }
      }





      void FA2SH (Coefs& SH, float FA, float ADC, float bvalue, int lmax, int precision)
      {
        float a = FA/sqrt(3.0 - 2.0*FA*FA);
        float ev1 = ADC*(1.0+2.0*a), ev2 = ADC*(1.0-a);

        Math::Vector sigs (precision);
        Math::Matrix SHT (precision, lmax/2+1);

        for (int i = 0; i < precision; i++) {
          float el = i*M_PI/(2.0*(precision-1));
          sigs[i] = exp(-bvalue*(ev1*cos(el)*cos(el) + ev2*sin(el)*sin(el)));
          for (int l = 0; l < lmax/2+1; l++)
            SHT(i,l) = value (0.0, el, 2*l, 0);
        }

        Math::Matrix SHinv (SHT.columns(), SHT.rows());
        Math::invert (SHinv, SHT);
        SH.V.multiply (SHinv, sigs);
      }




      void SH2RH (Math::Vector& RH, const Math::Vector& SH)
      {
        RH.allocate (SH.size());
        for (guint l = 0; l < SH.size(); l++)
          RH[l] = SH[l]/value (0.0, 0.0, 2*l, 0);
      }





      void precompute (int lmax, int num)
      {
        if (precomp_legendre) return;
        lmax_legendre = lmax;
        num_legendre_coefs = NforL_mpos (lmax_legendre);
        num_legendre_dirs = num;
        precomp_legendre = new float [num_legendre_coefs*num_legendre_dirs];
        inc_legendre = M_PI/(num_legendre_dirs-1);

        for (int n = 0; n < num_legendre_dirs; n++) {
          float* p = precomp_legendre + n*num_legendre_coefs;
          float cos_el = cos (n*inc_legendre);

          for (int l = 0; l <= lmax_legendre; l+=2) 
            for (int m = 0; m <= l; m++) 
              p[index_mpos(l,m)] = gsl_sf_legendre_sphPlm (l, m, cos_el);
        }
      }







      float value_precomputed (const float *values, const Point& unit_dir)
      {
        PrecomputedFraction f;
        calc_index_fractions (f, acos(unit_dir[2]));

        float az = atan2 (unit_dir[1], unit_dir[0]);
        float val = 0.0;
        for (int l = 0; l <= lmax_legendre; l+=2)
          val += values[index(l,0)]*legendre_precomputed (f, l, 0);

        for (int m = 1; m <= lmax_legendre; m++) {
          float caz = cos(m*az);
          float saz = sin(m*az);
          for (int l = 2*((m+1)/2); l <= lmax_legendre; l+=2) {
            float tmp = legendre_precomputed (f, l, m);
            val += values[index(l,m)] * tmp * caz;
            val += values[index(l,-m)] * tmp * saz;
          }
        }

        return (val);
      }





      float get_peak (const float* SH, int lmax, Point& unit_init_dir, bool precomputed)
      {
        float amplitude, dSH_del, dSH_daz, d2SH_del2, d2SH_deldaz, d2SH_daz2;
        float az = atan2 (unit_init_dir[1], unit_init_dir[0]);
        float el = acos (unit_init_dir[2]);
        float del, daz, dSH_dt, d2SH_dt2, dt;

        for (int i = 0; i < 50; i++) {
          az = atan2 (unit_init_dir[1], unit_init_dir[0]);
          el = acos (unit_init_dir[2]);
          derivatives (SH, lmax, el, az, amplitude, dSH_del, dSH_daz, d2SH_del2, d2SH_deldaz, d2SH_daz2, precomputed);

          del = sqrt (dSH_del*dSH_del + dSH_daz*dSH_daz);
          daz = dSH_daz/del;
          del = dSH_del/del;

          dSH_dt = daz*dSH_daz + del*dSH_del;
          d2SH_dt2 = daz*daz*d2SH_daz2 + 2.0*daz*del*d2SH_deldaz + del*del*d2SH_del2;
          dt = - dSH_dt / d2SH_dt2;

          if (dt < 0.0 || dt > MAX_DIR_CHANGE) dt = MAX_DIR_CHANGE;

          del *= dt;
          daz *= dt;

          unit_init_dir += Point (del*cos(az)*cos(el) - daz*sin(az), del*sin(az)*cos(el) + daz*cos(az), -del*sin(el));
          unit_init_dir.normalise();

          if (dt < ANGLE_TOLERANCE) return (amplitude);
        }

        unit_init_dir.invalidate();
        info ("failed to find SH peak!");
        return (GSL_NAN);
      }







      void derivatives (const float *SH, int lmax, float elevation, float azimuth, float &amplitude,
          float &dSH_del, float &dSH_daz, float &d2SH_del2, float &d2SH_deldaz, float &d2SH_daz2, bool  precomputed)
      {
        float sel = sin(elevation);
        bool atpole = sel < 1e-4;
        float legendre[NforL(lmax)];

        amplitude = dSH_del = dSH_daz = d2SH_del2 = d2SH_deldaz = d2SH_daz2 = 0.0;

        PrecomputedFraction f;
        if (precomputed) calc_index_fractions (f, elevation);
        elevation = cos (elevation);

        for (int l = 0; l <= (int) lmax; l+=2) {
          for (int m = 0; m <= l; m++)
            legendre[index_mpos(l,m)] = precomputed ? legendre_precomputed (f, l, m) : gsl_sf_legendre_sphPlm (l, m, elevation);

          amplitude += SH[index(l,0)] * legendre[index_mpos(l,0)];

          if (l) {
            dSH_del += SH[index(l,0)] * sqrt((float) l*(l+1)) * legendre[index_mpos(l,1)];
            d2SH_del2 += SH[index(l,0)] * (
                sqrt((float) l*(l+1)*(l-1)*(l+2)) * legendre[index_mpos(l,2)]
                - l*(l+1) * legendre[index_mpos(l,0)] )/2.0;
          }
        }

        for (int m = 1; m <= lmax; m++) {
          float caz = cos (m*azimuth);
          float saz = sin (m*azimuth);
          for (int l = 2*((m+1)/2); l <= lmax; l+=2) {
            amplitude += SH[index(l,m)] * legendre[index_mpos(l,m)] * caz;
            amplitude += SH[index(l,-m)] * legendre[index_mpos(l,m)] * saz;

            float tmp = sqrt((float) (l+m)*(l-m+1))*legendre[index_mpos(l,m-1)];
            if (l > m) tmp -= sqrt((float) (l-m)*(l+m+1))*legendre[index_mpos(l,m+1)];
            tmp /= -2.0;
            dSH_del += SH[index(l,m)] * tmp * caz;
            dSH_del += SH[index(l,-m)] * tmp * saz;

            float tmp2 = - ( (l+m)*(l-m+1) + (l-m)*(l+m+1) ) * legendre[index_mpos(l,m)];
            if (m == 1) tmp2 -= sqrt((float) gsl_pow_2((l+1)*l)) * legendre[index_mpos(l,1)];
            else tmp2 += sqrt((float) (l+m)*(l-m+1)*(l+m-1)*(l-m+2)) * legendre[index_mpos(l,m-2)];
            if (l > m+1) tmp2 += sqrt((float) (l-m)*(l+m+1)*(l-m-1)*(l+m+2)) * legendre[index_mpos(l,m+2)];
            tmp2 /= 4.0;
            d2SH_del2 += SH[index(l,m)] * tmp2 * caz;
            d2SH_del2 += SH[index(l,-m)] * tmp2 * saz;

            if (atpole) {
              dSH_daz -= SH[index(l,m)] * tmp * saz;
              dSH_daz += SH[index(l,-m)] * tmp * caz;
            }
            else {
              d2SH_deldaz -= m * SH[index(l,m)] * tmp * saz;
              d2SH_deldaz += m * SH[index(l,-m)] * tmp * caz;

              dSH_daz -= m * SH[index(l,m)] * legendre[index_mpos(l,m)] * saz;
              dSH_daz += m * SH[index(l,-m)] * legendre[index_mpos(l,m)] * caz;

              tmp =  m*m * legendre[index_mpos(l,m)];
              d2SH_daz2 -= SH[index(l,m)] * tmp * caz;
              d2SH_daz2 -= SH[index(l,-m)] * tmp * saz;
            }

          }
        }

        if (!atpole) {
          dSH_daz /= sel;
          d2SH_deldaz /= sel;
          d2SH_daz2 /= sel*sel;
        }
      }


    }
  }
}

