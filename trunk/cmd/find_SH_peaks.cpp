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


    03-10-2008 J-Donald Tournier <d.tournier@brain.org.au>
    * fix bug in looping structure to allow processing of whole data set.

*/

#include "app.h"
#include "image/position.h"
#include "dwi/SH.h"

#define DOT_THRESHOLD 0.99

using namespace std; 
using namespace MR; 

SET_VERSION_DEFAULT;

DESCRIPTION = {
  "identify the orientations of the N largest peaks of a SH profile",
  NULL
};

ARGUMENTS = {
  Argument ("SH", "SH coefficients image", "the input image of SH coefficients.").type_image_in (),
  Argument ("directions", "direction set", "the set of directions to use as seeds for the peak finding").type_file (),
  Argument ("ouput", "output image", "the output image. Each volume corresponds to the x, y & z component of each peak direction vector in turn.").type_image_out (),
  Argument::End
};

OPTIONS = {
  Option ("num", "number of peaks", "the number of peaks to extract (default is 3).")
    .append (Argument ("peaks", "number", "the number of peaks").type_integer (0, INT_MAX, 3)),

  Option ("direction", "specify direction", "the direction of a peak to estimate. The algorithm will attempt to find the same number of peaks as have been specified using this option.", false, true)
    .append (Argument ("phi", "azimuthal angle", "the azimuthal angle of the direction (in degrees).").type_float (GSL_NEGINF, GSL_POSINF, 0.0))
    .append (Argument ("theta", "elevation angle", "the elevation angle of the direction (in degrees, from the vertical z-axis).").type_float (GSL_NEGINF, GSL_POSINF, 0.0)),

  Option ("peaks", "true peaks image", "the program will try to find the peaks that most closely match those in the image provided.")
    .append (Argument ("image", "peaks image", "an image containing the true peaks to be estimated.").type_image_in ()),

  Option ("threshold", "amplitude threshold", "only peak amplitudes greater than the threshold will be considered.")
    .append (Argument ("value", "value", "the threshold value").type_float (GSL_NEGINF, GSL_POSINF, 0.0)),

  Option::End
};



class Direction {
  public:
    Direction () : a (GSL_NAN) { }
    Direction (const Direction& d) : a (d.a), v (d.v) { }
    Direction (float phi, float theta) : a (1.0), v (cos(phi)*sin(theta), sin(phi)*sin(theta), cos(theta)) { }
    float a;
    Point v;
    bool operator<(const Direction& d) const { return (a > d.a); }
};

EXECUTE {

  // Load direction set:
  Math::Matrix dirs;
  dirs.load (argument[1].get_string());
  if (dirs.columns() != 2) 
    throw ("expecting 2 columns for search directions matrix");
  
  std::vector<OptBase> opt = get_options (0); // num
  int npeaks = opt.size() ? opt[0][0].get_int() : 3;

  opt = get_options (1); // individual peaks 
  std::vector<Direction> true_peaks;
  for (guint n = 0; n < opt.size(); n++) {
    Direction p (M_PI*opt[n][0].get_float()/180.0, M_PI*opt[n][1].get_float()/180.0);
    true_peaks.push_back (p);
  }
  if (true_peaks.size()) npeaks = true_peaks.size();

  opt = get_options (3); // threshold
  float threshold = GSL_NEGINF;
  if (opt.size()) threshold = opt[0][0].get_float();

  Image::Object &SH_obj (*argument[0].get_image());
  Image::Header header (SH_obj);

  header.data_type = DataType::Float32;
  header.axes.set_ndim (4);

  Ptr<Image::Position> ipeaks;
  opt = get_options (2); // peaks image
  if (opt.size()) {
    if (true_peaks.size()) throw Exception ("you can't specify both a peaks file and orientations to be estimated at the same time");
    ipeaks = new Image::Position (*opt[0][0].get_image());
    if (ipeaks->dim(0) != header.dim(0) || ipeaks->dim(1) != header.dim(1) || ipeaks->dim(2) != header.dim(2))
      throw Exception ("dimensions of peaks image \"" + ipeaks->name() + "\" do not match that of SH coefficients image \"" + SH_obj.name() + "\"");
    npeaks = ipeaks->dim(3) / 3;
  }

  header.axes.dim[3] = 3 * npeaks;

  Image::Position SH (SH_obj);
  Image::Position out (*argument[2].get_image (header));

  std::vector<Direction> peaks_out (npeaks);
  std::vector<float> val (SH.dim(3));
  int lmax = DWI::SH::LforN (SH.dim(3));
  DWI::SH::precompute (lmax, 512);
  
 
  info ("using lmax = " + str (lmax));

  ProgressBar::init (SH.dim(0)*SH.dim(1)*SH.dim(2), "finding orientations of largest peaks...");

  for (SH.set(2,0), out.set(2,0); SH[2] < SH.dim(2); SH.inc(2), out.inc(2)) {
    if (ipeaks.get()) ipeaks->set(2,SH[2]);

    for (SH.set(1,0), out.set(1,0); SH[1] < SH.dim(1); SH.inc(1), out.inc(1)) {
      if (ipeaks.get()) ipeaks->set(1,SH[1]);

      for (SH.set(0,0), out.set(0,0); SH[0] < SH.dim(0); SH.inc(0), out.inc(0)) {

        bool skip = false;
        if (ipeaks.get()) {
          ipeaks->set(0, SH[0]);
          if (gsl_isnan (ipeaks->value())) skip = true;
        }
        
        if (!skip) {
          float min = GSL_POSINF, max = GSL_NEGINF;
          for (SH.set(3,0); SH[3] < SH.dim(3); SH.inc(3)) {
            val[SH[3]] = SH.value();
            if (gsl_isnan (val[SH[3]])) {
              skip = true;
              break;
            }
            if (val[SH[3]] < min) min = val[SH[3]];
            if (val[SH[3]] > max) max = val[SH[3]];
          }
          if (min == max) skip = true;
        }

        if (skip) for (out.set(3,0); out[3] < out.dim(3); out.inc(3)) out.value (GSL_NAN);
        else {
          std::vector<Direction> all_peaks;
          for (guint i = 0; i < dirs.rows(); i++) {
            Direction p (dirs(i,0), dirs(i,1)); 
            p.a = DWI::SH::get_peak (&val[0], lmax, p.v, true);
            
            if (gsl_finite (p.a)) {
              for (guint j = 0; j < all_peaks.size(); j++) {
                if (fabs (p.v.dot (all_peaks[j].v)) > DOT_THRESHOLD) {
                  p.a = NAN;
                  break;
                }
              }
            }
            if (gsl_finite (p.a) && p.a >= threshold) all_peaks.push_back (p);
          }
          
          if (ipeaks.get()) {
            for (int i = 0; i < npeaks; i++) {
              Point p;
              ipeaks->set(3, 3*i);
              for (int n = 0; n < 3; n++) { p[n] = ipeaks->value(); ipeaks->inc(3); }
              p.normalise();

              float mdot = 0.0;
              for (guint n = 0; n < all_peaks.size(); n++) {
                float f = fabs (p.dot (all_peaks[n].v));
                if (f > mdot) { 
                  mdot = f; 
                  peaks_out[i] = all_peaks[n];
                }
              }
            }
          }
          else if (true_peaks.size()) {
            for (int i = 0; i < npeaks; i++) {
              float mdot = 0.0;
              for (guint n = 0; n < all_peaks.size(); n++) {
                float f = fabs (all_peaks[n].v.dot (true_peaks[i].v));
                if (f > mdot) { 
                  mdot = f; 
                  peaks_out[i] = all_peaks[n];
                }
              }
            }
          }
          else std::partial_sort_copy (all_peaks.begin(), all_peaks.end(), peaks_out.begin(), peaks_out.end());



          
          int actual_npeaks = MIN (npeaks, (int) all_peaks.size());
          out.set (3, 0);
          for (int n = 0; n < actual_npeaks; n++) {
             out.value (peaks_out[n].a*peaks_out[n].v[0]); out.inc(3); 
             out.value (peaks_out[n].a*peaks_out[n].v[1]); out.inc(3); 
             out.value (peaks_out[n].a*peaks_out[n].v[2]); out.inc(3);
          }
          for (; out[3] < 3*npeaks; out.inc(3)) out.value (GSL_NAN);
        }

        ProgressBar::inc();
      }
    }
  }
  ProgressBar::done();
}

