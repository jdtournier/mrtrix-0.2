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

#include "histogram.h"
#include "min_max.h"
#include "image/position.h"

namespace MR {

  Histogram::Histogram (Image::Position& ima, guint num_buckets)
  {
    if (num_buckets < 10) throw Exception ("Error initialising Histogram: number of buckets must be greater than 10");

    info ("Initialising histogram with " + str(num_buckets) + " buckets...");
    list.resize (num_buckets);

    float min, max;
    get_min_max (ima, min, max);

    for (guint n = 0; n < list.size(); n++)
      list[n].value = min + ( (max-min) * (n + 0.5) / list.size() );

    ProgressBar::init (ima.voxel_count(), "building histogram...");

    do {
      float val = ima.re();

      if (gsl_finite(val)) { 
        guint pos = (guint) ( list.size() * ((val - min)/(max-min)) );
        if (pos >= list.size()) pos = list.size()-1;
        list[pos].frequency++;
      }
      if (ima.is_complex()) {
        val = ima.im();
        if (gsl_finite(val)) {
          guint pos = (guint) ( list.size() * (val - min)/(max-min) );
          if (pos >= list.size()) pos = list.size()-1;
          list[pos].frequency++;
        }
      }

      ProgressBar::inc();
    } while (ima++);

    ProgressBar::done();

  }







  float Histogram::first_min() const
  {
    guint first_peak = 0;
    guint first_peak_index = 0;
    guint second_peak = 0;
    guint second_peak_index = 0;
    guint first_minimum = 0;
    guint first_min_index = 0;
    guint range_step = list.size()/20;
    guint range = list.size()/20;
    guint index;

    for (index = 0; index < range; index++) {
      if (list[index].frequency > first_peak) {
        first_peak = list[index].frequency;
        first_minimum = first_peak;
        first_min_index = first_peak_index = index;
      }
    }

    range = first_peak_index + range_step;

    for (index = first_peak_index; index < range; index++) {
      if (list[index].frequency < first_minimum) {
        first_minimum = list[index].frequency;
        first_min_index = second_peak_index = index;
      }
    }

    range = first_min_index + range_step;

    for (index = first_min_index; index < range; index++) {
      if (list[index].frequency > second_peak) {
        second_peak = list[index].frequency;
        second_peak_index = index;
      }
    }

    return (list[first_min_index].value);
  }


}
