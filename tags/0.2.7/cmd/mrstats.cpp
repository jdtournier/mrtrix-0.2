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


    27-05-2009 J-Donald Tournier <d.tournier@brain.org.au>
    * added option to dump voxel intensities into a text file

    20-05-2009 J-Donald Tournier <d.tournier@brain.org.au>
    * added option to build intensity histogram

*/

#include <fstream>

#include "app.h"
#include "image/position.h"

using namespace std; 
using namespace MR; 

SET_VERSION_DEFAULT;

DESCRIPTION = {
  "compute images statistics.",
  NULL
};

ARGUMENTS = {
  Argument ("image", "image", "the input image from which statistics will be computed.").type_image_in (),
  Argument::End
};


OPTIONS = { 
  Option ("mask", "brain mask", "only perform computation within the specified binary brain mask image.")
    .append (Argument ("image", "image", "the mask image to use.").type_image_in ()),

  Option ("histogram", "generate histogram", "generate histogram of intensities and store in specified text file. Note that the first line of the histogram gives the centre of the bins.")
    .append (Argument ("file", "file", "the text file into which to store the histogram data.").type_file ()),

  Option ("bins", "number of histogram bins", "the number of bins to use to generate the histogram (default = 100).")
    .append (Argument ("num", "number", "the number of bins.").type_integer (2, G_MAXINT, 100)),

  Option ("dump", "dump voxel intensities", "dump the voxel intensities to a text file.")
    .append (Argument ("file", "file", "the text file to dump the values into.").type_file ()),

  Option::End 
};



template <class Function> inline void loop (Image::Position& pos, RefPtr<Image::Position>& mask, Function& func)
{
  if (mask) mask->set(2,0); 
  for (pos.set(2,0); pos[2] < pos.dim(2); pos.inc(2)) {
    if (mask) mask->set(1,0); 
    for (pos.set(1,0); pos[1] < pos.dim(1); pos.inc(1)) {
      if (mask) mask->set(0,0); 
      for (pos.set(0,0); pos[0] < pos.dim(0); pos.inc(0)) {

        bool skip = false;
        if (mask) if (mask->value() < 0.5) skip = true;
        if (!skip) {
          float val = pos.value();
          if (isfinite (val)) func (val);
        }
        if (mask) mask->inc (0);
      }
      if (mask) mask->inc (1);
    }
    if (mask) mask->inc (2);
  }
}


class GetStats {
  public:
    GetStats () : mean (0.0), std (0.0), min (INFINITY), max (-INFINITY), count (0) { }

    double mean, std;
    float min, max;
    size_t count;

    void operator() (float val) {
      mean += val;
      std += val*val;
      if (min > val) min = val;
      if (max < val) max = val;
      count++;
    }
    
    void finalise () {
      if (count == 0) throw Exception ("no voxels in mask - aborting");

      mean /= double(count);
      std = sqrt(std/double(count) - mean*mean);
    }
};


class GetMinMax {
  public:
    GetMinMax () : min (INFINITY), max (-INFINITY) { }

    float min, max;

    void operator() (float val) {
      if (val < min) min = val;
      if (val > max) max = val;
    }
};



class GetHistogram {
  public:
    GetHistogram (const GetMinMax& limits, int nbins) : N (nbins), min (limits.min), width ((limits.max - limits.min) / float (N+1)), data (new size_t [N]) { }
    ~GetHistogram () { delete [] data; }

    int N;
    float min, width;
    size_t* data;

    void clear () { memset (data, 0, N*sizeof(size_t)); }

    void operator() (float val) {
      int bin = int ((val-min) / width);
      if (bin < 0) bin = 0;
      else if (bin >= N) bin = N-1;
      data[bin]++;
    }
};





class DumpValues {
  public:
    DumpValues (const std::string& filename) : stream (filename.c_str()) { }

    std::ofstream stream;

    void operator() (float val) { stream << val << "\n"; }
};





EXECUTE {
  Image::Position ima (*argument[0].get_image());

  RefPtr<Image::Position> mask;
  std::vector<OptBase> opt = get_options (0); // mask
  if (opt.size()) {
    mask = new Image::Position (*opt[0][0].get_image ());
    if (mask->dim(0) != ima.dim(0) || mask->dim(1) != ima.dim(1) || mask->dim(2) != ima.dim(2)) 
      throw Exception ("dimensions of mask image do not match that of data image - aborting");
  }



  opt = get_options (1); // histogram
  if (opt.size()) {
    String filename = opt[0][0].get_string();

    opt = get_options (2); // bins
    size_t bins = 100;
    if (opt.size()) bins = opt[0][0].get_int();

    size_t nrows = 1;
    for (int i = 3; i < ima.ndim(); i++) nrows *= ima.dim(i);

    // calibrate:
    GetMinMax limits;
    ProgressBar::init (nrows, "calibrating...");
    do {
      loop (ima, mask, limits);
      ProgressBar::inc();
    } while (ima++);
    ProgressBar::done();

    // get histogram:
    GetHistogram hist (limits, bins);
    std::ofstream out (filename.c_str());

    // write out bin centres:
    for (int i = 0; i < hist.N; i++)
      out << (hist.min + hist.width/2.0) + i*hist.width << " ";
    out << "\n";

    ProgressBar::init (nrows, "building histogram...");
    do {
      hist.clear();
      loop (ima, mask, hist);
      for (int i = 0; i < hist.N; i++)
        out << hist.data[i] << " ";
      out << "\n";
      ProgressBar::inc();
    } while (ima++);
    ProgressBar::done();

    out.close();
  }
  else {
    bool header_shown = false;
    do {
      GetStats stats;
      loop (ima, mask, stats);
      stats.finalise();

      String s = "[ ";
      for (int n = 3; n < ima.ndim(); n++) s += str(ima[n]) + " ";
      s += "] ";

      if (!header_shown) print ("channel         mean        std. dev.   min         max         count\n");
      header_shown = true;
      print (MR::printf ("%-15s %-11g %-11g %-11g %-11g %-11d\n", s.c_str(), stats.mean, stats.std, stats.min, stats.max, stats.count));
    } while (ima++);
  }

  opt = get_options (3); // dump
  if (opt.size()) {
    DumpValues dump (opt[0][0].get_string());
    ProgressBar::init (ima.voxel_count(), "dumping values to file...");
    do {
      loop (ima, mask, dump);
      ProgressBar::inc();
    } while (ima++);
    ProgressBar::done();
  }
}

