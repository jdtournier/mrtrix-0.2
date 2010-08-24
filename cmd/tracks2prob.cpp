/*
    Copyright 2010 Brain Research Institute, Melbourne, Australia

    Written by Robert E. Smith and J-Donald Tournier, 14/07/10.

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

    04-08-2010 Robert E. Smith <r.smith@brain.org.au>
    * reverted to previous interpretation of voxel values as the number
      of tracks detected as entering each voxel, without normalisation

*/



#include "app.h"
#include "math/matrix.h"
#include "dwi/tractography/file.h"
#include "dwi/tractography/properties.h"

#include <set>



const gchar* data_type_choices[] = { "FLOAT32", "FLOAT32LE", "FLOAT32BE", "FLOAT64", "FLOAT64LE", "FLOAT64BE",
    "INT32", "UINT32", "INT32LE", "UINT32LE", "INT32BE", "UINT32BE",
    "INT16", "UINT16", "INT16LE", "UINT16LE", "INT16BE", "UINT16BE",
    "CFLOAT32", "CFLOAT32LE", "CFLOAT32BE", "CFLOAT64", "CFLOAT64LE", "CFLOAT64BE",
    "INT8", "UINT8", "BIT", NULL };




using namespace MR;



SET_VERSION_DEFAULT;

DESCRIPTION = {
  "convert a tracks file into a map of the fraction of tracks to enter each voxel.",
  NULL
};

ARGUMENTS = {
  Argument ("tracks", "track file", "the input track file.").type_file (),
  Argument ("output", "output image", "the output fraction image").type_image_out(),
  Argument::End
};


OPTIONS = {

  Option ("template", "template image", "an image file to be used as a template for the output (the output image wil have the same transform and field of view).")
    .append (Argument ("image", "template image", "the input template image").type_image_in()),

  Option ("vox", "voxel size", "provide either an isotropic voxel size, or comma-separated list of 3 voxel dimensions.")
    .append (Argument ("size", "voxel size", "the voxel size (in mm).").type_sequence_float()),

  Option ("colour", "coloured map", "add colour to the output image according to the direction of the tracks."),

  Option ("fraction", "output fibre fraction", "produce an image of the fraction of fibres through each voxel (as a proportion of the total number in the file), rather than the count."),

  Option ("datatype", "data type", "specify output image data type.")
    .append (Argument ("spec", "specifier", "the data type specifier.").type_choice (data_type_choices)),

  Option ("resample", "resample tracks", "resample the tracks at regular intervals using Hermite interpolation. If omitted, the program will select an appropriate interpolation factor automatically.")
    .append (Argument ("factor", "factor", "the factor by which to resample.").type_integer (1, INT_MAX, 1)),

  Option::End
};


// Value of tension for Hermite interpolation
#define HERMITE_TENSION 0.1

// When generating a custom image header from the data set, do not
// read more than this number of tracks to determine the real-space
// bounds of the data
#define MAX_TRACKS_READ_FOR_HEADER 1000000

// When automatically determining an appropriate interpolation factor,
// force the approximate step size of the interpolated track to be
// AT MOST this fraction of the minimum voxel dimension
#define INTERP_VOX_DIM_FRACTION 0.5



class Voxel
{

  public:
    Voxel () : x (0), y (0), z (0) { }

    Voxel (const size_t _x_, const size_t _y_, const size_t _z_) :
      x (_x_),
      y (_y_),
      z (_z_) { }

    Voxel (const Point& p) :
      x (round (p[0])),
      y (round (p[1])),
      z (round (p[2]))
    { 
      assert (finite (p[0]) && finite (p[1]) && finite (p[2]));
    }

    Voxel& operator+= (const Voxel& source) 
    {
      x += source.x;
      y += source.y;
      z += source.z;
      return (*this);
    }

    bool operator< (const Voxel& v) const
    {
      return (x == v.x ? (y == v.y ? (z < v.z ? 1 : 0) : (y < v.y)) : (x < v.x));
    }

    bool operator== (const Voxel& v) const
    {
      return (x == v.x && y == v.y && z == v.z);
    }

    bool test_bounds (const Image::Header& H) const 
    {
      return (x < size_t(H.dim(0)) && y < size_t(H.dim(1)) && z < size_t(H.dim(2)));
    }

    size_t x, y, z;

};


class VoxelDir : public Voxel
{

  public:
    VoxelDir () : Voxel () { }

    VoxelDir (const Point& p) :
      Voxel (p),
      dir (0.0, 0.0, 0.0) { }

    void set_dir (const Point& dp) 
    {
      dir[0] = fabs (dp[0]);
      dir[1] = fabs (dp[1]);
      dir[2] = fabs (dp[2]);
    }

    bool operator< (const VoxelDir& v) const
    {
      return (Voxel::operator< (v));
    }

    Point dir;

};







class Resampler
{

  public:
    Resampler (const Math::Matrix& interp_matrix, const size_t c) :
      M (interp_matrix),
      columns (c),
      data (4, c) { }

    ~Resampler() { }

    size_t get_columns () const { return (columns); }
    bool valid () const { return (M.is_valid()); }

    void init (const float* a, const float* b, const float* c)
    {
      for (size_t i = 0; i != columns; ++i) {
        data(0,i) = 0.0;
        data(1,i) = a[i];
        data(2,i) = b[i];
        data(3,i) = c[i];
      }
    }

    void increment (const float* a)
    {
      for (size_t i = 0; i != columns; ++i) {
        data(0,i) = data(1,i);
        data(1,i) = data(2,i);
        data(2,i) = data(3,i);
        data(3,i) = a[i];
      }
    }

    void interpolate (Math::Matrix& result) const
    {
      result.multiply (M, data);
    }

  private:
    const Math::Matrix& M;
    size_t columns;
    Math::Matrix data;

};





template <class T> class TrackMapper 
{

  public:
    TrackMapper (const Image::Position& pos, const Math::Matrix& interp_matrix) :
      H (pos.image.header()),
      interp (pos.image),
      resample_matrix (interp_matrix),
      R (resample_matrix, 3) { }

    void map (std::vector<Point>& tck, std::set<T>& output)
    {
      Math::Matrix data;
      if (R.valid()) {
        assert (resample_matrix.is_valid());
        assert (resample_matrix.rows());
        data.allocate (resample_matrix.rows(), 3);
        interp_track (tck, R, data);
      }
      voxelise (tck, output);
    }


  private:
    const Image::Header& H;
    Image::Interp interp;
    const Math::Matrix& resample_matrix;
    Resampler R;


    void tck_interp_prepare (std::vector<Point>& v)
    {
      const size_t s = v.size();
      v.insert    (v.begin(), v[ 0 ] + (2 * (v[ 0 ] - v[ 1 ])) - (v[ 1 ] - v[ 2 ]));
      v.push_back (           v[ s ] + (2 * (v[ s ] - v[s-1])) - (v[s-1] - v[s-2]));
    }

    void interp_track (
        std::vector<Point>& tck,
        Resampler& R,
        Math::Matrix& data)
    {
      std::vector<Point> out;
      tck_interp_prepare (tck);
      R.init (tck[0].get(), tck[1].get(), tck[2].get());
      for (size_t i = 3; i < tck.size(); ++i) {
        out.push_back (tck[i-2]);
        R.increment (tck[i].get());
        R.interpolate (data);
        for (size_t row = 0; row != data.rows(); ++row)
          out.push_back (Point (data(row,0), data(row,1), data(row,2)));
      }
      out.push_back (tck[tck.size() - 2]);
      out.swap (tck);
    }

    void voxelise (const std::vector<Point>& tck, std::set<T>& voxels) const;

};


template<> void TrackMapper<Voxel>::voxelise (const std::vector<Point>& tck, std::set<Voxel>& voxels) const
{
  for (std::vector<Point>::const_iterator i = tck.begin(); i != tck.end(); ++i) {
    Voxel vox (interp.R2P (*i));
    if (vox.test_bounds (H)) 
      voxels.insert (vox);
  }
}


template<> void TrackMapper<VoxelDir>::voxelise (const std::vector<Point>& tck, std::set<VoxelDir>& voxels) const
{
  std::vector<Point>::const_iterator prev = tck.begin();
  const std::vector<Point>::const_iterator last = tck.end() - 1;

  for (std::vector<Point>::const_iterator i = tck.begin(); i != last; ++i) {
    VoxelDir vox (interp.R2P (*i));
    if (vox.test_bounds (H)) {
      vox.set_dir ((*(i+1) - *prev).normalise());
      voxels.insert (vox);
    }
    prev = i;
  }
  VoxelDir vox (interp.R2P (*last));
  if (vox.test_bounds (H)) {
    vox.set_dir ((*last - *prev).normalise());
    voxels.insert (vox);
  }
}







template <class T> class MapWriterBase
{

  public:
    MapWriterBase (Image::Position& p, const float fraction_scaling_factor) :
      pos (p),
      H (p.image.header()),
      scale (fraction_scaling_factor)
    {
      buffer_size = H.dim(0) * H.dim(1) * H.dim(2);
    }

    virtual void write (const std::set<T>& voxels) = 0;

   protected:
    Image::Position& pos;
    const Image::Header& H;
    float scale;
    size_t buffer_size;

    size_t voxel_to_index(const T& vox)
    {
      return (vox.x + (vox.y * H.dim(0)) + (vox.z * H.dim(0) * H.dim(1)));
    }

};


class MapWriter : public MapWriterBase<Voxel> 
{

  public:
    MapWriter (Image::Position& p, const float fraction_scaling_factor) :
      MapWriterBase<Voxel> (p, fraction_scaling_factor),
      buffer (new size_t[buffer_size])
    {
      memset(buffer, 0, buffer_size * sizeof(size_t));
    }

    ~MapWriter ()
    {
      ProgressBar::init(H.dim(2), "writing image... ");
      size_t index = 0;
      for (pos.set(2,0); pos[2] < H.dim(2); pos.inc(2)) {
        for (pos.set(1,0); pos[1] < H.dim(1); pos.inc(1)) {
          for (pos.set(0,0); pos[0] < H.dim(0); pos.inc(0), ++index)
            pos.value (scale * buffer[index]);
        }
        ProgressBar::inc();
      }
      ProgressBar::done();
      delete[] buffer;
    }

    void write (const std::set<Voxel>& voxels)
    {
      for (std::set<Voxel>::const_iterator i = voxels.begin(); i != voxels.end(); ++i)
        buffer[voxel_to_index(*i)] += 1;
    }

  private:
    size_t* buffer;

};


class MapWriterColour : public MapWriterBase<VoxelDir>
{

  public:
    MapWriterColour(Image::Position& p, const float fraction_scaling_factor) :
      MapWriterBase<VoxelDir> (p, fraction_scaling_factor),
      buffer (new Point[buffer_size])
    {
      for (size_t i = 0; i != buffer_size; ++i)
        buffer[i] = Point (0.0, 0.0, 0.0);
    }

    ~MapWriterColour ()
    {
      ProgressBar::init (H.dim(2), "writing colour image... ");
      size_t index = 0;
      for (pos.set(2,0); pos[2] < H.dim(2); pos.inc(2)) {
        for (pos.set(1,0); pos[1] < H.dim(1); pos.inc(1)) {
          for (pos.set(0,0); pos[0] < H.dim(0); pos.inc(0), ++index) {
            pos.set(3, 0); pos.value ((buffer[index])[0]);
            pos.inc(3);    pos.value ((buffer[index])[1]);
            pos.inc(3);    pos.value ((buffer[index])[2]);
          }
        }
        ProgressBar::inc();
      }
      ProgressBar::done();
      delete[] buffer;
    }

    void write (const std::set<VoxelDir>& voxels)
    {
      std::set<VoxelDir>::const_iterator i = voxels.begin();
      while (i != voxels.end()) {
        const std::set<VoxelDir>::const_iterator this_voxel = i;
        Point this_voxel_dir(0.0, 0.0, 0.0);
        do {
          this_voxel_dir += i->dir;
        } while ((++i) != voxels.end() && *i == *this_voxel);
        buffer[voxel_to_index(*this_voxel)] += this_voxel_dir.normalise();
      }
    }


  private:
    Point* buffer;

};




class Hermite
{

  public:

    Hermite (float tension = 0.0) :
        t (0.5 * tension) { }

    void set (float position) {
      float p2 = position * position;
      float p3 = position * p2;
      w[0] = (0.5 - t) * ((2.0 * p2) - p3 - position);
      w[1] = 1.0 + ((1.5 + t) * p3) - ((2.5 + t) * p2);
      w[2] = ((2.0 + (2.0 * t)) * p2) + ((0.5 - t) * position) - ((1.5 + t) * p3);
      w[3] = (0.5- t) * (p3 - p2);
    }

    float coef (size_t i) const { return (w[i]); }

    float value (const float& a, const float& b, const float& c, const float& d) const
    {
      return (w[0]*a + w[1]*b + w[2]*c + w[3]*d);
    }

    float value (const float* vals) const
    {
      return (value (vals[0], vals[1], vals[2], vals[3]));
    }

  private:
    float w[4];
    float t;

};



void generate_header (Image::Header& header, DWI::Tractography::Reader& file, const std::vector<float>& voxel_size)
{

  std::vector<Point> tck;
  size_t track_counter = 0;

  Point min_values ( INFINITY,  INFINITY,  INFINITY);
  Point max_values (-INFINITY, -INFINITY, -INFINITY);

  ProgressBar::init (0, "creating new template image... ");

  while (file.next (tck) && track_counter++ < MAX_TRACKS_READ_FOR_HEADER) {
    for (std::vector<Point>::const_iterator i = tck.begin(); i != tck.end(); ++i) {
      min_values[0] = std::min (min_values[0], (*i)[0]);
      max_values[0] = std::max (max_values[0], (*i)[0]);
      min_values[1] = std::min (min_values[1], (*i)[1]);
      max_values[1] = std::max (max_values[1], (*i)[1]);
      min_values[2] = std::min (min_values[2], (*i)[2]);
      max_values[2] = std::max (max_values[2], (*i)[2]);
    }
    ProgressBar::inc();
  }

  min_values -= Point (3.0*voxel_size[0], 3.0*voxel_size[1], 3.0*voxel_size[2]);
  max_values += Point (3.0*voxel_size[0], 3.0*voxel_size[1], 3.0*voxel_size[2]);

  header.name = "tckmap image header";
  header.axes.set_ndim(3);

  for (size_t i = 0; i != 3; ++i) {
    header.axes.dim[i]    = ceil ((max_values[i] - min_values[i]) / voxel_size[i]);
    header.axes.vox[i]    = voxel_size[i];
    header.axes.axis[i]   = i;
    header.axes.units[i]  = Image::Axis::millimeters;
  }

  header.axes.desc[0] = Image::Axis::left_to_right;
  header.axes.desc[1] = Image::Axis::posterior_to_anterior;
  header.axes.desc[2] = Image::Axis::inferior_to_superior;

  Math::Matrix M (4, 4);
  M.identity();
  M(0,3) = min_values[0];
  M(1,3) = min_values[1];
  M(2,3) = min_values[2];
  header.set_transform (M);

  ProgressBar::done();

}



Math::Matrix gen_interp_matrix (const size_t os_factor)
{
  Math::Matrix M;
  if (os_factor > 1) {
    const size_t dim = os_factor - 1;
    Hermite interp (HERMITE_TENSION);
    M.allocate (dim, 4);
    for (size_t i = 0; i != dim; ++i) {
      interp.set ((i + 1.0) / float(os_factor));
      for (size_t j = 0; j != 4; ++j)
        M(i,j) = interp.coef (j);
    }
  }
  return (M);
}



void oversample_header (Image::Header& header, const std::vector<float>& voxel_size) 
{
  info ("oversampling header...");

  Math::Matrix M (header.transform());

  for (size_t i = 0; i != 3; ++i) {
    M(i, 3) += 0.5 * (voxel_size[i] - header.axes.vox[i]);
    header.axes.dim[i] = ceil (header.axes.dim[i] * header.axes.vox[i] / voxel_size[i]);
    header.axes.vox[i] = voxel_size[i];
  }

  header.set_transform(M);
}






EXECUTE {

  DWI::Tractography::Properties properties;
  DWI::Tractography::Reader file;
  file.open (argument[0].get_string(), properties);

  const size_t num_tracks = properties["count"]    .empty() ? 0   : to<size_t> (properties["count"]);
  const float  step_size  = properties["step_size"].empty() ? 0.0 : to<float>  (properties["step_size"]);

  const bool colour         = get_options (2).size();
  const bool fibre_fraction = get_options (3).size();

  std::vector<float> voxel_size;
  std::vector<OptBase> opt = get_options(1);
  if (opt.size())
    voxel_size = parse_floats (opt[0][0].get_string());

  if (voxel_size.size() == 1)
    voxel_size.assign (3, voxel_size.front());
  else if (!voxel_size.empty() && voxel_size.size() != 3)
    throw Exception ("voxel size must either be a single isotropic value, or a list of 3 comma-separated voxel dimensions");

  if (!voxel_size.empty())
    info("creating image with voxel dimensions [ " + str(voxel_size[0]) + " " + str(voxel_size[1]) + " " + str(voxel_size[2]) + " ]");

  Image::Header header;
  opt = get_options (0);
  if (opt.size()) {
    header = (*opt[0][0].get_image()).header();
    if (!voxel_size.empty())
      oversample_header (header, voxel_size);
  }
  else {
    if (voxel_size.empty())
      throw Exception ("please specify either a template image or the desired voxel size");
    generate_header (header, file, voxel_size);
    file.close();
    file.open (argument[0].get_string(), properties);
  }

  opt = get_options (4);
  if (opt.size())
    header.data_type.parse (data_type_choices[opt[0][0].get_int()]);
  else
    header.data_type = (fibre_fraction || colour) ? DataType::Float32 : DataType::UInt32;

  for (DWI::Tractography::Properties::const_iterator i = properties.begin(); i != properties.end(); ++i)
    header.comments.push_back (i->first + ": " + i->second);
  for (std::vector< RefPtr<DWI::Tractography::ROI> >::const_iterator i = properties.roi.begin(); i != properties.roi.end(); ++i)
    header.comments.push_back ("ROI: " + (*i)->specification());
  for (std::vector<std::string>::const_iterator i = properties.comments.begin(); i != properties.comments.end(); ++i)
    header.comments.push_back ("comment: " + *i);

  float scaling_factor = fibre_fraction ? 1.0 / float(num_tracks) : 1.0;
  header.comments.push_back("scaling_factor: " + str(scaling_factor));
  info ("intensity scaling factor set to " + str(scaling_factor));

  size_t resample_factor;
  opt = get_options (5);
  if (opt.size()) {
    resample_factor = opt[0][0].get_int();
    info ("track interpolation factor manually set to " + str(resample_factor));
  } 
  else if (step_size) {
    resample_factor = ceil (step_size / (minvalue (header.vox(0), header.vox(1), header.vox(2)) * INTERP_VOX_DIM_FRACTION));
    info ("track interpolation factor automatically set to " + str(resample_factor));
  } 
  else {
    resample_factor = 1;
    info ("track interpolation off; no track step size information in header");
  }

  Math::Matrix interp_matrix (gen_interp_matrix (resample_factor));
  std::vector<Point> tck;

  if (colour) {

    header.axes.set_ndim(4);
    header.axes.dim[3] = 3;
    header.axes.vox[3] = 0;
    header.axes.axis[0] = 1;
    header.axes.axis[1] = 2;
    header.axes.axis[2] = 3;
    header.axes.axis[3] = 0;
    header.axes.desc[3] = "direction";
    header.comments.push_back (std::string ("coloured track density map"));

    Image::Position       pos    (*argument[1].get_image (header));
    TrackMapper<VoxelDir> mapper (pos, interp_matrix);
    MapWriterColour       writer (pos, scaling_factor);

    ProgressBar::init (num_tracks, "mapping tracks to colour image... ");
    while (file.next (tck)) {
      std::set<VoxelDir> mapped_voxels;
      mapper.map (tck, mapped_voxels);
      writer.write (mapped_voxels);
      ProgressBar::inc();
    }
    ProgressBar::done();

  } 
  else {

    header.axes.set_ndim(3);
    header.comments.push_back (std::string (("track ") + str(fibre_fraction ? "fraction" : "count") + " map"));

    Image::Position    pos    (*argument[1].get_image(header));
    TrackMapper<Voxel> mapper (pos, interp_matrix);
    MapWriter          writer (pos, scaling_factor);

    ProgressBar::init (num_tracks, "mapping tracks to image... ");
    while (file.next (tck)) {
      std::set<Voxel> mapped_voxels;
      mapper.map (tck, mapped_voxels);
      writer.write (mapped_voxels);
      ProgressBar::inc();
    }
    ProgressBar::done();

  }

}

