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


    08-09-2008 J-Donald Tournier <d.tournier@brain.org.au>
    * fix handling of mosaic slice ordering (using SliceNormalVector entry in CSA header)

    03-03-2010 J-Donald Tournier <d.tournier@brain.org.au>
    * improved GE gradient information support

*/

#ifndef __file_dicom_image_h__
#define __file_dicom_image_h__

#include "ptr.h"
#include "math/vector.h"
#include "file/dicom/element.h"

namespace MR {
  namespace File {
    namespace Dicom {

      class Series;
      class Element;

      class Frame { 
        public:
          Frame () { 
            acq_dim[0] = acq_dim[1] = dim[0] = dim[1] = row_stride = instance = series_num = acq = sequence = UINT_MAX;
            position_vector[0] = position_vector[1] = position_vector[2] = GSL_NAN;
            orientation_x[0] = orientation_x[1] = orientation_x[2] = GSL_NAN;
            orientation_y[0] = orientation_y[1] = orientation_y[2] = GSL_NAN;
            orientation_z[0] = orientation_z[1] = orientation_z[2] = GSL_NAN;
            distance = GSL_NAN;
            pixel_size[0] = pixel_size[1] = slice_thickness = GSL_NAN; 
            scale_intercept = 0.0;
            scale_slope = 1.0;
            bvalue = G[0] = G[1] = G[2] = GSL_NAN;
            data = bits_alloc = data_size = frame_offset = 0;
            DW_scheme_wrt_image = false;
          }

          guint  acq_dim[2], dim[2], row_stride, series_num, instance, acq, sequence;
          gfloat position_vector[3], orientation_x[3], orientation_y[3], orientation_z[3], distance;
          gfloat pixel_size[2], slice_thickness, scale_slope, scale_intercept;
          gfloat bvalue, G[3];
          guint  data, bits_alloc, data_size, frame_offset;
          String filename;
          bool DW_scheme_wrt_image;
          std::vector<guint32> index;


          bool operator< (const Frame& frame) const {
            if (series_num != frame.series_num) return series_num < frame.series_num;
            if (acq != frame.acq) return acq < frame.acq;
            assert (!gsl_isnan(distance));
            assert (!gsl_isnan(frame.distance));
            if (distance != frame.distance) return distance < frame.distance;
            for (guint n = index.size(); n--;)
              if (index[n] != frame.index[n])
                return index[n] < frame.index[n];
            if (sequence != frame.sequence) return sequence < frame.sequence;
            if (instance != frame.instance) return instance < frame.instance;
            return false;
          }


          void calc_distance ()
          {
            if (gsl_isnan (orientation_z[0])) 
              Math::cross_product (orientation_z, orientation_x, orientation_y);
            else {
              float normal[3];
              Math::cross_product (normal, orientation_x, orientation_y);
              if (Math::dot_product (normal, orientation_z) < 0.0) {
                orientation_z[0] = -normal[0];
                orientation_z[1] = -normal[1];
                orientation_z[2] = -normal[2];
              }
              else {
                orientation_z[0] = normal[0];
                orientation_z[1] = normal[1];
                orientation_z[2] = normal[2];
              }
            }
            row_stride = dim[0];

            Math::normalise (orientation_z);
            distance = Math::dot_product (orientation_z, position_vector);
          }

          static std::vector<guint> count (const std::vector<Frame*>& frames);
          static gfloat get_slice_separation (const std::vector<Frame*>& frames, guint nslices);
          static Math::Matrix get_DW_scheme (const std::vector<Frame*>& frames, guint nslices, const Math::Matrix& image_transform);

          friend std::ostream& operator<< (std::ostream& stream, const Frame& item);
      };











      class Image : public Frame {

        public:
          Image (Series* parent = NULL) : 
            series (parent), 
            images_in_mosaic (0),
            is_BE (false),
            in_frames (false) { }

          Series* series;
          guint images_in_mosaic;
          String  sequence_name, manufacturer;
          bool   is_BE, in_frames;

          std::vector<guint32> frame_dim;
          std::vector< RefPtr<Frame> > frames;

          void read (bool print_DICOM_fields = false, bool print_CSA_fields = false);
          void parse_item (Element& item, bool print_DICOM_fields = false, bool print_CSA_fields = false);
          void decode_csa (const guint8* start, const guint8* end, bool print_fields = false);

          bool operator< (const Image& ima) const {
            return Frame::operator< (ima);
          }

          friend std::ostream& operator<< (std::ostream& stream, const Image& item);
      };






    }
  }
}


#endif


