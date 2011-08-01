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
#include "data_type.h"
#include "math/vector.h"
#include "file/dicom/element.h"

namespace MR {
  namespace File {
    namespace Dicom {

      class Series;
      class Element;

      class Image {

        public:
          Image (Series* parent = NULL);

          String            filename;
          String            sequence_name;
          String            manufacturer;
          Series*           series;

          guint             acq_dim[2], dim[2], instance, acq, sequence;
          gfloat            position_vector[3], orientation_x[3], orientation_y[3], orientation_z[3], distance;
          gfloat            pixel_size[2], slice_thickness, scale_slope, scale_intercept;
          gfloat            bvalue, G[3];
          guint             data, bits_alloc, images_in_mosaic, data_size;
          DataType          data_type;
          bool              is_BE;

          void              read ();
          void              parse_item (Element& item, const String& dirname = "");
          void              decode_csa (const guint8* start, const guint8* end);

          void              calc_distance ();

          bool              operator< (const Image& ima) const;

          void              print_fields (bool dcm, bool csa) const;
      };

      std::ostream& operator<< (std::ostream& stream, const Image& item);














      inline Image::Image (Series* parent) :
        series (parent)
      { 
        acq_dim[0] = acq_dim[1] = dim[0] = dim[1] = instance = acq = sequence = UINT_MAX;
        position_vector[0] = position_vector[1] = position_vector[2] = GSL_NAN;
        orientation_x[0] = orientation_x[1] = orientation_x[2] = GSL_NAN;
        orientation_y[0] = orientation_y[1] = orientation_y[2] = GSL_NAN;
        orientation_z[0] = orientation_z[1] = orientation_z[2] = GSL_NAN;
        distance = GSL_NAN;
        pixel_size[0] = pixel_size[1] = slice_thickness = GSL_NAN; 
        scale_intercept = 0.0;
        scale_slope = 1.0;
        bvalue = G[0] = G[1] = G[2] = GSL_NAN;
        data = bits_alloc = images_in_mosaic = data_size = 0;
        is_BE = false;
      }







      inline void Image::calc_distance ()
      {
        if (images_in_mosaic) {
          gfloat xinc = pixel_size[0] * (dim[0] - acq_dim[0]) / 2.0;
          gfloat yinc = pixel_size[1] * (dim[1] - acq_dim[1]) / 2.0;
          for (guint i = 0; i < 3; i++) 
            position_vector[i] += xinc * orientation_x[i] + yinc * orientation_y[i];

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
        else Math::cross_product (orientation_z, orientation_x, orientation_y);

        Math::normalise (orientation_z);
        distance = Math::dot_product (orientation_z, position_vector);
      }


    }
  }
}


#endif


