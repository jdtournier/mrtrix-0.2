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

    15-10-2008 J-Donald Tournier <d.tournier@brain.org.au>
    * remove MR::DICOM_DW_gradients_PRS flag

    31-10-2008 J-Donald Tournier <d.tournier@brain.org.au>
    * scale b-value by gradient magnitude and normalise gradient direction

    19-12-2008 J-Donald Tournier <d.tournier@brain.org.au>
    * handle cases where the data size is greater than expected, and interpret as multi-channel data

    03-03-2010 J-Donald Tournier <d.tournier@brain.org.au>
    * improved GE gradient information support

*/

#include "image/header.h"
#include "image/mapper.h"
#include "file/dicom/mapper.h"
#include "file/dicom/image.h"
#include "file/dicom/series.h"
#include "file/dicom/study.h"
#include "file/dicom/patient.h"
#include "file/dicom/tree.h"

namespace MR {
  namespace File {
    namespace Dicom {

      namespace {
        const gchar* FormatDICOM = "DICOM";

        inline void add_axis (const char* label, MR::Image::Axes& axes, guint& current_axis, guint& current_stride, guint size) 
        {
          if (size > 1) {
            axes.axis[current_axis] = current_stride++;
            axes.dim[current_axis] = size;
            axes.desc[current_axis] = label;
            ++current_axis;
          }
        }
      }


      void dicom_to_mapper (
          MR::Image::Mapper& dmap, 
          MR::Image::Header& H, 
          std::vector< RefPtr<Series> >& series_list)
      {
        assert (series_list.size() > 0);
        H.format = FormatDICOM;

        Patient* patient (series_list[0]->study->patient);
        String sbuf = ( patient->name.size() ? patient->name : "unnamed" );
        sbuf += " " + format_ID (patient->ID);
        if (series_list[0]->modality.size()) sbuf += String (" [") + series_list[0]->modality + "]";
        if (series_list[0]->name.size()) sbuf += String (" ") + series_list[0]->name;
        H.comments.push_back (sbuf);
        H.name = sbuf;

        // build up sorted list of frames:
        std::vector<Frame*> frames;

        // loop over series list:
        for (std::vector< RefPtr<Series> >::const_iterator series_it = series_list.begin(); series_it != series_list.end(); ++series_it) {
          Series& series (**series_it);

          try {
            series.read();
          }
          catch (Exception) { 
            throw Exception ("error reading series " + str (series.number) + " of DICOM image \"" + H.name + "\""); 
          }

          std::sort (series.begin(), series.end());

          // loop over images in each series:
          for (Series::const_iterator image_it = series.begin(); image_it != series.end(); ++image_it) {
            Image& image (**image_it);

            // if multi-frame, loop over frames in image:
            if (image.frames.size()) {
              std::sort (image.frames.begin(), image.frames.end());
              for (std::vector< RefPtr<Frame> >::const_iterator frame_it = image.frames.begin(); frame_it != image.frames.end(); ++frame_it) 
                frames.push_back (&**frame_it);
            }
            // otherwise add image frame:
            else 
              frames.push_back (&image);

          }
        }

        std::vector<guint> dim = Frame::count (frames);

        if (dim[0]*dim[1]*dim[2] > frames.size())
          throw Exception ("missing image frames for DICOM image \"" + H.name + "\"");

        if (dim[0] > 1) { // switch axes so slice dim is inner-most:
          std::vector<Frame*> list (frames);
          std::vector<Frame*>::iterator it = frames.begin();
          for (guint k = 0; k < dim[2]; ++k) 
            for (guint i = 0; i < dim[0]; ++i) 
              for (guint j = 0; j < dim[1]; ++j) 
                *(it++) = list[i+dim[0]*(j+dim[1]*k)];
        }

        gfloat slice_separation = Frame::get_slice_separation (frames, dim[1]);


        if (series_list[0]->study->name.size()) {
          sbuf = "study: " + series_list[0]->study->name;
          H.comments.push_back (sbuf);
        }

        if (patient->DOB.size()) {
          sbuf = "DOB: " + format_date (patient->DOB);
          H.comments.push_back (sbuf);
        }

        if (series_list[0]->date.size()) {
          sbuf = "DOS: " + format_date (series_list[0]->date);
          if (series_list[0]->time.size()) 
            sbuf += " " + format_time (series_list[0]->time);
          H.comments.push_back (sbuf);
        }





        MR::Image::Axes& axes (H.axes);
        const Image& image (*(*series_list[0])[0]);

        guint nchannels = image.frames.size() ? 1 : image.data_size / (image.dim[0] * image.dim[1] * (image.bits_alloc/8));
        if (nchannels > 1) 
          info ("data segment is larger than expected from image dimensions - interpreting as multi-channel data");

        axes.set_ndim (3 + (dim[0]*dim[2]>1) + (nchannels>1));

        guint current_axis = 3;
        guint current_stride = 0;

        add_axis ("channel", axes, current_axis, current_stride, nchannels);

        axes.axis[0] = current_stride++;
        axes.dim[0] = image.dim[0];
        axes.vox[0] = image.pixel_size[0];
        axes.desc[0] = MR::Image::Axis::left_to_right;
        axes.units[0] = MR::Image::Axis::millimeters;

        axes.axis[1] = current_stride++;
        axes.dim[1] = image.dim[1];
        axes.vox[1] = image.pixel_size[1];
        axes.desc[1] = MR::Image::Axis::posterior_to_anterior;
        axes.units[1] = MR::Image::Axis::millimeters;

        axes.axis[2] = current_stride++;
        axes.dim[2] = dim[1];
        axes.vox[2] = slice_separation;
        axes.desc[2] = MR::Image::Axis::inferior_to_superior;
        axes.units[2] = MR::Image::Axis::millimeters;

        add_axis ("volume", axes, current_axis, current_stride, dim[0]*dim[2]);

        if (image.bits_alloc == 8) H.data_type = DataType::UInt8;
        else if (image.bits_alloc == 16) {
          H.data_type = DataType::UInt16;
          if (image.is_BE) H.data_type.set_flag (DataType::BigEndian);
          else H.data_type.set_flag (DataType::LittleEndian);
        }
        else throw Exception ("unexpected number of allocated bits per pixel (" + str (image.bits_alloc) 
            + ") in file \"" + H.name + "\"");

        H.offset = image.scale_intercept;
        H.scale = image.scale_slope;

        Math::Matrix M(4,4);

        M(0,0) = -image.orientation_x[0];
        M(1,0) = -image.orientation_x[1];
        M(2,0) = image.orientation_x[2];
        M(3,0) = 0.0;

        M(0,1) = -image.orientation_y[0];
        M(1,1) = -image.orientation_y[1];
        M(2,1) = image.orientation_y[2];
        M(3,1) = 0.0;

        M(0,2) = -image.orientation_z[0];
        M(1,2) = -image.orientation_z[1];
        M(2,2) = image.orientation_z[2];
        M(3,2) = 0.0;

        M(0,3) = -image.position_vector[0];
        M(1,3) = -image.position_vector[1];
        M(2,3) = image.position_vector[2];
        M(3,3) = 1.0;

        H.set_transform (M);


        H.DW_scheme = Frame::get_DW_scheme (frames, dim[1], M);



        for (guint n = 1; n < frames.size(); ++n) // check consistency of data scaling:
          if (frames[n]->scale_intercept != frames[n-1]->scale_intercept ||
              frames[n]->scale_slope != frames[n-1]->scale_slope)
            throw Exception ("unable to load series due to inconsistent data scaling between DICOM images");


        if (image.frames.size()) { // need to preload and re-arrange:
          ProgressBar::init (frames.size(), "DICOM image contains multiple frames - reformating..."); 

          guint8* mem = NULL;
          try { 
            mem = new guint8 [H.memory_footprint()]; 
            if (!mem) throw;
          }
          catch (...) { 
            throw Exception ("failed to allocate memory for image data!"); 
          }

          guint8* dest = mem;
          const guint row_stride = nchannels * frames[0]->row_stride * (frames[0]->bits_alloc/8);
          const guint row_size = nchannels * frames[0]->dim[0] * (frames[0]->bits_alloc/8);
          File::MMap mmap;
          for (guint n = 0; n < frames.size(); ++n) {
            mmap.init (frames[n]->filename);
            mmap.map();
            const guint8* src = (guint8*) mmap.address() + frames[n]->data;
            for (guint row = 0; row < frames[n]->dim[1]; ++row) {
              memcpy (dest, src, row_size);
              dest += row_size;
              src += row_stride;
            }
            ProgressBar::inc();
          }
          ProgressBar::done();

          dmap.add (mem);

        }
        else { // use standard backend:
          
          for (guint n = 0; n < frames.size(); ++n) {
            const Image* image = static_cast<const Image*> (frames[n]);
            dmap.add (image->filename, image->data);
          }

        }


      
      }






    }
  }
}

