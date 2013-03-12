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

#include <zlib.h>
#include <glibmm/stringutils.h>

#include "file/nifti1.h"
#include "image/mapper.h"
#include "get_set.h"
#include "image/format/list.h"
#include "math/quaternion.h"


namespace MR {
  namespace Image {
    namespace Format {


      const gchar* FormatNIfTI = "NIfTI-1.1";

      bool NIfTI::read (Mapper& dmap, Header& H) const
      {
        if (!(Glib::str_has_suffix (H.name, ".nii") || Glib::str_has_suffix (H.name, ".nii.gz"))) 
          return (false);

        File::MMap fmap;
        String gzfilename;
        if (Glib::str_has_suffix (H.name, ".gz")) {
          gzfilename = H.name;
          fmap.init (Mapper::gunzip (H.name, "nii"));
        }
        else 
          fmap.init (H.name);

        fmap.map();
        H.format = FormatNIfTI;

        const nifti_1_header* NH = (const nifti_1_header*) fmap.address();

        bool is_BE = false;
        if (get<gint32> (&NH->sizeof_hdr, is_BE) != 348) {
          is_BE = true;
          if (get<gint32> (&NH->sizeof_hdr, is_BE) != 348) 
            throw Exception ("image \"" + H.name + "\" is not in NIfTI format (sizeof_hdr != 348)");
        }
        if (memcmp (NH->magic, "n+1\0", 4)) 
          throw Exception ("image \"" + H.name + "\" is not in NIfTI format (magic != \"n+1\\0\")");

        char db_name[19];
        strncpy (db_name, NH->db_name, 18);
        if (db_name[0]) {
          db_name[18] = '\0';
          H.comments.push_back (db_name);
        }

        // data set dimensions:
        int ndim = get<gint16> (&NH->dim, is_BE);
        if (ndim < 1) 
          throw Exception ("too few dimensions specified in NIfTI image \"" + H.name + "\"");
        if (ndim > 7) 
          throw Exception ("too many dimensions specified in NIfTI image \"" + H.name + "\"");
        H.axes.set_ndim (ndim);


        for (int i = 0; i < ndim; i++) {
          H.axes.dim[i] = get<gint16> (&NH->dim[i+1], is_BE);
          if (H.axes.dim[i] < 0) {
            info ("dimension along axis " + str(i) + " specified as negative in NIfTI image \"" + H.name + "\" - taking absolute value");
            H.axes.dim[i] = abs (H.axes.dim[i]);
          }
          if (!H.axes.dim[i]) H.axes.dim[i] = 1;
          H.axes.axis[i] = i;
          H.axes.forward[i] = true;
        }

        // data type:
        switch (get<gint16> (&NH->datatype, is_BE)) {
          case DT_BINARY:     H.data_type = DataType::Bit;        break;
          case DT_INT8:       H.data_type = DataType::Int8;       break;
          case DT_UINT8:      H.data_type = DataType::UInt8;      break;
          case DT_INT16:      H.data_type = DataType::Int16;      break;
          case DT_UINT16:     H.data_type = DataType::UInt16;     break;
          case DT_INT32:      H.data_type = DataType::Int32;      break;
          case DT_UINT32:     H.data_type = DataType::UInt32;     break;
          case DT_FLOAT32:    H.data_type = DataType::Float32;    break;
          case DT_FLOAT64:    H.data_type = DataType::Float64;    break;
          case DT_COMPLEX64:  H.data_type = DataType::CFloat32;   break;
          case DT_COMPLEX128: H.data_type = DataType::CFloat64;   break;
          default: throw Exception ("unknown data type for Analyse image \"" + H.name + "\"");
        }

        if ( !( H.data_type.is (DataType::Bit) || H.data_type.is (DataType::UInt8) || H.data_type.is (DataType::Int8) ) ) {
          if (is_BE) H.data_type.set_flag (DataType::BigEndian);
          else H.data_type.set_flag (DataType::LittleEndian);
        }

        if (get<gint16> (&NH->bitpix, is_BE) != (gint16) H.data_type.bits())
          error ("WARNING: bitpix field does not match data type in NIfTI image \"" + H.name + "\" - ignored");

        // voxel sizes:
        for (int i = 0; i < ndim; i++) {
          H.axes.vox[i] = get<float32> (&NH->pixdim[i+1], is_BE);
          if (H.axes.vox[i] < 0.0) {
            info ("voxel size along axis " + str(i) + " specified as negative in NIfTI image \"" + H.name + "\" - taking absolute value");
            H.axes.vox[i] = fabs (H.axes.vox[i]);
          }
        }


        // offset and scale:
        H.scale = get<float32> (&NH->scl_slope, is_BE);
        if (finite(H.scale) && H.scale != 0.0) {
          H.offset = get<float32> (&NH->scl_inter, is_BE);
          H.offset = finite (H.offset) ? H.offset : 0.0;
        }
        else {
          H.scale = 1.0;
          H.offset = 0.0;
        }

        gsize data_offset = (gsize) get<float32> (&NH->vox_offset, is_BE);

        char descrip[81];
        strncpy (descrip, NH->descrip, 80);
        if (descrip[0]) {
          descrip[80] = '\0';
          H.comments.push_back (descrip);
        }

        if (get<gint16> (&NH->sform_code, is_BE)) {
          Math::Matrix M (4,4);
          M(0,0) = get<float32> (&NH->srow_x[0], is_BE);
          M(0,1) = get<float32> (&NH->srow_x[1], is_BE);
          M(0,2) = get<float32> (&NH->srow_x[2], is_BE);
          M(0,3) = get<float32> (&NH->srow_x[3], is_BE);

          M(1,0) = get<float32> (&NH->srow_y[0], is_BE);
          M(1,1) = get<float32> (&NH->srow_y[1], is_BE);
          M(1,2) = get<float32> (&NH->srow_y[2], is_BE);
          M(1,3) = get<float32> (&NH->srow_y[3], is_BE);

          M(2,0) = get<float32> (&NH->srow_z[0], is_BE);
          M(2,1) = get<float32> (&NH->srow_z[1], is_BE);
          M(2,2) = get<float32> (&NH->srow_z[2], is_BE);
          M(2,3) = get<float32> (&NH->srow_z[3], is_BE);

          M(3,0) = M(3,1) = M(3,2) = 0.0;
          M(3,3) = 1.0;

          // get voxel sizes:
          H.axes.vox[0] = sqrt (gsl_pow_2 (M(0,0)) + gsl_pow_2 (M(1,0)) + gsl_pow_2 (M(2,0)));
          H.axes.vox[1] = sqrt (gsl_pow_2 (M(0,1)) + gsl_pow_2 (M(1,1)) + gsl_pow_2 (M(2,1)));
          H.axes.vox[2] = sqrt (gsl_pow_2 (M(0,2)) + gsl_pow_2 (M(1,2)) + gsl_pow_2 (M(2,2)));

          // normalize each transform axis:
          M(0,0) /= H.axes.vox[0];
          M(1,0) /= H.axes.vox[0];
          M(2,0) /= H.axes.vox[0];

          M(0,1) /= H.axes.vox[1];
          M(1,1) /= H.axes.vox[1];
          M(2,1) /= H.axes.vox[1];

          M(0,2) /= H.axes.vox[2];
          M(1,2) /= H.axes.vox[2];
          M(2,2) /= H.axes.vox[2];

          H.set_transform (M);
        }
        else if (get<gint16> (&NH->qform_code, is_BE)) {
          Math::Quaternion Q (
              get<float32> (&NH->quatern_b, is_BE), 
              get<float32> (&NH->quatern_c, is_BE),
              get<float32> (&NH->quatern_d, is_BE));
          float transform[9];
          Q.to_matrix (transform);
          Math::Matrix M (4,4);

          M(0,0) = transform[0];
          M(0,1) = transform[1];
          M(0,2) = transform[2];

          M(1,0) = transform[3];
          M(1,1) = transform[4];
          M(1,2) = transform[5];

          M(2,0) = transform[6];
          M(2,1) = transform[7];
          M(2,2) = transform[8];

          M(0,3) = get<float32> (&NH->qoffset_x, is_BE);
          M(1,3) = get<float32> (&NH->qoffset_y, is_BE);
          M(2,3) = get<float32> (&NH->qoffset_z, is_BE);

          M(3,0) = M(3,1) = M(3,2) = 0.0;
          M(3,3) = 1.0;

          // qfac:
          float qfac = get<float32> (&NH->pixdim[0], is_BE);
          if (qfac != 0.0) {
            M(0,2) = -M(0,2);
            M(1,2) = -M(1,2);
            M(2,2) = -M(2,2);
          }

          H.set_transform (M);
        }


        if (!H.axes.desc[0].size()) H.axes.desc[0] = Axis::left_to_right;
        if (!H.axes.units[0].size()) H.axes.units[0] = Axis::millimeters;
        if (H.axes.ndim() > 1) {
          if (!H.axes.desc[1].size()) H.axes.desc[1] = Axis::posterior_to_anterior;
          if (!H.axes.units[1].size()) H.axes.units[1] = Axis::millimeters;
          if (H.axes.ndim() > 2) {
            if (!H.axes.desc[2].size()) H.axes.desc[2] = Axis::inferior_to_superior;
            if (!H.axes.units[2].size()) H.axes.units[2] = Axis::millimeters;
          }
        }

        fmap.unmap();
        dmap.add_gz (fmap, gzfilename, data_offset);

        return (true);
      }





      bool NIfTI::check (Header& H, int num_axes) const
      {
        if (!(Glib::str_has_suffix (H.name, ".nii") || Glib::str_has_suffix (H.name, ".nii.gz"))) return (false);
        if (num_axes < 3) throw Exception ("cannot create NIfTI-1.1 image with less than 3 dimensions");
        if (num_axes > 8) throw Exception ("cannot create NIfTI-1.1 image with more than 8 dimensions");

        H.format = FormatNIfTI;

        H.axes.set_ndim (num_axes);
        for (int i = 0; i < H.axes.ndim(); i++) {
          if (H.axes.dim[i] < 1) H.axes.dim[i] = 1;
          H.axes.axis[i] = i;
          H.axes.forward[i] = true;
        }

        H.axes.desc[0] = Axis::left_to_right;
        H.axes.units[0] = Axis::millimeters;

        H.axes.desc[1] = Axis::posterior_to_anterior;
        H.axes.units[1] = Axis::millimeters;

        H.axes.desc[2] = Axis::inferior_to_superior;
        H.axes.units[2] = Axis::millimeters;

        return (true);
      }





      void NIfTI::create (Mapper& dmap, const Header& H) const
      {
        if (H.axes.ndim() > 7) 
          throw Exception ("NIfTI-1.1 format cannot support more than 7 dimensions for image \"" + H.name + "\"");

        guint msize = H.memory_footprint (H.ndim());

        File::MMap fmap;
        String gzfilename;
        if (Glib::str_has_suffix (H.name, ".gz")) {
          gzfilename = H.name;
          fmap.init (String(), 352 + msize, "nii");
        }
        else 
          fmap.init (H.name, 352 + msize);

        fmap.map();

        nifti_1_header* NH = (nifti_1_header*) fmap.address();

        bool is_BE = H.data_type.is_big_endian();

        // magic number:
        put<gint32> (348, &NH->sizeof_hdr, is_BE);

        strncpy ((gchar*) &NH->data_type, "dsr      \0", 10);
        strncpy ((gchar*) &NH->db_name, H.comments.size() ? H.comments[0].c_str() : "untitled\0\0\0\0\0\0\0\0\0\0\0", 18);
        put<gint32> (16384, &NH->extents, is_BE);
        NH->regular = 'r';
        NH->dim_info = 0;

        // data set dimensions:
        put<gint16> (H.ndim(), &NH->dim[0], is_BE);
        for (gint i = 0; i < H.ndim(); i++) 
          put<gint16> (H.dim(i), &NH->dim[i+1], is_BE);

        // pad out the other dimensions with 1, fix for fslview
        for (gint i = H.ndim() + 1; i < 8; i++) 
          put<gint16> (1, &NH->dim[i], is_BE);


        // data type:
        gint16 dt = 0;
        switch (H.data_type()) {
          case DataType::Bit:       dt = DT_BINARY; break;
          case DataType::Int8:      dt = DT_INT8;  break;
          case DataType::UInt8:     dt = DT_UINT8;  break;
          case DataType::Int16LE:   
          case DataType::Int16BE:   dt = DT_INT16; break;
          case DataType::UInt16LE: 
          case DataType::UInt16BE:  dt = DT_UINT16; break;
          case DataType::Int32LE: 
          case DataType::Int32BE:   dt = DT_INT32; break;
          case DataType::UInt32LE:
          case DataType::UInt32BE:  dt = DT_UINT32; break;
          case DataType::Float32LE:
          case DataType::Float32BE: dt = DT_FLOAT32; break;
          case DataType::Float64LE:
          case DataType::Float64BE: dt = DT_FLOAT64; break;
          case DataType::CFloat32LE:
          case DataType::CFloat32BE: dt = DT_COMPLEX64; break;
          case DataType::CFloat64LE:
          case DataType::CFloat64BE: dt = DT_COMPLEX128; break;
          default: throw Exception ("unknown data type for NIfTI-1.1 image \"" + H.name + "\"");
        }

        put<gint16> (dt, &NH->datatype, is_BE);  
        put<gint16> (H.data_type.bits(), &NH->bitpix, is_BE);

        // qfac:
        put<float32> (0.0, &NH->pixdim[0], is_BE);

        // voxel sizes:
        for (gint i = 0; i < H.ndim(); i++) 
          put<float32> (H.vox(i), &NH->pixdim[i+1], is_BE);

        put<float32> (352.0, &NH->vox_offset, is_BE);

        // offset and scale:
        put<float32> (H.scale, &NH->scl_slope, is_BE);
        put<float32> (H.offset, &NH->scl_inter, is_BE);

        NH->xyzt_units = SPACE_TIME_TO_XYZT (NIFTI_UNITS_MM, NIFTI_UNITS_SEC);

        gint pos = 0;
        char descrip[81];
        descrip[0] = '\0';
        for (guint i = 1; i < H.comments.size(); i++) {
          if (pos >= 75) break;
          if (i > 1) {
            descrip[pos++] = ';';
            descrip[pos++] = ' ';
          }
          strncpy (descrip + pos, H.comments[i].c_str(), 80-pos);
          pos += H.comments[i].size();
        }
        strncpy ((gchar*) &NH->descrip, descrip, 80);

        put<gint16> (NIFTI_XFORM_SCANNER_ANAT, &NH->qform_code, is_BE);
        put<gint16> (NIFTI_XFORM_SCANNER_ANAT, &NH->sform_code, is_BE);

        const Math::Matrix& M (H.transform());

        // qform:
        const float R [] = { 
          float(M(0,0)), float(M(0,1)), float(M(0,2)), 
          float(M(1,0)), float(M(1,1)), float(M(1,2)), 
          float(M(2,0)), float(M(2,1)), float(M(2,2)) };
        const Math::Quaternion Q (R);

        put<float32> (Q[1], &NH->quatern_b, is_BE);
        put<float32> (Q[2], &NH->quatern_c, is_BE);
        put<float32> (Q[3], &NH->quatern_d, is_BE);

        put<float32> (M(0,3), &NH->qoffset_x, is_BE);
        put<float32> (M(1,3), &NH->qoffset_y, is_BE);
        put<float32> (M(2,3), &NH->qoffset_z, is_BE);

        // sform:
        put<float32> (H.axes.vox[0]*M(0,0), &NH->srow_x[0], is_BE);
        put<float32> (H.axes.vox[1]*M(0,1), &NH->srow_x[1], is_BE);
        put<float32> (H.axes.vox[2]*M(0,2), &NH->srow_x[2], is_BE);
        put<float32> (M(0,3), &NH->srow_x[3], is_BE);

        put<float32> (H.axes.vox[0]*M(1,0), &NH->srow_y[0], is_BE);
        put<float32> (H.axes.vox[1]*M(1,1), &NH->srow_y[1], is_BE);
        put<float32> (H.axes.vox[2]*M(1,2), &NH->srow_y[2], is_BE);
        put<float32> (M(1,3), &NH->srow_y[3], is_BE);

        put<float32> (H.axes.vox[0]*M(2,0), &NH->srow_z[0], is_BE);
        put<float32> (H.axes.vox[1]*M(2,1), &NH->srow_z[1], is_BE);
        put<float32> (H.axes.vox[2]*M(2,2), &NH->srow_z[2], is_BE);
        put<float32> (M(2,3), &NH->srow_z[3], is_BE);


        strncpy ((gchar*) &NH->magic, "n+1\0", 4);
        fmap.unmap();

        dmap.add_gz (fmap, gzfilename, 352);
      }

    }
  }
}

