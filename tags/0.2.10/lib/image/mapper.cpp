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


    31-10-2008 J-Donald Tournier <d.tournier@brain.org.au>
    * use template get<T>() & put<T>() methods from lib/get_set.h
*/

#include "image/mapper.h"
#include "app.h"
#include "get_set.h"

#define DATAMAPPER_MAX_FILES 128

namespace MR {
  namespace Image {

    namespace {

      inline gsize calc_segsize (const Header& H, guint nfiles) 
      {
        gsize segsize = H.data_type.is_complex() ? 2 : 1;
        for (int i = 0; i < H.axes.ndim(); i++) segsize *= H.axes.dim[i]; 
        segsize /= nfiles;
        return (segsize);
      }

    }








    void Mapper::map (const Header& H)
    {
      debug ("mapping image \"" + H.name + "\"...");
      assert (list.size() || mem);
      assert (segment == NULL);

      if (list.size() > DATAMAPPER_MAX_FILES || 
          ( optimised && ( list.size() > 1 || H.data_type != DataType::Native )) ) {

        if (H.data_type == DataType::Bit) optimised = true;

        info (String ("loading ") + ( optimised ? "and optimising " : "" ) + "image \"" + H.name + "\"..."); 

        bool read_only = list[0].fmap.is_read_only();

        gsize bpp = optimised ? sizeof (float32) : H.data_type.bytes();
        mem = new guint8 [bpp*H.voxel_count()];
        if (!mem) throw Exception ("failed to allocate memory for image data!");

        if (files_new) memset (mem, 0, bpp*H.voxel_count());
        else {
          segsize = calc_segsize (H, list.size());

          for (guint n = 0; n < list.size(); n++) {
            list[n].fmap.map (); 

            if (optimised) {
              float32* data = (float32*) mem + n*segsize;
              guint8*   fdata = list[n].start();
              for (gsize i = 0; i < segsize; i++) 
                data[i] = get_func (fdata, i); 
            } 
            else memcpy (mem + n*segsize*bpp, list[n].start(), segsize*bpp);

            list[n].fmap.unmap();
          }
        }

        if (temporary || read_only) list.clear();
      }

      if (mem) {
        segment = new guint8* [1];
        segment[0] = mem;
        segsize = optimised ? sizeof (float32) : H.data_type.bytes();
        segsize *= H.voxel_count();
      }
      else {
        segment = new guint8* [list.size()];
        for (guint n = 0; n < list.size(); n++) {
          list[n].fmap.map();
          segment[n] = list[n].start();
        }
        segsize = calc_segsize (H, list.size());
      }


      debug ("data mapper for image \"" + H.name + "\" mapped with segment size = " + str (segsize) 
          + ( optimised ? " (optimised)" : "")); 
    }





    void Mapper::unmap (const Header& H)
    {
      if (mem && list.size()) {
        segsize = calc_segsize (H, list.size());
        if (!optimised) segsize *= H.data_type.bytes();

        info ("writing back data for image \"" + H.name + "\"...");
        for (guint n = 0; n < list.size(); n++) {
          bool err = false;
          try { list[n].fmap.map (); }
          catch (...) { err = true; error ("error writing data to file \"" + list[n].fmap.name() + "\""); }
          if (!err) {
            if (optimised) {
              const float32* data = (const float32*) mem + n*segsize;
              for (gsize i = 0; i < segsize; i++) 
                put_func (data[i], list[n].start(), i); 
            } 
            else memcpy (list[n].start(), mem + n*segsize, segsize);
            list[n].fmap.unmap();
          }
        }
      }

      delete [] mem;
      delete [] segment;
      mem = NULL;
      segment = NULL;
    }





    void Mapper::set_data_type (DataType dt)
    {
      switch (dt() & ~DataType::ComplexNumber) {
        case DataType::Bit:        get_func = getBit;        put_func = putBit;        return;
        case DataType::Int8:       get_func = getInt8;       put_func = putInt8;       return;
        case DataType::UInt8:      get_func = getUInt8;      put_func = putUInt8;      return;
        case DataType::Int16LE:    get_func = getInt16LE;    put_func = putInt16LE;    return;
        case DataType::UInt16LE:   get_func = getUInt16LE;   put_func = putUInt16LE;   return;
        case DataType::Int16BE:    get_func = getInt16BE;    put_func = putInt16BE;    return;
        case DataType::UInt16BE:   get_func = getUInt16BE;   put_func = putUInt16BE;   return;
        case DataType::Int32LE:    get_func = getInt32LE;    put_func = putInt32LE;    return;
        case DataType::UInt32LE:   get_func = getUInt32LE;   put_func = putUInt32LE;   return;
        case DataType::Int32BE:    get_func = getInt32BE;    put_func = putInt32BE;    return;
        case DataType::UInt32BE:   get_func = getUInt32BE;   put_func = putUInt32BE;   return;
        case DataType::Float32LE:  get_func = getFloat32LE;  put_func = putFloat32LE;  return;
        case DataType::Float32BE:  get_func = getFloat32BE;  put_func = putFloat32BE;  return;
        case DataType::Float64LE:  get_func = getFloat64LE;  put_func = putFloat64LE;  return;
        case DataType::Float64BE:  get_func = getFloat64BE;  put_func = putFloat64BE;  return;
        default: throw Exception ("invalid data type in image header");
      }
    }










    std::ostream& operator<< (std::ostream& stream, const Mapper& dmap)
    {
      stream << "mapper ";
      if (dmap.optimised) stream << " (optimised)";
      stream << ":\n  segment size = " << dmap.segsize << "\n  ";
      if (!dmap.segment) stream << "(unmapped)\n";
      else if (dmap.mem) stream << "in memory at " << (void*) dmap.mem << "\n";
      stream << "files:\n";
      for (guint i = 0; i < dmap.list.size(); i++) {
        stream << "    " << dmap.list[i].fmap.name() << ", offset " << dmap.list[i].offset << " (";
        if (dmap.list[i].fmap.is_mapped()) stream << "mapped at " << dmap.list[i].fmap.address();
        else stream << "unmapped";
        stream << ( dmap.list[i].fmap.is_read_only() ? ", read-only)\n" : ", read-write)\n" );
      }
      return (stream);
    }



    float32 Mapper::getBit       (const void* data, gsize i)  { return (get<bool>      (data, i)); }
    float32 Mapper::getInt8      (const void* data, gsize i)  { return (get<gint8>     (data, i)); }
    float32 Mapper::getUInt8     (const void* data, gsize i)  { return (get<guint8>    (data, i)); }
    float32 Mapper::getInt16LE   (const void* data, gsize i)  { return (getLE<gint16>  (data, i)); }
    float32 Mapper::getUInt16LE  (const void* data, gsize i)  { return (getLE<guint16> (data, i)); }
    float32 Mapper::getInt16BE   (const void* data, gsize i)  { return (getBE<gint16>  (data, i)); }
    float32 Mapper::getUInt16BE  (const void* data, gsize i)  { return (getBE<guint16> (data, i)); }
    float32 Mapper::getInt32LE   (const void* data, gsize i)  { return (getLE<gint32>  (data, i)); }
    float32 Mapper::getUInt32LE  (const void* data, gsize i)  { return (getLE<guint32> (data, i)); }
    float32 Mapper::getInt32BE   (const void* data, gsize i)  { return (getBE<gint32>  (data, i)); }
    float32 Mapper::getUInt32BE  (const void* data, gsize i)  { return (getBE<guint32> (data, i)); }
    float32 Mapper::getFloat32LE (const void* data, gsize i)  { return (getLE<float32> (data, i)); }
    float32 Mapper::getFloat32BE (const void* data, gsize i)  { return (getBE<float32> (data, i)); }
    float32 Mapper::getFloat64LE (const void* data, gsize i)  { return (getLE<float64> (data, i)); }
    float32 Mapper::getFloat64BE (const void* data, gsize i)  { return (getBE<float64> (data, i)); }

    void Mapper::putBit       (float32 val, void* data, gsize i) { put<bool>      (bool(val), data, i); }
    void Mapper::putInt8      (float32 val, void* data, gsize i) { put<gint8>     (gint8(val), data, i); }
    void Mapper::putUInt8     (float32 val, void* data, gsize i) { put<guint8>    (guint8(val), data, i); }
    void Mapper::putInt16LE   (float32 val, void* data, gsize i) { putLE<gint16>  (gint16(val), data, i); }
    void Mapper::putUInt16LE  (float32 val, void* data, gsize i) { putLE<guint16> (guint16(val), data, i); }
    void Mapper::putInt16BE   (float32 val, void* data, gsize i) { putBE<gint16>  (gint16(val), data, i); }
    void Mapper::putUInt16BE  (float32 val, void* data, gsize i) { putBE<guint16> (guint16(val), data, i); }
    void Mapper::putInt32LE   (float32 val, void* data, gsize i) { putLE<gint32>  (gint32(val), data, i); }
    void Mapper::putUInt32LE  (float32 val, void* data, gsize i) { putLE<guint32> (guint32(val), data, i); }
    void Mapper::putInt32BE   (float32 val, void* data, gsize i) { putBE<gint32>  (gint32(val), data, i); }
    void Mapper::putUInt32BE  (float32 val, void* data, gsize i) { putBE<guint32> (guint32(val), data, i); }
    void Mapper::putFloat32LE (float32 val, void* data, gsize i) { putLE<float32> (float32(val), data, i); }
    void Mapper::putFloat32BE (float32 val, void* data, gsize i) { putBE<float32> (float32(val), data, i); }
    void Mapper::putFloat64LE (float32 val, void* data, gsize i) { putLE<float64> (float64(val), data, i); }
    void Mapper::putFloat64BE (float32 val, void* data, gsize i) { putBE<float64> (float64(val), data, i); }

  }
}



