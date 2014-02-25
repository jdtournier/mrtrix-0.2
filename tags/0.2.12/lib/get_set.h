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
    * replace get::T and put::T() methods with template get<T>() & put<T>() methods
    * add get/put template specialisations for bool, int8 and uint8
    * remove obsolete ArrayXX classes
    * move MR::ByteOrder namespace & methods from lib/mrtrix.h to here

*/

#ifndef __get_set_h__
#define __get_set_h__

/** \defgroup Binary Binary access functions
 * \brief functions to provide easy access to binary data. */

#include "mrtrix.h"

#define BITMASK 0x01U << 7

#if G_BYTE_ORDER == G_BIG_ENDIAN
#define MRTRIX_IS_BIG_ENDIAN true
#define TO_LE(v) swap(v)
#define TO_BE(v) v
#else 
#define MRTRIX_IS_BIG_ENDIAN false
#define TO_LE(v) v
#define TO_BE(v) swap(v)
#endif 

namespace MR {

  /** \addtogroup Binary 
   * @{ */

  namespace ByteOrder {

    inline gint16    swap (gint16 v)          { return (GUINT16_SWAP_LE_BE  (v)); }
    inline guint16   swap (guint16 v)         { return (GUINT16_SWAP_LE_BE (v)); }
    inline gint32    swap (gint32 v)          { return (GUINT32_SWAP_LE_BE  (v)); }
    inline guint32   swap (guint32 v)         { return (GUINT32_SWAP_LE_BE (v)); }
    inline float32   swap (float32 v)
    {
      union { float32 f; guint32 i; } val = { v };
      val.i = GUINT32_SWAP_LE_BE (val.i);
      return (val.f);
    }

    inline float64 swap (float64 v)
    {
      union { float64 f; guint32 i[2]; } val = { v };
      val.i[1] = swap (val.i[0]);
      val.i[0] = swap (val.i[1]);
      return (val.f);
    }

    inline gint16 LE (gint16 v)     { return (GINT16_TO_LE (v)); }
    inline gint16 BE (gint16 v)     { return (GINT16_TO_BE (v)); }
    inline guint16 LE (guint16 v)   { return (GUINT16_TO_LE (v)); }
    inline guint16 BE (guint16 v)   { return (GUINT16_TO_BE (v)); }
    inline gint32 LE (gint32 v)     { return (GINT32_TO_LE (v)); }
    inline gint32 BE (gint32 v)     { return (GINT32_TO_BE (v)); }
    inline guint32 LE (guint32 v)   { return (GUINT32_TO_LE (v)); }
    inline guint32 BE (guint32 v)   { return (GUINT32_TO_BE (v)); }
    inline float32 LE (float32 v)   { return (TO_LE (v)); }
    inline float32 BE (float32 v)   { return (TO_BE (v)); }
    inline float64 LE (float64 v)   { return (TO_LE (v)); }
    inline float64 BE (float64 v)   { return (TO_BE (v)); }
  }


  template <typename T> inline T getLE (const void* address) { return (ByteOrder::LE (*((T*) address))); }
  template <typename T> inline T getBE (const void* address) { return (ByteOrder::BE (*((T*) address))); }
  template <typename T> inline T get (const void* address, bool is_big_endian = MRTRIX_IS_BIG_ENDIAN)
  { return (is_big_endian ? getBE<T> (address) : getLE<T> (address) ); } 

  template <typename T> inline void putLE (const T value, void* address) { *((T*) address) = ByteOrder::LE (value); }
  template <typename T> inline void putBE (const T value, void* address) { *((T*) address) = ByteOrder::BE (value); }
  template <typename T> inline void put (const T value, void* address, bool is_big_endian = MRTRIX_IS_BIG_ENDIAN)
  { is_big_endian ? putBE<T> (value, address) : putLE<T> (value, address); } 

  template <typename T> inline T getLE (const void* data, gsize i) { return (ByteOrder::LE (((T*) data)[i])); }
  template <typename T> inline T getBE (const void* data, gsize i) { return (ByteOrder::BE (((T*) data)[i])); }
  template <typename T> inline T get (const void* data, gsize i, bool is_big_endian = MRTRIX_IS_BIG_ENDIAN)
  { return (is_big_endian ? getBE<T> (data, i) : getLE<T> (data, i) ); } 

  template <typename T> inline void putLE (const T value, void* data, gsize i) { ((T*) data)[i] = ByteOrder::LE (value); }
  template <typename T> inline void putBE (const T value, void* data, gsize i) { ((T*) data)[i] = ByteOrder::BE (value); }
  template <typename T> inline void put (const T value, void* data, gsize i, bool is_big_endian = MRTRIX_IS_BIG_ENDIAN)
  { is_big_endian ? putBE<T> (value, data, i) : putLE<T> (value, data, i); } 


  template <> inline gint8 get<gint8> (const void* address, bool is_big_endian) { return (*((gint8*) address)); }
  template <> inline gint8 get<gint8> (const void* data, gsize i, bool is_big_endian) { return (((gint8 *) data)[i]); }

  template <> inline void put<gint8> (const gint8 value, void* address, bool is_big_endian) { *((gint8*) address) = value; }
  template <> inline void put<gint8> (const gint8 value, void* data, gsize i, bool is_big_endian) { ((gint8 *) data)[i] = value; }

  template <> inline guint8 get<guint8> (const void* address, bool is_big_endian) { return (*((guint8*) address)); }
  template <> inline guint8 get<guint8> (const void* data, gsize i, bool is_big_endian) { return (((guint8 *) data)[i]); }

  template <> inline void put<guint8> (const guint8 value, void* address, bool is_big_endian) { *((guint8*) address) = value; }
  template <> inline void put<guint8> (const guint8 value, void* data, gsize i, bool is_big_endian) { ((guint8*) data)[i] = value; }


  template <> inline bool get<bool> (const void* data, gsize i, bool is_big_endian)
  { return (((((guint8*) data)[i/8] << i%8) & BITMASK) ? true : false); } 

  template <> inline void put<bool> (const bool value, void* data, gsize i, bool is_big_endian)
  { 
    if (value) ((guint8*) data)[i/8] |= (BITMASK >> i%8); 
    else ((guint8*) data)[i/8] &= ~(BITMASK >> i%8); 
  }

  /** @} */

}


#endif
