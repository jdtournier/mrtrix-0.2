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

#ifndef __data_type_h__
#define __data_type_h__

#include "mrtrix.h"

namespace MR {

  class DataType {
    protected:
      guchar dt;

    public:
      DataType ();
      DataType (guchar type);
      DataType (const DataType& DT);

      guchar&                operator() ();
      const guchar&          operator() () const;
      bool                   operator== (guchar type) const;
      bool                   operator!= (guchar type) const;
      bool                   operator== (const DataType DT) const;
      bool                   operator!= (const DataType DT) const;
      const DataType&        operator= (const DataType DT);

      bool                   is (guchar type) const;
      bool                   is_complex () const;
      bool                   is_signed () const;
      bool                   is_little_endian () const;
      bool                   is_big_endian () const;
      void                   set_byte_order_native ();

      void                   parse (const String& spec);
      guint                  bits () const;
      guint                  bytes () const;
      const gchar*           description () const;
      const gchar*           specifier () const;

      void                   set_flag (guchar flag);
      void                   unset_flag (guchar flag);



      static const guchar     ComplexNumber = 0x10U;
      static const guchar     Signed        = 0x20U;
      static const guchar     LittleEndian  = 0x40U;
      static const guchar     BigEndian     = 0x80U;
      static const guchar     Text          = 0xFFU;
      static const guchar     GroupStart    = 0xFEU;
      static const guchar     GroupEnd      = 0xFDU;

      static const guchar     Undefined     = 0x00U;
      static const guchar     Bit           = 0x01U;
      static const guchar     UInt8         = 0x02U;
      static const guchar     UInt16        = 0x03U;
      static const guchar     UInt32        = 0x04U;
      static const guchar     Float32       = 0x05U;
      static const guchar     Float64       = 0x06U;

      static const guchar     Int8          = UInt8  | Signed;
      static const guchar     Int16         = UInt16 | Signed;
      static const guchar     Int16LE       = UInt16 | Signed | LittleEndian;
      static const guchar     UInt16LE      = UInt16 | LittleEndian;
      static const guchar     Int16BE       = UInt16 | Signed | BigEndian;
      static const guchar     UInt16BE      = UInt16 | BigEndian;
      static const guchar     Int32         = UInt32 | Signed;
      static const guchar     Int32LE       = UInt32 | Signed | LittleEndian;
      static const guchar     UInt32LE      = UInt32 | LittleEndian;
      static const guchar     Int32BE       = UInt32 | Signed | BigEndian;
      static const guchar     UInt32BE      = UInt32 | BigEndian;
      static const guchar     Float32LE     = Float32 | LittleEndian;
      static const guchar     Float32BE     = Float32 | BigEndian;
      static const guchar     Float64LE     = Float64 | LittleEndian;
      static const guchar     Float64BE     = Float64 | BigEndian;
      static const guchar     CFloat32      = ComplexNumber | Float32;
      static const guchar     CFloat32LE    = ComplexNumber | Float32 | LittleEndian;
      static const guchar     CFloat32BE    = ComplexNumber | Float32 | BigEndian;
      static const guchar     CFloat64      = ComplexNumber | Float64;
      static const guchar     CFloat64LE    = ComplexNumber | Float64 | LittleEndian;
      static const guchar     CFloat64BE    = ComplexNumber | Float64 | BigEndian;

      static const guchar     Native        = Float32 | 
#if G_BYTE_ORDER == G_BIG_ENDIAN
        BigEndian;
#else
        LittleEndian;
#endif
  };








  inline DataType::DataType () :
#ifdef BYTE_ORDER_BIG_ENDIAN
      dt (DataType::Float32BE)
#else
      dt (DataType::Float32LE)
#endif
  { }

  inline DataType::DataType (guchar type) : dt (type)               { }
  inline DataType::DataType (const DataType& DT) : dt (DT.dt)       { }
  inline guchar& DataType::operator() ()                            { return (dt); }
  inline const guchar& DataType::operator() () const                { return (dt); }
  inline bool DataType::operator== (guchar type) const              { return (dt == type); }
  inline bool DataType::operator!= (guchar type) const              { return (dt != type); }
  inline bool DataType::operator== (const DataType DT) const        { return (dt == DT.dt); }
  inline bool DataType::operator!= (const DataType DT) const        { return (dt != DT.dt); }
  inline const DataType& DataType::operator= (const DataType DT)    { dt = DT.dt; return (*this); }
  inline guint DataType::bytes () const                             { return ((bits()+7)/8); }
  inline bool DataType::is (guchar type) const                      { return (dt == type); }
  inline bool DataType::is_complex () const                         { return (dt & ComplexNumber); }
  inline bool DataType::is_signed () const                          { return (dt & Signed); }
  inline bool DataType::is_little_endian () const                   { return (dt & LittleEndian); }
  inline bool DataType::is_big_endian() const                       { return (dt & BigEndian); }
  inline void DataType::set_flag (guchar flag)                      { dt |= flag; }
  inline void DataType::unset_flag (guchar flag)                    { dt &= ~flag; }
  inline void DataType::set_byte_order_native ()
  {
    if ( dt != Bit && dt != Int8 && dt != UInt8 ) {
      if ( !is_little_endian() && !is_big_endian() ) {
#if G_BYTE_ORDER == G_BIG_ENDIAN
        dt |= BigEndian;
#else
        dt |= LittleEndian;
#endif
      }
    }
  }

}

#endif

