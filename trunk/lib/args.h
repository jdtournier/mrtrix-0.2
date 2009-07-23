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

#ifndef __args_h__
#define __args_h__

#include "ptr.h"
#include "image/object.h"

namespace MR {
  namespace Image { 
    class Object; 
    class Header;
  }

  namespace Dialog {
    class Window;
    class Argument;
  }

  typedef enum {
    Undefined,
    Integer,
    Float,
    Text,
    ArgFile,
    Choice,
    ImageIn,
    ImageOut,
    IntSeq,
    FloatSeq
  } ArgType;

  const gchar* argument_type_description (ArgType type);


  class ArgData {
    public:
      ArgData ();
      ArgType   type;
      union {
        int         i;
        float       f;
        const gchar* string;
      } data;

      RefPtr<Image::Object>  image;
      friend std::ostream& operator<< (std::ostream& stream, const ArgData& a)
      {
        stream << "{ arg type: " << a.type << " }";
        return (stream);
      }
  };





  class Argument;

  class ArgBase {
    protected:
      RefPtr<ArgData> data;

    public:
      ArgBase ();
      ArgBase (const Argument& arg, const gchar* string);

      int                 get_int () const;
      float               get_float () const;
      const gchar*         get_string () const;
      RefPtr<Image::Object>  get_image () const;
      RefPtr<Image::Object>  get_image (Image::Header& header) const;

      ArgType              type () const;

      friend std::ostream& operator<< (std::ostream& stream, const ArgBase& arg);
      friend class Dialog::Window;
      friend class Dialog::Argument;
      friend class Image::Object;
  };





  class OptBase : public std::vector<ArgBase> {
    public:
      guint index;

      friend std::ostream& operator<< (std::ostream& stream, const OptBase& opt);
      friend class Dialog::Window;
  };



  class Argument {
    friend class ArgBase;
    public:
      Argument ();
      Argument (const gchar* Short_Name, 
          const gchar* Long_Name, 
          const gchar* Description, 
          bool Mandatory = true, 
          bool Allow_Multiple = false);

      const gchar* sname;
      const gchar* lname;
      const gchar* desc;
      bool mandatory, allow_multiple;
      ArgType type;

      union {
        struct { int def, min, max; }     i;
        struct { float def, min, max; }   f;
        const gchar*                       string;
        const gchar**                      choice;
      } extra_info;

      const Argument& type_integer (int lowest, int highest, int default_value);
      const Argument& type_float   (float lowest, float highest, float default_value);
      const Argument& type_string  (const gchar* default_value = NULL);
      const Argument& type_file    ();
      const Argument& type_image_in ();
      const Argument& type_image_out ();
      const Argument& type_choice  (const gchar** choices);
      const Argument& type_sequence_int ();
      const Argument& type_sequence_float ();

      bool    is_valid () const;

      static const Argument End;
      friend std::ostream& operator<< (std::ostream& stream, const Argument& arg);
  };





  class Option : public std::vector<Argument> {
    public:
      Option ();
      Option (const gchar* Short_Name, 
          const gchar* Long_Name,
          const gchar* Description,
          bool Mandatory = false, 
          bool Allow_Multiple = false);

      const gchar* sname;
      const gchar* lname;
      const gchar* desc;
      bool mandatory, allow_multiple;

      Option& append (const Argument& argument);
      bool    is_valid () const;

      friend std::ostream& operator<< (std::ostream& stream, const Option& opt);
      friend std::ostream& operator<< (std::ostream& stream, const OptBase& opt);

      static const Option End;
  };









  inline ArgData::ArgData () : type (Undefined) { data.string = NULL; }

  inline ArgBase::ArgBase () { }

  inline int         ArgBase::get_int () const     { return (data->data.i); }
  inline float       ArgBase::get_float () const   { return (data->data.f); }
  inline const gchar* ArgBase::get_string () const  { return (data->data.string); }
  inline RefPtr<Image::Object> ArgBase::get_image () const   { assert (type() == ImageIn); return (data->image); }
  inline RefPtr<Image::Object> ArgBase::get_image (Image::Header& header) const
  {
    assert (type() == ImageOut);
    data->image->create (get_string(), header);
    return (data->image); 
  }


  inline ArgType ArgBase::type () const { return (!data ? Undefined : data->type ); }







  inline Argument::Argument () : type (Undefined) { sname = lname = desc = NULL; } 


  inline Argument::Argument (
      const gchar* Short_Name,
      const gchar* Long_Name,
      const gchar* Description,
      bool Mandatory,
      bool Allow_Multiple) :
    sname (Short_Name),
    lname (Long_Name),
    desc (Description),
    mandatory (Mandatory),
    allow_multiple (Allow_Multiple)
  {
  }



  inline bool Argument::is_valid () const { return (sname); }




  inline const Argument& Argument::type_integer (int lowest, int highest, int default_value)
  {
    type = Integer;
    extra_info.i.def = default_value;
    extra_info.i.min = lowest;
    extra_info.i.max = highest;
    return (*this);
  }

  inline const Argument& Argument::type_float (float lowest, float highest, float default_value)
  {
    type = Float;
    extra_info.f.def = default_value;
    extra_info.f.min = lowest;
    extra_info.f.max = highest;
    return (*this);
  }

  inline const Argument& Argument::type_string (const gchar* default_value)
  {
    type = Text;
    extra_info.string = default_value;
    return (*this);
  }

  inline const Argument& Argument::type_choice (const gchar** choices)
  {
    type = Choice;
    extra_info.choice = choices;
    return (*this);
  }

  inline const Argument& Argument::type_file ()
  {
    type = ArgFile;
    extra_info.string = NULL;
    return (*this);
  }

  inline const Argument& Argument::type_image_in ()
  {
    type = ImageIn;
    extra_info.string = NULL;
    return (*this);
  }

  inline const Argument& Argument::type_image_out ()
  {
    type = ImageOut;
    extra_info.string = NULL;
    return (*this);
  }

  inline const Argument& Argument::type_sequence_int ()
  {
    type = IntSeq;
    extra_info.string = NULL;
    return (*this);
  }

  inline const Argument& Argument::type_sequence_float ()
  {
    type = FloatSeq;
    extra_info.string = NULL;
    return (*this);
  }






  inline Option::Option () { sname = lname = desc = NULL; } 

  inline Option::Option (
      const gchar* Short_Name, 
      const gchar* Long_Name,
      const gchar* Description,
      bool Mandatory,
      bool Allow_Multiple) :
    sname (Short_Name),
    lname (Long_Name),
    desc (Description),
    mandatory (Mandatory),
    allow_multiple (Allow_Multiple)
  {
  }

  inline Option& Option::append (const Argument& argument)
  {
    std::vector<Argument>::push_back (argument);
    return (*this);
  }

  inline bool Option::is_valid () const { return (sname); }
  
}

#endif

