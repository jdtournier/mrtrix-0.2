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

#ifndef __image_interp_h__
#define __image_interp_h__

#include "image/position.h"
#include "image/object.h"
#include "point.h"

namespace MR {
  namespace Image {

    //! \addtogroup Image
    // @{

    //! This class provides access to the voxel intensities of an image, using tri-linear interpolation.
    /*! Interpolation is only performed along the first 3 (spatial) axes. 
     * The (integer) position along the remaining axes can be set using the MR::Image::Position functions. 
     * The spatial coordinates can be set using the functions P(), I(), and R(). 
     * For example:
     * \code
     * Image::Interp interp (image_object); // create an Interp object using image_object as the parent Image::Object
     * interp.R (10.2, 3.59, 54.1);         // set the real-space position to [ 10.2 3.59 54.1 ]
     * interp.set (3, 0);                   // set the fourth dimension (e.g. time) to zero
     * float value = interp.value();        // get the value at this position
     * \endcode */

    class Interp : public Position {
      public:
        //! construct an Interp object to point to the data contained in the MR::Image::Object \p parent
        /*! All non-spatial coordinates (i.e. axes > 3) will be initialised to zero. 
         * The spatial coordinates are left undefined. */
        Interp (Object& parent);
        ~Interp() { }

        Interp&   operator= (const int i[MRTRIX_MAX_NDIMS]) { Position::operator= (i); return (*this); }

        //! test whether current position is within bounds.
        /*! \return true if the current position is out of bounds, false otherwise */
        bool  operator! () const { return (out_of_bounds); }

        //! Transform the position \p r from real-space to pixel-space
        Point R2P (const Point& r) const;
        //! Transform the position \p r from pixel-space to real-space
        Point P2R (const Point& r) const;
        //! Transform the position \p r from image-space to pixel-space
        Point I2P (const Point& r) const { return (Point (r[0]/vox(0), r[1]/vox(1), r[2]/vox(2))); }
        //! Transform the position \p r from pixel-space to image-space
        Point P2I (const Point& r) const { return (Point (r[0]*vox(0), r[1]*vox(1), r[2]*vox(2))); }

        //! Transform the orientation \p r from real-space to pixel-space
        Point vec_R2P (const Point& r) const;
        //! Transform the orientation \p r from pixel-space to real-space
        Point vec_P2R (const Point& r) const;

        bool P (const Point& pos);
        bool I (const Point& pos) { return (P (I2P (pos))); }
        bool R (const Point& pos) { return (P (R2P (pos))); }

        float value () const { return (Interp::re()); }
        float re () const;
        float im () const;

        float re_abs () const;
        float im_abs () const;

        void  get (OutputType format, float& val, float& val_im);
        void  abs (OutputType format, float& val, float& val_im);

        class Map {
          public:
            Map () { }
            Map (const Image::Object& dest, const Image::Object& source) { set (dest, source); }
            void set (const Image::Object& dest, const Image::Object& source) 
            {
              Math::Matrix M;
              M.multiply (dest.header().R2P(), source.header().P2R());
              R[0] = M(0,0); R[1] = M(0,1); R[2] = M(0,2); R[3] = M(0,3); 
              R[4] = M(1,0); R[5] = M(1,1); R[6] = M(1,2); R[7] = M(1,3);
              R[8] = M(2,0); R[9] = M(2,1); R[10]= M(2,2); R[11]= M(2,3);
            }
            Point operator() (const Point& P) const 
            {
              return (Point (
                    R[0]*P[0] + R[1]*P[1] + R[2]*P[2] + R[3],
                    R[4]*P[0] + R[5]*P[1] + R[6]*P[2] + R[7],
                    R[8]*P[0] + R[9]*P[1] + R[10]*P[2] + R[11] ));
            }
            Point operator() (const Position& P) const { return (operator() (Point (P[0], P[1], P[2]))); }
          protected:
            float R[12];

        };

      protected:
        float         PR[3][4], RP[3][4];
        float         bounds[3];
        bool          out_of_bounds;
        float         faaa, faab, faba, fabb, fbaa, fbab, fbba, fbbb;

        Point         set_fractions (const Point& pos);
    };

    //! @}










    inline Point Interp::R2P (const Point& r) const
    {
      return (Point (
        RP[0][0]*r[0] + RP[0][1]*r[1] + RP[0][2]*r[2] + RP[0][3],
        RP[1][0]*r[0] + RP[1][1]*r[1] + RP[1][2]*r[2] + RP[1][3],
        RP[2][0]*r[0] + RP[2][1]*r[1] + RP[2][2]*r[2] + RP[2][3] ));
    }



    
    inline Point Interp::P2R (const Point& r) const
    {
      return (Point (
        PR[0][0]*r[0] + PR[0][1]*r[1] + PR[0][2]*r[2] + PR[0][3],
        PR[1][0]*r[0] + PR[1][1]*r[1] + PR[1][2]*r[2] + PR[1][3],
        PR[2][0]*r[0] + PR[2][1]*r[1] + PR[2][2]*r[2] + PR[2][3] ));
    }


    inline Point Interp::vec_R2P (const Point& r) const
    {
      return (Point (
        RP[0][0]*r[0] + RP[0][1]*r[1] + RP[0][2]*r[2],
        RP[1][0]*r[0] + RP[1][1]*r[1] + RP[1][2]*r[2],
        RP[2][0]*r[0] + RP[2][1]*r[1] + RP[2][2]*r[2] ));
    }


    inline Point Interp::vec_P2R (const Point& r) const
    {
      return (Point (
        PR[0][0]*r[0] + PR[0][1]*r[1] + PR[0][2]*r[2],
        PR[1][0]*r[0] + PR[1][1]*r[1] + PR[1][2]*r[2],
        PR[2][0]*r[0] + PR[2][1]*r[1] + PR[2][2]*r[2] ));
    }






    inline Point Interp::set_fractions (const Point& pos)
    {
      if (pos[0] < -0.5 || pos[0] > bounds[0] || 
          pos[1] < -0.5 || pos[1] > bounds[1] || 
          pos[2] < -0.5 || pos[2] > bounds[2]) {
        out_of_bounds = true;
        return (Point (GSL_NAN, GSL_NAN, GSL_NAN));
      }

      out_of_bounds = false;
      set (0, int (pos[0]));
      set (1, int (pos[1]));
      set (2, int (pos[2]));

      return (Point (pos[0]-x[0], pos[1]-x[1], pos[2]-x[2]));
    }






    inline bool Interp::P (const Point& pos)
    {
      Point f = set_fractions (pos);
      if (out_of_bounds) return (true);

      if (pos[0] < 0.0) { f[0] = 0.0; set (0,0); }
      else if (pos[0] > bounds[0]-0.5) f[0] = 0.0;

      if (pos[1] < 0.0) { f[1] = 0.0; set (1,0); }
      else if (pos[1] > bounds[1]-0.5) f[1] = 0.0;

      if (pos[2] < 0.0) { f[2] = 0.0; set (2,0); }
      else if (pos[2] > bounds[2]-0.5) f[2] = 0.0;

      faaa = (1.0-f[0]) * (1.0-f[1]) * (1.0-f[2]); if (faaa < 1e-6) faaa = 0.0;
      faab = (1.0-f[0]) * (1.0-f[1]) *      f[2];  if (faab < 1e-6) faab = 0.0;
      faba = (1.0-f[0]) *      f[1]  * (1.0-f[2]); if (faba < 1e-6) faba = 0.0;
      fabb = (1.0-f[0]) *      f[1]  *      f[2];  if (fabb < 1e-6) fabb = 0.0;
      fbaa =      f[0]  * (1.0-f[1]) * (1.0-f[2]); if (fbaa < 1e-6) fbaa = 0.0;
      fbab =      f[0]  * (1.0-f[1])      * f[2];  if (fbab < 1e-6) fbab = 0.0;
      fbba =      f[0]  *      f[1]  * (1.0-f[2]); if (fbba < 1e-6) fbba = 0.0;
      fbbb =      f[0]  *      f[1]  *      f[2];  if (fbbb < 1e-6) fbbb = 0.0;

      return (false);
    }









    inline float Interp::re () const
    {
      if (out_of_bounds) return (GSL_NAN);
      float val = 0.0;
      gsize os (offset);
      if (faaa) val  = faaa * image.re (os); os += stride[2];
      if (faab) val += faab * image.re (os); os += stride[1];
      if (fabb) val += fabb * image.re (os); os -= stride[2];
      if (faba) val += faba * image.re (os); os += stride[0];
      if (fbba) val += fbba * image.re (os); os -= stride[1];
      if (fbaa) val += fbaa * image.re (os); os += stride[2];
      if (fbab) val += fbab * image.re (os); os += stride[1];
      if (fbbb) val += fbbb * image.re (os);
      return (val);
    }





    inline float Interp::im () const
    {
      if (out_of_bounds) return (GSL_NAN);
      float val = 0.0;
      gsize os (offset);
      if (faaa) val  = faaa * image.im (os); os += stride[2];
      if (faab) val += faab * image.im (os); os += stride[1];
      if (fabb) val += fabb * image.im (os); os -= stride[2];
      if (faba) val += faba * image.im (os); os += stride[0];
      if (fbba) val += fbba * image.im (os); os -= stride[1];
      if (fbaa) val += fbaa * image.im (os); os += stride[2];
      if (fbab) val += fbab * image.im (os); os += stride[1];
      if (fbbb) val += fbbb * image.im (os);
      return (val);
    }




    inline float Interp::re_abs () const
    {
      if (out_of_bounds) return (GSL_NAN);
      float val = 0.0;
      gsize os (offset);
      if (faaa) val  = faaa * fabs (image.re (os)); os += stride[2];
      if (faab) val += faab * fabs (image.re (os)); os += stride[1];
      if (fabb) val += fabb * fabs (image.re (os)); os -= stride[2];
      if (faba) val += faba * fabs (image.re (os)); os += stride[0];
      if (fbba) val += fbba * fabs (image.re (os)); os -= stride[1];
      if (fbaa) val += fbaa * fabs (image.re (os)); os += stride[2];
      if (fbab) val += fbab * fabs (image.re (os)); os += stride[1];
      if (fbbb) val += fbbb * fabs (image.re (os));
      return (val);
    }





    inline float Interp::im_abs () const
    {
      if (out_of_bounds) return (GSL_NAN);
      float val = 0.0;
      gsize os (offset);
      if (faaa) val  = faaa * fabs (image.im (os)); os += stride[2];
      if (faab) val += faab * fabs (image.im (os)); os += stride[1];
      if (fabb) val += fabb * fabs (image.im (os)); os -= stride[2];
      if (faba) val += faba * fabs (image.im (os)); os += stride[0];
      if (fbba) val += fbba * fabs (image.im (os)); os -= stride[1];
      if (fbaa) val += fbaa * fabs (image.im (os)); os += stride[2];
      if (fbab) val += fbab * fabs (image.im (os)); os += stride[1];
      if (fbbb) val += fbbb * fabs (image.im (os));
      return (val);
    }




    inline void Interp::get (OutputType format, float& val, float& val_im)
    {
      if (out_of_bounds) { val = val_im = GSL_NAN; return; }
      switch (format) {
        case Default:   val = re(); return;
        case Real:      val = re(); return;
        case Imaginary: val = im(); return;
        case Magnitude: val = re(); val_im = im(); val = sqrt (val*val + val_im*val_im); return;
        case Phase:     val = re(); val_im = im(); val = atan2 (val_im, val); return;
        case RealImag:  val = re(); val_im = im(); return;
      }
      assert (false);
    }





    inline void Interp::abs (OutputType format, float& val, float& val_im)
    {
      if (out_of_bounds) { val = val_im = GSL_NAN; return; }
      switch (format) {
        case Default:   val = re_abs(); return;
        case Real:      val = re_abs(); return;
        case Imaginary: val = im_abs(); return;
        case Magnitude: val = re_abs(); val_im = im_abs(); val = sqrt (val*val + val_im*val_im); return;
        case Phase:     val = re_abs(); val_im = im_abs(); val = atan2 (val_im, val); return;
        case RealImag:  val = re_abs(); val_im = im_abs(); return;
      }
      assert (false);
    }


  }
}

#endif

