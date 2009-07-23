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

#ifndef __image_position_h__
#define __image_position_h__

#include "image/object.h"
#include "math/complex_number.h"

namespace MR {
  namespace Image {

    class Object;

    //! \addtogroup Image
    // @{

    //! This class provides access to the voxel intensities of an image.
    class Position {
      public:
        //! construct a Position object to point to the data contained in the MR::Image::Object \p parent
        /*! All coordinates will be initialised to zero. */
        explicit Position (Object& parent) : image (parent), offset (image.start), stride (image.stride) { memset (x, 0, ndim()*sizeof(int)); image.map(); }

        Object&     image; //!< The MR::Image::Object containing the image data.

        //! test whether current position is within bounds.
        /*! \return true if the current position is out of bounds, false otherwise */
        bool        operator! () const { for (int n = 0; n < ndim(); n++) if (x[n] < 0 || x[n] >= dim(n)) return (true); return (false); }

        Position&   operator= (const Position& pos) { return (operator= (pos.x)); }
        Position&   operator= (const int pos[MRTRIX_MAX_NDIMS]);

        //! used to loop over all image coordinates.
        /*! This operator increments the current position to the next voxel, incrementing to the next axis as required.
         * It is used to process all voxels in turn. For example:
         * \code
         * MR::Image::Position position (image_object);
         * do {
         *   process (position.value());
         * } while (position++);
         * \endcode
         * \return true once the last voxel has been reached (i.e. the next increment would bring the current position out of bounds), false otherwise. */
        bool        operator++ (int notused);

        //! reset all coordinates to zero. 
        void        zero () { offset = image.start; memset (x, 0, ndim()*sizeof(int)); }


        //! %set individual coordinate components.
        /*! Sets the coordinate along the axis specified to the position specified. */
        void        set (guint axis, int position)    { offset += stride[axis] * gssize(position - x[axis]); x[axis] = position; }

        //! increment individual coordinate components.
        /*! Increment the coordinate along the axis specified by one. */
        void        inc (guint axis)                  { offset += stride[axis]; x[axis]++; }
        
        //! modify individual coordinate components.
        /*! Move the position along the axis specified by the amount specified in \p increment. */
        void        move (guint axis, int increment)  { offset += stride[axis]*gssize(increment); x[axis] += increment; } 

        //! return the coordinate along the specified axis.
        int         operator[] (guint axis) const     { return (x[axis]); }

        //! %get the image name.
        /*! \return the name of the parent image */
        const String& name () const         { return (image.name()); }

        //! %get the image dimension along the axis specified.
        /*! \return the dimension of the parent image along axis \p axis */
        int         dim (guint index) const { return (image.dim (index)); }

        //! %get the number of dimensions of the image.
        /*! \return the number of dimensions of the parent image */
        int         ndim () const           { return (image.ndim()); }

        //! %get the voxel size along the axis specified.
        /*! \return the voxel size of the parent image along the axis specified */
        float       vox (guint index) const { return (image.vox (index)); }

        //! %get whether the image data is complex
        /*! \return true if the image data are complex */
        bool        is_complex () const      { return (image.is_complex()); }

        //! %get the number of voxels in the image
        /*! By default, the total number of voxels is returned.
         * if \p up_to_dim is specified, only dimensions up to the value specified are included in the calculation.
         * \return the number of voxels contained in a volume of dimensionality \p up_to_dim */
        gsize       voxel_count (guint up_to_dim = MRTRIX_MAX_NDIMS) const { return (image.voxel_count (up_to_dim)); }

        //! %get the value of the voxel at the current position
        float       value () const           { return (re()); }

        //! %set the value of the voxel at the current position
        void        value (float val)        { re (val); }

        //! %get the real component of the voxel at the current position
        float       re () const              { return (image.re (offset)); }

        //! %set the real component of the voxel at the current position
        void        re (float val)           { image.re (offset, val); }

        //! %get the imaginary component of the voxel at the current position
        /*! \note No check is performed to ensure the image is actually complex. Calling this function on real-valued data will produce undefined results */
        float       im () const              { return (image.im (offset)); }

        //! %set the imaginary component of the voxel at the current position
        /*! \note No check is performed to ensure the image is actually complex. Calling this function on real-valued data will produce undefined results */
        void        im (float val)           { image.im (offset, val); }

        //! %get the complex value stored at the current position
        /*! \note No check is performed to ensure the image is actually complex. Calling this function on real-valued data will produce undefined results */
        Math::ComplexNumber<float> Z () const          { return (Math::ComplexNumber<float> (re(),im())); }

        //! %set the complex value stored at the current position
        /*! \note No check is performed to ensure the image is actually complex. Calling this function on real-valued data will produce undefined results */
        void        Z (const Math::ComplexNumber<float>& val) { re (val.re()); im (val.im()); }
        //
        //! %set the complex value stored at the current position
        /*! \note No check is performed to ensure the image is actually complex. Calling this function on real-valued data will produce undefined results */
        void        Z (float val_re, float val_im)     { re (val_re); im (val_im); }

        //! %get the voxel data stored at the current position
        /*! This sets the parameters \p val and \p val_im using the voxel data at the current image position, according the \p format speficier.
         * \note If \p format refers to a complex data type, no check is performed to ensure the image is complex. 
         * In this case, calling this function on real-valued data will produce undefined results */
        void        get (OutputType format, float& val, float& val_im);


      protected:
        int       x[MRTRIX_MAX_NDIMS]; //!< the current image coordinates
        gsize     offset; //!< the offset in memory to the current voxel
        const gssize* stride; //!< the offset in memory between adjacent image voxels along each axis

        friend class Entry;
        friend std::ostream& operator<< (std::ostream& stream, const Position& pos);
    };

    //! @}





















    inline Position& Position::operator= (const int pos[MRTRIX_MAX_NDIMS])
    {
      memcpy (x, pos, image.ndim()*sizeof(int));
      gssize shift = 0;
      for (int n = 0; n < image.ndim(); n++) shift += image.stride[n] * gssize(x[n]);
      offset = image.start + shift;
      return (*this);
    }



    inline bool Position::operator++ (int notused)
    {
      int axis = 0;
      do {
        inc (axis);
        if (x[axis] < image.dim(axis)) return (true);
        set (axis, 0);
        axis++;
      } while (axis < image.ndim());
      return (false);
    }



    inline void   Position::get (OutputType format, float& val, float& val_im)
    {
      switch (format) {
        case Default:   val = value(); return;
        case Real:      val = re(); return;
        case Imaginary: val = im(); return;
        case Magnitude: val = Z().mod(); return;
        case Phase:     val = Z().phase(); return;
        case RealImag:  val = re(); val_im = im(); return;
      }
      assert (false);
    }



    inline std::ostream& operator<< (std::ostream& stream, const Position& pos)
    {
      stream << "position for image \"" << pos.image.name() << "\" = [ ";
      for (int n = 0; n < pos.image.ndim(); n++) stream << pos.x[n] << " ";
      stream << "]\n  current offset = " << pos.offset;
      return (stream);
    }


  }
}

#endif


