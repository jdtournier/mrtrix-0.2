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

#ifndef __math_vector_h__
#define __math_vector_h__

#include <gsl/gsl_vector.h>
#include "math/matrix.h"

namespace MR {
  namespace Math {

    class Vector {
      protected:
        gsl_vector*   V;

      public:
        Vector ();
        Vector (guint size);
        Vector (const Vector& vec);
        ~Vector ();

        Vector&       operator= (const Vector& vec);

        guint         size () const;
        bool          is_valid () const;
        void          allocate (const Vector& vec);
        void          allocate (guint num_elements);
        void          copy (const Vector& vec);
        void          copy (const gsl_vector* vec);
        void          set_all (double value);
        void          reset ();
        void          zero ();
        void          zero (guint num_elements);
        double&       operator[] (guint i) const;
        gsl_vector*   get_gsl_vector () const;
        void          set_gsl_vector (gsl_vector* vec);
        gsl_vector*   disown_gsl_vector ();

        double        min () const;
        guint         min_index () const;
        double        max () const;
        guint         max_index () const;

        void          load (const String& filename);
        void          save (const String& filename) const;

        void          normalise ();
        double        magnitude () const;
        double        norm2 () const;
        double        mean () const;
        double        sum () const;
        double        dot (const Vector& vec) const;
        double        dot (const Matrix& M, guint row) const;

        void          add (const Vector& vec);
        void          add (double d, const Vector& vec); // x' = x + d*vec
        void          sub (const Vector& vec);
        void          div (const Vector& vec);

        void          multiply (double val);
        void          multiply (const Matrix& M, const Vector& vec);
        void          multiply (const Matrix& M, const gsl_vector* vec);
        void          multiply_trans (const Matrix& M, const Vector& vec);
        void          multiply_trans (const Matrix& M, const gsl_vector* vec);

        void          print () const;

        friend std::ostream& operator<< (std::ostream& stream, const Vector& vec);
    };

    void  normalise (float *v);
    void  normalise (float &v1, float &v2, float &v3);
    float magnitude (const float *v);
    float magnitude (float v1, float v2, float v3);
    void  cross_product (float *a, const float *b, const float *c);
    float dot_product (const float *a, const float *b);















    inline void  normalise (float *v)    { float a; a = sqrt(v[0]*v[0]+v[1]*v[1]+v[2]*v[2]); v[0] /= a; v[1] /= a; v[2] /= a; }
    inline void  normalise (float &v1, float &v2, float &v3) { float a; a = sqrt(v1*v1+v2*v2+v3*v3); v1 /= a; v2 /= a; v3 /= a; }
    inline float magnitude (const float *v)                   { return (sqrt(v[0]*v[0]+v[1]*v[1]+v[2]*v[2])); }
    inline float magnitude (float v1, float v2, float v3)     { return (sqrt(v1*v1+v2*v2+v3*v3)); }
    inline void  cross_product (float *a, const float *b, const float *c) { a[0] = b[1]*c[2]-b[2]*c[1]; a[1] = b[2]*c[0]-b[0]*c[2]; a[2] = b[0]*c[1]-b[1]*c[0]; }
    inline float dot_product (const float *a, const float *b)  { return (a[0]*b[0]+a[1]*b[1]+a[2]*b[2]); }

    /** @} */



    inline Vector::Vector ()                               { V = NULL; }
    inline Vector::Vector (guint num_elements)             { V = gsl_vector_alloc (num_elements); }
    inline Vector::Vector (const Vector& vec)              { V = NULL; copy (vec); }
    inline Vector::~Vector ()                              { if (V) gsl_vector_free (V); }
    inline Vector& Vector::operator= (const Vector& vec)   { copy (vec); return (*this); }
    inline void Vector::allocate (const Vector& vec)       { allocate (vec.size()); }
    inline void Vector::copy (const Vector &vec)           { copy (vec.V); }
    inline void Vector::copy (const gsl_vector* vec)       { allocate (vec->size); gsl_vector_memcpy (V, vec); }
    inline void Vector::reset ()                           { if (V) gsl_vector_free (V); V = NULL; }
    inline void Vector::set_all (double value)             { if (V) gsl_vector_set_all (V, value); }

    inline guint Vector::size () const                            { if (V) return (V->size); else return (0); }
    inline bool Vector::is_valid () const                         { return (V); }
    inline void Vector::zero ()                                   { if (V) gsl_vector_set_zero (V); }
    inline void Vector::zero (guint num_elements)                 { allocate (num_elements); zero(); }
    inline double &Vector::operator[] (guint i) const             { return (V->data[i * V->stride]); }
    inline gsl_vector *Vector::get_gsl_vector () const            { return (V); }
    inline void Vector::set_gsl_vector (gsl_vector* vec)          { if (V) gsl_vector_free (V); V = vec; }
    inline gsl_vector* Vector::disown_gsl_vector ()               { gsl_vector* vec = V; V = NULL; return(vec); }
    inline double Vector::min () const                            { if (V) return (gsl_vector_min (V)); else return (NAN); }
    inline guint   Vector::min_index () const                     { if (V) return (gsl_vector_min_index (V)); else return (0); }
    inline double Vector::max () const                            { if (V) return (gsl_vector_max (V)); else return (NAN); }
    inline guint   Vector::max_index () const                     { if (V) return (gsl_vector_max_index (V)); else return (0); }
    inline double Vector::magnitude () const                      { return (sqrt(norm2())); }

    inline void Vector::normalise ()                              { multiply (1.0/magnitude()); }

    inline double Vector::norm2 () const
    {
      double val = 0.0;
      for (guint i = 0; i < size(); i++) val += (*this)[i] * (*this)[i];
      return (val);
    }

    inline double Vector::mean () const                           { return (sum()/size()); }
    inline double Vector::sum () const
    {
      double val = 0.0;
      for (guint i = 0; i < size(); i++) val += (*this)[i];
      return (val);
    }

    inline double Vector::dot (const Vector& vec) const
    {
      double val = 0.0;
      for (guint i = 0; i < size(); i++) val += (*this)[i] * vec[i];
      return (val);
    }
    inline double Vector::dot (const Matrix& M, guint row) const
    {
      double val = 0.0;
      for (guint i = 0; i < size(); i++) val += (*this)[i] * M(row, i);
      return (val);
    }

    inline void Vector::add (const Vector& vec)                   { if (gsl_vector_add (V, vec.V)) throw Exception("vector"); }
    inline void Vector::add (double d, const Vector& vec)
    {
      if (!V) throw Exception ("vector");
      if (size() != vec.size()) throw Exception ("vector");
      for (guint n = 0; n < size(); n++)
        (*this)[n] += d*vec[n];
    }

    inline void Vector::sub (const Vector& vec)                   { if (gsl_vector_sub (V, vec.V)) throw Exception ("vector"); }
    inline void Vector::div (const Vector& vec)                   { if (gsl_vector_div (V, vec.V)) throw Exception ("vector"); }
    inline void Vector::multiply (double val)                     { if (gsl_vector_scale (V, val)) throw Exception ("vector"); }
    inline void Vector::multiply (const Matrix& M, const Vector& vec) { multiply (M, vec.V); }
    inline void Vector::multiply (const Matrix& M, const gsl_vector* vec)
    {
      allocate (M.rows());
      if (gsl_blas_dgemv (CblasNoTrans, 1.0, M.get_gsl_matrix(), vec, 0.0, V)) throw Exception ("vector");
    }

    inline void Vector::multiply_trans (const Matrix& M, const Vector& vec) { multiply_trans (M, vec.V); }
    inline void Vector::multiply_trans (const Matrix& M, const gsl_vector* vec)
    {
      allocate (M.columns());
      if (gsl_blas_dgemv (CblasTrans, 1.0, M.get_gsl_matrix(), vec, 0.0, V)) throw Exception ("vector");
    }

    inline void Vector::allocate (guint num_elements)
    {
      if (V) {
        if (size() == num_elements) return;
        gsl_vector_free (V);
      }
      V = gsl_vector_alloc (num_elements);
    }

}
}


#endif

