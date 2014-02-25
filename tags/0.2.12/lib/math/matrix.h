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

#ifndef __math_matrix_h__
#define __math_matrix_h__

#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>

#include "mrtrix.h"

namespace MR {
  namespace Math {

    class Matrix
    {
      protected:
        gsl_matrix*   M;

      public:
        Matrix ();
        Matrix (const Matrix& mat);
        Matrix (guint n_rows, guint n_columns);
        ~Matrix ();

        Matrix&       operator= (const Matrix& A);
        guint         rows () const;
        guint         columns () const;

        bool          is_valid () const;
        void          allocate (const Matrix& mat);
        void          allocate (guint n_rows, guint n_columns);
        void          copy (const Matrix &mat);
        void          copy (guint n_rows, guint n_columns, const double* data);
        void          reset ();
        void          zero ();
        void          identity ();
        double&       operator() (guint i, guint j) const;

        double        min () const;
        double        max () const;

        void           load (const String& filename);
        void           save (const String& filename) const;

        gsl_matrix*   get_gsl_matrix () const;

        void           add (const Matrix& A);
        void           add (double v);
        void           subtract (const Matrix& A);

        void           multiply (const Matrix& A, const Matrix& B);
        void           multiply (double v);

        void           multiply_elements (const Matrix& A);
        void           divide_elements (const Matrix& A);

        void           transpose (const Matrix& A);
        void           transpose ();

        void          print () const;


        friend std::ostream& operator<< (std::ostream& stream, const Matrix& M);
    };



    bool operator!= (const Matrix& left, const Matrix& right);
    bool operator== (const Matrix& left, const Matrix& right);


















    inline Matrix::Matrix ()                               { M = NULL; }
    inline Matrix::Matrix (const Matrix& mat)              { M = NULL; copy (mat); }
    inline Matrix::Matrix (guint n_rows, guint n_columns)  { M = NULL; allocate (n_rows, n_columns); }
    inline Matrix::~Matrix ()                              { if (M) gsl_matrix_free (M); }
    inline Matrix& Matrix::operator= (const Matrix& A)     { copy (A); return (*this); }
    inline void Matrix::copy (const Matrix &A)             { allocate(A); if (M) gsl_matrix_memcpy (M, A.M); }
    inline void Matrix::copy (guint n_rows, guint n_columns, const double* data) 
    {
      allocate (n_rows, n_columns); 
      for (guint row = 0; row < n_rows; row++) 
        for (guint col = 0; col < n_columns; col++)
          (*this)(row, col) = data [col + n_columns*row];
    }
    inline void Matrix::reset ()                           { if (M) gsl_matrix_free (M); M = NULL; }
    inline void Matrix::allocate (const Matrix& mat)       { allocate (mat.rows(), mat.columns()); }

    inline guint Matrix::rows () const                           { if (!M) return(0); return (M->size1); }
    inline guint Matrix::columns () const                        { if (!M) return(0); return (M->size2); }
    inline bool Matrix::is_valid () const                       { return (M); }
    inline void Matrix::zero ()                                 { gsl_matrix_set_zero (M); }
    inline void Matrix::identity ()                             { gsl_matrix_set_identity (M); }
    inline double &Matrix::operator() (guint i, guint j) const    { return (M->data[i * M->tda + j]); }
    inline gsl_matrix *Matrix::get_gsl_matrix () const          { return (M); }
    inline double Matrix::min () const                          { return (gsl_matrix_min(M)); }
    inline double Matrix::max () const                          { return (gsl_matrix_max(M)); }
    inline void Matrix::transpose (const Matrix& A)
    {
      allocate(A.columns(), A.rows()); 
      if (gsl_matrix_transpose_memcpy (M, A.M)) throw Exception ("matrix"); 
    }
    inline void Matrix::transpose ()                             { if (gsl_matrix_transpose(M)) throw Exception ("matrix"); }

    inline void Matrix::add (const Matrix& A)                    { if (gsl_matrix_add (M, A.M)) throw Exception ("matrix"); }
    inline void Matrix::add (double v)                           { if (gsl_matrix_add_constant (M, v)) throw Exception ("matrix"); }
    inline void Matrix::subtract (const Matrix& A   )            { if (gsl_matrix_sub (M, A.M)) throw Exception ("matrix"); }
    inline void Matrix::multiply (double v)                      { if (gsl_matrix_scale (M, v)) throw Exception ("matrix"); }
    inline void Matrix::multiply_elements (const Matrix& A)      { if (gsl_matrix_mul_elements (M, A.M)) throw Exception ("matrix"); }
    inline void Matrix::divide_elements (const Matrix& A)        { if (gsl_matrix_div_elements (M, A.M)) throw Exception ("matrix"); }

    inline void Matrix::allocate (guint n_rows, guint n_columns)
    {
      if (M) {
        if (rows() == n_rows && columns() == n_columns) return;
        gsl_matrix_free (M);
      }
      M = ( n_rows && n_columns ? gsl_matrix_alloc (n_rows, n_columns) : NULL );
    }

    inline void Matrix::multiply (const Matrix& A, const Matrix& B)
    {
      allocate (A.rows(), B.columns());
      if (gsl_blas_dgemm (CblasNoTrans, CblasNoTrans, 1.0, A.M, B.M, 0.0, M)) throw Exception ("matrix");
    }

    inline bool operator!= (const Matrix& left, const Matrix& right)
    {
      if (left.rows() != right.rows() || left.columns() != right.columns()) return (true);
      for (guint j = 0; j < left.columns(); j++)
        for (guint i = 0; i < left.rows(); i++)
          if (left(i,j) != right(i,j))
            return (true);
      return (false);
    }


    inline bool operator== (const Matrix& left, const Matrix& right) { return (! (left != right)); }




  }
}

#endif

