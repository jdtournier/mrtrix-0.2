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

#ifndef __math_linalg_h__
#define __math_linalg_h__

#include <gsl/gsl_linalg.h>

#include "math/vector.h"
#include "math/matrix.h"

namespace MR {
  namespace Math {

    class PseudoInverter {
      protected:
        gsl_vector*   singular_values;
        gsl_vector*   SVD_work;
        Matrix*       SVD_V;
        Matrix*       SVD_U;
        Matrix*       SVD_Ut;
        Matrix*       SVD_S;
        Matrix*       SVD_D;

      public:
        PseudoInverter ();
        PseudoInverter (const Matrix& inv, const Matrix& src);
        ~PseudoInverter ();

        void      init (const Matrix& inv, const Matrix& src);
        void      invert (Matrix& inv, const Matrix& src, double threshold = 0.0);
    };



    inline void invert (Matrix& inv, const Matrix& src, double threshold = 0.0)
    {
      PseudoInverter I (inv, src);
      I.invert (inv, src, threshold);
    }


    void QR (Matrix& A, Vector& tau);
    void QR (Matrix& A, Matrix& Q, Matrix& R);
    void QR_solve (Matrix& A, Vector& b, Vector& x);
    void QR_solve (Matrix& A, Vector& tau, Vector& b, Vector& x);
    void QR_LS_solve (Matrix& A, Vector& b, Vector& x, Vector& residuals);

    void eig_init (Matrix& src, bool compute_eigenvectors);
    void eig (Matrix& src, Vector& evals);
    void eig (Matrix& src, double* evals);
    void eig (Matrix& src, double* evals, Matrix& evec);
    void eig_end ();









    inline void QR(Matrix& A, Vector& tau)
    {
      tau.allocate(A.columns());
      if (gsl_linalg_QR_decomp(A.get_gsl_matrix(), tau.get_gsl_vector())) throw Exception ("matrix");
    }




    inline void QR(Matrix& A, Matrix& Q, Matrix& R)
    {
      Vector tau;
      QR (A, tau);
      Q.allocate (A.rows(), A.rows());
      R.allocate (A);
      if (gsl_linalg_QR_unpack(A.get_gsl_matrix(), tau.get_gsl_vector(), Q.get_gsl_matrix(), R.get_gsl_matrix()))
        throw Exception ("matrix");
    }



    inline void QR_solve(Matrix& A, Vector& b, Vector& x)
    {
      Vector tau;
      QR (A, tau);
      x.allocate(A.rows());
      if (gsl_linalg_QR_solve (A.get_gsl_matrix(), tau.get_gsl_vector(), b.get_gsl_vector(), x.get_gsl_vector()))
        throw Exception ("matrix");
    }


    inline void QR_LS_solve(Matrix& A, Vector& b, Vector& x, Vector& residuals)
    {
      Vector tau;
      x.allocate (A.columns());
      residuals.allocate (A.rows());
      QR (A, tau);
      if (gsl_linalg_QR_lssolve (
            A.get_gsl_matrix(), tau.get_gsl_vector(), b.get_gsl_vector(), 
            x.get_gsl_vector(), residuals.get_gsl_vector())) 
        throw Exception ("matrix");
    }

  }
}

#endif

