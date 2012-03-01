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

#include <gsl/gsl_eigen.h>
#include <gsl/gsl_sort_vector.h>

#include "math/linalg.h"

namespace MR {
  namespace Math {

    namespace {

      gsl_vector*                 eigen_values;
      gsl_eigen_symm_workspace*   eig_work;
      gsl_eigen_symmv_workspace*  eigv_work;

    }




    PseudoInverter::PseudoInverter ()
    {
      singular_values = SVD_work = NULL;
      SVD_V = SVD_U = SVD_Ut = SVD_S = SVD_D = NULL;
    }




    PseudoInverter::PseudoInverter (const Matrix& inv, const Matrix& src)
    {
      singular_values = SVD_work = NULL;
      SVD_V = SVD_U = SVD_Ut = SVD_S = SVD_D = NULL;
      init(inv, src);
    }



    void PseudoInverter::init (const Matrix& inv, const Matrix& src)
    {
      if (singular_values) { delete singular_values; singular_values = NULL; }
      if (SVD_work) { delete SVD_work; SVD_work = NULL; }
      if (SVD_V) { delete SVD_V; SVD_V = NULL; }
      if (SVD_U) { delete SVD_U; SVD_U = NULL; }
      if (SVD_Ut) { delete SVD_Ut; SVD_Ut = NULL; }
      if (SVD_S) { delete SVD_S; SVD_S = NULL; }
      if (SVD_D) { delete SVD_D; SVD_D = NULL; }

      SVD_V = SVD_U = SVD_Ut = SVD_S = SVD_D = NULL;
      if (src.rows() < src.columns()) 
        throw Exception ("Cannot invert MxN matrix when M < N");

      singular_values = gsl_vector_alloc(src.columns());
      SVD_work = gsl_vector_alloc(src.columns());

      SVD_U = new Matrix(src.rows(), src.columns());
      SVD_Ut = new Matrix(src.columns(), src.rows());
      SVD_V = new Matrix(src.columns(), src.columns());
      SVD_S = new Matrix(src.columns(), src.columns());
      SVD_D = new Matrix(src.columns(), src.rows());

      SVD_S->zero();
    }




    void PseudoInverter::invert (Matrix& inv, const Matrix& src, double threshold)
    {
      SVD_S->zero();
      SVD_U->copy(src);

      if (gsl_linalg_SV_decomp(SVD_U->get_gsl_matrix(), SVD_V->get_gsl_matrix(), singular_values, SVD_work)) 
        throw Exception ("error computing SVD for pseudo-inverse");

      for (guint x = 0; x < src.columns(); x++)
        (*SVD_S)(x,x) = gsl_vector_get(singular_values, x) > threshold ?
          1.0/gsl_vector_get(singular_values, x) : 0.0;

      //     VAR(gsl_vector_get(singular_values, 0)/ gsl_vector_get(singular_values,src.columns()-1));
      SVD_Ut->transpose(*SVD_U);

      SVD_D->multiply(*SVD_S, *SVD_Ut);
      inv.multiply(*SVD_V, *SVD_D);
    }





    PseudoInverter::~PseudoInverter()
    {
      gsl_vector_free(singular_values);
      gsl_vector_free(SVD_work);
      delete SVD_U;
      delete SVD_Ut;
      delete SVD_V;
      delete SVD_S;
      delete SVD_D;
    }





    void eig_init (Matrix& src, bool compute_eigenvectors)
    {
      if (src.rows() != src.columns()) 
        throw Exception ("can't calculate eigenvalues for non-square matrices");

      eigen_values = gsl_vector_alloc(src.rows());
      eig_work = NULL;
      eigv_work = NULL;
      if (compute_eigenvectors) eigv_work = gsl_eigen_symmv_alloc (src.rows());
      else eig_work = gsl_eigen_symm_alloc (src.rows());
    }


    void eig (Matrix& src, double* evals)
    {
      gsl_eigen_symm (src.get_gsl_matrix(), eigen_values, eig_work);
      gsl_sort_vector (eigen_values);
      for (guint i = 0; i < src.rows(); i++)
        evals[i] = gsl_vector_get(eigen_values, i);
    }

    void eig (Matrix& src, Vector& evals)
    {
      evals.allocate (src.rows());
      gsl_eigen_symm (src.get_gsl_matrix(), evals.get_gsl_vector(), eig_work);
      gsl_sort_vector (evals.get_gsl_vector());
    }


    void eig (Matrix& src, double* evals, Matrix& evec)
    {
      gsl_eigen_symmv (src.get_gsl_matrix(), eigen_values, evec.get_gsl_matrix(), eigv_work);
      gsl_eigen_symmv_sort (eigen_values, evec.get_gsl_matrix(), GSL_EIGEN_SORT_VAL_ASC);
      for (guint i = 0; i < src.rows(); i++)
        evals[i] = gsl_vector_get (eigen_values, i);
    }



    void eig_end ()
    {
      if (eig_work) gsl_eigen_symm_free (eig_work);
      if (eigv_work) gsl_eigen_symmv_free (eigv_work);
      gsl_vector_free (eigen_values);
    }


  }
}

