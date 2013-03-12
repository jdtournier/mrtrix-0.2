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

#include "app.h"
#include "image/position.h"
#include "math/linalg.h"
#include "dwi/tensor.h"

using namespace std; 
using namespace MR; 

SET_VERSION_DEFAULT;

DESCRIPTION = {
  "generate maps of tensor-derived parameters.",
  NULL
};

ARGUMENTS = {
  Argument ("tensor", "input tensor image", "the input diffusion tensor image.").type_image_in (),
  Argument::End
};


OPTIONS = { 
  Option ("adc", "mean ADC", "compute the mean apparent diffusion coefficient (ADC) of the diffusion tensor.")
    .append (Argument ("image", "image", "the output mean ADC image.").type_image_out ()),

  Option ("fa", "fractional anisotropy", "compute the fractional anisotropy of the diffusion tensor.")
    .append (Argument ("image", "image", "the output fractional anisotropy image.").type_image_out ()),

  Option ("num", "numbers", "specify the desired eigenvalue/eigenvector(s). Note that several eigenvalues can be specified as a number sequence. For example, '1,3' specifies the major (1) and minor (3) eigenvalues/eigenvectors (default = 1).")
    .append (Argument ("image", "image", "the mask image to use.").type_string ()),

  Option ("vector", "eigenvector", "compute the selected eigenvector(s) of the diffusion tensor.")
    .append (Argument ("image", "image", "the output eigenvector image.").type_image_out ()),

  Option ("value", "eigenvalue", "compute the selected eigenvalue(s) of the diffusion tensor.")
    .append (Argument ("image", "image", "the output eigenvalue image.").type_image_out ()),

  Option ("mask", "brain mask", "only perform computation within the specified binary brain mask image.")
    .append (Argument ("image", "image", "the mask image to use.").type_image_in ()),

  Option::End 
};


EXECUTE {
  Image::Position dt (*argument[0].get_image());
  Image::Header header (dt.image);

  if (header.ndim() != 4) 
    throw Exception ("base image should contain 4 dimensions");

  if (header.dim(3) != 6) 
    throw Exception ("expecting dimension 3 of image \"" + header.name + "\" to be 6");

  header.data_type = DataType::Float32;
  std::vector<OptBase> opt;

  std::vector<int> vals(1);
  vals[0] = 1;
  opt = get_options (2); // num
  if (opt.size()) {
    vals = parse_ints (opt[0][0].get_string());
    if (vals.empty()) throw Exception ("invalid eigenvalue/eigenvector number specifier");
    for (size_t i = 0; i < vals.size(); i++)
      if (vals[i] < 1 || vals[i] > 3) throw Exception ("eigenvalue/eigenvector number is out of bounds");
  }

  RefPtr<Image::Position> adc, fa, eval, evec, mask;

  opt = get_options (3); // vector
  if (opt.size()) {
    header.axes.dim[3] = 3*vals.size();
    evec = new Image::Position (*opt[0][0].get_image (header));
  }

  opt = get_options (4); // value
  if (opt.size()) {
    header.axes.dim[3] = vals.size();
    eval = new Image::Position (*opt[0][0].get_image (header));
  }

  header.axes.set_ndim (3);
  opt = get_options (0); // adc
  if (opt.size()) adc = new Image::Position (*opt[0][0].get_image (header));

  opt = get_options (1); // FA
  if (opt.size()) fa = new Image::Position (*opt[0][0].get_image (header));

  opt = get_options (5); // mask
  if (opt.size()) {
    mask = new Image::Position (*opt[0][0].get_image ());
    if (mask->dim(0) != dt.dim(0) || mask->dim(1) != dt.dim(1) || mask->dim(2) != dt.dim(2)) 
      throw Exception ("dimensions of mask image do not match that of tensor image - aborting");
  }


  if ( ! (adc || fa || eval || evec))
    throw Exception ("no output metric specified - aborting");


  for (size_t i = 0; i < vals.size(); i++)
    vals[i] = 3-vals[i];
 

  Math::Matrix V(3,3), M(3,3);
  double ev[3];
  float el[6];

  Math::eig_init (M, evec);

  ProgressBar::init (dt.dim(0)*dt.dim(1)*dt.dim(2), "computing tensor metrics...");

  for (dt.set(2,0); dt[2] < dt.dim(2); dt.inc(2)) {
    if (mask) mask->set(1,0); 
    if (fa) fa->set(1,0); 
    if (adc) adc->set(1,0); 
    if (eval) eval->set(1,0); 
    if (evec) evec->set(1,0);
    for (dt.set(1,0); dt[1] < dt.dim(1); dt.inc(1)) {
      if (mask) mask->set(0,0); 
      if (fa) fa->set(0,0); 
      if (adc) adc->set(0,0); 
      if (eval) eval->set(0,0); 
      if (evec) evec->set(0,0);
      for (dt.set(0,0); dt[0] < dt.dim(0); dt.inc(0)) {

        bool skip = false;
        if (mask) if (mask->value() < 0.5) skip = true;

        if (!skip) {

          for (dt.set(3,0); dt[3] < dt.dim(3); dt.inc(3)) 
            el[dt[3]] = dt.value();

          if (adc) adc->value (DWI::tensor2ADC (el));
          if (fa) fa->value (DWI::tensor2FA (el));

          if (eval || evec) {
            M(0,0) = el[0];
            M(1,1) = el[1];
            M(2,2) = el[2];
            M(0,1) = M(1,0) = el[3];
            M(0,2) = M(2,0) = el[4];
            M(1,2) = M(2,1) = el[5];

            if (evec) {
              Math::eig (M, ev, V);
              evec->set(3,0);
              for (size_t i = 0; i < vals.size(); i++) {
                evec->value (V(0,vals[i])); evec->inc(3);
                evec->value (V(1,vals[i])); evec->inc(3);
                evec->value (V(2,vals[i])); evec->inc(3);
              }
            }
            else Math::eig (M, ev);

            if (eval) {
              for (eval->set(3,0); (*eval)[3] < (int) vals.size(); eval->inc(3))
                eval->value (ev[vals[(*eval)[3]]]); 
            }
          }
        }

        ProgressBar::inc();

        if (mask) mask->inc(0);
        if (fa) fa->inc(0);
        if (adc) adc->inc(0);
        if (eval) eval->inc(0);
        if (evec) evec->inc(0);
      }
      if (mask) mask->inc(1);
      if (fa) fa->inc(1);
      if (adc) adc->inc(1);
      if (eval) eval->inc(1);
      if (evec) evec->inc(1);
    }
    if (mask) mask->inc(2);
    if (fa) fa->inc(2);
    if (adc) adc->inc(2);
    if (eval) eval->inc(2);
    if (evec) evec->inc(2);
  }

  ProgressBar::done();
  Math::eig_end();
}

