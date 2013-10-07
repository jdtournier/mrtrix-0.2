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

#include <fstream>

#include "app.h"
#include "image/position.h"
#include "dwi/gradient.h"
#include "dwi/tensor.h"
#include "dwi/SH.h"

using namespace std; 
using namespace MR; 

SET_VERSION_DEFAULT;

DESCRIPTION = {
  "estimate the fibre response function for use in spherical deconvolution.",
  NULL
};

ARGUMENTS = {
  Argument ("dwi", "input DW image", "the input diffusion-weighted images.").type_image_in (),
  Argument ("mask", "single-fibre mask image", "the mask image of the voxels assumed to contain a single fibre population.").type_image_in (),
  Argument ("response", "response file", "the output text file where the even l, m=0 SH coefficients of the response function will be stored.").type_file (),
  Argument::End
};


OPTIONS = {
  Option ("grad", "supply gradient encoding", "specify the diffusion-weighted gradient scheme used in the acquisition. The program will normally attempt to use the encoding stored in image header.", false, true)
    .append (Argument ("encoding", "gradient encoding", "the gradient encoding, supplied as a 4xN text file with each line is in the format [ X Y Z b ], where [ X Y Z ] describe the direction of the applied gradient, and b gives the b-value in units (1000 s/mm^2).").type_file ()),

  Option ("lmax", "maximum harmonic order", "set the maximum harmonic order to be estimated. By default, the program will use the highest possible lmax given the number of diffusion-weighted images.")
    .append (Argument ("order", "order", "the maximum harmonic order to use.").type_integer (2, 30, 8)),

  Option ("normalise", "normalise to b=0", "normalise the DW signal to the b=0 image"),

  Option::End
};


EXECUTE {

  Image::Object &dwi_obj (*argument[0].get_image());
  Image::Header header (dwi_obj);

  if (header.ndim() != 4) 
    throw Exception ("dwi image should contain 4 dimensions");

  Math::Matrix grad;

  std::vector<OptBase> opt = get_options (0);
  if (opt.size()) grad.load (opt[0][0].get_string());
  else {
    if (!header.DW_scheme.is_valid()) 
      throw Exception ("no diffusion encoding found in image \"" + header.name + "\"");
    grad.copy (header.DW_scheme);
  }

  if (grad.rows() < 7 || grad.columns() != 4) 
    throw Exception ("unexpected diffusion encoding matrix dimensions");

  info ("found " + str(grad.rows()) + "x" + str(grad.columns()) + " diffusion-weighted encoding");

  if (header.dim(3) != (int) grad.rows()) 
    throw Exception ("number of studies in base image does not match that in encoding file");

  DWI::normalise_grad (grad);

  std::vector<int> bzeros, dwis;
  DWI::guess_DW_directions (dwis, bzeros, grad);

  {
    String msg ("found b=0 images in studies [ ");
    for (guint n = 0; n < bzeros.size(); n++) msg += str(bzeros[n]) + " ";
    msg += "]";
    info (msg);
  }

  info ("found " + str(dwis.size()) + " diffusion-weighted directions");

  opt = get_options (1);
  int lmax = opt.size() ? opt[0][0].get_int() : DWI::SH::LforN (dwis.size());
  if (lmax & 1) {
    lmax--;
    error ("odd number supplied for lmax - using lmax = " + str (lmax));
  }
  if (lmax > DWI::SH::LforN (dwis.size())) {
    info ("warning: not enough data to estimate spherical harmonic components up to order " + str(lmax));
    lmax = DWI::SH::LforN (dwis.size());
  }
  info ("calculating even spherical harmonic components up to order " + str(lmax));

  bool normalise = get_options(2).size();

  Math::Matrix bmat, binv;
  DWI::grad2bmatrix (bmat, grad);
  Math::invert (binv, bmat);

  Image::Object &mask_obj (*argument[1].get_image());
  if (mask_obj.dim(0) != header.dim(0) ||
      mask_obj.dim(1) != header.dim(1) ||
      mask_obj.dim(2) != header.dim(2)) 
    throw Exception ("mask & DWI image dimensions do not match");


  Image::Position dwi (dwi_obj);
  Image::Position mask (mask_obj);


  Math::Matrix dirs, T, E(3,3), V(3,3), D(3,3), rotated_grad (grad), SHT, iSHT;
  Math::Vector vec(3), rot(3), sig(dwis.size());
  Math::Vector response(lmax/2+1);
  DWI::SH::Coefs SH;


  DWI::gen_direction_matrix (dirs, grad, dwis);
  DWI::SH::init_transform (SHT, dirs, lmax);
  Math::PseudoInverter inverter (iSHT, SHT);
  Math::eig_init (D, true);

  int count = 0;
  for (mask.set(2,0); mask[2] < mask.dim(2); mask.inc(2)) 
    for (mask.set(1,0); mask[1] < mask.dim(1); mask.inc(1)) 
      for (mask.set(0,0); mask[0] < mask.dim(0); mask.inc(0)) 
	if (mask.value() > 0.5) count++;

  info ("found " + str(count) + " single-fibre voxels");

  response.zero();

  std::vector<float> val (dwi.dim(3));
  ProgressBar::init (count, "estimating response function...");

  for (dwi.set(2,0), mask.set(2,0); dwi[2] < dwi.dim(2); dwi.inc(2), mask.inc(2)) {
    for (dwi.set(1,0), mask.set(1,0); dwi[1] < dwi.dim(1); dwi.inc(1), mask.inc(1)) {
      for (dwi.set(0,0), mask.set(0,0); dwi[0] < dwi.dim(0); dwi.inc(0), mask.inc(0)) {

	if (mask.value() > 0.5) {

          for (dwi.set(3,0); dwi[3] < dwi.dim(3); dwi.inc(3)) {
	    val[dwi[3]] = dwi.value();
            val[dwi[3]] = val[dwi[3]] > 0.0 ? -log (val[dwi[3]]) : -1e-12;
          }

          float dt[6];
          for (int n = 0; n < 6; n++) {
            dt[n] = 0.0;
            for (int i = 0; i < dwi.dim(3); i++)
              dt[n] += (float) (binv(n,i) * val[i]);
          }

	  D(0,0) = dt[0];
	  D(1,1) = dt[1];
	  D(2,2) = dt[2];
	  D(0,1) = D(1,0) = dt[3];
	  D(0,2) = D(2,0) = dt[4];
	  D(1,2) = D(2,1) = dt[5];

          double evals[3];
          Math::eig (D, evals, V);
	  V.transpose();

	  for (guint i = 0; i < grad.rows(); i++) {
	    vec[0] = grad(i,0);
	    vec[1] = grad(i,1);
	    vec[2] = grad(i,2);
	    rot.multiply (V, vec);
	    rotated_grad(i,0) = rot[0];
	    rotated_grad(i,1) = rot[1];
	    rotated_grad(i,2) = rot[2];
	  }

	  DWI::gen_direction_matrix (dirs, rotated_grad, dwis);

          DWI::SH::init_transform (SHT, dirs, lmax);
	  inverter.invert (iSHT, SHT);

          float norm = 1.0;
	  if (normalise) {
	    norm = 0.0;
	    for (guint i = 0; i < bzeros.size(); i++) { dwi.set(3, bzeros[i]); norm += dwi.value(); }
	    norm /= bzeros.size();
	  }

	  for (guint i = 0; i < dwis.size(); i++) {
            dwi.set (3, dwis[i]);
	    sig[i] = dwi.value()/norm;
          }

	  SH.V.multiply (iSHT, sig);
/*
	  for (guint i = 0; i < SH.V.size(); i++) 
            std::cout << SH.V[i] << " ";
          std::cout << "\n";
*/
	  for (int l = 0; l <= lmax; l+=2) 
	    response[l/2] += SH(l,0);

          ProgressBar::inc();
	}
      }
    }
  }
  ProgressBar::done();
  Math::eig_end ();

  std::ofstream response_file (argument[2].get_string());

  for (guint i = 0; i < response.size(); i++) {
    response[i] /= count;
    response_file << response[i] << " ";
  }
  response_file << "\n";

  info ("estimated response SH coefficients: " + str (response));
}

