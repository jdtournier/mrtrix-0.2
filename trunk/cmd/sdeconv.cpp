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


    31-10-2008 J-Donald Tournier <d.tournier@brain.org.au>
    * use MR::Ptr instead of std::auto_ptr
   
    25-09-2009 J-Donald Tournier <d.tournier@brain.org.au>
    * fix documentation of SH coefficient storage convention

*/

#include "app.h"
#include "ptr.h"
#include "image/position.h"
#include "math/linalg.h"
#include "dwi/gradient.h"
#include "dwi/SH.h"

using namespace std; 
using namespace MR; 

SET_VERSION_DEFAULT;

DESCRIPTION = {
  "perform linear spherical deconvolution.",

  "Note that this program makes use of implied symmetries in the diffusion profile. First, the fact the signal attenuation profile is real implies that it has conjugate symmetry, i.e. Y(l,-m) = Y(l,m)* (where * denotes the complex conjugate). Second, the diffusion profile should be antipodally symmetric (i.e. S(x) = S(-x)), implying that all odd l components should be zero. Therefore, this program only computes the even elements.",

  "Note that the spherical harmonics equations used here differ slightly from those conventionally used, in that the (-1)^m factor has been omitted. This should be taken into account in all subsequent calculations.",

  "Each volume in the output image corresponds to a different spherical harmonic component, according to the following convention:\n"
    "[0]     Y(0,0)\n"
    "\n"
    "[1] Im {Y(2,2)}\n"
    "[2] Im {Y(2,1)}\n"
    "[3]     Y(2,0)\n"
    "[4] Re {Y(2,1)}\n"
    "[5] Re {Y(2,2)}\n"
    "\n"
    "[6] Im {Y(4,4)}\n"
    "[7] Im {Y(4,3)}\n"
    "etc...\n",

  NULL
};

ARGUMENTS = {
  Argument ("dwi", "input DW image", "the input diffusion-weighted image.").type_image_in (),
  Argument ("response", "response function", "the diffusion-weighted signal response function for a single fibre population.").type_file (),
  Argument ("SH", "output SH image", "the output spherical harmonics coefficients image.").type_image_out (),
  Argument::End
};


OPTIONS = { 
  Option ("grad", "supply gradient encoding", "specify the diffusion-weighted gradient scheme used in the acquisition. The program will normally attempt to use the encoding stored in image header.", false, true)
    .append (Argument ("encoding", "gradient encoding", "the gradient encoding, supplied as a 4xN text file with each line is in the format [ X Y Z b ], where [ X Y Z ] describe the direction of the applied gradient, and b gives the b-value in units (1000 s/mm^2).").type_file ()),

  Option ("lmax", "maximum harmonic order", "set the maximum harmonic order for the output series. By default, the program will use the highest possible lmax given the number of diffusion-weighted images.")
    .append (Argument ("order", "order", "the maximum harmonic order to use.").type_integer (2, 30, 8)),

  Option ("mask", "brain mask", "only perform computation within the specified binary brain mask image.")
    .append (Argument ("image", "image", "the mask image to use.").type_image_in ()),

  Option ("filter", "angular frequency filter", "apply the linear frequency filtering parameters specified to the spherical deconvolution.")
    .append (Argument ("spec", "specification", "a text file containing the filtering coefficients for each even harmonic order.").type_file ()),

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
    String msg ("found b-zeros in studies [ ");
    for (guint n = 0; n < bzeros.size(); n++) msg += str(bzeros[n]) + " ";
    msg += "]";
    info (msg);
  }

  info ("found " + str(dwis.size()) + " diffusion-weighted directions");

  opt = get_options (1);
  int lmax = opt.size() ? opt[0][0].get_int() : DWI::SH::LforN (dwis.size());
  if (lmax > DWI::SH::LforN (dwis.size())) {
    info ("warning: not enough data to estimate spherical harmonic components up to order " + str(lmax));
    lmax = DWI::SH::LforN (dwis.size());
  }
  info ("calculating even spherical harmonic components up to order " + str(lmax));


  info (String ("setting response function from file \"") + argument[1].get_string() + "\"");
  Math::Vector response;
  response.load (argument[1].get_string());
  info ("setting response function using even SH coefficients: " + str(response));


  opt = get_options (3);
  Math::Vector filter;
  if (opt.size()) {
    filter.load (opt[0][0].get_string());
    if (filter.size() < (guint) ((lmax/2)+1)) 
      throw Exception ("not enough filter coefficients supplied for lmax = " + str(lmax));
  }
  else {
    filter.allocate (lmax/2+1);
    filter.set_all (1.0);
  }
  info ("using filtering coefficients: " + str(filter));




  {
    Math::Vector tmp;
    DWI::SH::SH2RH (tmp, response);
    for (int i = 0; i <= lmax/2; i++)
      if (response[i]) response[i] = filter[i]/tmp[i];
  }
  info ("setting inverse response function using even RH coefficients: " + str(response));



  Math::Matrix dirs;
  DWI::gen_direction_matrix (dirs, grad, dwis);
  DWI::SH::Transform SHT (dirs, lmax);
  SHT.set_filter (response);


  header.axes.dim[3] = SHT.n_SH();
  header.data_type = DataType::Float32;
  header.axes.axis[0] = 1; header.axes.forward[0] = true;
  header.axes.axis[1] = 2; header.axes.forward[1] = true;
  header.axes.axis[2] = 3; header.axes.forward[2] = true;
  header.axes.axis[3] = 0; header.axes.forward[3] = true;

  Image::Position dwi (dwi_obj);
  Image::Position sh (*argument[2].get_image (header));

  opt = get_options (2);
  Ptr<Image::Position> mask;
  if (opt.size()) 
    mask = new Image::Position (*opt[0][0].get_image());

  DWI::SH::Coefs res (lmax);
  Math::Vector sigs (dwis.size());


  bool normalise = get_options(4).size();

  ProgressBar::init (dwi.dim(0)*dwi.dim(1)*dwi.dim(2), "performing spherical deconvolution...");

  for (dwi.set(2,0), sh.set(2,0); dwi[2] < dwi.dim(2); dwi.inc(2), sh.inc(2)) {
    if (mask.get()) mask->set(2, dwi[2]);

    for (dwi.set(1,0), sh.set(1,0); dwi[1] < dwi.dim(1); dwi.inc(1), sh.inc(1)) {
      if (mask.get()) mask->set(1, dwi[1]);

      for (dwi.set(0,0), sh.set(0,0); dwi[0] < dwi.dim(0); dwi.inc(0), sh.inc(0)) {
        bool proceed = true;
        if (mask.get()) {
          mask->set(0, dwi[0]);
          if (mask->value() < 0.5) proceed = false;
        }

        if (proceed) {
          double norm = 0.0;
          if (normalise) {
            for (guint n = 0; n < bzeros.size(); n++) {
              dwi.set (3, bzeros[n]);
              norm += dwi.value ();
            }
            norm /= bzeros.size();
          }

          for (guint n = 0; n < dwis.size(); n++) {
            dwi.set (3, dwis[n]);
            sigs[n] = dwi.value(); 
            if (sigs[n] < 0.0) sigs[n] = 0.0;
            if (normalise) sigs[n] /= norm;
          }

          SHT.A2SH (res, sigs);

          for (sh.set(3,0); sh[3] < sh.dim(3); sh.inc(3))
            sh.value (res.V[sh[3]]);
        }

        ProgressBar::inc();
      }
    }
  }

  ProgressBar::done();
}

