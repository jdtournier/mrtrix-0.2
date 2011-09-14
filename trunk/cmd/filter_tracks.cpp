/*
    Copyright 2008 Brain Research Institute, Melbourne, Australia

    Written by Robert E. Smith and J-Donald Tournier, 14/09/11.

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
#include "image/interp.h"
#include "math/vector.h"
#include "point.h"
#include "dwi/tractography/file.h"
#include "dwi/tractography/roi.h"
#include "dwi/tractography/tracker/base.h"


using namespace std; 
using namespace MR; 
using namespace MR::DWI; 
using namespace MR::DWI::Tractography; 

SET_VERSION_DEFAULT;

DESCRIPTION = {
  "Use regions-of-interest to select a sub-set of tracks from a given track file.\n ",
  "Each region-of-interest should be either the path to a binary mask image, or a comma-separated list of 4 floating-point values, specifying the [x,y,z] coordinates of the centre and radius of a spherical ROI.",
  NULL
};

ARGUMENTS = {
  Argument ("input",  "input tracks file",  "the input file containing the tracks to be filtered.").type_file(),
  Argument ("output", "output tracks file", "the output file containing the tracks selected.")     .type_file(),
  Argument::End
};


OPTIONS = {

  Option ("include", "inclusion ROI", "specify an inclusion region of interest. Tracks that enter ALL such inclusion ROI will be included.", false, true)
    .append (Argument ("spec", "ROI specification", "specifies the parameters necessary to define the ROI.").type_string()),

  Option ("exclude", "exclusion ROI", "specify an exclusion region of interest. Tracks that enter ANY such exclusion ROI will NOT be included.", false, true)
    .append (Argument ("spec", "ROI specification", "specifies the parameters necessary to define the ROI.").type_string()),

  Option ("nomaskinterp", "no interpolation of mask regions", "do NOT perform tri-linear interpolation of mask images."),

  Option::End
};



class ROI_filter {

  public:
    ROI_filter (Properties& properties)
    {

      bool no_mask_interp;
      if (properties["no_mask_interp"].empty()) {
        no_mask_interp = false;
        properties["no_mask_interp"] = "0";
      } else {
        no_mask_interp = to<bool> (properties["no_mask_interp"]);
      }

      for (std::vector< RefPtr<ROI> >::iterator i = properties.roi.begin(); i != properties.roi.end(); ++i) {
        ROI& roi (**i);
        if (!roi.mask.empty() && !roi.mask_object) {
          roi.mask_object = new MR::Image::Object;
          roi.mask_object->open (roi.mask);
        }
        switch (roi.type) {

          case ROI::Seed:
            throw Exception ("filter_tracks should not receive any seed regions");

          case ROI::Include:
            if (roi.mask.empty())
              spheres.include.push_back (Tracker::Base::Sphere (roi.position, roi.radius));
            else
              masks  .include.push_back (Tracker::Base::Mask   (*roi.mask_object, no_mask_interp));
            break;

          case ROI::Exclude:
            if (roi.mask.empty())
              spheres.exclude.push_back (Tracker::Base::Sphere (roi.position, roi.radius));
            else
              masks  .exclude.push_back (Tracker::Base::Mask   (*roi.mask_object, no_mask_interp));
            break;

          case ROI::Mask:
            throw Exception ("filter_tracks should not receive any mask regions");

          default:
            assert (0);

        }
      }

    }


    bool accept_track (const std::vector<Point>& tck)
    {

      for (std::vector<Tracker::Base::Sphere>::iterator i = spheres.include.begin(); i != spheres.include.end(); ++i)
        i->included = false;
      for (std::vector<Tracker::Base::Mask  >::iterator i = masks  .include.begin(); i != masks  .include.end(); ++i)
        i->included = false;

      for (std::vector<Point>::const_iterator pos = tck.begin(); pos != tck.end(); ++pos) {

        for (std::vector<Tracker::Base::Sphere>::iterator i = spheres.exclude.begin(); i != spheres.exclude.end(); ++i)
          if (i->contains (*pos)) return (false);

        for (std::vector<Tracker::Base::Mask  >::iterator i = masks  .exclude.begin(); i != masks  .exclude.end(); ++i)
          if (i->contains (*pos)) return (false);

        for (std::vector<Tracker::Base::Sphere>::iterator i = spheres.include.begin(); i != spheres.include.end(); ++i)
          if (!i->included && i->contains (*pos)) i->included = true;

        for (std::vector<Tracker::Base::Mask  >::iterator i = masks  .include.begin(); i != masks  .include.end(); ++i)
          if (!i->included && i->contains (*pos)) i->included = true;

      }

      for (std::vector<Tracker::Base::Sphere>::iterator i = spheres.include.begin(); i != spheres.include.end(); ++i)
        if (!i->included) return false;

      for (std::vector<Tracker::Base::Mask  >::iterator i = masks  .include.begin(); i != masks  .include.end(); ++i)
        if (!i->included) return false;

      return true;

    }


  private:
    Tracker::Base::ROISphere spheres;
    Tracker::Base::ROIMask   masks;

};




EXECUTE
{

  Reader reader;
  Properties properties;
  Writer writer;

  reader.open (argument[0].get_string(), properties);
  const float progress_multiplier = properties["count"].empty() ? 0.0 : 100.0 / to<float> (properties["count"]);

  properties.roi.clear(); // remove those used to generate the input track file
  properties.erase ("count");
  properties.erase ("total_count");
  properties["source"] = argument[0].get_string();

  std::vector<OptBase> opt = get_options (0); // include
  for (std::vector<OptBase>::iterator i = opt.begin(); i != opt.end(); ++i)
    properties.roi.push_back (RefPtr<ROI> (new ROI (ROI::Include, (*i)[0].get_string())));

  opt = get_options (1); // exclude
  for (std::vector<OptBase>::iterator i = opt.begin(); i != opt.end(); ++i)
    properties.roi.push_back (RefPtr<ROI> (new ROI (ROI::Exclude, (*i)[0].get_string())));

  opt = get_options (2); // nomaskinterp
  properties["no_mask_interp"] = opt.size() ? "1" : "0"; // need to override any existing property entry

  ROI_filter filter (properties);
  writer.create (argument[1].get_string(), properties);

  std::vector<Point> tck;
  while (reader.next (tck)) {
    ++writer.total_count;
    if (filter.accept_track (tck))
      writer.append (tck);
    fprintf (stderr, "\r%8u read, %8u selected    [%3d%%]",
      writer.total_count, writer.count, int(progress_multiplier * writer.total_count));
  }
  reader.close();
  writer.close();
  fprintf (stderr, "\r%8u read, %8u selected    [100%%]\n",
    writer.total_count, writer.count);

}
