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
  "Concatenate two track files",
  NULL
};

ARGUMENTS = {
  Argument ("output", "output tracks file", "the output file containing the tracks.")     .type_file(),
  Argument ("track1", "first input tracks file",  "the first input tracks file.").type_file(),
  Argument ("track2", "second tracks file", "the second input tracks file.", true, true)     .type_file(),
  Argument::End
};


OPTIONS = {

  Option::End
};

EXECUTE
{
  	int num_tracks = argument.size()-1;
	Reader reader;
	Properties properties;
	Writer writer;
	Properties outputproperties;

	reader.open (argument[1].get_string(), properties);

	//const float progress_multiplier = properties["count"].empty() ? 0.0 : 100.0 / to<float> (properties["count"]);
	
	outputproperties = Properties(properties);

	outputproperties["count"] = properties["count"];
	outputproperties["total_count"] = properties["total_count"];

	//properties.erase ("count");
	//properties.erase ("total_count");
	outputproperties["source"] = argument[0].get_string();
	
	writer.create (argument[0].get_string(), outputproperties);

	reader.close();
	for(int i = 1; i <= num_tracks; i++)
	{
		reader.open (argument[i].get_string(), properties);
		
		std::vector<Point> tck;
		while (reader.next (tck))
		{
			++writer.total_count;
			writer.append (tck);
		}
		reader.close();
	}	
	writer.close();
}
