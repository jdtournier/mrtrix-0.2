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

#ifndef __file_dicom_series_h__
#define __file_dicom_series_h__

#include "ptr.h"
#include "file/dicom/image.h"

namespace MR {
  namespace File {
    namespace Dicom {

      class Study;
      class Image;

      class Series : public std::vector< RefPtr <Image> > {
        public:
          Series (
              Study* parent, 
              const String& series_name, 
              guint         series_number,
              const String& series_modality = "", 
              const String& series_date = "", 
              const String& series_time = "");

          Study*   study;
          String   name;
          guint    number;
          String   modality;
          String   date;
          String   time;

          void               read ();
          std::vector<int>   count () const;
          void               print_fields (bool dcm, bool csa) const;

          bool               operator< (const Series& s) const;

          friend std::ostream& operator<< (std::ostream& stream, const Series& item);
      };













      inline Series::Series (
          Study*  parent, 
          const String& series_name, 
          guint         series_number,
          const String& series_modality, 
          const String& series_date, 
          const String& series_time) :
        study (parent),
        name (series_name),
        modality (series_modality),
        date (series_date),
        time (series_time)
      { 
        number = series_number; 
      }




      inline void Series::read ()
      {
        ProgressBar::init (size(), "reading DICOM series \"" + name + "\"...");
        for (guint i = 0; i < size(); i++) {
          (*this)[i]->read(); 
          ProgressBar::inc();
        }
        ProgressBar::done();
      }




      inline void Series::print_fields (bool dcm, bool csa) const 
      {
        for (guint i = 0; i < size(); i++) 
          (*this)[i]->print_fields (dcm, csa);
      }


    }
  }
}


#endif



