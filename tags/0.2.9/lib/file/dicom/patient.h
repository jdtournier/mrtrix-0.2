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

#ifndef __file_dicom_patient_h__
#define __file_dicom_patient_h__

#include "ptr.h"
#include "file/dicom/study.h"

namespace MR {
  namespace File {
    namespace Dicom {

      class Study;

      class Patient : public std::vector< RefPtr<Study> > {
        public:
          Patient (
              const String& patient_name, 
              const String& patient_ID = "", 
              const String& patient_DOB = "");
          String   name;
          String   ID;
          String   DOB;

          RefPtr<Study>     find (
              const String& study_name, 
              const String& study_ID = "", 
              const String& study_date = "", 
              const String& study_time = "");

      };


      std::ostream& operator<< (std::ostream& stream, const Patient& item);








      inline Patient::Patient (
          const String& patient_name, 
          const String& patient_ID, 
          const String& patient_DOB) :
        name (patient_name),
        ID (patient_ID),
        DOB (patient_DOB)
      {
      }


    }
  }
}


#endif



