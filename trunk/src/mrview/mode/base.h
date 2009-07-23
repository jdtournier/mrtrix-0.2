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

#ifndef __mrview_mode_base_h__
#define __mrview_mode_base_h__

#include <gdk/gdkevents.h>

#include "ptr.h"
#include "image/object.h"


#define NUM_VIEW_MODES 1

namespace MR {

  namespace Math { class Vector; }

  namespace Viewer {

    class Window;
    class Pane;
    
    namespace Mode {

      class Base 
      {
        public:
          Base (Pane& segment, guint type_identifier) : pane (segment), type_id (type_identifier), OS (1), OS_x (0), OS_y (0) { }   
          virtual ~Base () { }

          virtual void      draw () = 0;
          virtual void      realize ();
          virtual void      configure () = 0;
          virtual bool      on_key_press (GdkEventKey* event);
          virtual bool      on_button_press (GdkEventButton* event);
          virtual bool      on_button_release (GdkEventButton* event);
          virtual bool      on_motion (GdkEventMotion* event);
          virtual bool      on_scroll (GdkEventScroll* event);
          guint             type () const { return (type_id); }
          Pane&             pane;

          void              set_oversampling (guint factor = 1, guint x = 0, guint y = 0) { OS = factor; OS_x = x; OS_y = y; }
          guint             get_oversampling () const { return (OS); }

        protected:
          guint             type_id, OS, OS_x, OS_y;
      };


      RefPtr<Base> create (Pane& parent, guint index);

    }
  }
}

#endif

