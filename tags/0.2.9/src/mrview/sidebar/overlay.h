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

#ifndef __mrview_sidebar_overlay_h__
#define __mrview_sidebar_overlay_h__

#include <gtkmm/checkbutton.h>
#include <gtkmm/scrolledwindow.h>
#include <gtkmm/treeview.h>
#include <gtkmm/treerowreference.h>
#include <gtkmm/frame.h>
#include <gtkmm/scale.h>
#include <gtkmm/liststore.h>
#include <gtkmm/spinbutton.h>
#include <gtkmm/label.h>
#include <gtkmm/table.h>

#include "mrview/sidebar/base.h"
#include "mrview/slice.h"

namespace MR {
  namespace Viewer {
    namespace SideBar {

      class OverlayList;

      class Overlay : public Base
      {
        public:
          Overlay ();
          ~Overlay ();

          void draw();

          void set_scaling (float min, float max);

        protected:
          class OverlayList : public Gtk::TreeView {
            public:
              OverlayList (Overlay& sidebar);
              virtual ~OverlayList();
              void draw (int transparency);

              void set_scaling (float min, float max);

            protected:
              Overlay& parent;
              Gtk::TreeModel::Row row;
              bool min_max_changed;
              int colourmap;
              float min_from_parent, max_from_parent;

              class OL {
                public:
                  OL (RefPtr<MR::Image::Object> obj) : image (new Image (obj)) { }
                  RefPtr<Image> image;
                  Slice::Renderer render;
              };

              class Columns : public Gtk::TreeModel::ColumnRecord {
                public:
                  Columns() { add (show); add (name); add (colourmap), add (min), add (max), add (overlay); }

                  Gtk::TreeModelColumn<bool> show;
                  Gtk::TreeModelColumn<String> name;
                  Gtk::TreeModelColumn<int> colourmap;
                  Gtk::TreeModelColumn<float> min;
                  Gtk::TreeModelColumn<float> max;
                  Gtk::TreeModelColumn<RefPtr<OL> > overlay;
              };

              Columns columns;
              Glib::RefPtr<Gtk::ListStore> model;
              Gtk::Menu popup_menu, colourmap_menu;

              bool on_button_press_event(GdkEventButton *event);
              void on_open ();
              void on_close ();
              void on_colourmap (int mode);
              void colourmap_callback (const Gtk::TreeModel::iterator& iter);
              void on_clear ();
              void on_selection ();
              void selected_row_callback (const Gtk::TreeModel::iterator& iter);
              void set_scaling_callback (const Gtk::TreeModel::iterator& iter);
              void on_tick (const String& path);
              void load (RefPtr<MR::Image::Object> image);
          };

          Gtk::CheckButton     show_overlays;
          Gtk::Frame           overlay_frame, transparency_frame, scaling_frame;
          Gtk::Label           min_label, max_label;
          Gtk::HScale          transparency;
          Gtk::Table           scaling_table;
          Gtk::SpinButton      min_value, max_value;
          Gtk::ScrolledWindow  overlay_scrolled_window;
          OverlayList          overlay_list;

          void on_change ();
          void on_scaling ();
      };

    }
  }
}

#endif



