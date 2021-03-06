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

#ifndef __mrview_sidebar_roi_analysis_roi_list_h__
#define __mrview_sidebar_roi_analysis_roi_list_h__

#include <gtkmm/treeview.h>
#include <gtkmm/liststore.h>
#include <gtkmm/menu.h>
//#include <gtkmm/tooltip.h>

#include "mrview/slice.h"
#include <deque>

namespace MR {
  namespace Viewer {
    namespace SideBar {

      class ROIAnalysis;

      class DP_ROIList : public Gtk::TreeView
      {
        public:
          DP_ROIList (const ROIAnalysis& sidebar);
          virtual ~DP_ROIList();

          void draw (int transparency);
	  bool on_button_press (GdkEventButton* event, float brush, bool brush3d, bool isobrush);
  	  bool on_motion (GdkEventMotion* event, float brush, bool brush3d, bool isobrush) 
	  { 
	    if (editing) 
	      { 
	      process (event->x, event->y, brush, brush3d, isobrush);
	      return (true); 
	      } 
	    return (false); 
	  };

          bool on_button_release (GdkEventButton* event) 
	  { 
	    if (editing) 
	      { 
	      editing = false;
	      // push undo buffer onto queue, then clear it
	      AddToUndo(processUndoBuff);
	      processUndoBuff.clear();
	      return (true); 
	      } 
	    return (false); 
	  }

	  bool on_key_press (GdkEventKey* event);

          bool undo();
          bool redo();
          bool floodfill();
	  bool copyslice(gint offset);

        protected:
          const ROIAnalysis& parent;
          bool  set, editing;
          Gtk::TreeModel::Row row;

	  // not the most efficient way of doing this, but will do for starters.
	  unsigned MaxUndoSize;
          class EdVox {
            public:
              EdVox (bool value, gsize offset) : value (value), offset (offset) { }
              bool value;
              gsize offset;
          };

          typedef std::vector<EdVox> EdVecType;
          typedef std::deque< EdVecType > EditQueueType;

          // undo/redo queues
          EditQueueType UndoQueue, RedoQueue;
	  // a global undo buffer so we don't need to change the 
          // mouse buttonpress interface
	  EdVecType processUndoBuff;
	  // applies the undo, and modifies EV to become a suitable redo list.
	  void ApplyUndo(EdVecType &EV);
	  void AddToUndo(EdVecType EV);

          void AddVox(MR::Image::Position ima, EdVecType &EV) { EV.push_back (EdVox (ima.value(), ima.getoffset())); }

          void AddVox(MR::Image::Position ima, EdVecType &EV, const float value)
          {
            float val = ima.value();
              // only push if it has changed.
              // specifically for drawing
            if (val != value) 
              EV.push_back (EdVox (val, ima.getoffset()));
          }

          class ROI {
            public:
              ROI (RefPtr<MR::Image::Object> image, guint32 C) : mask (new Image (image)), render (false), colour (C) { mask->image->set_read_only (false);}
              RefPtr<Image> mask;
              Slice::Renderer render;
              guint32 colour;
          };

          class Columns : public Gtk::TreeModel::ColumnRecord {
            public:
              Columns() { add (show); add (pix); add (name); add (roi); }

              Gtk::TreeModelColumn<bool> show;
              Gtk::TreeModelColumn<Glib::RefPtr<Gdk::Pixbuf> > pix;
              Gtk::TreeModelColumn<String> name;
              Gtk::TreeModelColumn<RefPtr<ROI> > roi;
          };

          Columns columns;
          Glib::RefPtr<Gtk::ListStore> model;
          Gtk::Menu popup_menu;

          bool on_button_press_event(GdkEventButton *event);
          void on_open ();
          void on_new ();
          void on_close ();
          void on_set_colour ();
          void on_clear ();
          void on_tick (const String& path);

          bool set_selected_row ();

	  void process (gdouble x, gdouble y, float brush, bool brush3d, bool isobrush);
          Point position (gdouble x, gdouble y);

          void load (RefPtr<MR::Image::Object> image);
      };

    }
  }
}

#endif

