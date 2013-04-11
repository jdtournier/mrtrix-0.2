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

#include "mrview/sidebar/roi_analysis.h"
#include "mrview/window.h"

namespace MR {
  namespace Viewer {
    namespace SideBar {

      ROIAnalysis::ROIAnalysis () : 
        Base (1),
        show_ROIs ("show ROIs"),
        roi_frame ("ROIs"), 
        transparency_frame ("opacity"), 
        brush_size_frame ("brush size"), 
        transparency (0.0, 256, 1.0),
        brush_size (1.0, 20.0, 1.0),
	brush3d("3D brush"),
	isobrush("Isotropic brush"),
        fill_button ("fill"),
        undo_button ("undo"),
        redo_button ("redo"),
        duplicate_from_previous_button ("copy previous"),
        duplicate_from_next_button ("copy next"),
        roi_list (*this)
      { 
        show_ROIs.set_active (true);

        transparency.set_draw_value (false);
        transparency.set_value (255.0);
        transparency.set_update_policy (Gtk::UPDATE_DELAYED);

        brush_size.set_draw_value (true);
        brush_size.set_value (1.0);
	brush3d.set_active(false);
        roi_scrolled_window.add (roi_list);
        roi_scrolled_window.set_policy (Gtk::POLICY_AUTOMATIC, Gtk::POLICY_AUTOMATIC);
        roi_scrolled_window.set_shadow_type (Gtk::SHADOW_IN);
        roi_scrolled_window.set_border_width (3);
        roi_frame.add (roi_scrolled_window);

        fill_button.set_tooltip_text ("fill region from crosshairs location\n\nShortcut: F");
        undo_button.set_tooltip_text ("undo last edit\n\nShortcut: Ctrl+Z");
        redo_button.set_tooltip_text ("re-apply last undo\n\nShortcut: Ctrl+Shift+Z");
        duplicate_from_previous_button.set_tooltip_text ("copy contents of previous slice onto current slice\n\nShortcut: P");
        duplicate_from_next_button.set_tooltip_text ("copy contents of next slice onto current slice\n\nShortcut: N");

        transparency_frame.add (transparency);
        brush_size_frame.add (brush_size);
        pack_start (show_ROIs, Gtk::PACK_SHRINK);
        pack_start (roi_frame);
        pack_start (transparency_frame, Gtk::PACK_SHRINK);
        pack_start (brush_size_frame, Gtk::PACK_SHRINK);
        pack_start (brush3d, Gtk::PACK_SHRINK);
        pack_start (isobrush, Gtk::PACK_SHRINK);
        pack_start (undo_button, Gtk::PACK_SHRINK);
        pack_start (redo_button, Gtk::PACK_SHRINK);
        pack_start (fill_button, Gtk::PACK_SHRINK);
        pack_start (duplicate_from_previous_button, Gtk::PACK_SHRINK);
        pack_start (duplicate_from_next_button, Gtk::PACK_SHRINK);

	isobrush.set_active();
        show_all();

        Window::Main->pane().activate (this);

        transparency.signal_value_changed().connect (sigc::mem_fun (*this, &ROIAnalysis::on_change));
        show_ROIs.signal_toggled().connect (sigc::mem_fun (*this, &ROIAnalysis::on_change));
        undo_button.signal_clicked().connect (sigc::mem_fun (*this, &ROIAnalysis::on_undo));
        redo_button.signal_clicked().connect (sigc::mem_fun (*this, &ROIAnalysis::on_redo));
        fill_button.signal_clicked().connect (sigc::mem_fun (*this, &ROIAnalysis::on_fill));
        duplicate_from_previous_button.signal_clicked().connect (sigc::mem_fun (*this, &ROIAnalysis::on_copy_previous));
        duplicate_from_next_button.signal_clicked().connect (sigc::mem_fun (*this, &ROIAnalysis::on_copy_next));
      }


      

      ROIAnalysis::~ROIAnalysis () {  }

      void ROIAnalysis::draw () { if (show_ROIs.get_active()) roi_list.draw ((int) transparency.get_value()); }
      void ROIAnalysis::on_change () { Window::Main->update (this); }

      bool ROIAnalysis::on_button_press (GdkEventButton* event) { return (roi_list.on_button_press (event, brush_size.get_value(), brush3d.get_active(), isobrush.get_active())); }
      bool ROIAnalysis::on_motion (GdkEventMotion* event) { return (roi_list.on_motion (event, brush_size.get_value(), brush3d.get_active(), isobrush.get_active())); }
      bool ROIAnalysis::on_button_release (GdkEventButton* event) { return (roi_list.on_button_release (event)); }

      bool ROIAnalysis::on_key_press (GdkEventKey* event)	
      {
	return(roi_list.on_key_press (event));
      }

      void ROIAnalysis::on_undo () { roi_list.undo(); }
      void ROIAnalysis::on_redo () { roi_list.redo(); }
      void ROIAnalysis::on_fill () { roi_list.floodfill(); }
      void ROIAnalysis::on_copy_previous () { roi_list.copyslice (-1); }
      void ROIAnalysis::on_copy_next () { roi_list.copyslice (1); }

    }
  }
}


