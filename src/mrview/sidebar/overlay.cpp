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

#include <gtkmm/stock.h>

#include "mrview/sidebar/overlay.h"
#include "mrview/window.h"
#include "mrview/colourmap.h"
#include "dialog/file.h"

namespace MR {
  namespace Viewer {
    namespace SideBar {


      Overlay::OverlayList::OverlayList (Overlay& sidebar) :
        parent (sidebar) {
        using namespace Gtk::Menu_Helpers;
        popup_menu.items().push_back (StockMenuElem(Gtk::Stock::OPEN, sigc::mem_fun(*this, &OverlayList::on_open) ) );
        popup_menu.items().push_back (StockMenuElem(Gtk::Stock::CLOSE, sigc::mem_fun(*this, &OverlayList::on_close) ) );
        popup_menu.items().push_back (SeparatorElem());
        popup_menu.items().push_back (MenuElem("_Colourmap"));
        popup_menu.items().back().set_submenu (colourmap_menu);
        popup_menu.items().push_back (SeparatorElem());
        popup_menu.items().push_back (StockMenuElem(Gtk::Stock::CLEAR, sigc::mem_fun(*this, &OverlayList::on_clear) ) );
        popup_menu.accelerate (*this);

        Gtk::RadioMenuItem::Group colourmap_group;
        colourmap_menu.items().push_back (RadioMenuElem (colourmap_group, "_Gray", sigc::bind<int> (sigc::mem_fun (*this, &OverlayList::on_colourmap), 0)));
        colourmap_menu.items().push_back (RadioMenuElem (colourmap_group, "_Hot", sigc::bind<int> (sigc::mem_fun (*this, &OverlayList::on_colourmap), 1)));
        colourmap_menu.items().push_back (RadioMenuElem (colourmap_group, "_Cool", sigc::bind<int> (sigc::mem_fun (*this, &OverlayList::on_colourmap), 2)));
        colourmap_menu.items().push_back (RadioMenuElem (colourmap_group, "_Jet", sigc::bind<int> (sigc::mem_fun (*this, &OverlayList::on_colourmap), 3)));
        colourmap_menu.items().push_back (SeparatorElem());
        colourmap_menu.items().push_back (RadioMenuElem (colourmap_group, "_RGB", sigc::bind<int> (sigc::mem_fun (*this, &OverlayList::on_colourmap), COLOURMAP_RGB)));
        colourmap_menu.items().push_back (RadioMenuElem (colourmap_group, "Comple_x", sigc::bind<int> (sigc::mem_fun (*this, &OverlayList::on_colourmap), COLOURMAP_COMPLEX)));

        model = Gtk::ListStore::create (columns);
        set_model (model);

        int tick_column_index = append_column_editable ("", columns.show) - 1;
        append_column ("file", columns.name);

        set_tooltip_text ("right-click for more options");

        set_headers_visible (false);

        Gtk::CellRendererToggle* tick = dynamic_cast<Gtk::CellRendererToggle*> (get_column_cell_renderer (tick_column_index));
        tick->signal_toggled().connect (sigc::mem_fun (*this, &OverlayList::on_tick));
        get_selection()->signal_changed().connect (sigc::mem_fun (*this, &OverlayList::on_selection));
      }




      Overlay::OverlayList::~OverlayList() { }




      void Overlay::OverlayList::draw (int transparency)
      {
        Gtk::TreeModel::Children overlays = model->children();
        if (overlays.size() == 0) return;

        Pane& pane (Window::Main->pane());
        Slice::Current slice (pane);

        glEnable(GL_BLEND);
        glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
        glEnable (GL_TEXTURE_2D);
        glDisable (GL_DEPTH_TEST);
        glTexEnvf (GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_MODULATE);
        glDepthMask (GL_FALSE);
        glColorMask (GL_TRUE, GL_TRUE, GL_TRUE, GL_TRUE);

        for (Gtk::TreeModel::Children::iterator iter = overlays.begin(); iter != overlays.end(); ++iter) {
          bool show = (*iter)[columns.show];
          if (show) {
            RefPtr<OL> overlay = (*iter)[columns.overlay];

            Slice::Info S;
            S.image = overlay->image;
            S.focus = slice.focus;
            S.colourmap = (*iter)[columns.colourmap];
            float min = (*iter)[columns.min];
            float max = (*iter)[columns.max];
            if (gsl_finite (min) && gsl_finite (max)) {
              S.scaling.multiplier = 255.0/(max-min); 
              S.scaling.offset = -S.scaling.multiplier*min; 
            }

            if (!slice.orientation && slice.image->image->header().transform() == overlay->image->image->header().transform()) { 
              S.projection = slice.projection;
            }
            else {
              const GLdouble* M (pane.get_modelview());
              float matrix[] = { 
                M[0], M[1], M[2],
                M[4], M[5], M[6],
                M[8], M[9], M[10]
              };
              S.orientation.from_matrix (matrix);
              S.projection = 2;
            }
            S.interpolate = true;

            Slice::Current current (S);

            overlay->render.update (current);
            glColor4ub (255, 255, 255, transparency);
            overlay->render.draw();

            if (!gsl_finite (min) || !gsl_finite (max)) {
              (*iter)[columns.min] = - S.scaling.offset / S.scaling.multiplier;
              (*iter)[columns.max] = 255.0 / S.scaling.multiplier + (*iter)[columns.min];
              if (get_selection()->is_selected (iter)) 
                parent.set_scaling ((*iter)[columns.min], (*iter)[columns.max]);
            }
          }
        }
        glDepthMask (GL_TRUE);
        glDisable (GL_TEXTURE_2D);
      }




      void Overlay::OverlayList::set_scaling (float min, float max)
      {
        min_max_changed = false;
        min_from_parent = min;
        max_from_parent = max;
        get_selection()->selected_foreach_iter (sigc::mem_fun (*this, &OverlayList::set_scaling_callback) );
        if (min_max_changed) 
          Window::Main->update (&parent);
      }



      void Overlay::OverlayList::set_scaling_callback (const Gtk::TreeModel::iterator& iter)
      {
        Gtk::TreeModel::Row row = *iter;
        if (row[columns.min] != min_from_parent) { min_max_changed = true; row[columns.min] = min_from_parent; }
        if (row[columns.max] != max_from_parent) { min_max_changed = true; row[columns.max] = max_from_parent; }
      }





      bool Overlay::OverlayList::on_button_press_event (GdkEventButton *event)
      {
        if ((event->type == GDK_BUTTON_PRESS) && (event->button == 3)) {
          Gtk::TreeModel::Path path;
          Gtk::TreeViewColumn* col;
          int x, y;
          bool is_row = get_path_at_pos ((int) event->x, (int) event->y, path, col, x, y);
          if (is_row) {
            if (!get_selection()->is_selected (path)) 
              TreeView::on_button_press_event(event);
          }
          else get_selection()->unselect_all();

          popup_menu.items()[2].set_sensitive (is_row);
          popup_menu.items()[4].set_sensitive (is_row);
          popup_menu.popup (event->button, event->time);
          return (true);
        }
        return (TreeView::on_button_press_event(event));
      }




      void Overlay::OverlayList::on_colourmap (int mode) 
      {
        colourmap = mode;
        get_selection()->selected_foreach_iter (sigc::mem_fun (*this, &OverlayList::colourmap_callback) );
        Window::Main->update (&parent);
      }



      void Overlay::OverlayList::colourmap_callback (const Gtk::TreeModel::iterator& iter) 
      {
        Gtk::TreeModel::Row row = *iter;
        row[columns.colourmap] = colourmap;
      }


      void Overlay::OverlayList::on_open ()
      {
        Dialog::File dialog ("Open overlay image", true, true);

        if (dialog.run() == Gtk::RESPONSE_OK) {
          std::vector<RefPtr<MR::Image::Object> > selection = dialog.get_images();
          if (selection.size()) {
            for (guint n = 0; n < selection.size(); n++) 
              load (selection[n]);
            Window::Main->update (&parent);
          }
        }
      }






      void Overlay::OverlayList::on_close ()
      {
        std::list<Gtk::TreeModel::Path> paths = get_selection()->get_selected_rows();
        std::list<Gtk::TreeModel::RowReference> rows;

        for (std::list<Gtk::TreeModel::Path>::iterator row = paths.begin(); row != paths.end(); row++) 
          rows.push_back (Gtk::TreeModel::RowReference (model, *row));

        for (std::list<Gtk::TreeModel::RowReference>::iterator row = rows.begin(); row != rows.end(); row++) 
          model->erase (model->get_iter (row->get_path()));

        Window::Main->update (&parent);
      }








      void Overlay::OverlayList::on_selection ()
      {
        get_selection()->selected_foreach_iter (sigc::mem_fun (*this, &OverlayList::selected_row_callback) );
      }



      void Overlay::OverlayList::selected_row_callback (const Gtk::TreeModel::iterator& iter)
      {
        Gtk::TreeModel::Row row = *iter;
        parent.set_scaling ((*iter)[columns.min], (*iter)[columns.max]);
      }



      void Overlay::OverlayList::on_clear ()
      {
        model->clear(); 
        Window::Main->update (&parent);
      }




      void Overlay::OverlayList::on_tick (const String& path) 
      {
        get_selection()->select (Gtk::TreePath (path));
        Window::Main->update (&parent); 
      }





      void Overlay::OverlayList::load (RefPtr<MR::Image::Object> image)
      { 
        Gtk::TreeModel::Row row = *(model->append());
        row[columns.show] = true;
        row[columns.name] = image->name();
        row[columns.colourmap] = 1;
        row[columns.min] = GSL_NAN;
        row[columns.max] = GSL_NAN;
        row[columns.overlay] = RefPtr<OL> (new OL (image));;
        get_selection()->select (row);
      }










      Overlay::Overlay () : 
        Base (1),
        show_overlays ("show overlays"),
        overlay_frame ("Images"),
        transparency_frame ("opacity"), 
        scaling_frame ("scaling"),
        min_label ("min: "),
        max_label ("max: "),
        transparency (0.0, 256, 1.0),
        scaling_table (2, 2),
        overlay_list (*this)
      { 
        show_overlays.set_active (true);

        min_value.set_range (-INFINITY, INFINITY);
        min_value.set_digits (4);
        max_value.set_range (-INFINITY, INFINITY);
        max_value.set_digits (4);

        transparency.set_draw_value (false);
        transparency.set_value (128.0);
        transparency.set_update_policy (Gtk::UPDATE_CONTINUOUS);

        overlay_scrolled_window.add (overlay_list);
        overlay_scrolled_window.set_policy (Gtk::POLICY_AUTOMATIC, Gtk::POLICY_AUTOMATIC);
        overlay_scrolled_window.set_shadow_type (Gtk::SHADOW_IN);
        overlay_scrolled_window.set_border_width (3);
        overlay_frame.add (overlay_scrolled_window);

        transparency_frame.add (transparency);
        
        scaling_table.attach (min_label, 0, 1, 0, 1);
        scaling_table.attach (min_value, 1, 2, 0, 1);
        scaling_table.attach (max_label, 0, 1, 1, 2);
        scaling_table.attach (max_value, 1, 2, 1, 2);

        scaling_frame.add (scaling_table);

        pack_start (show_overlays, Gtk::PACK_SHRINK);
        pack_start (overlay_frame);
        pack_start (transparency_frame, Gtk::PACK_SHRINK);
        pack_start (scaling_frame, Gtk::PACK_SHRINK);
        show_all();

        Window::Main->pane().activate (this);

        transparency.signal_value_changed().connect (sigc::mem_fun (*this, &Overlay::on_change));
        show_overlays.signal_toggled().connect (sigc::mem_fun (*this, &Overlay::on_change));
        min_value.signal_value_changed().connect (sigc::mem_fun (*this, &Overlay::on_scaling));
        max_value.signal_value_changed().connect (sigc::mem_fun (*this, &Overlay::on_scaling));
      }




      Overlay::~Overlay () { }



      void Overlay::draw () 
      { 
        if (show_overlays.get_active())
          overlay_list.draw ((int) transparency.get_value()); 
      }



      void Overlay::set_scaling (float min, float max) 
      {
        min_value.set_value (min);
        max_value.set_value (max);
      }



      void Overlay::on_scaling () 
      {
        double min = min_value.get_value(), max = max_value.get_value();
        double mean = 0.5*(min+max);
        int digits = floor (log (mean) / log(10.0));
        double inc = pow(10.0, digits-2), linc = 10.0*inc;

        min_value.set_increments (inc, linc);
        max_value.set_increments (inc, linc);

        digits = MAX (-digits+2, 0);
        min_value.set_digits (digits);
        max_value.set_digits (digits);

        overlay_list.set_scaling (min, max);
      }


      void Overlay::on_change () { Window::Main->update (this); }



    }
  }
}



