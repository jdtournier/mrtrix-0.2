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

#include <gtkmm/treerowreference.h>
#include <gtkmm/colorselection.h>
#include <gtkmm/stock.h>
#include <gdkmm/pixbuf.h>

#include "mrview/sidebar/roi_analysis/roi_list.h"
#include "mrview/sidebar/roi_analysis.h"
#include "mrview/window.h"
#include "dialog/file.h"

#include <queue>

namespace MR {
  namespace Viewer {
    namespace SideBar {

      DP_ROIList::DP_ROIList (const ROIAnalysis& sidebar) : parent (sidebar), set (true), editing (false)
      {
        using namespace Gtk::Menu_Helpers;
        popup_menu.items().push_back (StockMenuElem(Gtk::Stock::OPEN, sigc::mem_fun(*this, &DP_ROIList::on_open) ) );
        popup_menu.items().push_back (StockMenuElem(Gtk::Stock::NEW, sigc::mem_fun(*this, &DP_ROIList::on_new) ) );
        popup_menu.items().push_back (StockMenuElem(Gtk::Stock::CLOSE, sigc::mem_fun(*this, &DP_ROIList::on_close) ) );
        popup_menu.items().push_back (SeparatorElem());
        popup_menu.items().push_back (MenuElem("_Change colour", sigc::mem_fun(*this, &DP_ROIList::on_set_colour)));
        popup_menu.items().push_back (SeparatorElem());
        popup_menu.items().push_back (StockMenuElem(Gtk::Stock::CLEAR, sigc::mem_fun(*this, &DP_ROIList::on_clear) ) );
        popup_menu.accelerate (*this);

        model = Gtk::ListStore::create (columns);
        set_model (model);

        int tick_column_index = append_column_editable ("", columns.show) - 1;
        append_column ("", columns.pix);
        append_column ("file", columns.name);

        set_tooltip_text ("right-click for more options");

        set_headers_visible (false);

        Gtk::CellRendererToggle* tick = dynamic_cast<Gtk::CellRendererToggle*> (get_column_cell_renderer (tick_column_index));
        tick->signal_toggled().connect (sigc::mem_fun (*this, &DP_ROIList::on_tick));

        MaxUndoSize=10;

      }






      DP_ROIList::~DP_ROIList()  { }





      void DP_ROIList::draw (int transparency)
      {
        Gtk::TreeModel::Children rois = model->children();
        if (rois.size() == 0) return;

        Pane& pane (Window::Main->pane());
        Slice::Current slice (pane);

        glEnable(GL_BLEND);
        glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
        glEnable (GL_TEXTURE_2D);
        glDisable (GL_DEPTH_TEST);
        glTexEnvf (GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_MODULATE);
        glDepthMask (GL_FALSE);
        glColorMask (GL_TRUE, GL_TRUE, GL_TRUE, GL_TRUE);

        for (Gtk::TreeModel::Children::iterator iter = rois.begin(); iter != rois.end(); ++iter) {
          bool show = (*iter)[columns.show];
          if (show) {
            RefPtr<ROI> roi = (*iter)[columns.roi];

            Slice::Info S;
            S.image = roi->mask;
            S.focus = slice.focus;

            bool transforms_differ = false;
            if (slice.orientation) 
              transforms_differ = true;
            else {
              for (guint i = 0; i < 4; ++i) {
                for (guint j = 0; j < 3; ++j) {
                  if (fabs (slice.image->image->header().transform()(i,j) - roi->mask->image->header().transform()(i,j)) > 1e-5) {
                    transforms_differ = true;
                    break;
                  }
                }
              }
            }

            if (transforms_differ) {
              const GLdouble* M (pane.get_modelview());
              float matrix[] = { 
                float(M[0]), float(M[1]), float(M[2]),
                float(M[4]), float(M[5]), float(M[6]),
                float(M[8]), float(M[9]), float(M[10])
              };
              S.orientation.from_matrix (matrix);
              S.projection = 2;
            }
            else 
              S.projection = slice.projection;

            S.interpolate = false;

            Slice::Current current (S);

            roi->render.update (current);
            glColor4ub ((roi->colour & 0xFF000000) >> 24, (roi->colour & 0x00FF0000) >> 16, (roi->colour & 0x0000FF00) >> 8, transparency);
            roi->render.draw();
          }
        }
        glDepthMask (GL_TRUE);
        glDisable (GL_TEXTURE_2D);
      }







      void DP_ROIList::load (RefPtr<MR::Image::Object> image)
      { 
        guint32 colour = 0xFFFF00FF;
        Glib::RefPtr<Gdk::Pixbuf> pix = Gdk::Pixbuf::create  (Gdk::COLORSPACE_RGB, false, 8, 16, 16);
        pix->fill (colour);
        Gtk::TreeModel::Row row = *(model->append());
        row[columns.show] = true;
        row[columns.pix] = pix;
        row[columns.name] = image->name();
        row[columns.roi] = RefPtr<ROI> (new ROI (image, colour));;
        get_selection()->select (row);
      }







      bool DP_ROIList::on_button_press_event (GdkEventButton* event)
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




      void DP_ROIList::on_set_colour () 
      {
        Gtk::ColorSelectionDialog dialog ("Choose colour for ROI");
        if (dialog.run() == Gtk::RESPONSE_OK) {
          Gtk::TreeModel::iterator iter = get_selection()->get_selected();
          RefPtr<ROI> roi = (*iter)[columns.roi];
          Gdk::Color colour (dialog.get_colorsel()->get_current_color());
          GLubyte C[] = { 
            GLubyte (colour.get_red() >> 8),
            GLubyte (colour.get_green() >> 8),
            GLubyte (colour.get_blue() >> 8)
          };
          roi->colour = C[0] << 24 | C[1] << 16 | C[2] << 8 | 255;
          Glib::RefPtr<Gdk::Pixbuf> pix = (*iter)[columns.pix];
          pix->fill (roi->colour);
          Window::Main->update (&parent);
        }
      }









      void DP_ROIList::on_open () 
      {
        Dialog::File dialog ("Open mask image", true, true);

        if (dialog.run() == Gtk::RESPONSE_OK) {
          std::vector<RefPtr<MR::Image::Object> > selection = dialog.get_images();
          if (selection.size()) {
            for (guint n = 0; n < selection.size(); n++) 
              load (selection[n]);
            Window::Main->update (&parent);
          }
        }
      }







      void DP_ROIList::on_new () 
      {
        Dialog::File dialog ("Create mask image", false, false);

        if (dialog.run() == Gtk::RESPONSE_OK) {
          std::vector<String> selection = dialog.get_selection();
          if (selection.size()) {
            MR::Image::Header header (*Window::Main->image->image);
            header.data_type = DataType::Bit;
            header.axes.set_ndim (3);
            RefPtr<MR::Image::Object> obj (new MR::Image::Object);
            try {
              obj->create (selection[0], header); 
              load (obj);
            }
            catch (...) { error ("error creating mask image \"" + selection[0] + "\": " + Glib::strerror (errno)); }
          }
        }
      }








      void DP_ROIList::on_close ()
      {
        std::list<Gtk::TreeModel::Path> paths = get_selection()->get_selected_rows();
        std::list<Gtk::TreeModel::RowReference> rows;

        for (std::list<Gtk::TreeModel::Path>::iterator row = paths.begin(); row != paths.end(); row++) 
          rows.push_back (Gtk::TreeModel::RowReference (model, *row));

        for (std::list<Gtk::TreeModel::RowReference>::iterator row = rows.begin(); row != rows.end(); row++) 
          model->erase (model->get_iter (row->get_path()));

        Window::Main->update (&parent);
      }





      void DP_ROIList::on_clear () 
      {
        model->clear(); 
        Window::Main->update (&parent);
      }




      void DP_ROIList::on_tick (const String& path) { Window::Main->update (&parent); }


      bool DP_ROIList::on_button_press (GdkEventButton* event, float brush, bool brush3d, bool isobrush) 
      {
        Gtk::TreeModel::iterator iter = get_selection()->get_selected();
        if (!iter) return (false);
        guint state = (event->state & MODIFIERS) & (~GDK_BUTTON1_MASK);
        if (state != GDK_SHIFT_MASK && 
            state != ( GDK_SHIFT_MASK | CTRL_CMD_MASK ))
          return (false);

        row = *iter;
        bool show = row[columns.show];
        if (!show) return (false);

        set = state == GDK_SHIFT_MASK;
        editing = true;
        process (event->x, event->y, brush, brush3d, isobrush);
        return (true);
      }

      bool DP_ROIList::set_selected_row ()
      {
        Gtk::TreeModel::iterator iter = get_selection()->get_selected();
        if (!iter) return (false);
        row = *iter;
        bool show = row[columns.show];
        if (!show) return (false);
        return (true);
      }


      bool DP_ROIList::on_key_press (GdkEventKey* event)	
      {
        // ignore releases
        if (event->type == GDK_KEY_RELEASE) return (false);
        switch (event->keyval) {
          case GDK_KEY_F:
          case GDK_KEY_f: 
            return floodfill();
	    break;
          case GDK_KEY_N:
          case GDK_KEY_n: 
            return copyslice (1);
	    break;
          case GDK_KEY_P:
          case GDK_KEY_p: 
            return copyslice (-1);
	    break;
          case GDK_KEY_Z:
            return (event->state & CTRL_CMD_MASK) ? redo() : false;
	    break;
          case GDK_KEY_z:
            return (event->state & CTRL_CMD_MASK) ? undo() : false;
	    break;
          default:
            break;
        }
        return false;
      }






      bool DP_ROIList::undo () 
      {
        if (!set_selected_row())
          return false;

        if (UndoQueue.size() > 0) {
          EdVecType ThisUndo = UndoQueue.front();
          UndoQueue.pop_front();
          ApplyUndo(ThisUndo);
          RedoQueue.push_front(ThisUndo);
        }

        return true;
      }


      bool DP_ROIList::redo () 
      {
        if (!set_selected_row())
          return false;

        if (RedoQueue.size() > 0) {
          EdVecType ThisUndo = RedoQueue.front();
          RedoQueue.pop_front();
          ApplyUndo(ThisUndo);
          UndoQueue.push_front(ThisUndo);
        }

        return true;
      }


      void DP_ROIList::ApplyUndo(EdVecType &EV)
      {
        editing=true;
        RefPtr<ROI> roi = row[columns.roi];
        MR::Image::Position ima (*roi->mask->image);
        for (unsigned k=0; k < EV.size(); k++)
        {
          bool newval=EV[k].value;
          ima.setoffset(EV[k].offset);
          EV[k].value = ima.value();
          ima.value(newval);
        }
        editing=false;
        Window::Main->update (&parent);
      }


      void DP_ROIList::AddToUndo(EdVecType EV)
      {
        UndoQueue.push_front(EV);
        // as soon as something is edited, we clear the redo list
        RedoQueue.clear();
        if (UndoQueue.size() > MaxUndoSize)
          UndoQueue.pop_back();
      }


      bool DP_ROIList::copyslice (gint offset)
      {
        if (!set_selected_row())
          return false;

        editing=true;

        RefPtr<ROI> roi = row[columns.roi];
        // figure out our current slice
        Point pos (roi->mask->interp->R2P (position (0, 0)));
        MR::Image::Position ima (*roi->mask->image);
        MR::Image::Position imb (*roi->mask->image);
        int p[] = { round (pos[0]), round(pos[1]), round(pos[2]) };

        Pane& pane (Window::Main->pane());
        const Slice::Current S (pane);
        unsigned projection(S.projection);

        int sourceslice =  p[projection] + offset;
        // don't do anything if the source slice is invalid
        if ((sourceslice < 0) || (sourceslice >= ima.dim(projection))) {
          std::cout << "Slice copy : out of range" << std::endl;
	  editing=false;
          return false;
        }
        unsigned int ax1 = 0, ax2 = 1;
        switch (projection)
        {
          case 0:
            ax1 = 1; ax2 = 2;
            break;
          case 1:
            ax1 = 0; ax2 = 2;
            break;
          case 2:
            ax1 = 0; ax2 = 1;
            break;
          default:
            break;
        }
        // ima points to current slice
        ima.set (projection, p[projection]);
        ima.set (ax1, 0);
        ima.set (ax2, 0);

        // imb points to the source slice
        imb.set (projection, sourceslice);
        imb.set (ax1, 0);
        imb.set (ax2, 0);
        EdVecType UndoVec;

        for (imb.set(ax1, 0), ima.set(ax1, 0); ima[ax1] < ima.dim(ax1); ima.inc(ax1), imb.inc(ax1)) {
          for (imb.set(ax2, 0), ima.set(ax2, 0); ima[ax2] < ima.dim(ax2); ima.inc(ax2), imb.inc(ax2)) {
            AddVox(ima, UndoVec);
            ima.value(imb.value());
          }
        }
        AddToUndo(UndoVec);
        Window::Main->update (&parent);

        editing=false;
        return true;
      }




      bool DP_ROIList::floodfill ()
      {
        if (!set_selected_row())
          return false;

        editing=true;

        Pane& pane (Window::Main->pane());
        const Slice::Current S (pane);

        RefPtr<ROI> roi = row[columns.roi];
        Point pos (roi->mask->interp->R2P (S.focus));
        MR::Image::Position ima (*roi->mask->image);
        int p[] = { round (pos[0]), round(pos[1]), round(pos[2]) };

        EdVecType UndoInfo;

        // set the position we are starting from
        ima.set(0, p[0]);
        ima.set(1, p[1]);
        ima.set(2, p[2]);
        // should we adapt to be flood delete too??
        float bgvalue = ima.value();
        float fillvalue = !bgvalue;
        {

          // we are only flooding in 2D, so pick which axes this
          // corresponds to.
          unsigned ax1=0, ax2=1;
          switch (S.projection)
          {
            case 0:
              ax1=1; ax2=2;
              break;
            case 1:
              ax1=0; ax2=2;
              break;
            case 2:
              ax1=0; ax2=1;
              break;
            default:
              break;
          }
          // set the other values
          // standard flooding algorithm :
          // put voxel on queue
          // Pop top of queue: visit neighbours : label unlabelled
          // neighbours and return place them on queue
          std::queue<MR::Image::Position> fqueue;
          AddVox(ima, UndoInfo);
          ima.value(fillvalue);
          fqueue.push(ima);

          int offset1[8]={-1, -1, -1, 0, 0, 1, 1, 1};
          int offset2[8]={-1, 0, 1, -1, 1, -1, 0, 1};

          while (!fqueue.empty()) {
            MR::Image::Position idx = fqueue.front();
            fqueue.pop();
            for (unsigned P=0;P<8;P++) {
              MR::Image::Position nidx = idx;	
              // set up the neighbour
              int a1 = idx[ax1]+offset1[P];
              int a2 = idx[ax2]+offset2[P];
              if (!((a1 < 0) || (a2 < 0) || (a1 >= idx.dim(ax1)) || (a2 >= idx.dim(ax2)) ) ) {
                nidx.set(ax1, a1);
                nidx.set(ax2, a2);
                if (nidx.value()==bgvalue) {
                  AddVox(nidx, UndoInfo);
                  nidx.value(fillvalue);
                  fqueue.push(nidx);
                }
              }
            }
          }
        }
        AddToUndo(UndoInfo);
        Window::Main->update (&parent);
        editing = false;

        return true;
      }





      void DP_ROIList::process (gdouble x, gdouble y, float brush, bool brush3d, bool isobrush)
      {
        RefPtr<ROI> roi = row[columns.roi];
        Point pos (roi->mask->interp->R2P (position (x, y)));

        MR::Image::Position ima (*roi->mask->image);
        int p[] = { round (pos[0]), round(pos[1]), round(pos[2]) };
        int e = ceil(brush/2.0);
        int e0(e), e1(e), e2(e);
        float dist = (brush*brush)/4.0;
        float Sc0(1), Sc1(1), Sc2(1);
        if (isobrush)
        {
          // interpret size as mm
          float v0 = ima.vox(0);
          float v1 = ima.vox(1);
          float v2 = ima.vox(2);

          e0 = ceil(e/v0);
          e1 = ceil(e/v1);
          e2 = ceil(e/v2);
          Sc0=(v0*v0);
          Sc1=(v1*v1);
          Sc2=(v2*v2);
        }

        Pane& pane (Window::Main->pane());
        const Slice::Current S (pane);

        if (brush3d) {
          for (ima.set (2, p[2]-e2); ima[2] <= p[2]+e2; ima.inc(2)) {
            if (ima[2] < 0 || ima[2] >= ima.dim(2)) continue;
            for (ima.set (1, p[1]-e1); ima[1] <= p[1]+e1; ima.inc(1)) {
              if (ima[1] < 0 || ima[1] >= ima.dim(1)) continue;
              for (ima.set (0, p[0]-e0); ima[0] <= p[0]+e0; ima.inc(0)) {
                if (ima[0] < 0 || ima[0] >= ima.dim(0)) continue;
                if ((ima[0]-p[0])*(ima[0]-p[0])*Sc0 + (ima[1]-p[1])*(ima[1]-p[1])*Sc1 + (ima[2]-p[2])*(ima[2]-p[2])*Sc2 < dist) {
                  AddVox(ima, processUndoBuff, set ? 1.0 : 0.0);
                  ima.value (set ? 1.0 : 0.0);
                }
              }
            }
          }
        }
        else {
          const unsigned projection(S.projection);
          switch (projection) {
            case 0:
              // sagittal
              ima.set (0, p[0]);
              for (ima.set (2, p[2]-e2); ima[2] <= p[2]+e2; ima.inc(2)) {
                if (ima[2] < 0 || ima[2] >= ima.dim(2)) continue;
                for (ima.set (1, p[1]-e1); ima[1] <= p[1]+e1; ima.inc(1)) {
                  if (ima[1] < 0 || ima[1] >= ima.dim(1)) continue;
                  if (Sc1*(ima[1]-p[1])*(ima[1]-p[1]) + Sc2*(ima[2]-p[2])*(ima[2]-p[2]) < dist) {
                    AddVox(ima, processUndoBuff, set ? 1.0 : 0.0);
                    ima.value (set ? 1.0 : 0.0);
                  }
                }
              }
              break;
            case 1:
              // coronal
              ima.set (1, p[1]);
              for (ima.set (2, p[2]-e2); ima[2] <= p[2]+e2; ima.inc(2)) {
                if (ima[2] < 0 || ima[2] >= ima.dim(2)) continue;
                for (ima.set (0, p[0]-e0); ima[0] <= p[0]+e0; ima.inc(0)) {
                  if (ima[0] < 0 || ima[0] >= ima.dim(0)) continue;
                  if (Sc0*(ima[0]-p[0])*(ima[0]-p[0]) + Sc2*(ima[2]-p[2])*(ima[2]-p[2]) < dist) {
                    AddVox(ima, processUndoBuff, set ? 1.0 : 0.0);
                    ima.value (set ? 1.0 : 0.0);
                  }
                }
              }
              break;
            case 2:
              // axial
              ima.set (2, p[2]);
              for (ima.set (1, p[1]-e1); ima[1] <= p[1]+e1; ima.inc(1)) {
                if (ima[1] < 0 || ima[1] >= ima.dim(1)) continue;
                for (ima.set (0, p[0]-e0); ima[0] <= p[0]+e0; ima.inc(0)) {
                  if (ima[0] < 0 || ima[0] >= ima.dim(0)) continue;
                  if (Sc0*(ima[0]-p[0])*(ima[0]-p[0]) + Sc1*(ima[1]-p[1])*(ima[1]-p[1]) < dist) {
                    AddVox(ima, processUndoBuff, set ? 1.0 : 0.0);
                    ima.value (set ? 1.0 : 0.0);
                  }
                }
              }
              break;
            default:
              assert (0);
          }
        }
        Window::Main->update (&parent);
      }



      Point DP_ROIList::position (gdouble x, gdouble y)
      {
        Pane& pane (Window::Main->pane());
        const Slice::Current S (pane);
        Point f = pane.model_to_screen (S.focus);
        f[0] = x; 
        f[1] = pane.height() - y;
        return (pane.screen_to_model (f));
      }

    }
  }
}


