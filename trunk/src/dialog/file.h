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

#ifndef __dialog_file_decl_h__
#define __dialog_file_decl_h__

#include <gtkmm/dialog.h>
#include <gtkmm/treeview.h>
#include <gtkmm/liststore.h>
#include <gtkmm/button.h>
#include <gtkmm/box.h>
#include <gtkmm/label.h>
#include <gtkmm/paned.h>
#include <gtkmm/entry.h>
#include <gtkmm/scrolledwindow.h>
#include <glibmm/fileutils.h>

#include "file/dicom/tree.h"


namespace MR {
  namespace Image {
    class Object;
  }

  namespace File {
    namespace Dicom {
      class Series;
      class Tree;
    }
  }


  namespace Dialog {

    class File : public Gtk::Dialog
    {
      public:
        File (const String& message, bool multiselection, bool expand_DICOM);
        ~File ();

        void                               set_selection (const String& filename) { selection_entry.set_text (filename); }
        std::vector<String>                get_selection ();
        std::vector<RefPtr<Image::Object> >  get_images ();

        static const String& get_cwd () { if (!cwd.size()) cwd = Glib::get_current_dir(); return (cwd); }

      protected:

        class FolderColumns : public Gtk::TreeModel::ColumnRecord {
          public:
            FolderColumns () { add (name); }
            Gtk::TreeModelColumn<String> name;
        };

        class FileColumns : public Gtk::TreeModel::ColumnRecord {
          public:
            FileColumns () { add (name); add (series); }
            Gtk::TreeModelColumn<String> name;
            Gtk::TreeModelColumn<RefPtr<MR::File::Dicom::Series> > series;
        };


        FolderColumns          folder_columns;
        FileColumns            file_columns;
 
        Gtk::Entry             path_entry, selection_entry;
        Gtk::HPaned            splitter;
        Gtk::Button            up_button, home_button, refresh_button;
        Gtk::HBox              button_box, path_box, selection_box;
        Gtk::VBox              main_box;
        Gtk::Label             path_label, selection_label;
        Gtk::ScrolledWindow    files_viewport, folders_viewport;
        Glib::Dir*             dir;
        String                 next_file;
        bool                   filter_images, folders_read, updating_selection;

        Glib::RefPtr<Gtk::ListStore> file_list;
        Glib::RefPtr<Gtk::ListStore> folder_list;
        Gtk::TreeView          files, folders;

        sigc::connection       idle_connection;

        void                   on_up ();
        void                   on_home ();
        void                   on_refresh ();
        void                   on_folder_selected (const Gtk::TreeModel::Path& path, Gtk::TreeViewColumn* column);
        void                   on_file_selected (const Gtk::TreeModel::Path& path, Gtk::TreeViewColumn* column);
        void                   on_file_highlighted ();
        void                   on_selection_entry ();
        void                   on_selection_activated ();
        void                   on_path ();
        bool                   on_idle ();
        bool                   on_file_tooltip (int x, int y, bool keyboard_tooltip, const Glib::RefPtr<Gtk::Tooltip>& tooltip);

        static String          cwd;
        static int             window_position_x, window_position_y, window_size_x, window_size_y;

        void                   update ();
        void                   check_image (const String& path, const String& base);
        void                   check_dicom (const String& path, const String& base);
        void                   update_dicom ();
        MR::File::Dicom::Tree  dicom_tree;
    };

  }
}

#endif

