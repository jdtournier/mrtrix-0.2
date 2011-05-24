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

    
    14-07-2008 J-Donald Tournier <d.tournier@brain.org.au>
    * fixed Track::index() for use on 64 systems.
      now uses gsize rather than guint in pointer arithmetic

    15-12-2008 J-Donald Tournier <d.tournier@brain.org.au>
    * a few bug fixes + memory performance improvements for the depth blend option
    
    12-03-2010 J-Donald Tournier <d.tournier@brain.org.au>
    * a few bug fixes for the colour handling and support for depth blend on
    * 64 bit systems
    
    28-07-2010 J-Donald Tournier <d.tournier@brain.org.au>
    * increased allocator slab size to 32MB to force release of RAM to OS on close()

*/

#ifndef __mrview_sidebar_tractography_tracklist_item_h__
#define __mrview_sidebar_tractography_tracklist_item_h__

#include <list>
#include <glib/gstdio.h>

#include "use_gl.h"
#include "ptr.h"
#include "point.h"
#include "math/vector.h"
#include "dwi/tractography/properties.h"

#define TRACK_ALLOCATOR_SLAB_SIZE 0x2000000U
#define TRANSPARENCY_EXPONENT 4.9 

namespace MR {
  namespace Viewer {
    namespace SideBar {

      class TrackListItem 
      {
        public:
          class Allocator;

          class Track {
            public:

              class Point {
                public:
                  float   pos[3];
                  GLubyte C[4];

                  void set_pos (const MR::Point& p) { pos[0] = p[0]; pos[1] = p[1]; pos[2] = p[2]; }
                  void set_colour (const MR::Point& dir, GLubyte alpha) { 
                    C[0] = abs(int(255*dir[0]));
                    C[1] = abs(int(255*dir[1])); 
                    C[2] = abs(int(255*dir[2])); 
                    C[3] = alpha;
                  }
                  void set_colour (const GLubyte c[3], GLubyte alpha) {
                    C[0] = c[0]; 
                    C[1] = c[1]; 
                    C[2] = c[2]; 
                    C[3] = alpha;
                  }
                  gsize index () const { return (this-root); }

                  bool operator< (const Point& b) const { 
                    return ((pos[0]-b.pos[0])*normal[0] + (pos[1]-b.pos[1])*normal[1] + (pos[2]-b.pos[2])*normal[2] < 0.0);
                  }
                  static MR::Point normal;
                  static const Point* root;
              };


              Track (Allocator& alloc, guint num);

              gsize          size () const               { return (num_p); }
              const Point&   operator[] (int n) const    { return (data[n]); }
              Point&         operator[] (int n)          { return (data[n]); }

              static const int stride = sizeof (Point);

            protected:
              Point* data;
              guint  num_p;
          };



          class Allocator {
            public:
              Allocator () : next (NULL), end (NULL), min ((guint8*) (-1)) { }
              ~Allocator () { clear(); }

              Track::Point* operator () (guint size) { 
                size *= sizeof(Track::Point);
                assert (size <= TRACK_ALLOCATOR_SLAB_SIZE);
                if (next + size > end) new_block();
                Track::Point* p = (Track::Point*) next;
                next += size;
                return (p);
              }

              void clear () { 
                for (std::list<guint8*>::iterator i = blocks.begin(); i != blocks.end(); ++i) 
                  delete [] *i; 
                blocks.clear(); 
                next = end = NULL;
                min = (guint8*) (-1);
              }

              const guint8* root() const { return (min); }

            private:
              std::list<guint8*> blocks;
              guint8* next;
              guint8* end;
              guint8* min;

              void new_block () {
                next = new guint8 [TRACK_ALLOCATOR_SLAB_SIZE];
                end = next + TRACK_ALLOCATOR_SLAB_SIZE;
                blocks.push_back (next); 
                if (guint rem = (next - (guint8*) NULL ) % sizeof(Track::Point)) 
                  next += sizeof(Track::Point)-rem;
                if (next < min) min = next;
              }
          };


          TrackListItem () { colour_by_dir = true; alpha = 1.0; }


          String file;
          time_t mtime;
          std::list<Track>  tracks;
          DWI::Tractography::Properties properties;
          GLubyte colour[3];
          bool  colour_by_dir, colour_by_dir_previous;
          float alpha, alpha_previous;

          void update_RGBA (bool force = false) 
          {
            if (force || colour_by_dir != colour_by_dir_previous || alpha != alpha_previous) {
              GLubyte A = (GLubyte) (255.0*get_alpha());
              GLubyte default_colour[] = { 255, 255, 255 };
              for (std::list<Track>::iterator i = tracks.begin(); i != tracks.end(); ++i) {
                Track& tck (*i);

                if (colour_by_dir) {
                  if (tck.size() > 1) {
                    gsize n;
                    Point dir (tck[1].pos[0]- tck[0].pos[0], tck[1].pos[1]- tck[0].pos[1], tck[1].pos[2]- tck[0].pos[2]);
                    tck[0].set_colour (dir.normalise(), A);
                    for (n = 1; n < tck.size()-1; n++) {
                      dir.set (tck[n+1].pos[0]- tck[n-1].pos[0], tck[n+1].pos[1]- tck[n-1].pos[1], tck[n+1].pos[2]- tck[n-1].pos[2]);
                      tck[n].set_colour (dir.normalise(), A);
                    }
                    dir.set (tck[n].pos[0]- tck[n-1].pos[0], tck[n].pos[1]- tck[n-1].pos[1], tck[n].pos[2]- tck[n-1].pos[2]);
                    tck[n].set_colour (dir.normalise(), A);
                  }
                  else tck[0].set_colour (default_colour, A);
                }
                else {
                  for (gsize n = 0; n < tck.size(); n++) 
                    tck[n].set_colour (colour, A);
                }
              }
              colour_by_dir_previous = colour_by_dir;
              alpha_previous = alpha;
            }
          }

          GLclampf get_alpha () const { return (exp (TRANSPARENCY_EXPONENT * ( alpha - 1.0))); }

          void load (const String& filename);
          bool refresh ()
          {
            struct_stat64 S;
            if (STAT64 (file.c_str(), &S))
              throw Exception ("error accessing tracks file \"" + file + "\": " + Glib::strerror (errno));
            if (mtime != S.st_mtime) { 
              load (file);
              return (true);
            }
            return (false);
          }

          void draw ();

          guint count () const
          {
            guint n = 0;
            for (std::list<Track>::const_iterator i = tracks.begin(); i != tracks.end(); ++i) n += i->size();
            return (n);
          }

          const Track::Point* root () { return ((const Track::Point*) alloc.root()); }

          void add (std::vector<guint>& vertices, float min_dist, float max_dist)
          {
            for (std::list<Track>::iterator i = tracks.begin(); i != tracks.end(); ++i) {
              for (gsize n = 0; n < i->size(); n++) {
                Track::Point& P ((*i)[n]);
                float Z = Track::Point::normal.dot (Point (P.pos));
                if (Z > min_dist && Z < max_dist) 
                  vertices.push_back (P.index());
              }
            }
          }

        protected:
          Allocator alloc;
      };






      inline TrackListItem::Track::Track (Allocator& alloc, guint num) : num_p (num)       { data = alloc (num_p); }

    }
  }
}

#endif

