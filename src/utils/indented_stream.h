// utils.h --- A stream which can maintain a indentation state
// Copyright (C) 2006-2008 Wouter Boomsma
//
// This file is part of Phaistos
//
// Phaistos is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.

// Phaistos is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.

// You should have received a copy of the GNU General Public License
// along with Phaistos.  If not, see <http://www.gnu.org/licenses/>.

#ifndef INDENTED_STREAM_H
#define INDENTED_STREAM_H

#include <boost/iostreams/filtering_stream.hpp>

namespace phaistos {

//! Indented stream filter. Maintains fixed indentation state during output
class IndentedStreamFilter: public boost::iostreams::output_filter {
private:

     //! Indentation level
     int indentation;

     //! Whether current position is at beginning of line
     bool start_of_line;
     
public:

     //! Constructor
     //!
     //! \param indentation Indentation level
     explicit IndentedStreamFilter(int indentation)
          : indentation(indentation), 
            start_of_line(true) { 
     }

     //! Override put function
     template<typename Sink>
     bool put(Sink& dest, int c) {
          if (c == '\n') {
              start_of_line = true;
          } else {
               if (start_of_line) {
                    for(int n=0; n<indentation; ++n) {
                         boost::iostreams::put(dest, ' ');
                    }
                    start_of_line = false;
               }
          }
          return boost::iostreams::put(dest, c);
     }
    
     //! Override close function
     template<typename Sink>
     void close(Sink&) {
          indentation = 0;
          start_of_line = true;
     }
};

}

#endif

