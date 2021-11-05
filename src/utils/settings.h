// settings.h --- Base class for settings objects
// Copyright (C) 2006-2011 Wouter Boomsma
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

#ifndef SETTINGS_H
#define SETTINGS_H

//! Settings base class
class Settings {
public:
     //! Destructor
     virtual ~Settings(){}

     //! virtual output method - derived classes can override this to define output functionality
     virtual void output(std::ostream &o) const {};

     // Overload output operator - calls virtual output method
     friend std::ostream &operator<<(std::ostream &o, const Settings &settings) {
          settings.output(o);
          return o;
     } 
};

//! Base class for Settings call back objects.
//! Useful for specifying callbacks to member functions in setting objects
class SettingsCallbackBase {
public:

     //! Make callback call
     virtual void execute() =0;
};

//! Settings call back: Wrapper for a settings object and a member function pointer.
//! Useful for specifying callbacks to member functions in setting objects
template <typename SETTINGS_TYPE, typename CALLBACK_TYPE>
class SettingsCallback: public SettingsCallbackBase {
public:

     //! Constructor
     //! \param settings Settings object
     //! \param callback function
     SettingsCallback(SETTINGS_TYPE *settings,
                      const CALLBACK_TYPE callback)
          : settings(settings),
            callback(callback) {}

     virtual ~SettingsCallback() {}

     //! Make callback call
     virtual void execute () {
          ((*settings).*(callback))();          
     }

private:
     //! Settings object
     SETTINGS_TYPE *settings;

     //! Call back function pointer
     CALLBACK_TYPE callback;

};


#endif
