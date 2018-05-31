#include "MDRExcept.h"

//     Copyright 2008 Marylyn Ritchie

//     This file is part of MDR.
// 
//     MDR is free software: you can redistribute it and/or modify
//     it under the terms of the GNU General Public License as published by
//     the Free Software Foundation, either version 3 of the License, or
//     (at your option) any later version.
// 
//     MDR is distributed in the hope that it will be useful,
//     but WITHOUT ANY WARRANTY; without even the implied warranty of
//     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//     GNU General Public License for more details.
// 
//     You should have received a copy of the GNU General Public License
//     along with MDR.  If not, see <http://www.gnu.org/licenses/>.

namespace mdr{

// constructor
// Arg: none
// Ret: none
MDRExcept::MDRExcept(){
        error = "No message set by code.";
}

// constructor
// Arg: Error message
// Ret: none
MDRExcept::MDRExcept(string message){
        error = message;
}

// returns error message
// Arg: none
// Ret: error message in object
string MDRExcept::get_error(){
        return error;
}

// overloads << operator to print error message
// Arg:  output stream
//       object to print
// Ret:  output stream
ostream & operator << (ostream & os, const MDRExcept & ke){
        os << ke.error;
        return os;
}

}

