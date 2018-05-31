#ifndef __MDRExcept_H__
#define __MDRExcept_H__

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

///
/// Thrown by functions that need to indicate errors
/// Error message is stored within the object
///

#include <string>
#include <iostream>
using namespace std;

namespace mdr{

class MDRExcept{
        public:
                MDRExcept();
                MDRExcept(string message);
                string get_error();
                friend ostream & operator << (ostream & os, const MDRExcept & ke);
        private:
                string error;
};

}

#endif

