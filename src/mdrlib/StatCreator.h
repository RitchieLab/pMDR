// StatCreator.h

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

#ifndef __STATCREATOR_H__
#define __STATCREATOR_H__

using namespace std;

#include "Stat.h"
#include "fitness.h"

namespace mdr{
///
/// Abstract base class that calculates stat to be used in ranking models
///
class StatCreator{

  public:
  
    Stat* create_stat(fitnessTypes type);
    
  private:

};



}

#endif
