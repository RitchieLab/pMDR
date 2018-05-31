// StatCreator.cpp

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

#include "StatCreator.h"
#include "OddsRatioStat.h"
#include "BalancedAccStat.h"
#include "ErrorStat.h"
#include <iostream>
using namespace std;

namespace mdr{

///
/// Dynamically creates new Stat matching fitnessType
/// @param type fitnessType
/// @return Stat*
///
Stat* StatCreator::create_stat(fitnessTypes type){

  Stat* calculator;

  switch(type){
    case BalancedAccuracy:
      calculator = new BalancedAccStat;
      break;
    case Error:
      calculator = new ErrorStat;
      break;
    case OddsRatio:
      calculator = new OddsRatioStat;
      break;
    default:
      calculator = new BalancedAccStat;
  };
  
  return calculator;
  
}


}
