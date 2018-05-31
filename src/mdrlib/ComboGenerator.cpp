// ComboGenerator.cpp

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


#include<iostream>
using namespace std;
#include "ComboGenerator.h"

namespace mdr{

///
/// Constructor initializes parameters
///
ComboGenerator::ComboGenerator(){
  initialize();
}

///
/// Initializes variables for generator
/// @return none
///
void ComboGenerator::initialize(){
  ComboInterval = 1000;
  AlreadyStarted = false;
  NumLoci = 0;
  ComboStart=0;
  ComboEnd=0;
}


///
///Sets the interval for generating combinations.  Used so that
///don't excessively use memory 
//@param cmbInterval interval to use in producing combinations
//@return none
///
void ComboGenerator::SetComboInterval(int cmbInterval){
  ComboInterval = cmbInterval;
}


///
///Sets the start combination number and ending combination
///@param  combStart minimum size of combination
///@param  combEnd maximum size of combination
///@return  none
///
void ComboGenerator::ComboEnds(int combStart, int combEnd){
  ComboStart = combStart;
  ComboEnd = combEnd;
//  initialize();
}

///
///Sets total number of loci for generator
///@param nLoci number of loci
///@return none
///
void ComboGenerator::SetLoci(int nLoci){
  NumLoci = nLoci;
}

///
/// Empty default implementation
///
void ComboGenerator::PassInitConfig(std::map<std::string, std::string> config){
}

}

