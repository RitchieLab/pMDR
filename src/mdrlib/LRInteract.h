// LRInteract

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

#ifndef __LRINTERACT_H__
#define __LRINTERACT_H__

#include "FlatIndex.h"
#include "LogistReg.h"
#include "MDR.h"

namespace mdr{

///
/// Takes an mdr model and creates a two dimensional data
/// vector where individuals are coded at each locus.  Each locus
/// is converted into a high or low.
///


/// Runs logistic regression
class LRInteract{
  
  public:
    LRInteract();
    
    void calculate_lr_stat(Dataset& set, Model& model);
    
  private:
    
    void recode_data(Model& mod, vector<unsigned int>& data, unsigned int curr_loc, Dataset& set);
    
    Model analyze_single(Dataset& set, unsigned int curr_loc);
    
    void initialize();
    
    MDR mdr_calc;
    LogistReg lr_calc;

};

}
#endif
