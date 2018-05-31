//KnuthComboGenerator.h

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



#ifndef __KNUTHCOMBOGENERATOR_H__
#define __KNUTHCOMBOGENERATOR_H__

#include "ComboGenerator.h"

namespace mdr{

///
/// Generates combinations of numbers <br>
/// For example, will generate all 2 number combinations from
/// 1 to 100: <br> <i> 1-2, 1-3, 1-4, 1-5, 1-6, 1-7, etc.</i> <br>
/// Adapted from Will Bush's sMDR code
///

/// Generates numerical combinations
class KnuthComboGenerator:public ComboGenerator{
  
  public:
    KnuthComboGenerator();
    bool GenerateCombinations();  
    bool GenerateCombinations(int new_interval);
    
    bool AdvanceParameters();
    
    double calc_combinations();
    
    void initialize_state();
    
    int get_count_generated(){return counter+1;}
    
    int param_j(){return j;}
    void param_j(int new_j){j=new_j;}
    
    int param_x(){return x;}
    void param_x(int new_x){x=new_x;}
    
    int param_kdec(){return kdec;}
    void param_kdec(int new_kdec){kdec=new_kdec;}
    
    int param_c(int index){return c[index];}
    void set_param_c(int index, int value){c[index] = value;}
   
    /// Reads checkpoint information into generator
    virtual void ReadCheckpoint(ifstream& instream);
    
    /// Writes checkpoint information from generator
    virtual void WriteCheckpoint(ofstream& outstream);    
  
  private:
    void initialize(); 
    int j, x, kdec, counter;
    int * c;
    
};

}

#endif
