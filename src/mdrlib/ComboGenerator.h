//ComboGenerator.h

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



#ifndef __COMBOGENERATOR_H__
#define __COMBOGENERATOR_H__

#include <vector>
#include "MDRExcept.h"
#include "Dataset.h"
#include <map>

namespace mdr{

///
/// Generates combinations of numbers <br>
/// For example, will generate all 2 number combinations from
/// 1 to 100: <br> <i> 1-2, 1-3, 1-4, 1-5, 1-6, 1-7, etc.</i> <br>
/// Adapted from Will Bush's sMDR code
///

/// Generates numerical combinations
class ComboGenerator{
  
  public:
    ComboGenerator();
    virtual ~ComboGenerator(){}
    
    virtual bool GenerateCombinations()=0;  
    virtual bool GenerateCombinations(int new_interval)=0;
    
    virtual bool AdvanceParameters()=0;
    virtual void initialize_state()=0;   
    
    virtual double calc_combinations()=0;
    
    /// Pass a map with parameters needed by a generator (filename for example)
    virtual void PassInitConfig(std::map<std::string, std::string> config);
    
    void SetComboInterval(int combInterval);
    int getComboInterval(){return ComboInterval;}
    void ComboEnds(int combStart, int combEnd);
    void SetLoci(int nLoci);
    int get_combo_max(){return ComboEnd;}
    int get_combo_min(){return ComboStart;}
    std::vector < std::vector<unsigned int> > ComboList;
        
    virtual int get_count_generated(){return counter+1;}
    
    bool param_AlreadyStarted(){return AlreadyStarted;}
    void param_AlreadyStarted(bool new_AlreadyStarted){AlreadyStarted=new_AlreadyStarted;}
    
    void set_dataset(Dataset* set){dataset=set;}
    
    /// Reads checkpoint information into generator
    virtual void ReadCheckpoint(ifstream& instream)=0;
    
    /// Writes checkpoint information from generator
    virtual void WriteCheckpoint(ofstream& outstream)=0;
  
  protected:
    void initialize(); 
    int ComboInterval, ComboStart, ComboEnd, NumLoci;
    int counter;
    bool AlreadyStarted;
    Dataset* dataset;
    
};

}

#endif

