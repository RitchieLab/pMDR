//BioComboGenerator.h

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

#ifndef __BIOCOMBOGENERATOR_H__
#define __BIOCOMBOGENERATOR_H__

#include "ComboGenerator.h"
#include "BioFileAccess.h"

namespace mdr{

class BioComboGenerator:public ComboGenerator{

  public:
    BioComboGenerator();
    
    ~BioComboGenerator();
    
    void PassInitConfig(map<string, string> config);
    
    bool GenerateCombinations();  
    bool GenerateCombinations(int new_interval);
    
    bool AdvanceParameters();
    
    double calc_combinations();
    
    int get_count_generated();
    
    void initialize_state();

    /// Reads checkpoint information into generator
    virtual void ReadCheckpoint(ifstream& instream);
    
    /// Writes checkpoint information from generator
    virtual void WriteCheckpoint(ofstream& outstream);  

  protected:
    void FreeMemory();
    void initialize();
    BioFileAccess* _fileReader;
    double _current_model;
    unsigned int _current_index, _total_models_in_file;

};

}

#endif
