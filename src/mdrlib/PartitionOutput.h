// PartitionOutput.h

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

#ifndef __PARTITIONOUTPUT_H__
#define __PARTITIONOUTPUT_H__

#include "Model.h"
#include "Dataset.h"
#include <iostream>
#include <iomanip>

using namespace std;

namespace mdr{

class PartitionOutput{

  public:
  
    PartitionOutput();
  
    void add_header(ostream& os, int max_locus);
    
    void output_model(ostream& os, Model& model, int curr_cv, Dataset& set);
    
    void output_crossval_header(ostream& os, int curr_cv);
    
    void output_best(ostream& os, vector<Model>& best_models, Dataset& set);
    
  private:

    int maximum_model_size, last_cv;
    vector<vector<int> > includedIndexes;
    
};



}

#endif

