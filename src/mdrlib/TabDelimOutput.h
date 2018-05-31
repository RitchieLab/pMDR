// TabDelimOutput.h

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

#ifndef __TABDELIMOUTPUT_H__
#define __TABDELIMOUTPUT_H__

#include "Model.h"
#include "Dataset.h"
#include "Stat.h"
#include <iostream>

using namespace std;

namespace mdr{

class TabDelimOutput{

  public:
  
    void add_header(ostream& os, int max_locus, string fitness_used);
    
    void output_model(ostream& os, Model& model, int curr_cv, Dataset& set, Stat* calculator);
        
  private:

    unsigned int maximum_model_size;
};


}

#endif

