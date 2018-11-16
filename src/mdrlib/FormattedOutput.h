//FormattedOutput.h

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

#ifndef __FORMATTEDOUTPUT_H__
#define __FORMATTEDOUTPUT_H__

#include "Model.h"
#include "Dataset.h"
#include "Stat.h"
#include <iostream>


using namespace std;

namespace mdr{

class FormattedOutput{

  public:
    
    FormattedOutput();
  
    void add_header(ostream& os, int max_locus, string fitness_used);
    
    void output_model(ostream& os, Model& model, int curr_cv, int total_cv, Dataset& set,
      Stat* calculator);
    
    void output_crossval_header(ostream& os, int cv, int total_cv);
    
    void missing_data(bool missing);
    
    void output_best(ostream& os, vector<Model>& best_models, int total_cv, Dataset& set,
      Stat* calculator);
    
    void output_single_best(ostream& os, Model& model, Dataset& set);
    
    void output_p_values(ostream& os, vector<Model>& models, Dataset& set, int p_tests);
    
    void output_p_values(ostream& os, vector<vector<Model> >& models, Dataset& set, int p_tests);
    
    void output_lr_p_values(ostream& os, vector<Model>& models, Dataset& set, int p_tests);
    
  private:
    bool data_missing;

    string fitness_name;

    int maximum_model_size, last_cv;
};



}

#endif

