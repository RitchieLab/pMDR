//OutputWriter.h

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

// outputs models to all the different files 
// as requested

#ifndef __OUTPUTWRITER_H__
#define __OUTPUTWRITER_H__

#include <fstream>
#include "Model.h"
#include "FormattedOutput.h"
#include "TabDelimOutput.h"
#include "RulesOutput.h"
#include "PartitionOutput.h"
#include "MatrixOutput.h"
#include "AllModelOutput.h"
#include "MDRExcept.h"

using namespace std;

namespace mdr{


class OutputWriter{

  public:
  
    OutputWriter();
    ~OutputWriter();
    
    void initialize();
    
    void open_form(string formfilename, int max_locus, string fitness_used, bool any_missing);
    void open_tab(string tabfilename, int maxlocus, string fitness_used);
    void open_rules(string rulesfilename, int maxlocus);
    void open_part(string partfilename, int maxlocus);
    void open_matrix(string matrixfilename, int maxlocus);
    void open_all_out(string alloutfilename, int cv);
    
    void output_all(vector<Model>& models,std::vector<unsigned int>& combination);
    
    void output_model(int curr_cv, Model& model, int total_cv, Dataset& set, Stat* calculator);
    
    void output_best(vector<Model>& best_models, int total_cv, Dataset& set, Stat* calculator);
    
    void output_single_best(Model& best, Dataset& set);
    
    void output_pvalues(vector<Model>& best_models, Dataset& set, int p_tests);
    
    void output_pvalues(vector<vector<Model> >& best_models, Dataset& set, int p_tests);
    
    void output_lr_pvalues(vector<Model>& best_models, Dataset& set, int p_tests);
    
    void output_all_models();
    
    void close_files();
        
  private:

    ofstream form_os, tab_os, rules_os, part_os, matrix_os, all_os;

    FormattedOutput* form_out;
    TabDelimOutput* tab_out;
    RulesOutput* rules_out;
    PartitionOutput* part_out;
    MatrixOutput* matrix_out;
    AllModelOutput* all_out;
    
    
};

}

#endif

