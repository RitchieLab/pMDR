// OutputWriter.cpp

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

#include "OutputWriter.h"

namespace mdr{

///
/// constructor
///
OutputWriter::OutputWriter(){
  initialize();
}
  
///
/// Initialize
///
void OutputWriter::initialize(){
  form_out = NULL;
  tab_out = NULL;
  rules_out = NULL;
  part_out = NULL;
  matrix_out = NULL;
  all_out = NULL;
}


///
/// Destructor
///
OutputWriter::~OutputWriter(){
  close_files();
}


///
/// Closes files
///
void OutputWriter::close_files(){
  
  if(form_out != NULL){
    form_os.close();
    delete form_out;
  }
  
  if(tab_out != NULL){
    tab_os.close();
    delete tab_out;
  }
 
  if(rules_out != NULL){
    rules_os.close();
    delete rules_out;
  } 

  if(part_out != NULL){
    part_os.close();
    delete part_out;
  } 
 
  if(matrix_out != NULL){
    matrix_os.close();
    delete matrix_out;
  } 
  
  if(all_out != NULL){
  	all_os.close();
  	delete all_out;
  }
 
}


///
/// Opens stream for full model output and
/// constructs class to work on it
/// @param alloutfilename
/// @param min_locus Minimum number of loci that can be in model
/// @param max_locus Maximum number of loci that can be in model
/// @param cv Number of cross-validations used
///
void OutputWriter::open_all_out(string alloutfilename, int cv){

  all_os.open(alloutfilename.c_str(), ios::out);
  if(!all_os.is_open()){
    throw MDRExcept("Cannot open " + alloutfilename + " for writing\n");
  }
  
  all_out = new AllModelOutput;

  all_out->add_header(all_os, cv);

}


///
/// Opens stream for tab delimited output and
/// constructs class to work on it
/// @param tabfilename
/// @param max_locus Maximum number of loci that can be in model
/// @param fitness_used
///
void OutputWriter::open_tab(string tabfilename, int max_locus, string fitness_used){

  tab_os.open(tabfilename.c_str(), ios::out);
  if(!tab_os.is_open()){
    throw MDRExcept("Cannot open " + tabfilename + " for writing\n");
  }
  
  tab_out = new TabDelimOutput;

  tab_out->add_header(tab_os, max_locus, fitness_used);

}


///
/// Opens stream for tab formatted output and
/// constructs class to work on it
/// @param formfilename 
/// @param max_locus Maximum number of loci that can be in model
/// @param any_missing Any missing data in the set
///
void OutputWriter::open_form(string formfilename, int max_locus, string fitness_used,
  bool any_missing){
 
  form_os.open(formfilename.c_str(), ios::out);
  if(!form_os.is_open()){
    throw MDRExcept("Cannot open " + formfilename + " for writing\n");
  }
  
  form_out = new FormattedOutput;
  form_out->missing_data(any_missing);
  form_out->add_header(form_os, max_locus, fitness_used);

}


///
/// Opens stream for rules output and
/// constructs class to work on it
/// @param rulesfilename 
///
void OutputWriter::open_rules(string rulesfilename, int maxlocus){
  
  rules_os.open(rulesfilename.c_str(), ios::out);
  if(!rules_os.is_open()){
    throw MDRExcept("Cannot open " + rulesfilename + " for writing\n");
  }
  
  rules_out = new RulesOutput;
  rules_out->add_header(rules_os, maxlocus);

}


///
/// Opens stream for partition output and
/// constructs class to work on it
/// @param partfilename 
///
void OutputWriter::open_matrix(string matrixfilename, int maxlocus){
  
  matrix_os.open(matrixfilename.c_str(), ios::out);
  if(!matrix_os.is_open()){
    throw MDRExcept("Cannot open " + matrixfilename + " for writing\n");
  }
  
  matrix_out = new MatrixOutput;
  matrix_out->add_header(matrix_os, maxlocus);

}


///
/// Opens stream for matrix output and
/// constructs class to work on it
/// @param matrixfilename 
/// @param maxlocus
///
void OutputWriter::open_part(string partfilename, int maxlocus){
  
  part_os.open(partfilename.c_str(), ios::out);
  if(!part_os.is_open()){
    throw MDRExcept("Cannot open " + partfilename + " for writing\n");
  }
  
  part_out = new PartitionOutput;
  part_out->add_header(part_os, maxlocus);

}



///
/// Outputs model to all the types specified
/// @param curr_cv current crossvalidation interval being displayed
/// @param model Model to output
/// @param total_cv total number of crossvalidations performed in analysis
/// @param set
/// @param calculator Stat for proper display of values in model
///
void OutputWriter::output_model(int curr_cv, Model& model, int total_cv, Dataset& set, 
  Stat* calculator){
  
  if(form_out != NULL)
    form_out->output_model(form_os,model, curr_cv, total_cv, set, calculator);
  
  if(rules_out != NULL)
    rules_out->output_model(rules_os,model, curr_cv, set);
  
  if(tab_out != NULL)
    tab_out->output_model(tab_os,model, curr_cv, set, calculator);
  
  if(part_out != NULL)
    part_out->output_model(part_os,model, curr_cv, set);
    
  if(matrix_out != NULL)
    matrix_out->output_model(matrix_os, model, curr_cv, set);


}


///
/// Outputs models 
/// @param modelMap
///
void OutputWriter::output_all(vector<Model>& models,std::vector<unsigned int>& combination){
	all_out->output_model(all_os, models, combination);
}

///
/// Outputs best model to the formatted output
/// @param best Model
/// @param set Dataset
///
void OutputWriter::output_single_best(Model& best, Dataset& set){
  
  if(form_out != NULL)
    form_out->output_single_best(form_os, best, set);

}

///
/// Output p values to formatted output
/// @param best_models
/// @param set Dataset
/// 
void OutputWriter::output_pvalues(vector<Model>& best_models, Dataset& set, int p_tests){
  if(form_out != NULL)
    form_out->output_p_values(form_os, best_models, set, p_tests);
}

///
/// Output LR p values to formatted output
/// @param best_models
/// @param set Dataset
/// 
void OutputWriter::output_lr_pvalues(vector<Model>& best_models, Dataset& set, int p_tests){
  if(form_out != NULL)
    form_out->output_lr_p_values(form_os, best_models, set, p_tests);
}

///
/// Outputs best model information
/// @param best_models vector<Model> best model at each size
/// @param total_cv number of CV performed in analysis
/// @param calculator Stat*
///
void OutputWriter::output_best(vector<Model>& best_models, int total_cv, Dataset& set,
  Stat* calculator){
  if(form_out != NULL)
    form_out->output_best(form_os, best_models, total_cv, set, calculator);
    
  if(total_cv > 1)
    if(part_out != NULL)
      part_out->output_best(part_os, best_models, set);
}

}

