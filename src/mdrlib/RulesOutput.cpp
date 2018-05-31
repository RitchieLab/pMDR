//RulesOutput.cpp

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

#include "RulesOutput.h"

namespace mdr{


///
/// Constructor
///
RulesOutput::RulesOutput(){
  last_cv = -1;
}

///
/// Adds header 
/// @param os ostream for writing
/// @param max_locus maximum model size
/// @return 
///
void RulesOutput::add_header(ostream& os, int max_locus){
  
  maximum_model_size = max_locus;
  
  os << "*****************************************************" << endl;
  
}


///
/// Outputs crossvalidation header
/// @param os ostream
/// @param curr_cv Crossvalidation interval to output
/// @return
/// 
void RulesOutput::output_crossval_header(ostream& os, int curr_cv){

  os << "\n-----------------------------------------------------" << endl;
  os << "Model Partitions - Crossvalidation Interval " << curr_cv << endl;
  os << "-----------------------------------------------------" << endl;
  
}


/// 
/// Output rules for the model
/// @param os ostream
/// @param model Model
/// @param curr_cv Crossvalidation interval for this model
/// @param set Dataset
/// @return
///
void RulesOutput::output_model(ostream& os, Model& model, int curr_cv, Dataset& set){
  
  unsigned int model_size = model.combination.size();

  if(includedIndexes.size() ==0)
    includedIndexes = includedIndexes = set.converter().get_included_indexes();

  vector<int> genotypes, indexes;
  unsigned int i;
  
  os << "Model ";
  for(i=0; i< model_size; i++)
    os << set.snp_names[model.combination[i]] << " ";
  os << endl; 
  indexes = includedIndexes[model_size];
  unsigned int curr_geno;
  
  for(i=0; i<indexes.size(); i++){
    genotypes = set.converter().decode_index(indexes[i], model_size);
    os << "IF " << set.snp_names[model.combination[0]] << " = " << genotypes[0];
    for(curr_geno=1; curr_geno < genotypes.size(); curr_geno++){
      os << " AND " << set.snp_names[model.combination[curr_geno]] << " = " << genotypes[curr_geno];
    }
    
    os << " THEN STATUS = " << model.cell_status[indexes[i]] << endl;
  }
  
  os << endl;
  
}

}

