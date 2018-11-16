//PartitionOutput.cpp

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

#include "PartitionOutput.h"
#include <sstream>

namespace mdr{


PartitionOutput::PartitionOutput(){
  maximum_model_size = 5;
  last_cv = -1;
}


///
/// Adds header 
/// @param os ostream for writing
/// @param max_locus maximum model size
/// @return 
///
void PartitionOutput::add_header(ostream& os, int max_locus){
  
  maximum_model_size = max_locus;
  
  os << "************************************************************" << endl;
  
}


///
/// Outputs crossvalidation header
/// @param os ostream
/// @param curr_cv Crossvalidation interval to output
/// @return
/// 
void PartitionOutput::output_crossval_header(ostream& os, int curr_cv){
  os << "Model Partitions - Crossvalidation Interval " << curr_cv << endl;
  os << "------------------------------------------------------------" << endl;
}


///
/// Adds model partitions to the output
/// @param os
/// @param model
/// @param curr_cv
/// @param set Dataset
/// @return
///
void PartitionOutput::output_model(ostream& os, Model& model, int curr_cv,
  Dataset& set){

  if(includedIndexes.size() ==0)
    includedIndexes = includedIndexes = set.converter().get_included_indexes();

  if(curr_cv != last_cv){
    last_cv = curr_cv;
    output_crossval_header(os, curr_cv);
  }
    
  vector<int> indexes;
  unsigned int i;
  
  os << "Model [ ";
  unsigned int model_size = model.combination.size();
  for(i=0; i<model_size; i++)
    os << set.snp_names[model.combination[i]] << " ";
  os << "]" << endl;
  
  os << "Genotype" << setw(5) << "A" << setw(7) << "U" << setw(7) << "TOT" << setw(18) 
     << "Ratio/Threshold" << setw(15) << "Risk" << endl;
  os << "------------------------------------------------------------" << endl;
  
  indexes = includedIndexes[model_size];
  
  vector<unsigned int> combo = model.combination;
  
  for(unsigned int currIndex=0; currIndex<indexes.size(); currIndex++){
  
    vector<int> genotype = set.converter().decode_index(indexes[currIndex], combo.size());
    stringstream ss;
    for(i=0; i<genotype.size(); i++)
      ss << genotype[i] << " ";
      
    // output the stats
    os << left << setw(6) << ss.str() << right << setprecision(5) << setw(7)
      << model.training.affected[indexes[currIndex]]
      << setw(7) << model.training.unaffected[indexes[currIndex]]
      << setw(7) << model.training.unaffected[indexes[currIndex]] + 
        model.training.affected[indexes[currIndex]]
      << setw(9) << setprecision(3) << model.training.affected[indexes[currIndex]] /
        float(model.training.unaffected[indexes[currIndex]])
      << " / "
      << model.get_threshold()
      << setw(15);
      
    switch(model.cell_status[indexes[currIndex]]){
      case 0:      // IF LOW RISK
        os << "Low Risk"  << endl;
        break;
      case 1:      // IF HIGH RISK
        os << "High Risk" << endl;
        break;
      default:
        os << model.cell_status[indexes[currIndex]] << endl;
        break;  
    };
  
  }
  os << endl;
  
}
  
///
/// for combined result just generate a model that has all the stats
/// with no crossvalidation and then run through the function above
/// @param os ostream to write to
/// @param best_models
/// @param set Datset
///
void PartitionOutput::output_best(ostream& os, vector<Model>& best_models,
  Dataset& set){
  
  last_cv = 0;
  
  os << "************************************************************" << endl;
  os << "------- Partitions For Best Models Across All Data ---------" << endl;
  
  
  for(unsigned int i=0; i<best_models.size(); i++){
    if(best_models[i].combination.size() > 0)
      output_model(os, best_models[i], 0, set);
  }
  
  
}


}

