//TabDelimOutput.cpp

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

#include "TabDelimOutput.h"

namespace mdr{

///
/// Adds header to the tab output file
/// @param os ostream for writing
/// @param max_locus int with maximum model size
/// @param fitness_used
/// @return
///
void TabDelimOutput::add_header(ostream& os, int max_locus, string fitness_used){
  
  maximum_model_size = max_locus;

  os << "CV\tModel Size\t";
	for(int i=1; i<=max_locus; i++){
		os << "Locus " << i << "\tChrom\tPos\t";
	}
	os << "Training " << fitness_used << "\t";
	os << "Testing " << fitness_used << "\t";
	
	os << "Prediction Error";
	
	os << "\tNum Rules";
	
	os << "\tTraining True Pos\tTraining False Pos\tTraining True Neg\tTraining False Neg";
	
	os << "\tTesting True Pos\tTesting False Pos\tTesting True Neg\tTesting False Neg";
	
	os << "\n";

}

    
///
/// Outputs model to tab delimited file
///
void TabDelimOutput::output_model(ostream& os, Model& model, int curr_cv, Dataset& set,
  Stat* calculator){ 
  os << curr_cv + 1 << "\t";
  os << model.combination.size() << "\t";
  unsigned int loc;
  for(loc=0; loc < model.combination.size(); loc++){
    os << set.snp_names[model.combination[loc]] << "\t" << set.chromosomes[model.combination[loc]] <<
      "\t" << set.positions[model.combination[loc]] << "\t";
  }
  for(;loc < maximum_model_size; loc++)
    os << "\t\t\t";
  
  os << calculator->convert_to_display(model.get_train_stats().balanced_error) << "\t";
  os << calculator->convert_to_display(model.get_test_stats().balanced_error) << "\t";
  os << model.get_test_stats().error << "\t";

  int num_rules = 0;
  for(unsigned int i=0; i<model.cell_status.size(); i++){
    if(model.cell_status[i] != -1)
      num_rules++;
  }

  os << num_rules << "\t";
  os << model.get_train_stats().classhigh << "\t";
  os << model.get_train_stats().misclasshigh << "\t";
  os << model.get_train_stats().classlow << "\t";
  os << model.get_train_stats().misclasslow << "\t";
  os << model.get_test_stats().classhigh << "\t";
  os << model.get_test_stats().misclasshigh << "\t";
  os << model.get_test_stats().classlow << "\t";
  os << model.get_test_stats().misclasslow << "\t";
  os << "\n";

}


}

