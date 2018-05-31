//MatrixOutput.cpp

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

#include "MatrixOutput.h"

#include <iomanip>
#include <sstream>

namespace mdr{


///
/// constructor
///
MatrixOutput::MatrixOutput(){
 last_cv = -1;
}

///
/// Outputs header 
/// @param os ostream
/// @param max_locus
/// @param fitness_used 
///
void MatrixOutput::add_header(ostream& os, int max_locus){
  os << "Confustion matrices " << endl;
}


///
/// Outputs crossvalidation header
/// @param os ostream
/// @param cv crossvalidation interval to output
///
void MatrixOutput::output_crossval_header(ostream& os, int cv){
  os << endl << "-----------------------------------------------------" << endl;  
  os << "Confusion Matrix - Crossvalidation Interval " << cv << endl;   
	os << "-----------------------------------------------------" << endl;  
	os << right << setw(27) << "Training" << right << setw(21) << "Testing"<< endl;
	os << left<< setw(14) <<  "Model" << setw(5) << "T+" << setw(5) << "F+"
	  << setw(5) << "T-" << setw(5) << "F-" << "| " << setw(5) << "T+" << setw(5) << "F+"
	  << setw(5) << "T-" << setw(5) << "F-" << endl;
}



///
/// Outputs model in list of models that were kept from run
/// @param os
/// @param model Model to output
/// @param curr_cv current crossvalidation interval for this model
/// @param set Dataset
///
void MatrixOutput::output_model(ostream& os, Model& model, int curr_cv, Dataset& set){
  vector<unsigned int> combo = model.combination;
  
  if(curr_cv != last_cv){
    last_cv = curr_cv;
    output_crossval_header(os, curr_cv+1);
  }
  
  stringstream ss;
	for(unsigned int i = 0; i < combo.size(); i++){
	  ss << set.snp_names[combo[i]] << " ";
	}
  os << left << setw(14) << ss.str();        
  os << setw(5) << model.training.classhigh;
  os << setw(5) << model.training.misclasshigh;
  os << setw(5) << model.training.classlow;
  os << setw(5) << model.training.misclasslow;
  os << "| ";
  os << setw(5) << model.testing.classhigh;
  os << setw(5) << model.testing.misclasshigh;
  os << setw(5) << model.testing.classlow;
  os << setw(5) << model.testing.misclasslow;       
  os << endl;  

}




}

