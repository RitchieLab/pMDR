//AllModelOutput.cpp

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

#include "AllModelOutput.h"

#include <iomanip>
#include <sstream>

namespace mdr{


///
/// constructor
///
AllModelOutput::AllModelOutput(){
	total_cv = 0;
}

///
/// Outputs header 
/// @param os ostream
/// @param max_locus
/// @param fitness_used 
///

void AllModelOutput::add_header(ostream& os, int cv){
  total_cv = cv;  
  os << "All models\n";
  os << "Model";
  for(int i=1; i<=total_cv; i++){
  	os << "\tCV " << i << " Train\tCV " << i << " Test";
  }
  os << "\n";
}

///
/// Outputs model in list of models that were kept from run
/// @param os
/// @param model Model to output
/// @param curr_cv current crossvalidation interval for this model
/// @param set Dataset
///
void AllModelOutput::output_model(ostream& os, vector<Model>& models,std::vector<unsigned int>& combination){
  
  for(vector<unsigned int>::iterator iter=combination.begin(); iter != combination.end();
  	++iter){
    os << *iter << " ";
  }
  
  for(int curr_cv=0; curr_cv < total_cv; curr_cv++){
  	os << "\t" << models[curr_cv].get_combined_fitness() << "\t" << (100-models[curr_cv].testing.balanced_error);
  }
  os << "\n";
  
}


}

