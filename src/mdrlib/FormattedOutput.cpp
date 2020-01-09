//FormattedOutput.cpp

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

#include "FormattedOutput.h"
#include <iomanip>
#include <sstream>

namespace mdr{


///
/// constructor
///
FormattedOutput::FormattedOutput(){
  last_cv = -1;
  data_missing = false;
}


///
/// Outputs header 
/// @param os ostream
/// @param max_locus
/// @param fitness_used 
///
void FormattedOutput::add_header(ostream& os, int max_locus, string fitness_used){
  fitness_name = fitness_used;
  os << "Fitness used: " << fitness_used << "\n";
}


///
/// Outputs crossvalidation header
/// @param os ostream
/// @param cv crossvalidation interval to output
///
void FormattedOutput::output_crossval_header(ostream& os, int cv, int total_cv){

  os << "\n" << "------------------------Subgroup # " << cv << " of " << total_cv << "-------------------------";
  if(data_missing){
    os << "---------------------";
  }
  os << "\n";

  if(data_missing)
    os << setw(61) << fitness_name << "\n";
  else
    os << setw(50) << fitness_name<< "\n";
      
  os << setw(22) << left << "Loci";

  os << setw(12) << "Training";
    
  if(data_missing){
    os << setw(11) << left << "TrainMiss";
  }

  os << setw(10) << "Testing"; 
  if(data_missing){
    os << setw(9) << left << "TestMiss";
  }
   
  os << setw(11) << "Pred Error";
    
  os << "\n";
      
  os << "------------------------------------------------------------------";
  if(data_missing){
    os << "---------------------";
  }
  os << "\n";  

}


///
/// Sets whether any genotypes are missing from the set
/// @param missing true for missing
/// 
void FormattedOutput::missing_data(bool missing){
  data_missing = missing;
}


///
/// Outputs model in list of models that were kept from run
/// @param os
/// @param model Model to output
/// @param curr_cv current crossvalidation interval for this model
/// @param set Dataset
/// @param calculator Stat* 
///
void FormattedOutput::output_model(ostream& os, Model& model, int curr_cv, int total_cv,
  Dataset& set, Stat* calculator){

  if(curr_cv != last_cv){
    last_cv = curr_cv;
    output_crossval_header(os, curr_cv+1, total_cv);
  }
 
  stringstream ss;
  for(unsigned int i=0; i<model.combination.size(); i++)
    ss << set.snp_names[model.combination[i]] << " ";
  
  os << setiosflags(ios::fixed | ios::showpoint) << setprecision(2);
  os << setw(22) << left << ss.str();

  os << setw(12) << calculator->convert_to_display(model.training.balanced_error);
  if(data_missing)
    os << setw(11) << model.training.totalmissing;
  os << setprecision(2);
  
  os << setw(10) << calculator->convert_to_display(model.testing.balanced_error);
  if(data_missing)
    os << setw(9) << model.testing.totalmissing;
  os << model.testing.error;

  os << "\n" << right << fixed;

}

///
/// Outputs information about best model for each order
/// @param os ostream
/// @param best_models vector containing best models
/// @param total_cv number of CV in dataset
/// @param calculator Stat*
///
void FormattedOutput::output_best(ostream& os, vector<Model>& best_models, int total_cv,
  Dataset& set, Stat* calculator){
  
  os << "\n" << "\n"<< "--------------- Best Models for Each Order -------------";
  os << "-----------" << "\n";

  os << left <<  setw(18) << left << "loci";


  if(total_cv> 1){
    os << setw(23) << ("Avg Pred " + fitness_name) << setw(13) << "Avg Pred Err";
  }
  else{
    os << setw(28) << fitness_name;
  }

  os << setw(10) << "TotalMiss";
  os << setw(3) << "CV" << "\n";
  os << "--------------------------------------------------------";
  os << "-----------" << "\n";
 
  for(unsigned int i=1; i<best_models.size(); i++){
    if(best_models[i].combination.size() == 0)
      continue;
    stringstream ss;
    
    for(unsigned int j=0; j<best_models[i].combination.size(); j++){
      ss << set.snp_names[best_models[i].combination[j]] << " ";
    }
    
    os << setiosflags(ios::fixed | ios::showpoint) << setprecision(2);
    os << setw(20) << left << ss.str();

    if(total_cv > 1){
      os << setw(21) << calculator->convert_to_display(best_models[i].get_balpredavg());
      os << setprecision(2);
      os << setw(13) << best_models[i].get_predictavg();
    }
    else{
      os << setw(28) << calculator->convert_to_display(best_models[i].training.balanced_error);
    }

    os << setprecision(9);
    os << fixed;

    os << setw(10) << int(best_models[i].get_totalmissdata());
    
    os << setprecision(2);
    os << setw(3) << best_models[i].get_cvc();

    os << "\n" << fixed << right;
  }
 
}


///
/// Output p values
/// @param os ostream 
/// @param models
/// @param set
/// @param p_tests Total number of tests run
///
void FormattedOutput::output_p_values(ostream& os, vector<Model>& models, Dataset& set,
  int p_tests){
  
  float min = 1/float(p_tests);
  
  // precision determined by number of p_tests done
  int precision = 0;
  int result = p_tests;
  
  do{
    result = result/10;  
    precision++;
  }while (result > 1);
  
  os << "Permutation p values:" << "\n";
  os << setprecision(precision);
  for(unsigned int i=0; i<models.size(); i++){
    if(models[i].combination.size() > 0){
      stringstream ss;
      for(unsigned int j=0; j<models[i].combination.size(); j++){
       ss << set.snp_names[models[i].combination[j]] << " ";
      }
      os << setw(20) << left << ss.str();
      if(models[i].get_pvalue() < min)
        os << "< " << min;
      else
        os << models[i].get_pvalue();
      
      os << "\n";
    }
  }

}


///
/// Output p values
/// @param os ostream 
/// @param models
/// @param set
/// @param p_tests Total number of tests run
///
void FormattedOutput::output_p_values(ostream& os, vector<vector<Model> >& models, Dataset& set,
  int p_tests){
  
  float min = 1/float(p_tests);
  
  // precision determined by number of p_tests done
  int precision = 0;
  int result = p_tests;
  
  do{
    result = result/10;  
    precision++;
  }while (result > 1);
  
  os << "Permutation p values:" << "\n";
  os << setprecision(precision);
  for(unsigned int i=0; i<models.size(); i++){
  	for(unsigned j=0; j<models[i].size(); j++){
		if(models[i][j].combination.size() > 0){
		  stringstream ss;
		  for(unsigned int k=0; k<models[i][j].combination.size(); k++){
		   ss << set.snp_names[models[i][j].combination[k]] << " ";
		  }
		  os << setw(20) << left << ss.str();
		  if(models[i][j].get_pvalue() < min)
			os << "< " << min;
		  else
			os << models[i][j].get_pvalue();
	  
		  os << "\n";
		}  		
  	}
  }
  
}

///
/// Output LR p values
/// @param os ostream 
/// @param models
/// @param set
/// @param p_tests Total number of tests run
///
void FormattedOutput::output_lr_p_values(ostream& os, vector<Model>& models, Dataset& set,
  int p_tests){
  
  os << "\n";
  
  float min = 1/float(p_tests);
  
  // precision determined by number of p_tests done
  int precision = 0;
  int result = p_tests;
  
  do{
    result = result/10;  
    precision++;
  }while (result > 1);
  
  os << "LR p values:" << "\n";
  os << setprecision(precision);
  for(unsigned int i=0; i<models.size(); i++){
    if(models[i].combination.size() > 0){
      stringstream ss;
      for(unsigned int j=0; j<models[i].combination.size(); j++){
       ss << set.snp_names[models[i].combination[j]] << " ";
      }
      os << setw(20) << left << ss.str();
      if(models[i].get_lr_pvalue() < min)
        os << "< " << min;
      else
        os << models[i].get_lr_pvalue();
      
      os << "\n";
    }
  }

}




///
/// Outputs best model information
/// @param os ostream to write to
/// @param model Model
/// @param set Dataset
///
void FormattedOutput::output_single_best(ostream& os, Model& model, Dataset& set){

  os << "\n" << "******* Single Best Model **********" << "\n";
  os << " Model [ ";
  for(unsigned int i=0; i<model.combination.size(); i++)
    os << set.snp_names[model.combination[i]] << " ";
  os << "]" << "\n" << "\n";

}

}

