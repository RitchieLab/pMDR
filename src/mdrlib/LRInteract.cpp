//LRInteract.cpp

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

#include "LRInteract.h"

namespace mdr{

///
/// Constructor
///
LRInteract::LRInteract(){
  initialize();
}


///
/// Establishes object 
/// @param
///
void LRInteract::initialize(){
}


///
/// Passed model is used to create single locus models
/// These models are used in assigning new values to indivduals
/// in the set.  An individual gets a 1 if it is high risk and
/// a zero if low risk.  Missing data is set to 2.
/// @param set Dataset
/// @param model Model
/// @return logistical regression statistic
///
void LRInteract::calculate_lr_stat(Dataset& set, Model& model){

  // for each locus, calculate model
  Model single_model;
  
  // initialize set for recoding data
  vector<vector<unsigned int> > recoded_genos(model.combination.size(), vector<unsigned int>(0,0));
  
  float full_lr, reduced_lr;
  
  for(unsigned int curr_loc=0; curr_loc < model.combination.size(); curr_loc++){
  
    // get model calculated on all data
    single_model = analyze_single(set, model.combination[curr_loc]);
    
    // reassign genotype values to new values
    recode_data(single_model, recoded_genos[curr_loc], model.combination[curr_loc], set);
  }
  
  vector<unsigned int> status;
  unsigned int num_inds = set.get_num_inds();
  // individual status is first index of each row of the 2-d data array
  for(unsigned int curr_ind=0; curr_ind < num_inds; curr_ind++){
    status.push_back(set.data[curr_ind][0]);
  }
  
  // pass vector to LR routine for full interaction
  full_lr = lr_calc.run_lr(recoded_genos, status, set.any_missing_data(), 1, true); 
  // pass vector to LR routine removing last interaction
  reduced_lr = lr_calc.run_lr(recoded_genos, status, set.any_missing_data(), 1, false);

  // calculate difference and store in model
  model.set_interact_llr(reduced_lr - full_lr);
}


///
/// Recodes data in order of the dataset to match the single locus
/// high or low cell 
/// @param mod ModelP1
/// @param data vector to store new values
/// @param curr_loc index of current locus
/// @param set Dataset
///
void LRInteract::recode_data(Model& mod, vector<unsigned int>& data, unsigned int curr_loc, Dataset& set){
  
  unsigned int num_inds = set.get_num_inds();
  
  vector<unsigned int> recode(4,0);
  recode[3] = 2; // missing value is recoded as a 2
  // others will be either 0 or 1
  // will be 1 when it is high risk cell in model
  // will be 0 when it is low risk cell in model
  for(unsigned int cell=0; cell < 3; cell++){
    if(mod.training.affected[cell] / float(mod.training.unaffected[cell]) >= mod.get_threshold())
      recode[cell] = 1;
    else
      recode[cell] = 0;
  }
  
  // recode data as 0, 1, or 2
  for(unsigned int curr_ind=0; curr_ind < num_inds; curr_ind++){
      data.push_back(recode[(unsigned int)(set.data[curr_ind][curr_loc])]);
  }
  
}


///
/// Returns model for locus passed
/// @param set Dataset
/// @param curr_loc locus to test in model
///
Model LRInteract::analyze_single(Dataset& set, unsigned int curr_loc){
  Model new_model;
  new_model.combination.assign(1, curr_loc); 
  mdr_calc.set_cv(set.get_crossval_interval());
  // calculate mdr on model using entire dataset 
  mdr_calc.full_model(new_model, set);

  return new_model;
}

}
