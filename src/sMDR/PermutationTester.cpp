//PermutationTester.cpp

//     Copyright 2008 Marylyn Ritchie

//     This file is part of sMDR.
// 
//     sMDR is free software: you can redistribute it and/or modify
//     it under the terms of the GNU General Public License as published by
//     the Free Software Foundation, either version 3 of the License, or
//     (at your option) any later version.
// 
//     sMDR is distributed in the hope that it will be useful,
//     but WITHOUT ANY WARRANTY; without even the implied warranty of
//     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//     GNU General Public License for more details.
// 
//     You should have received a copy of the GNU General Public License
//     along with sMDR.  If not, see <http://www.gnu.org/licenses/>.

#include "PermutationTester.h"
#include <algorithm>

#include <iostream>

///
/// Runs permutations and creates tree from best results
/// @param number of permutations
/// @param set Dataset
/// @param config Config
/// @param log_out LogOutput
///
void PermutationTester::run_permutations(unsigned int num_permutations, Dataset& set, Config& config,
  LogOutput& log_out){
  
// cout << "in run_permutations" << endl;  
  
  Config perm_config = config;
  perm_config.models_to_keep(1);  
  vector<string> none_always_included;
//  perm_config.loci_always_included(none_always_included);
  perm_config.force_loci_included(none_always_included);
  perm_config.check_point_file("");
// cout << "perm_config.check_point_file()=" << perm_config.check_point_file() << endl;
  perm_config.check_interval(0);
  int total_cv = perm_config.num_crossval();
  total_perms = num_permutations;
  
  LRInteract lr_calc;
 
  for(unsigned int curr_perm = 1; curr_perm <= num_permutations; curr_perm++){
// cout << "curr_perm=" << curr_perm << endl;  
//    srand(curr_perm);

    set.permute_set(curr_perm);
  
    // create analysis object
    sAnalysis permuted_analysis;
    permuted_analysis.train_models(set, perm_config, log_out, false);

    // get best models for each order    
    vector<Model> best_models =  permuted_analysis.get_best_models(set, perm_config);
    
    // get single best model   
    Model single_best = permuted_analysis.get_single_best(best_models, perm_config);
    
    float fitness;
    
    // add this best one to the tree
    if(total_cv == 1)
      fitness = single_best.training.balanced_error;
    else
      fitness = single_best.get_balpredavg();

    if(config.regress_test()){
// cout << "PermutationTester calculate lr stat" << endl;
      lr_calc.calculate_lr_stat(set, single_best);
      // store lr result in vector
      lr_results.push_back(single_best.get_interact_llr());
    }

// cout << "fitness added=" << fitness << endl;
    permuted_models.Add(fitness, single_best.combination);
    log_out.add_perm_model(single_best, curr_perm, set, permuted_analysis.get_calculator(), 
      config.regress_test());
  
  }
  
  // sort the lr results
  sort(lr_results.begin(), lr_results.end(), sort_ascending());
  
}


///
/// Returns p value corresponding to the value passed
/// @param value to compare distribution against
/// @return p value 
///
float PermutationTester::get_p_value(float value){
  
  int counter=0;
  ModelTreeNode* curr_node;
  
  for(curr_node = permuted_models.GetLast(); curr_node != NULL;
    curr_node = curr_node->GetPrev()){

    if(curr_node->GetKey() > value)
      break;
    counter++;
  }
  return float(counter)/total_perms;
}



///
/// Returns p value corresponding to the value passed for an LR score 
/// @param value to compare distribution against
/// @return p value 
///
float PermutationTester::get_lr_p_value(float value){

  int counter=0;
  for(; counter < total_perms; counter++){
    if(value > lr_results[counter])
      break;
  }
  
//   ModelTreeNode* curr_node;
//   
//   for(curr_node = permuted_models.GetLast(); curr_node != NULL;
//     curr_node = curr_node->GetPrev()){
//     
//     if(curr_node->GetKey() > value)
//       break;
//     counter++;
//   }

  return float(counter)/total_perms;
}



