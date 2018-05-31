//ParallelPerm.cpp

//     Copyright 2008 Marylyn Ritchie

//     This file is part of pMDR.
// 
//     pMDR is free software: you can redistribute it and/or modify
//     it under the terms of the GNU General Public License as published by
//     the Free Software Foundation, either version 3 of the License, or
//     (at your option) any later version.
// 
//     pMDR is distributed in the hope that it will be useful,
//     but WITHOUT ANY WARRANTY; without even the implied warranty of
//     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//     GNU General Public License for more details.
// 
//     You should have received a copy of the GNU General Public License
//     along with pMDR.  If not, see <http://www.gnu.org/licenses/>.


#include "ParallelPerm.h"
#include <algorithm>


///
/// Returns p value corresponding to the value passed
/// @param value to compare distribution against
/// @return p value 
///
float ParallelPerm::get_p_value(float value){

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
float ParallelPerm::get_lr_p_value(float value){

  int counter=0;
  for(; counter < total_perms; counter++){
    if(value > lr_results[counter])
      break;
  }

  return float(counter)/total_perms;
}



///
/// Runs permutations and creates tree from best results
/// @param number of permutations
/// @param set Dataset
/// @param config Config
/// @param nproc Total number of processors in run
///
void ParallelPerm::run_master_permutations(int num_permutations, Dataset& set, Config& config,
  int nproc, LogOutput& log_out){

  Config perm_config = config;
  vector<string> none_always_included;

  perm_config.force_loci_included(none_always_included);
  perm_config.models_to_keep(1);
  perm_config.check_point_file("");
  perm_config.check_interval(0);

  int total_cv = perm_config.num_crossval();
  total_perms = num_permutations;
  int curr_perm;
  
  LRInteract lr_calc;
  
  for(curr_perm = 1; curr_perm <= num_permutations; curr_perm++){
    // send seed to the slave nodes
    send_seed(curr_perm);
    set.permute_set(curr_perm);

    // create analysis object
    MasterAnalysis permuted_analysis;
    permuted_analysis.run_master_analysis(set, perm_config, nproc, log_out, false);

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
      lr_calc.calculate_lr_stat(set, single_best);
      // store lr result in vector
      lr_results.push_back(single_best.get_interact_llr());
    }

    permuted_models.Add(fitness, single_best.combination);
    log_out.add_perm_model(single_best, curr_perm, set,  permuted_analysis.get_calculator(),
      config.regress_test());

  }
  
  // send seed which exceeds number of permutations to perform
  send_seed(curr_perm);
  
  // sort the lr results
  sort(lr_results.begin(), lr_results.end(), sort_ascending());

}


///
/// Run permutations for slave
/// @param number of permutations
/// @param set Dataset
/// @param config Config
///
void ParallelPerm::run_slave_permutations(int num_permutations, Dataset& set, Config& config){

  Config perm_config = config;
  perm_config.models_to_keep(1);

  total_perms = num_permutations;
  
  int seed = receive_seed();
  while(seed <= num_permutations){
    set.permute_set(seed);

    SlaveAnalysis analyzer;
    analyzer.run_slave_analysis(set, perm_config);
    
    seed = receive_seed();
  }

}


///
/// Broadcasts seed to all slaves for use in permutation tests
/// @param seed 
///
void ParallelPerm::send_seed(int seed){
  MPI_Bcast(&seed, 1, MPI_INT, 0, MPI_COMM_WORLD);
}


///
/// Receives seed from master
/// @return seed
///
int ParallelPerm::receive_seed(){

  int seed;
  MPI_Bcast(&seed, 1, MPI_INT, 0, MPI_COMM_WORLD);
  return seed;
  
}




