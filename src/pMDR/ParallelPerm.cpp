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


float ParallelPerm::get_p_value(float value, int rank){
  
  unsigned int counter;
   
  // use top rank tree if non ranked is selected
  if(permtype != Ranked){
  	rank=0;
  }

  
  // find nearest max to the value
   PermModelTreeNode* curr_node=permuted_models[rank].FindNearestMin(value);
    if(curr_node != NULL)
  	 counter=curr_node->GetData();
  	else
  	  counter=0;


  if(permtype == Combined){
  	return float(counter)/(permuted_models[rank].GetCount());
  }
  else{
	 return float(counter)/total_perms;
  }
    
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
  permtype = config.perm_method();
//   pval_sorter.SetMaxSize(models_kept);
  vector<string> none_always_included;

  perm_config.force_loci_included(none_always_included);
//   perm_config.models_to_keep(1);
  perm_config.check_point_file("");
  perm_config.check_interval(0);

  int total_cv = perm_config.num_crossval();
  total_perms = num_permutations;
  int curr_perm;
  
  LRInteract lr_calc;
  
  int models_kept = config.models_to_keep();
  PermModelTree pval_sorter;
 
  if(permtype==Combined || permtype == TopRank){
  	permuted_models.assign(1, pval_sorter);
  	permuted_models[0].SetMaxSize(models_kept*total_perms);
  	pval_sorter.SetMaxSize(models_kept);
  }
  else{ // Ranked
  	permuted_models.assign(config.models_to_keep(), pval_sorter);
  } 
  
  pval_sorter.SetMaxSize(models_kept);
  
  for(curr_perm = 1; curr_perm <= num_permutations; curr_perm++){
    // send seed to the slave nodes
    send_seed(curr_perm);
    set.permute_set(curr_perm);

    // create analysis object
    MasterAnalysis permuted_analysis;
    permuted_analysis.run_master_analysis(set, perm_config, nproc, log_out, false);

    // get best models for each order    
    vector<Model> best_models =  permuted_analysis.get_best_models(set, perm_config);

    vector<vector<Model> > pvalmodels;
    permuted_analysis.get_pval_models(set, perm_config, pvalmodels);

    // get single best model   
    Model single_best = permuted_analysis.get_single_best(best_models, perm_config);

    float fitness;

    // iterate through the pval models and insert the fitness values 
    // into tree 
	for(unsigned int i=perm_config.model_size_start(); i<=perm_config.model_size_end(); i++){	  	  
	  	for(unsigned int j=0; j<pvalmodels[i].size(); j++){
	  		if(total_cv == 1)
	  			fitness=pvalmodels[i][j].training.balanced_error;
	  		else{
	  			fitness=pvalmodels[i][j].get_balpredavg();
	  		}
	  		pval_sorter.Add(fitness, 0);
	  	}
	}

	PermModelTreeNode* curr_node=pval_sorter.GetLast(); 
	// add the model values into correct tree for ranking
	int max_models=models_kept;
	if(config.perm_method()==TopRank){
		max_models=1;
	}

	// for TopRank only the best model is included in the first rank
	// for Combined all the models are inserted in the first rank
	// for Ranked put models in appropriate ranks
	for(int i=0; i<max_models; i++){
		if(config.perm_method()==Ranked)
			permuted_models[i].Add(curr_node->GetKey(), 0);
		else
			permuted_models[0].Add(curr_node->GetKey(), 0);
		curr_node=curr_node->GetPrev();	
		if(curr_node == NULL)
			break;
	}


    if(config.regress_test()){
      lr_calc.calculate_lr_stat(set, single_best);
      // store lr result in vector
      lr_results.push_back(single_best.get_interact_llr());
    }

    log_out.add_perm_model(single_best, curr_perm, set,  permuted_analysis.get_calculator(),
      config.regress_test());
    pval_sorter.Clear();

  }
  
  // send seed which exceeds number of permutations to perform
  send_seed(curr_perm);
  
  // sort the lr results
  sort(lr_results.begin(), lr_results.end(), sort_ascending());
  
  PermModelTreeNode* perm_node;
  // assign the ranking to the tree nodes
  for(unsigned int i=0; i<permuted_models.size(); i++){	
  	unsigned int counter=0;
  	for(perm_node=permuted_models[i].GetLast(); perm_node != NULL;
	  	perm_node=perm_node->GetPrev()){
  		perm_node->SetData(++counter);
  	}
  }  

}


///
/// Run permutations for slave
/// @param number of permutations
/// @param set Dataset
/// @param config Config
///
void ParallelPerm::run_slave_permutations(int num_permutations, Dataset& set, Config& config){

  Config perm_config = config;
//   perm_config.models_to_keep(1);

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




