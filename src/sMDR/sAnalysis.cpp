//Analysis.cpp

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

#include "sAnalysis.h"

///
/// Sets configuration parameters for run and
/// then analyzes models
/// @param set Dataset
/// @param config Config containing parameters
/// @param log_out LogOutput
/// @param log_progress bool specifying whether progress should be recorded in log
///
void sAnalysis::train_models(Dataset& set, Config& config, LogOutput& log_out, bool log_all){
  
  set_generator_type(config.get_biofilter_filename());
  set_parameters(config, set, log_all);
  train_models(set, log_out);
}


///
/// Runs through combinations and stores best models in tree
/// @param set Dataset 
/// @param log_out LogOutput
/// @return
///
void sAnalysis::train_models(Dataset& set, LogOutput& log_out){
  double completed_mods=0.0;
  double total_mods = generator->calc_combinations();

  int numIntervals = 20;
  if(total_mods > 10000000)
    numIntervals = 100;
  double reportInterval = total_mods / numIntervals;
  double nextInterval = reportInterval;
  
  // establish nextCheck as too big a number to ever reach
  double nextCheck = total_mods * 2;
  double checkModels = 0.0;
  
  if(log_training){  
    log_out.add_time();
    log_out.total_models(total_mods);
    log_out.add_line("\n");
  }
  
  vector<Model> models;

  models.resize(total_cv);
  unsigned int  curr_combo, num_combos;
  current_model_size = 0;

  training_results.initialize_list(total_cv, generator->get_combo_max()+always_included.size(), max_train);

  // if a checkpoint file is being used to restart an analysis
  // read it in 
  if(checkpointfile.length() > 0){
    generator->initialize_state();
    generator->param_AlreadyStarted(true);
    completed_mods = read_checkpoint(checkpointfile);
    nextInterval = completed_mods + reportInterval;
  }


// cout << "check_interval=" << check_interval << endl;
  // establish a checkpoint interval if needed
  if(check_interval > 0){
    // checkpoints are indicated as a percentage of the 
    // full run -- reset nextCheck for that
    checkModels = total_mods * (check_interval/100);
    nextCheck = checkModels + completed_mods;
  }
// cout << "nextCheck=" << nextCheck << endl;
// cout << "completed_mods=" << completed_mods << endl;

  bool combos_done;
// 

  if(output_all)
	  out_writer.open_all_out(flush_fn, total_cv);

//   vector<vector<uint> >::iterator combo_iter;
  if(forced_loci.size() > 0){
//   while(!generator.GenerateCombinations()){
    test_single_model(forced_loci, models, set);
  }
  else{
  do{
    combos_done = generator->GenerateCombinations();

    if(always_included.size() > 0)
      add_included_snps(generator->ComboList);

    num_combos = generator->ComboList.size();
//cout << "num_combos=" << num_combos << endl;
    for(curr_combo=0; curr_combo < num_combos; curr_combo++){
      test_single_model(generator->ComboList[curr_combo], models, set, output_all);
      if(output_all)
      	output_all_models(models, generator->ComboList[curr_combo]);
    }
    
    completed_mods += num_combos;

    if(log_training){  
      if(completed_mods >= nextInterval){
        while(nextInterval < completed_mods){
          nextInterval += reportInterval;
        }
        log_out.report_completed(completed_mods, total_mods);
      }
    }
    if(completed_mods > nextCheck && !combos_done){
      float curr_check_interval = completed_mods / total_mods * 100;
      write_checkpoint(curr_check_interval, basecheckpoint, completed_mods);
      while(nextCheck < completed_mods)
        nextCheck += checkModels;
// cout << "called write_checkpoint nextCheck=" << nextCheck << endl;
    }
    
  }while(!combos_done);
  }

}


