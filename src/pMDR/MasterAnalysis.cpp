//MasterAnalysis.cpp

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

#include "MasterAnalysis.h"
#include "pBioComboGenerator.h"
#include "pKnuthComboGenerator.h"

///
/// Runs analysis and sends parameters to slave nodes as needed
///
void MasterAnalysis::run_master_analysis(Dataset& set, Config& config, int nproc,
  LogOutput& log_out, bool log_all){

  set_generator_type(config.get_biofilter_filename());

  set_parameters(config, set, log_all);

  bool done = false;

  double completed_mods = 0.0;
  double total_mods = generator->calc_combinations();

  int numIntervals = 20;
  if(total_mods > 10000000)
    numIntervals = 100;
  double reportInterval = total_mods / numIntervals;
  double nextInterval = reportInterval;

  double nextCheck = total_mods * 2;
  double checkModels = 0.0;

  // For use with master processing
  vector<Model> models;
  models.resize(total_cv);
  current_model_size = 0;

  training_results.initialize_list(total_cv, generator->get_combo_max()+always_included.size(),
    max_train);

  // if a checkpoint file is being used to restart an analysis
  // read it in
  if(checkpointfile.length() > 0){
    generator->initialize_state();
    generator->param_AlreadyStarted(true);
    completed_mods = read_checkpoint(checkpointfile);
    nextInterval = completed_mods + reportInterval;
  }

  // establish a checkpoint interval if needed
  if(check_interval > 0){
    // checkpoints are indicated as a percentage of the
    // full run -- reset nextCheck for that
    checkModels = total_mods * (check_interval/100);
    nextCheck = checkModels + completed_mods;
  }


  if(forced_loci.size() > 0){
    total_mods = 1;
    done = true;
    test_single_model(forced_loci, models, set);
  }

  if(log_training){
    log_out.add_time();
    log_out.total_models(total_mods);
    log_out.add_line("\n");
  }

  unsigned int combo_interval = generator->getComboInterval();
  unsigned int master_model_interval = combo_interval / (4 * nproc);

  if(nproc > 50)
    master_model_interval = 0;

  // Request Variables
  MPI_Request RequestParameterspending;

  // Status Variables
  MPI_Status status;

  int RequestReceived = 0;

  // MPI Tags
  int srRequestParameters = 62;

  // Number of nodes that have terminated
  int termcount = 1;
  int recvjunk = 0;
  unsigned int num_combos = 0;
  unsigned int original_num_combos = 0;
  unsigned int curr_combo = 0;
  
  if(output_all)
    out_writer.open_all_out(flush_fn, total_cv);

  // Listen for Nodes Requesting Parameters
  MPI_Irecv(&recvjunk, 1, MPI_INT, MPI_ANY_SOURCE, srRequestParameters,
    MPI_COMM_WORLD, &RequestParameterspending);

  // Until all nodes have finished
  while(termcount < nproc){
    // Check to see if a Node has requested Parameters
    MPI_Test(&RequestParameterspending, &RequestReceived , &status);

    if(RequestReceived){
			// if all combinations are not completed send parameters
      if(!done){
        paramMover->SendParameters(status.MPI_SOURCE, srParameters);

        // run through to increment the parameters in combo generator
        done = generator->AdvanceParameters();

        if(log_training){
          completed_mods += generator->get_count_generated();
          if(completed_mods >= nextInterval){
            while(nextInterval < completed_mods){
              nextInterval += reportInterval;
            }
            log_out.report_completed(completed_mods, total_mods);
          }
        }

        // listen again
        MPI_Irecv(&recvjunk, 1, MPI_INT, MPI_ANY_SOURCE, srRequestParameters,
          MPI_COMM_WORLD, &RequestParameterspending);

        if(completed_mods > nextCheck && !done){
          receive_checkpoint_models(&RequestParameterspending, nproc);
          float curr_check_interval = completed_mods / total_mods * 100;
          write_checkpoint(curr_check_interval, basecheckpoint, completed_mods);
          while(nextCheck < completed_mods)
            nextCheck += checkModels;
          // start receiving requests again
          MPI_Irecv(&recvjunk, 1, MPI_INT, MPI_ANY_SOURCE, srRequestParameters,
            MPI_COMM_WORLD, &RequestParameterspending);
        }

      }
      else{ // send termination notice to slave
        paramMover->SendStop(status.MPI_SOURCE, srParameters);
        // listen again if not all processors done
        if(++termcount < nproc)
          MPI_Irecv(&recvjunk, 1, MPI_INT, MPI_ANY_SOURCE, srRequestParameters,
            MPI_COMM_WORLD, &RequestParameterspending);
      }
    }
    else{ // no node has requested parameters, so master can do a little work
      if(!done && master_model_interval){
        done = generator->GenerateCombinations(master_model_interval);

        original_num_combos = generator->ComboList.size();
        if(always_included.size() > 0)
          add_included_snps(generator->ComboList);

        num_combos = generator->ComboList.size();
        for(curr_combo=0; curr_combo < num_combos; curr_combo++){
          test_single_model(generator->ComboList[curr_combo], models, set, output_all);
          if(output_all)
            for(unsigned int i=0; i<total_cv; i++){
              if(models[i].training.balanced_error <= error_threshold){
                 output_all_models(models, generator->ComboList[curr_combo]);
                 break;
              }
            }
//      	    output_all_models(models, generator->ComboList[curr_combo]);
        }
        if(log_training){
          completed_mods += original_num_combos;
          if(completed_mods >= nextInterval){
            while(nextInterval < completed_mods){
              nextInterval += reportInterval;
            }
            log_out.report_completed(completed_mods, total_mods);
          }
        }
        // write check point if needed
        if(completed_mods > nextCheck && !done){
          receive_checkpoint_models(&RequestParameterspending, nproc);
          float curr_check_interval = completed_mods / total_mods * 100;
          write_checkpoint(curr_check_interval, basecheckpoint, completed_mods);
          while(nextCheck < completed_mods)
            nextCheck += checkModels;
            // start receiving requests again
          MPI_Irecv(&recvjunk, 1, MPI_INT, MPI_ANY_SOURCE, srRequestParameters,
            MPI_COMM_WORLD, &RequestParameterspending);
        }
      }
    }
  } // all analysis finished

  // gather models and store in master
  for(int i=1; i < nproc; i++){
    receive_models(i);
  }

}


///
/// Creates the correct generator type based on whether a file
/// is passed for biofilter.
/// @param biofilename string
///
void MasterAnalysis::set_generator_type(string biofilename){
  if(biofilename.empty()){
    pKnuthComboGenerator *kgen = new pKnuthComboGenerator();
    generator = kgen;
    paramMover = kgen;
  }
  else{
    pBioComboGenerator *bgen = new pBioComboGenerator();
    generator = bgen;
    paramMover = bgen;
    map<string, string> temp_map;
    temp_map["BIOFILENAME"] = biofilename;
    generator->PassInitConfig(temp_map);
  }

}

///
/// Receives models from all slaves
/// As each slave requests parameters tell slave to
/// send models
/// @param nproc total number of nodes
///
void MasterAnalysis::receive_checkpoint_models(MPI_Request* RequestParameterspending,
  int nproc){

  MPI_Status status;

  int RequestReceived = 0;

  int recvjunk = 0;
  int nodes_sent = 1;

  // MPI Tags
  int srRequestParameters = 62;

  while(nodes_sent < nproc){
    MPI_Test(RequestParameterspending, &RequestReceived , &status);
    // send indicator to slave node to send models
    if(RequestReceived){
      paramMover->SendModelSignal(status.MPI_SOURCE, srParameters);

      receive_models(status.MPI_SOURCE);

      ++nodes_sent;

      // check again
      if(nodes_sent < nproc)
        MPI_Irecv(&recvjunk, 1, MPI_INT, MPI_ANY_SOURCE, srRequestParameters,
          MPI_COMM_WORLD, RequestParameterspending);
    }
  }

  // Broadcast continue signal to nodes
  send_continue_signal();

}


///
/// Send continue signal to nodes
/// @param node Node to signal
///
void MasterAnalysis::send_continue_signal(){
  int continueAnalysis = 1;

  MPI_Bcast(&continueAnalysis, 1, MPI_INT, 0, MPI_COMM_WORLD);

}


///
/// receive models from slaves
///
void MasterAnalysis::receive_models(int node){

  MPI_Status status;
  int index = 0;
  int order_min = generator->get_combo_min() + always_included.size();
  int order_max = generator->get_combo_max() + always_included.size();

  int models_kept = training_results.max_models_kept();
  int temp_array_size = (order_max+1) * total_cv * models_kept * (order_max-order_min+1);

  int curr_order, curr_loc, curr_model;
  unsigned int curr_cv;
  int total_models;

  vector<unsigned int> loci;
  float model_fitness;

  float * temp = new float[temp_array_size];
  MPI_Recv(temp, temp_array_size, MPI_FLOAT, node, 22, MPI_COMM_WORLD, &status);

  for(curr_cv = 0; curr_cv < total_cv; curr_cv++){
    total_models=0;
    for(curr_order = order_min; curr_order <= order_max; curr_order++){
      // check to see if model is empty
      for(curr_model = 0; curr_model < models_kept; curr_model++){
        if(temp[index] == -1){
          index += order_max+1;
        }
        else{
          loci.clear();
          for(curr_loc =0; curr_loc < order_max; curr_loc++){
            if(temp[index] != -1)
              loci.push_back(int(temp[index]));
            index++;
          }
          model_fitness = temp[index++];

          training_results.Insert(model_fitness, loci, curr_cv);
        }
      }
    }
  }

  delete [] temp;
}


///
/// Set parameters for doing combination production
/// @param total_loci
/// @param min_size
/// @param max_size
/// @return
///
void MasterAnalysis::set_combination_limits(int total_loci, int min_size, int max_size){

  Analysis::set_combination_limits(total_loci, min_size, max_size);

  // initializes current state of the combination generator
  generator->initialize_state();

}



