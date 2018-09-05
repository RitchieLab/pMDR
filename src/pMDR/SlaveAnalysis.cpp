//SlaveAnalysis.cpp

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


#include "SlaveAnalysis.h"
#include "pBioComboGenerator.h"
#include "pKnuthComboGenerator.h"


SlaveAnalysis::~SlaveAnalysis(){
}

///
/// Runs analysis
///
void SlaveAnalysis::run_slave_analysis(Dataset& set, Config& config){

  set_generator_type(config.get_biofilter_filename(), config.get_id_size());
  set_parameters(config, set, false);

  bool done = false;
  unsigned int num_combos, curr_combo;
  current_model_size = 0;
  vector<Model> models;
  models.resize(total_cv);
  training_results.initialize_list(total_cv, generator->get_combo_max()+always_included.size(), 
    max_train);

    
  if(output_all)
    out_writer.open_all_out(flush_fn, total_cv);    
  
  int parameter_signal=0;
  int request_number = 1;
  
  while(1){
   parameter_signal = receive_slave_parameters(request_number++);

    if(parameter_signal==0){
      done = generator->GenerateCombinations();

      if(always_included.size() > 0)
        add_included_snps(generator->ComboList);
      num_combos = generator->ComboList.size();
//cout << "slave got " << num_combos << " from generator" << endl;
      for(curr_combo=0; curr_combo < num_combos; curr_combo++){
    	 test_single_model(generator->ComboList[curr_combo], models, set, output_all);
         if(output_all){
            for(unsigned int i=0; i<total_cv; i++){
              if(models[i].training.balanced_error <= error_threshold){
                 output_all_models(models, generator->ComboList[curr_combo]);
                 break;
              }
            }
         }
//       if(output_all)
//     	 output_all_models(models, generator->ComboList[curr_combo]);
      }
    }
    else if(parameter_signal == -1){
      send_models();
      training_results.clear(); // clear trees before proceeding
      receive_continue_signal(); // wait for continue signal 
    }
    else{ // received terminate notice from master
      send_models();
      break;
    }
  }
}


///
/// Creates the correct generator type based on whether a file
/// is passed for biofilter.
/// @param biofilename string
///
void SlaveAnalysis::set_generator_type(string biofilename, int idSize){
  if(biofilename.empty()){
    pKnuthComboGenerator *kgen = new pKnuthComboGenerator();
    generator = kgen;
    paramMover = kgen;
  }
  else{
    pBioComboGenerator *bgen = new pBioComboGenerator();
    generator = bgen;
    generator->set_id_size(idSize);
    paramMover = bgen;
    map<string, string> temp_map;
    temp_map["BIOFILENAME"] = biofilename;
    generator->PassInitConfig(temp_map);
  }
}


///
/// Receive combination generator parameters from master
///
int SlaveAnalysis::receive_slave_parameters(int request_number){
  int srRequestParameters = 62;

int sendjunk = request_number;
  MPI_Send(&sendjunk, 1, MPI_INT, 0, srRequestParameters, MPI_COMM_WORLD);

  int parameter_signal = paramMover->ReceiveParameters(srParameters);
  
  return parameter_signal;
}


///
/// Set parameters for doing combination production
/// @param total_loci
/// @param min_size
/// @param max_size
/// @return
///
void SlaveAnalysis::set_combination_limits(int total_loci, int min_size, int max_size){

  Analysis::set_combination_limits(total_loci, min_size, max_size);
  
  // initializes current state of the combination generator
  generator->initialize_state();
  
}
 
///
/// Send models to master
///
void SlaveAnalysis::send_models(){
  int index = 0;
  int order_min = generator->get_combo_min() + always_included.size();
  int order_max = generator->get_combo_max() + always_included.size();

  int models_kept = training_results.max_models_kept();
  int temp_array_size = (order_max+1) * total_cv * models_kept * (order_max-order_min+1);

  int curr_order, curr_loc;
  unsigned int curr_cv;
  int total_models;
  ModelTreeNode* curr_node;
  vector<unsigned int> loci;
  
  float * temp = new float[temp_array_size];

  for(curr_cv = 0; curr_cv < total_cv; curr_cv++){
    for(curr_order = order_min; curr_order <= order_max; curr_order++){
      total_models=0;
      for(curr_node = training_results.result_trees[curr_cv][curr_order]->GetLast();
        curr_node != NULL; curr_node = curr_node->GetPrev()){
        
        // first send loci
        loci = curr_node->GetData();
        int num_loci = loci.size();
        for(curr_loc=0; curr_loc < num_loci; curr_loc++)
          temp[index++] = loci[curr_loc];
        for(; curr_loc < order_max; curr_loc++)
          temp[index++] = -1;
        
        temp[index++] = curr_node->GetKey();
        total_models++;
      }

      while(total_models < models_kept){
        for(curr_loc=0; curr_loc < order_max; curr_loc++)
          temp[index++] = -1;
        temp[index++] = -1;
        total_models++;
      }
    }
  }

  // Send the temp array
  MPI_Send(temp, temp_array_size, MPI_FLOAT, 0, 22, MPI_COMM_WORLD);
  
  delete [] temp;
}


///
/// Receive broadcast for continuing 
///
void SlaveAnalysis::receive_continue_signal(){
  int continueAnalysis ;
  MPI_Bcast(&continueAnalysis, 1,MPI_INT, 0, MPI_COMM_WORLD);

}
