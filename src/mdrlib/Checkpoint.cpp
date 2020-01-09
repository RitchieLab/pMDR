//Checkpoint.cpp

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

#include "Checkpoint.h"


namespace mdr{

///
/// Writes checkpoint file
///
void Checkpoint::write_checkpoint(TrainList& models, ComboGenerator* generator, string filename,
  float interval, double mods_complete){

  char format_interval[30];
  sprintf(format_interval, "%.3f", interval);
  string checkpointname = filename + "_" + format_interval;

  checkpoint_name(checkpointname);

  ofstream checkstream;
  checkstream.open(checkpointname.c_str(), ofstream::out);
  
  if(!checkstream.is_open())
    throw MDRExcept("Unable to open " + checkpointname + " for writing");
  
  checkstream << mods_complete << endl;
  
  generator->WriteCheckpoint(checkstream);

  // output information listing numbers of models, cv, etc
  
  // output total CV
  int total_cv = models.num_crossvalidation_intervals();
  checkstream << total_cv << endl;
  
  int curr_cv;
  unsigned int curr_model_size, total_models, curr_loc;
  ModelTreeNode* curr_node;
 
  // output models 
  for(curr_cv = 0; curr_cv < total_cv; ++curr_cv){
    checkstream << curr_cv << endl;
  
    // only output model sizes where count is greater than 0
    for(curr_model_size=0; curr_model_size < models.result_trees[curr_cv].size();
      ++curr_model_size){
     
      total_models = models.result_trees[curr_cv][curr_model_size]->GetCount();
      // skip if no models
      if(total_models == 0)
        continue;
      
      // output size and total models for that size
      checkstream << curr_model_size << " " << total_models << endl;
      
      for(curr_node = models.result_trees[curr_cv][curr_model_size]->GetLast();
        curr_node != NULL; curr_node = curr_node->GetPrev()){
        
        for(curr_loc=0; curr_loc < curr_model_size; ++curr_loc){
          checkstream << curr_node->GetData().combination[curr_loc] << " ";
        }
        checkstream << curr_node->GetKey() << endl;
      }
     
    }
    // indicate end of cv interval
    checkstream << "-1 -1" << endl;
  }
  checkstream << curr_cv << endl;
  checkstream.close();

}


///
/// Reads checkpoint
///
double Checkpoint::read_checkpoint(TrainList& models, ComboGenerator* generator, string checkfilename){
  checkpoint_name(checkfilename);
  double already_completed = read_checkpoint(models, generator);
  checkpoint_name("");
  return already_completed;
}


///
/// Reads checkpoint file and fills objects with correct values
/// @returns number of models completed
///
double Checkpoint::read_checkpoint(TrainList& models, ComboGenerator* generator){
  
  ifstream checkstream;
  
  checkstream.open(filename.c_str());
  
  if(!checkstream.is_open())
    throw MDRExcept("Unable to open checkpoint file: " + filename);
  
  double completed_mods;
  
  // get parameters for combination generator
  checkstream >> completed_mods;
  
  generator->ReadCheckpoint(checkstream);
  
  int total_cv, curr_model_size, loc, total_models,curr_model;

  vector<unsigned int> model_loci;
  float fitness;
  
  checkstream >> total_cv;
  int curr_cv;
  checkstream >> curr_cv;
  
  Model mod;

  while(curr_cv < total_cv){
    
    checkstream >> curr_model_size >> total_models;
    
    while(curr_model_size != -1){
      model_loci.assign(curr_model_size, 0);
      for(curr_model=0; curr_model < total_models; ++curr_model){
        for(loc=0; loc < curr_model_size; ++loc)
          checkstream >> model_loci[loc];
        checkstream >> fitness;
        mod.combination = model_loci;
        models.Insert(fitness, mod, curr_cv);
      }
      
      checkstream >> curr_model_size >> total_models;

    }
    checkstream >> curr_cv;
  }
  
  checkstream.close();

  return completed_mods;
}

}
