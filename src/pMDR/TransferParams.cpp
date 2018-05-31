//TransferParams.cpp

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


#include "TransferParams.h"

TransferParams::TransferParams(){
  colsplit = 10;
}

///
/// send parameters to slaves
/// @param config Config
/// @param set Dataset
///
void TransferParams::send_params(Config& config, Dataset& set){
  int NumCfgOptions = 20;
  int * cfg = new int[NumCfgOptions];
  
  cfg[0] = set.get_num_inds();
  cfg[1] = set.get_num_loci();
  cfg[2] = set.get_max_locus_value();
  cfg[3] = set.get_missing_value();
  cfg[4] = config.model_size_start();
  cfg[5] = config.model_size_end();
  cfg[6] = config.num_crossval();
  cfg[7] = config.build_interval();
  cfg[8] = config.random_seed();
  cfg[9] = config.models_to_keep();
  cfg[10] = config.shuffle_data();
  cfg[11] = config.not_recognized_response();
  cfg[12] = config.tie_cell_value();
  cfg[13] = set.any_missing_data();
  cfg[14] = config.num_ptests();
  cfg[15] = config.num_ptests_completed();
  cfg[16] = config.calc_threshold_each_model();
  unsigned int num_included_loci = config.loci_always_included().size();
  cfg[17] = num_included_loci;
  cfg[18] = config.fitness();
  cfg[19] = config.get_biofilter_filename().length();
  
   
  MPI_Bcast(cfg, NumCfgOptions, MPI_INT, 0, MPI_COMM_WORLD);
  
  delete [] cfg;
  
  // when have loci that are always included pass to slaves
  if(num_included_loci > 0){
    int * loci = new int[num_included_loci];
    vector<string> included_snps = config.loci_always_included();
    for(unsigned int i=0; i<num_included_loci; i++)
      loci[i] = set.locus_name_index(included_snps[i]);
    MPI_Bcast(loci, num_included_loci, MPI_INT, 0, MPI_COMM_WORLD);
    delete [] loci;
  }
  
  
  // If have to pass name of bio filter file do it here
  if(config.get_biofilter_filename().length() > 0){
    char bfilename[config.get_biofilter_filename().length() +1];
    strcpy(bfilename, config.get_biofilter_filename().c_str());
    MPI_Bcast(bfilename,config.get_biofilter_filename().length()+1,
      MPI_CHAR, 0, MPI_COMM_WORLD);
  }
  
}
  
///
/// receives parameters from master
/// @param config Config holding parameters
/// @param set Dataset
/// @returns vector containing locus indexes of any loci that are
/// to be included in all models analyzed
///
vector<unsigned int> TransferParams::receive_params(Config& config, Dataset& set){
  int NumCfgOptions = 20;
  int * cfg = new int[NumCfgOptions];

  MPI_Bcast(cfg, NumCfgOptions, MPI_INT, 0, MPI_COMM_WORLD);

  set.set_num_loci(cfg[1]);
  set.set_num_inds(cfg[0]);
  set.set_max_genotype_value(cfg[2]);
  config.max_locus(cfg[2]);

  set.set_missing_value(cfg[3]);

  config.model_size_start(cfg[4]);
  config.model_size_end(cfg[5]);
  config.num_crossval(cfg[6]);
  config.build_interval(cfg[7]);
  config.random_seed(cfg[8]);
  set.set_shuffle_seed(config.random_seed());
  config.models_to_keep(cfg[9]);
  config.shuffle_data(bool(cfg[10]));
  config.not_recognized_response(cfg[11]);
  config.tie_cell_value(cfg[12]);
  
  set.set_contains_missing(bool(cfg[13]));

  config.num_ptests(cfg[14]);
  config.num_ptests_completed(cfg[15]);
  config.calc_threshold_each_model(cfg[16]);
  int num_included_loci = cfg[17];
  config.fitness(static_cast<fitnessTypes>(cfg[18]));
  int biofilename_len = cfg[19];
  
  
  delete [] cfg;
  
  vector<unsigned int> always_included_loci;
  
  if(num_included_loci > 0){
    int * loci = new int [num_included_loci];
    MPI_Bcast(loci, num_included_loci, MPI_INT, 0, MPI_COMM_WORLD);    
    for(int i=0; i<num_included_loci; i++){
      always_included_loci.push_back(loci[i]);
    }
    delete [] loci;
  }

  if(biofilename_len > 0){
    char bfile[biofilename_len+1];
    MPI_Bcase(bfile, biofilename_len+1, MPI_CHAR, 0, MPI_COMM_WORLD);
    config.set_biofilter_filename(bfile);
  }

  return always_included_loci;
  
}

///
/// send data
/// @param set Dataset to transmit
///
void TransferParams::send_data(Dataset& set){

  float threshold = set.get_set_threshold();
  // send set threshold
  MPI_Bcast(&threshold, 1, MPI_FLOAT, 0, MPI_COMM_WORLD);
  int cols = set.get_num_loci()+1;
  int rows = set.get_num_inds();

  
  int numofsplits = cols/colsplit; // additional column is for status
  int leftovers = cols % colsplit;
  
  int colcounter = 0;
  
  int * senddata = new int[rows * colsplit];
  
  // send each split
  for(int s=0; s<numofsplits; s++){
    for(int i=0; i<rows; i++){
      for(int j=0; j<colsplit; j++){
        senddata[i*colsplit+j] = set.data[i][j+colcounter];
      }
    }
    MPI_Bcast(senddata, rows*colsplit, MPI_INT, 0, MPI_COMM_WORLD);
    colcounter += colsplit;
  }
  
  delete [] senddata;
  
  // send remainder
  int * leftdata = new int[rows*leftovers];
  
  for(int i=0; i<rows; i++)
    for(int j=0; j<leftovers; j++)
      leftdata[i*leftovers+j] = set.data[i][j+colcounter];
  
  MPI_Bcast(leftdata, rows*leftovers, MPI_INT, 0, MPI_COMM_WORLD);
  
  delete [] leftdata;
  
}

/// 
/// receive data
/// @param set Dataset to hold data
///
void TransferParams::receive_data(Dataset& set){
  // get threshold for entire set
  float threshold;
  MPI_Bcast(&threshold, 1, MPI_FLOAT, 0, MPI_COMM_WORLD);
  set.set_entire_threshold(threshold);
  
  // create array to hold entire set
  set.initialize_data();
  
  int cols = set.get_num_loci()+1;
  int rows = set.get_num_inds();
  
  int numofsplits = cols/colsplit; // additional column is for status
  int leftovers = cols % colsplit;
  
  int * receivedata = new int[rows * colsplit];
  
  int colcounter = 0;
  
  for(int s=0; s<numofsplits; s++){
    MPI_Bcast(receivedata, rows*colsplit, MPI_INT, 0, MPI_COMM_WORLD);
    for(int i=0; i < rows; i++)
      for(int j=0; j < colsplit; j++)
        set.data[i][j+colcounter] = receivedata[i*colsplit+j];
    colcounter = colcounter + colsplit;
  }
  
  delete [] receivedata;
  
  int * leftdata = new int[rows*leftovers];

  // Receive the leftovers
  MPI_Bcast(leftdata, rows*leftovers, MPI_INT, 0, MPI_COMM_WORLD);

  // Store them in the data array
  for(int i=0; i < rows; i++)
    for(int j=0; j < leftovers; j++)
      set.data[i][j+colcounter] = leftdata[i*leftovers+j];

  delete [] leftdata;
  
}
