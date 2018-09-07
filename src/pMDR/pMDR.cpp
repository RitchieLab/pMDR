//pMDR.cpp

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


// Parallel version of MDR
#include "MasterAnalysis.h"
#include "SlaveAnalysis.h"
#include "ParallelPerm.h"
#include <Stringmanip.h>
#include <ConfigFileReader.h>
#include <FileManager.h>
#include <FlatIndex.h>
#include <MapFileReader.h>
#include <ModelFileReader.h>
#include <stdio.h>

using namespace std;
using namespace mdr;


void stop_program(MDRExcept& me);
void set_dataset_params(Config& config_info, Dataset& set);
string get_basename(string configfilename);
void send_id_size(int idSize);
int receive_id_size();


int main(int argc, char* argv[]){
  string version_num = "1.0.1";

  int myrank;
  int nproc;

  // set up MPI
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
  MPI_Comm_size(MPI_COMM_WORLD, &nproc);

  if(myrank ==0){
    cout << "\n\tParallel MDR " << version_num << "\t09/07/18" << endl;
  }

  if(argc < 2){
    cout << "\n\tUsage:  pMDR <config file> [output]\n\n";
    return 1;
  }

  string configfilename(argv[1]);

  ConfigFileReader reader;
  Config config_info;
  Dataset set;
  string base_out;

  try{
    config_info = reader.read_config(configfilename);
    // set random seed
    srand(config_info.random_seed());
    set_dataset_params(config_info, set);

    FileManager file_reader;
    file_reader.parse_input_file(config_info.datafile(), set, config_info.input_format());
    set.initialize_converter(config_info.model_size_start(), config_info.model_size_end());

    // read map file
    if(config_info.mapfile().length() > 0){
      MapFileReader map_reader;
      map_reader.parse_map_file(config_info.mapfile(), set);
    }
    else{
      set.reference_snp_names();
    }

  }
  catch(MDRExcept& me){
    stop_program(me);
    return(1);
  }

      if(argc==3){
        base_out = argv[2];
      }
      else{
        base_out = get_basename(configfilename);
      }
      config_info.basename(base_out);

  string tempModelfn=base_out + "mdrmodels.tmp.txt";

  // need to split here between master and slave analysis
  if(myrank ==0){ // master
    try{

      LogOutput log_out(base_out + ".mdr.log");
      if(config_info.getModelFileName().length() > 0){
        int idSize =  ModelFileReader::convert_model_file(config_info.getModelFileName(),tempModelfn);
        config_info.set_biofilter_filename(tempModelfn);
        send_id_size(idSize);
        config_info.set_id_size(idSize);
      }
      MasterAnalysis analyzer;
      
      if(config_info.getOutputAllModels()){
	    analyzer.set_flushfn(base_out + ".1.allmodels.txt");
	  }
      
      cout << "Analyzing models...";
      cout.flush();
      analyzer.run_master_analysis(set, config_info, nproc, log_out, true);
      cout << endl;
      
     if(config_info.getOutputAllModels()){
//    	cout << "\nReported all models to " << base_out << ".1.allmodels.txt\n\n";
        cout << "\nReported all models to:\n";
        for(int i=0; i<nproc; i++){
          string fn = base_out + "." + Stringmanip::itos(i+1) + ".allmodels.txt";
          cout << fn << "\n";
        }
    	MPI_Finalize();
    	return 0;
     } 
      
       // output results
      analyzer.output_models(set, config_info, base_out);

      vector<Model> best_models =  analyzer.get_best_models(set, config_info);

      LRInteract lr_calc;
      if(config_info.regress_test()){
        for(unsigned int i=config_info.model_size_start(); i<best_models.size(); i++)
          lr_calc.calculate_lr_stat(set, best_models[i]);
      }

      analyzer.output_best(best_models, set);

      Model single_best = analyzer.get_single_best(best_models, config_info);
      analyzer.output_single_best(single_best, set);
      if(config_info.num_ptests() >= 1){
        cout << "Performing permutation tests...";
        cout.flush();
        ParallelPerm perms;
        perms.run_master_permutations(config_info.num_ptests(), set, config_info, nproc, log_out);

        float value;
        // assign p values
        for(unsigned int i=config_info.model_size_start(); i<best_models.size(); i++){
          if(config_info.num_crossval() == 1)
            value = best_models[i].training.balanced_error;
          else
            value = best_models[i].get_balpredavg();
          best_models[i].set_pvalue(perms.get_p_value(value));

          if(config_info.regress_test()){
            best_models[i].set_lr_pvalue(perms.get_lr_p_value(best_models[i].get_interact_llr()));
          }
        }
        analyzer.output_p_values(best_models, set, config_info.num_ptests());
        if(config_info.regress_test()){
          analyzer.output_lr_p_values(best_models, set, config_info.num_ptests());
        }
      }
   // delete temporary model file
       if(config_info.getModelFileName().length() > 0){
             remove(tempModelfn.c_str());
        }
      cout << endl;
    }
    catch(MDRExcept& me){
      stop_program(me);
      return(1);
    }
  }
  else{ //slave

    if(config_info.getModelFileName().length() > 0){
        config_info.set_biofilter_filename(tempModelfn);
        int idSize=receive_id_size();
        config_info.set_id_size(idSize);
    }

    SlaveAnalysis analyzer;

    analyzer.set_rank(myrank);
    
    if(config_info.getOutputAllModels()){
//       analyzer.set_flushfn(base_out + ".1.allmodels.txt");
       string fn = base_out + "." + Stringmanip::itos(myrank+1) + ".allmodels.txt";
       analyzer.set_flushfn(fn);
//       analyzer.set_flushfn(base_out + "." + Stringmanip::itos(myrank+1) + ".allmodels.txt");
//               cout << fn << endl;
	 }    

    analyzer.run_slave_analysis(set, config_info);
    
	if(config_info.getOutputAllModels()){
    	MPI_Finalize();
    	return 0;
     }    
     
    if(config_info.num_ptests() >= 1){
      ParallelPerm perms;
      perms.run_slave_permutations(config_info.num_ptests(), set, config_info);
    }
  }

  MPI_Finalize();
  return 0;

}


///
/// Stop program
/// @param me MDRExcept with error message
///
void stop_program(MDRExcept& me){
  int continueAnalysis=0;
  MPI_Bcast(&continueAnalysis, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Finalize();
  cout << me.get_error() << endl;
}

///
/// Send model file ID size
///
void send_id_size(int idSize){
  MPI_Bcast(&idSize, 1, MPI_INT, 0, MPI_COMM_WORLD);
}

///
/// Receive model file ID size
///
int receive_id_size(){
  int idSize;
  MPI_Bcast(&idSize, 1, MPI_INT, 0, MPI_COMM_WORLD);
  return idSize;
}


///
/// Sets parameters in dataset
/// @param config_info Config
/// @param set Dataset
///
void set_dataset_params(Config& config_info, Dataset& set){
  // prepare dataset
  set.set_max_genotype_value(config_info.max_locus());
  set.set_missing_value(config_info.missing_data_value());
  set.set_random_shuffle(config_info.shuffle_data());
  set.set_shuffle_seed(config_info.random_seed());
  set.set_permodelthreshold(config_info.calc_threshold_each_model());
}

///
/// Returns base name from configuration file
/// @param string configuration file name
/// @returns string with base name to use in creating new files
///
string get_basename(string configfilename){
  string base_out;
  string::size_type directorySep, period;
  directorySep = configfilename.rfind("/");
  if(directorySep==string::npos)
    directorySep=0;
  else
    directorySep++;
  period = configfilename.rfind(".");
  if(period==string::npos)
    period = configfilename.length();
  base_out = configfilename.substr(0,period);
  return base_out;
}


