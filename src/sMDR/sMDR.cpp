// sMDR.cpp

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
#include "PermutationTester.h"
#include <ConfigFileReader.h>
#include <FileManager.h>
#include <FlatIndex.h>
#include <MapFileReader.h>
#include <ModelFileReader.h>
#include <stdio.h>

using namespace std;
using namespace mdr;

void exit_app(MDRExcept& me);
void usage();
void version(string version_date, string version_num);


int main(int argc, char* argv[]){

  string version_date = "11/05/18";
  string version_num = "1.0.1";

  version(version_date, version_num);
  // check for configuration file
  if(argc < 2){
    usage();
    MDRExcept ex("");
    exit_app(ex);
  }

  string configfilename(argv[1]);

  try{
    ConfigFileReader reader;
    Config config_info = reader.read_config(configfilename);

    // set random seed
    srand(config_info.random_seed());

    // read data
    Dataset set;
    set.set_max_genotype_value(config_info.max_locus());
    set.set_missing_value(config_info.missing_data_value());
    set.set_random_shuffle(config_info.shuffle_data());
    set.set_shuffle_seed(config_info.random_seed());
    set.set_permodelthreshold(config_info.calc_threshold_each_model());

    FileManager file_reader;
    file_reader.parse_input_file(config_info.datafile(), set, config_info.input_format());
    set.initialize_converter(config_info.model_size_start(), config_info.model_size_end());

    if(config_info.mapfile().length() > 0){
      MapFileReader map_reader;
      map_reader.parse_map_file(config_info.mapfile(), set);
    }
    else{
      set.reference_snp_names();
    }

//     if(set.any_missing_data())
//       FlatIndex::set_genos_per_locus(set.get_max_locus_value()+2);
//     else
//       FlatIndex::set_genos_per_locus(set.get_max_locus_value()+1);
//
//     FlatIndex::set_included_indexes(config_info.model_size_start(), config_info.model_size_end(),
//       set.any_missing_data(), set.get_max_locus_value()+1);

    // set up log file
//     stringstream logstream;
//     logstream <<  "Number of models: "<< total_mods << "\n";
//     logfile.add_line(logstream.str());
//     logfile.add_time();
//     logfile.add_line("\n");

    std::string base_out;
    if(argc==3){
      base_out = argv[2];
    }
    else{
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
    }

    config_info.basename(base_out);
    LogOutput log_out(base_out + ".mdr.log");

    string tempModelfn=base_out + "mdrmodels.tmp.txt";

    // convert model file for use in sMDR
    if(config_info.getModelFileName().length() > 0){
      int idSize =  ModelFileReader::convert_model_file(config_info.getModelFileName(),tempModelfn);
      config_info.set_biofilter_filename(tempModelfn);
      config_info.set_id_size(idSize);
    }

    sAnalysis analyzer;

    cout << "\n\tAnalyzing models...";
    cout.flush();
    
    if(config_info.getOutputAllModels()){
	    analyzer.set_flushfn(base_out + ".allmodels.txt");
	}

    analyzer.train_models(set, config_info, log_out, true);
    cout << endl;
    
    if(config_info.getOutputAllModels()){
    	cout << "\nReported all models to " << base_out << ".allmodels.txt\n\n";
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
      cout << "\n\tPerforming permutation tests...";
      cout.flush();
      
      vector<vector<Model> > pvalmodels;
      analyzer.get_pval_models(set, config_info, pvalmodels);

      PermutationTester perms;
      perms.run_permutations(config_info.num_ptests(), set, config_info, log_out);
      
	  // assign p values
	  float value;


	  for(unsigned int i=config_info.model_size_start(); i<=config_info.model_size_end(); i++){	  	  
	  	for(unsigned int j=0; j<pvalmodels[i].size(); j++){
	  		if(config_info.num_crossval() == 1)
	  			value=pvalmodels[i][j].training.balanced_error;
	  		else{
	  			value=pvalmodels[i][j].get_balpredavg();
	  		}
	  			
	  		pvalmodels[i][j].set_pvalue(perms.get_p_value(value));
	  		if(config_info.regress_test()){
	  			pvalmodels[i][j].set_lr_pvalue(perms.get_lr_p_value(pvalmodels[i][j].get_interact_llr()));
		  	}
	  	}
	  }      
            
//       float value;
//       // assign p values
//       for(unsigned int i=config_info.model_size_start(); i<best_models.size(); i++){
//         if(config_info.num_crossval() == 1)
//           value = best_models[i].training.balanced_error;
//         else
//           value = best_models[i].get_balpredavg();
// 
//         best_models[i].set_pvalue(perms.get_p_value(value));
// 
//         if(config_info.regress_test()){
//           best_models[i].set_lr_pvalue(perms.get_lr_p_value(best_models[i].get_interact_llr()));
//         }
// 
//       }
      analyzer.output_p_values(pvalmodels, set, config_info.num_ptests());
      if(config_info.regress_test()){
        analyzer.output_lr_p_values(best_models, set, config_info.num_ptests());
      }
    }

   // delete temporary model file
   if(config_info.getModelFileName().length() > 0){
      remove(tempModelfn.c_str());
    }

  }catch(MDRExcept ex){
    exit_app(ex);
  }

  cout << endl << endl;

}


///
/// Output version  information
/// @param version_date
///
void version(string version_date, string version_num){
  cout << endl << "\t\tsMDR version: " << version_num << endl
    << "\t\t\t" << version_date << endl << endl;
}



///
/// Output usage information
/// @param version_date string with sMDR version
///
void usage(){
  cout << "\tusage:\t\tsMDR  <configuration file> [output name]" << endl;
  cout << "\texample:\tsMDR mdr.config sample" << endl << endl;
}

///
/// Outputs message in MDRExcept and exits application
/// @param  me MDRExcept error object
///
void exit_app(MDRExcept & me){
  std::cout << me.get_error() << std::endl;
  exit(1);
}

