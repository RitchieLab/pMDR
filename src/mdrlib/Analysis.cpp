//Analysis.cpp

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

#include "Analysis.h"
#include "KnuthComboGenerator.h"
#include "BioComboGenerator.h"
#include <map>
#include <algorithm>  

namespace mdr{

///
/// Set parameters for doing combination production
/// @param total_loci
/// @param min_size
/// @param max_size
/// @return
///
void Analysis::set_combination_limits(int total_loci, int min_size, int max_size){

  // set size of models to generate
  generator->ComboEnds(min_size, max_size);

  // set number of loci to use in generating the combinations
  generator->SetLoci(total_loci);

  generator->SetComboInterval(10000);

}


Analysis::~Analysis(){
if(generator != NULL) delete generator;
}


///
/// Creates the correct generator type based on whether a file
/// is passed for biofilter
/// @param biofilename string
///
void Analysis::set_generator_type(string biofilename, int idSize){
  if(biofilename.empty()){
    generator = new KnuthComboGenerator();
  }
  else{
    generator = new BioComboGenerator();
    generator->set_id_size(idSize);
    map<string, string> temp_map;
    temp_map["BIOFILENAME"] = biofilename;
    generator->PassInitConfig(temp_map);
  }

}



///
/// Sets parameters for analysis
/// @param config Config containing parameters
/// @param set
/// @param log_all
///
void Analysis::set_parameters(Config& config, Dataset& set, bool log_all){

  // establish the indexes of any loci that are always included
  // in models
  vector<string> included_snps = config.loci_always_included();
  always_included.clear();
  for(unsigned int i=0; i<included_snps.size(); i++){
    always_included.push_back(set.locus_name_index(included_snps[i]));
  }

  error_threshold = config.get_balacc_thresh();
  vector<string> force_loci_snps = config.force_loci_included();
  forced_loci.clear();

  for(unsigned int i=0; i<force_loci_snps.size(); i++){
    forced_loci.push_back(set.locus_name_index(force_loci_snps[i]));
  }

  unsigned int num_included = always_included.size();

  // set parameters on combination generator
  set_combination_limits(set.get_num_loci(), config.model_size_start()-num_included, config.model_size_end()-num_included);

  generator->SetComboInterval(config.build_interval());
  generator->set_dataset(&set);

  // set parameters for analysis
  set_max_models(config.models_to_keep());
  mdr.set_tiecell_value(config.tie_cell_value());
  mdr.set_not_in_training(config.not_recognized_response());
  mdr.initialize(config.num_crossval(), config.fitness());

  total_cv = config.num_crossval();
  set.set_crossval_interval(total_cv);

  log_all_training(log_all);

  check_interval = config.check_interval();
  checkpointfile = config.check_point_file();
  basecheckpoint = config.basename();

}

///
/// Adds any included snps to the combinations for testing
/// @param combinations
///
void Analysis::add_included_snps(vector<vector<unsigned int> >& combinations){
  unsigned int combo, loc, curr_combo;
  unsigned int always_included_size = always_included.size();
  unsigned int num_combos = combinations.size();
  bool match;

  vector<vector<unsigned int> > new_combinations;


  for(combo = 0; combo < num_combos; combo++){
      match = false;
      unsigned int combo_size = combinations[combo].size();
      for(loc = 0; loc < always_included_size; loc++){
        for(curr_combo=0; curr_combo < combo_size; curr_combo++){
          if(always_included[loc] == combinations[combo][curr_combo]){
            match = true;
            break;
          }
        }
      }
      if(!match){
        combinations[combo].insert(combinations[combo].end(), always_included.begin(), always_included.end());
        new_combinations.push_back(combinations[combo]);
      }
  }
  combinations = new_combinations;
}


///
/// Tests single model
/// @param loci vector<int>
/// @param models
/// @param set
///
void Analysis::test_single_model(vector<unsigned int>& loci, vector<Model>& models, Dataset& set,
	bool calc_prediction){
   unsigned int size_combo = loci.size();
   unsigned int i, curr_cv;

   if(size_combo != current_model_size){
      for(i=0; i<total_cv; i++){
        models[i].resize_cell_status_vector(set.converter().get_size_array(size_combo));
      }
      current_model_size = size_combo;
    }

    for(i=0; i<total_cv; i++){
      models[i].reset_stats();
      models[i].combination = loci;
    }

	if(!calc_prediction)
	    mdr.test_model(models, set);
	else
		mdr.test_model_complete(models, set);
    
    for(curr_cv=0; curr_cv < total_cv; curr_cv++){
    	if(models[curr_cv].training.balanced_error <= error_threshold)
	      training_results.Insert(models[curr_cv].training.balanced_error, models[curr_cv].combination, curr_cv);
    }
}


///
/// Writes all models currently held to output and flushes the current list
/// @param config Config with output details
/// 


///
/// Outputs models
/// @param set Dataset
/// @param config Config with details on output
/// @param base_out basename for output
///
void Analysis::output_models(Dataset& set, Config& config, string base_out){

  out_writer.open_form(base_out + ".mdr.out", config.model_size_end(), mdr.fitness_name(),
    set.any_missing_data());
  out_writer.open_tab(base_out + ".mdr.tab.txt", config.model_size_end(),
    mdr.fitness_name());

  if(config.display_partitions())
    out_writer.open_part(base_out + ".mdr.partitions.txt", config.model_size_end());

  if(config.display_totals())
    out_writer.open_matrix(base_out + ".mdr.totals.txt", config.model_size_end());

  if(config.display_rules())
    out_writer.open_rules(base_out + ".mdr.rules.txt", config.model_size_end());

  int model_size_start = config.model_size_start();
  int model_size_end = config.model_size_end();
  int curr_size;
  ModelTreeNode* curr_node;
  vector<int> loci;
  Model curr_model;

  // loop through each model in each crossvalidation and
  // output all information for each
  for(unsigned int curr_cv=0; curr_cv < total_cv; curr_cv++){
    for(curr_size=model_size_start; curr_size <= model_size_end; curr_size++){
      // resize model to match current size
      curr_model.resize_cell_status_vector(set.converter().get_size_array(curr_size));

      for(curr_node = training_results.result_trees[curr_cv][curr_size]->GetLast();
        curr_node != NULL; curr_node = curr_node->GetPrev()){
          // blank stats in model
          curr_model.reset_stats();
          curr_model.combination = curr_node->GetData();
          mdr.test_model(curr_model, set, curr_cv);
          out_writer.output_model(curr_cv, curr_model, total_cv, set,
            mdr.get_calculator());
      }
    }
  }

}


///
/// Returns best model for each order based on cross-validation
/// consistency.  Ties are broken with average prediction error
/// for the models that tie
/// @param set Dataset
/// @param config Config
///
vector<Model> Analysis::get_best_models(Dataset& set, Config& config){

  ModelTreeNode* best_node;
  int model_size_start = config.model_size_start();
  int model_size_end = config.model_size_end();
  int curr_size;
  vector<unsigned int> best_combo;

  Model curr_model;
  vector<Model> best_models(model_size_end+1, curr_model);

  map<vector<unsigned int>, Model>::iterator map_iter;

  Stat* calculator = mdr.get_calculator();

  // for each model order determine best model
  // only look at best models in each tree when calculating
  for(curr_size = model_size_start; curr_size <= model_size_end; curr_size++){

    curr_model.resize_cell_status_vector(set.converter().get_size_array(curr_size));

    float avg_bal_error = 0.0;
    float avg_error = 0.0;

    map<vector<unsigned int>, Model> model_map;

    for(unsigned int curr_cv=0; curr_cv < total_cv; curr_cv++){
      best_node=training_results.result_trees[curr_cv][curr_size]->GetLast();
      curr_model.reset_stats();
      curr_model.combination = best_node->GetData();

      mdr.test_model(curr_model, set, curr_cv);
      if(model_map.find(curr_model.combination) == model_map.end()){
        curr_model.set_cvc(1);
        curr_model.set_predictavg(curr_model.testing.error);
        curr_model.set_balpredavg(curr_model.testing.balanced_error);
        model_map[curr_model.combination] = curr_model;
      }
      else{
        model_map[curr_model.combination].set_cvc(1 + model_map[curr_model.combination].get_cvc());
        model_map[curr_model.combination].set_predictavg(curr_model.testing.error +
           model_map[curr_model.combination].get_predictavg());
        model_map[curr_model.combination].set_balpredavg(curr_model.testing.balanced_error +
          model_map[curr_model.combination].get_balpredavg());
      }

      avg_bal_error += calculator->convert_to_display(curr_model.testing.balanced_error);
      avg_error += curr_model.testing.error;
    }

    // find best cvc, in case of tie use better balanced prediction error
    int highest_cvc = 0;
    for(map_iter = model_map.begin(); map_iter != model_map.end(); map_iter++){
      map_iter->second.set_balpredavg(map_iter->second.get_balpredavg() /
        map_iter->second.get_cvc());
      map_iter->second.set_predictavg(map_iter->second.get_predictavg() /
        map_iter->second.get_cvc());

      if(map_iter->second.get_cvc() > highest_cvc){
        highest_cvc = map_iter->second.get_cvc();
        best_combo =  map_iter->first;
      }
      else if(map_iter->second.get_cvc() == highest_cvc){
        if(map_iter->second.get_balpredavg() < model_map[best_combo].get_balpredavg())
          best_combo = map_iter->first;
      }
    }

    // put best model on to vector
    best_models[curr_size] = model_map[best_combo];
    best_models[curr_size].set_totalmissdata(best_models[curr_size].training.totalmissing + best_models[curr_size].testing.totalmissing);

    best_models[curr_size].set_balpredavg(calculator->convert_to_display(avg_bal_error / total_cv));
    best_models[curr_size].set_predictavg(avg_error / total_cv);
  }

  return best_models;

}


///
/// Returns best model for each order based on cross-validation
/// consistency.  Ties are broken with average prediction error
/// for the models that tie
/// @param set Dataset
/// @param config Config
///
void Analysis::get_pval_models(Dataset& set, Config& config, 
	vector<vector<Model > >& pvalmodels){

  ModelTreeNode* best_node;
  int model_size_start = config.model_size_start();
  int model_size_end = config.model_size_end();
  int curr_size=0;
  vector<unsigned int> best_combo;
  
  vector<Model> tmpModels;
  pvalmodels.assign(model_size_end+1, tmpModels);

  Model curr_model;
  map<vector<unsigned int>, Model>::iterator map_iter;

  Stat* calculator = mdr.get_calculator();

  // for each model order determine best model
  // only look at best models in each tree when calculating
  for(int curr_size = model_size_start; curr_size <= model_size_end; curr_size++){
	unsigned int models_kept = training_results.result_trees[0][curr_size]->GetCount();
    curr_model.resize_cell_status_vector(set.converter().get_size_array(curr_size));

    map<vector<unsigned int>, Model> model_map;
    map<vector<unsigned int>, int> cvc_count;

	vector<float> avgbalerror(training_results.result_trees[0][curr_size]->GetCount(), 0.0);
	vector<float> avgerror(training_results.result_trees[0][curr_size]->GetCount(), 0.0);

	unsigned int modrank=0;
	// need a vector of bal error totals for each position 
	// need a vector of models so that you can use the bal error average to set model worth for position
	for(unsigned int curr_cv=0; curr_cv < total_cv; curr_cv++){
		modrank=0;
		for(best_node=training_results.result_trees[curr_cv][curr_size]->GetLast(); best_node != NULL;
			best_node=best_node->GetPrev()){
    	  curr_model.reset_stats();
	      curr_model.combination = best_node->GetData();
	      mdr.test_model(curr_model, set, curr_cv);	
      	  if(model_map.find(curr_model.combination) == model_map.end()){
		        curr_model.set_cvc(1);
    	    curr_model.set_predictavg(curr_model.testing.error);
        	curr_model.set_balpredavg(curr_model.testing.balanced_error);
	        model_map[curr_model.combination] = curr_model;
	        cvc_count[curr_model.combination] = 1;
    	  }	      
    	  else{
	        model_map[curr_model.combination].set_cvc(1 + model_map[curr_model.combination].get_cvc());
	        model_map[curr_model.combination].set_predictavg(curr_model.testing.error +
           	model_map[curr_model.combination].get_predictavg());
    	    model_map[curr_model.combination].set_balpredavg(curr_model.testing.balanced_error +
          	model_map[curr_model.combination].get_balpredavg());
          	cvc_count[curr_model.combination]++;  	
    	  }
    	  
    	  avgbalerror[modrank] += calculator->convert_to_display(curr_model.testing.balanced_error);
    	  avgerror[modrank] +=  curr_model.testing.error;
      		modrank++;
		}
		
	}


    for(map_iter=model_map.begin(); map_iter != model_map.end(); ++map_iter){
    	pvalmodels[curr_size].push_back(map_iter->second);
    	pvalmodels[curr_size].back().set_balpredavg(pvalmodels[curr_size].back().get_balpredavg() /
    		cvc_count[pvalmodels[curr_size].back().combination]);
    	pvalmodels[curr_size].back().set_predictavg(pvalmodels[curr_size].back().get_predictavg() /
    		cvc_count[pvalmodels[curr_size].back().combination]);
    } 
    
    std::sort(pvalmodels[curr_size].begin(), pvalmodels[curr_size].end(),sortModelsByCVC);
    
    // need to limit this calculation to the number of MODELSTOKEEP listed in config
	// so that we can calculate the balanced accuracy correctly
    // vector at this size now holds the ordered list of models
    // update values in the vector


    for(unsigned int i=0; i < models_kept; i++){
    	pvalmodels[curr_size][i].set_totalmissdata(pvalmodels[curr_size][i].training.totalmissing + 
    		pvalmodels[curr_size][i].testing.totalmissing);
    	pvalmodels[curr_size][i].set_balpredavg(calculator->convert_to_display(avgbalerror[i]/total_cv));
    	pvalmodels[curr_size][i].set_predictavg(avgerror[i]/total_cv);
    }
    if(pvalmodels[curr_size].size() > models_kept)
	    pvalmodels[curr_size].erase(pvalmodels[curr_size].begin()+models_kept, pvalmodels[curr_size].end());
  }

}


///
/// Outputs best models for each size
/// @param best_models Model
/// @param set Dataset
///
void Analysis::output_best(vector<Model> best_models, Dataset& set){

  // convert the models to the full model
  for(unsigned int i=1; i<best_models.size(); i++){
    // skip blank models
    if(best_models[i].combination.size() == 0)
      continue;
    best_models[i].reset_stats();
    mdr.full_model(best_models[i], set);
  }

  out_writer.output_best(best_models, total_cv, set, mdr.get_calculator());
}


///
/// Outputs single best model
/// @param best_model Single best model across all orders
/// @param set Dataset used in determining single best model
///
void Analysis::output_single_best(Model& best_model, Dataset& set){
  out_writer.output_single_best(best_model, set);
}


///
/// Selects single best model across all model sizes
/// @param best_models vector of Model where each model is the best for the order
/// @param config Config
///
Model Analysis::get_single_best(vector<Model> best_models, Config& config){

  int model_size_start = config.model_size_start();
  int model_size_end = config.model_size_end();

  int bestCV = model_size_start;
  int bestAvgPredError = model_size_start;

  for(int size=model_size_start+1; size <= model_size_end; size++){
    if(best_models[size].get_cvc() > best_models[bestCV].get_cvc())
      bestCV = size;
    if(best_models[size].get_balpredavg() < best_models[bestAvgPredError].get_balpredavg())
      bestAvgPredError = size;
  }

  if(bestCV == bestAvgPredError)
    return best_models[bestCV];
  // when have same cross-validation consistency
  // select one with lower best average balanced prediction error
  else if(best_models[bestCV].get_cvc() == best_models[bestAvgPredError].get_cvc())
    return best_models[bestAvgPredError];
  // otherwise return one with fewer loci
  else if(best_models[bestCV].combination.size() < best_models[bestAvgPredError].combination.size()){
    return best_models[bestCV];
  }
  else{
    return best_models[bestAvgPredError];
  }
}


///
/// Output p values
/// @param models
/// @param set
/// @param p_tests
///
void Analysis::output_p_values(vector<Model>& models, Dataset& set, int p_tests){
  out_writer.output_pvalues(models, set, p_tests);
}


///
///
/// Output p values with multiple models for each size
/// @param models
/// @param set
/// @param p_tests
///
void Analysis::output_p_values(vector<vector<Model> >& models, Dataset& set, int p_tests){
  out_writer.output_pvalues(models, set, p_tests);
}


///
/// Output LR p values
/// @param models
/// @param set
/// @param p_tests
///
void Analysis::output_lr_p_values(vector<Model>& models, Dataset& set, int p_tests){
  out_writer.output_lr_pvalues(models, set, p_tests);
}


///
/// Reads checkpoint file and sets the generator and
/// the training list
/// @param
/// @param
/// @return number of models completed
///
double Analysis::read_checkpoint(string checkpointfile){
  return check_pointer.read_checkpoint(training_results, generator, checkpointfile);
}


///
/// Produces checkpoint file
/// @param
/// @throws MDRExcept when can't write checkpoint
///
void Analysis::write_checkpoint(float interval, string basefile, double completed_mods){

  // get previous filename
  string last_file = check_pointer.checkpoint_name();

  check_pointer.write_checkpoint(training_results, generator, basefile, interval, completed_mods);
  if(last_file.length() > 0)
    remove(last_file.c_str());

}

///
/// Output all models
///
void Analysis::output_all_models(vector<Model>& models,std::vector<unsigned int>& combination){

  out_writer.output_all(models, combination);

	// clear to save memory
// 	training_results.clear();
}
bool sortModelsByCVC(mdr::Model& left, mdr::Model& right){
	if(left.get_cvc() != right.get_cvc())
    	return left.get_cvc() > right.get_cvc();
    return left.get_balpredavg() < right.get_balpredavg();
}
}
