//MDR.cpp

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

#include "MDR.h"
#include "stdlib.h"


namespace mdr{

///
/// Constructor
///
MDR::MDR(){
  initialize(5, BalancedAccuracy);
}


///
/// Constructor
///
MDR::MDR(int crossval){
  initialize(crossval, BalancedAccuracy);
}

///
/// Sets cv and intitializes parameters
///
void MDR::set_cv(int CrossVal){
  initialize(CrossVal, BalancedAccuracy);
}

///
/// Destructor
///
MDR::~MDR(){
  delete fitness_calculator;
  delete error_calculator;
}


///
/// Initializes vectors for storing totals during model testing
/// @param crossval number of crossvalidations to perform
/// @param type fitnessTypes
/// @return
///
void MDR::initialize(int crossval, fitnessTypes type){

  fitness_calculator = NULL;
  error_calculator = NULL;

  total_cv = crossval;

  last_model_size = 0;
  
  vector<int> temp;
  zero_cv.assign(2,temp);
  
  for(int z=0; z<crossval; z++){
    zero_totals.assign(crossval, zero_cv);
  }
  
  set_calculators(type);
  
}


///
/// Creates and sets the calculators for model
/// @param type fitnessType
///
void MDR::set_calculators(fitnessTypes type){

  if(fitness_calculator != NULL)
    delete fitness_calculator;
  if(error_calculator != NULL)
    delete error_calculator;

  error_calculator = stat_creator.create_stat(Error);
  fitness_calculator = stat_creator.create_stat(type);
  
}



///
/// Checks to see if includedIndexes has been initialized 
/// and if not, does so
/// @param set Dataset
///
void MDR::check_included(Dataset& set){
  if(includedIndexes.size() ==0){
    includedIndexes = set.converter().get_included_indexes();
  }
}



///
/// Runs MDR on indicated model for all crossvalidation intervals
/// and returns results in mdr_results vector
/// @param models input/output vector for results from analysis with each index corresponding
/// to a crossvalidation interval
/// @param set Dataset to use in testing
/// @return
///
void MDR::test_model(vector<Model>& models, Dataset& set){

  check_included(set);

  int array_size = set.converter().get_size_array(models[0].combination.size());
  int total_inds;
  if(last_model_size != models[0].combination.size()){
    last_model_size = models[0].combination.size();
    included_cells = includedIndexes[last_model_size];
    // establish correct totals size
    for(int current_cv=0; current_cv < total_cv; current_cv++){
      zero_totals[current_cv][0].assign(array_size,0);
      zero_totals[current_cv][1].assign(array_size,0);
    }
    zero_cv[0].assign(array_size, 0);
    zero_cv[1].assign(array_size, 0);
  }
  // store the totals across all crossvalidation intervals
  fill_cv_totals(totals, set, models[0].combination);
  
  float threshold = set.determine_threshold(models[0].combination);

  for(int current_cv=0; current_cv < total_cv; current_cv++){
    models[current_cv].set_threshold(threshold);
 
    
    if(total_cv > 1){
      total_inds = set_cv_totals(totals, cv_totals, current_cv, total_cv);
      calculate_err(cv_totals, cv_totals, models[current_cv].get_threshold(), included_cells, 
        models[current_cv].training, models[current_cv]);
    }
    else{
      calculate_err(totals[0], totals[0], models[current_cv].get_threshold(), included_cells, 
        models[current_cv].training, models[current_cv]);
      total_inds = set.get_num_inds();
    }
    
    models[current_cv].set_totalmissdata(total_inds-(models[current_cv].training.classhigh + models[current_cv].training.classlow
      + models[current_cv].training.misclasshigh + models[current_cv].training.misclasslow));
    
    models[current_cv].set_combined_fitness(100-models[current_cv].training.balanced_error);
  }
  
}

///
/// Runs MDR on indicated model for all crossvalidation intervals 
/// and returns results in mdr_results vector.  Calculates testing in addition to training
/// results.
/// @param models input/output vector for results from analysis with each index corresponding
/// to a crossvalidation interval
/// @param set Dataset to use in testing
/// @return
///
void MDR::test_model_complete(vector<Model>& models, Dataset& set){

  check_included(set);

  int array_size = set.converter().get_size_array(models[0].combination.size());
  int total_inds;
  if(last_model_size != models[0].combination.size()){
    last_model_size = models[0].combination.size();
    included_cells = includedIndexes[last_model_size];
    // establish correct totals size
    for(int current_cv=0; current_cv < total_cv; current_cv++){
      zero_totals[current_cv][0].assign(array_size,0);
      zero_totals[current_cv][1].assign(array_size,0);
    }
    zero_cv[0].assign(array_size, 0);
    zero_cv[1].assign(array_size, 0);
  }
  // store the totals across all crossvalidation intervals
  fill_cv_totals(totals, set, models[0].combination);
  
  float threshold = set.determine_threshold(models[0].combination);

  for(int current_cv=0; current_cv < total_cv; current_cv++){
    models[current_cv].set_threshold(threshold);
 
    
    if(total_cv > 1){
      total_inds = set_cv_totals(totals, cv_totals, current_cv, total_cv);
      calculate_err(cv_totals, cv_totals, models[current_cv].get_threshold(), included_cells, 
        models[current_cv].training, models[current_cv]);
      // calculate testing 
      calculate_err(cv_totals, totals[current_cv], models[current_cv].get_threshold(), included_cells,  
        models[current_cv].testing, models[current_cv]);
      models[current_cv].testing.totalmissing = set.get_testing_crossval_total(current_cv) - (models[current_cv].testing.classhigh 
        + models[current_cv].testing.classlow + models[current_cv].testing.misclasshigh 
        + models[current_cv].testing.misclasslow);
    }
    else{
      calculate_err(totals[0], totals[0], models[current_cv].get_threshold(), included_cells, 
        models[current_cv].training, models[current_cv]);
      total_inds = set.get_num_inds();
    }
    
    models[current_cv].set_totalmissdata(total_inds-(models[current_cv].training.classhigh + models[current_cv].training.classlow
      + models[current_cv].training.misclasshigh + models[current_cv].training.misclasslow));
    
    models[current_cv].set_combined_fitness(100-models[current_cv].training.balanced_error);
  }
  
}


///
/// Runs mdr and fills information on the model into 
/// vectors contained within the Model passed into it
/// @param curr_model Model to store results
/// @param set Dataset to test model on
/// @param cv Crossvalidation to test
/// @return
///
void MDR::test_model(Model& curr_model, Dataset& set, int current_cv){
  check_included(set);


  if(last_model_size != curr_model.combination.size()){
    int array_size = set.converter().get_size_array(curr_model.combination.size());
    last_model_size = curr_model.combination.size();
    included_cells = includedIndexes[last_model_size];
    for(int current_cv=0; current_cv < total_cv; current_cv++){
      zero_totals[current_cv][0].assign(array_size,0);
      zero_totals[current_cv][1].assign(array_size,0);
    }
    zero_cv[0].assign(array_size, 0);
    zero_cv[1].assign(array_size, 0);    
  }
  // store the totals across all crossvalidation intervals
  fill_cv_totals(totals, set, curr_model.combination);
  curr_model.set_threshold(set.determine_threshold(curr_model.combination));

  if(total_cv > 1){
    set_cv_totals(totals, cv_totals, current_cv, total_cv);
    // fill training
    calculate_err(cv_totals, cv_totals, curr_model.get_threshold(), included_cells,
      curr_model.training, curr_model);

    // fill testing
    calculate_err(cv_totals, totals[current_cv], curr_model.get_threshold(), included_cells,  
      curr_model.testing, curr_model);
    curr_model.testing.totalmissing = set.get_testing_crossval_total(current_cv) - (curr_model.testing.classhigh 
      + curr_model.testing.classlow + curr_model.testing.misclasshigh 
      + curr_model.testing.misclasslow);
  }
  else{
    calculate_err(totals[0], totals[0], curr_model.get_threshold(), included_cells,
      curr_model.training, curr_model);
    
  }


  curr_model.training.totalmissing = set.get_training_crossval_total(current_cv) - (curr_model.training.classhigh + curr_model.training.classlow
      + curr_model.training.misclasshigh + curr_model.training.misclasslow);
  curr_model.set_totalmissdata(set.get_num_inds()-curr_model.training.totalmissing
    -curr_model.testing.totalmissing);

  curr_model.set_combined_fitness(100-curr_model.training.balanced_error);

}


///
/// Calculates model on all data within the set and
/// fills model with this information
/// @param curr_model Model
/// @param set Dataset
///
void MDR::full_model(Model& curr_model, Dataset& set){
  check_included(set);
  if(last_model_size != curr_model.combination.size()){
    int array_size = set.converter().get_size_array(curr_model.combination.size());
    last_model_size = curr_model.combination.size();
    included_cells = includedIndexes[last_model_size];
    for(int current_cv=0; current_cv < total_cv; current_cv++){
      zero_totals[current_cv][0].assign(array_size,0);
      zero_totals[current_cv][1].assign(array_size,0);
    }
    zero_cv[0].assign(array_size, 0);
    zero_cv[1].assign(array_size, 0);    
  }
  fill_cv_totals(totals, set, curr_model.combination);

  curr_model.set_threshold(set.determine_threshold(curr_model.combination));
  // passing -1 means all cv will be included in the calculations
  int total_inds = set_cv_totals(totals, cv_totals, -1, total_cv);
  // fill training
  calculate_err(cv_totals, cv_totals, curr_model.get_threshold(), included_cells,
    curr_model.training, curr_model);
  curr_model.set_totalmissdata(total_inds-(curr_model.training.classhigh + curr_model.training.classlow
      + curr_model.training.misclasshigh + curr_model.training.misclasslow)-(curr_model.testing.classhigh + curr_model.testing.classlow
      + curr_model.testing.misclasshigh + curr_model.testing.misclasslow));

}


///
/// Fills the totals vector with totals for each testing cv
/// these will be used in each cv to calculate the fitness
/// @param
/// @param
///
void MDR::fill_cv_totals(vector<vector<vector<int> > >& totals, Dataset& MasterData,
  vector<unsigned int>& loci){

  totals = zero_totals;
  
  int crossvals = MasterData.get_crossval_interval();
  if(crossvals > 1)
    for(int cv=0; cv<crossvals; cv++){
      LinearMDR(MasterData.testing[cv], MasterData, loci, totals[cv]);
    }
  else
    LinearMDR(MasterData.training[0], MasterData, loci, totals[0]);
}

///
/// sums the totals for the training cross-validations for use in 
/// the particular cross-validation being tested
/// @param totals vector containing totals for all crossvalidations
/// @param cv_totals vector to store training totals for the current crossvalidation
/// @param cv current crossvalidation
/// @param total_cv total crossvalidations
/// @return total number of individuals included
///
int MDR::set_cv_totals(vector<vector<vector<int> > >& totals, 
  vector<vector<int> >& cv_totals, int cv, int total_cv){

    unsigned int curr_cell;
    unsigned int total_cells = totals[0][0].size();
    cv_totals = zero_cv;
    int total_inds = 0;

    for(int curr_cv=0; curr_cv < total_cv; curr_cv++){
      if(curr_cv==cv)
        continue;
      for(curr_cell=0; curr_cell<total_cells; curr_cell++){
        cv_totals[0][curr_cell] += totals[curr_cv][0][curr_cell];
        cv_totals[1][curr_cell] += totals[curr_cv][1][curr_cell];
        total_inds += totals[curr_cv][0][curr_cell] + totals[curr_cv][1][curr_cell];
      }
      
    }
  return total_inds;
  
}


void MDR::LinearMDR(vector<int>& current_set, Dataset& MasterData,
  vector<unsigned int>& loci, vector<vector<int> >& totals){
  
  // establish a vector of correct size that can be used with
  vector<unsigned int> genotype = loci;

  unsigned int currLoc;
  unsigned int total_inds = current_set.size();
  unsigned int total_loci = loci.size();
  
  // place each individual in proper location in model linear arrays
  // based on location and genoype
  for(unsigned int currInd=0; currInd < total_inds; currInd++){
    // create genotype vector
    for(currLoc=0; currLoc < total_loci;
      currLoc++){
      genotype[currLoc] = MasterData.data[current_set[currInd]][loci[currLoc]];  
    } 

    // unaffected will be in vector 0, affected in vector 1  
    totals[MasterData.data[current_set[currInd]][0]][MasterData.converter().flatten_indexes(genotype)]++;
  }  
  
}

/// 
/// Calculates and returns error and balanced error
/// @param thresh_totals totals in each cell of the model of unaffected 
/// and affected individuals
/// @param totals contains totals to be tested in determining the error
/// @param threshold ratio for determining whether a cell is high or low risk
/// @param includedCells contains indexes into thresh_totals and totals vectors
/// for cells to be included in this calculation
/// @param stats model_stats to store results in
/// @param model vector recording cell status
/// @return
///
void MDR::calculate_err(vector<vector<int> >& thresh_totals, 
        vector<vector<int> >& totals, float threshold, vector<int>& includedCells,
        model_stats& stats, Model& model){

  float calculatedRisk;
  unsigned int numCells=includedCells.size();
  
  // set totals in the model to match calculated totals
  stats.affected = totals[1];
  stats.unaffected = totals[0];
  
  for(unsigned int currCell=0; currCell < numCells; currCell++){  
    // check that unaffected total in cell is greater than zero
    if(thresh_totals[0][includedCells[currCell]]>0)
      calculatedRisk = float(thresh_totals[1][includedCells[currCell]]) /
        thresh_totals[0][includedCells[currCell]];
    else if(thresh_totals[1][includedCells[currCell]] > 0)
      calculatedRisk = threshold + 1;
    else{
      if(notInTraining == 1){
        stats.incorrectclass += totals[0][includedCells[currCell]];
        stats.correctclass += totals[1][includedCells[currCell]];
        stats.classhigh += totals[1][includedCells[currCell]];
        stats.misclasshigh += totals[0][includedCells[currCell]];
      }
      else if(notInTraining == 0){
        stats.incorrectclass += totals[1][includedCells[currCell]];
        stats.correctclass += totals[0][includedCells[currCell]];
        stats.classlow += totals[0][includedCells[currCell]];
        stats.misclasslow += totals[1][includedCells[currCell]];
      }
      continue; //skip to next cell
    }
   
    if(calculatedRisk > threshold){
      model.cell_status[includedCells[currCell]] = 1;
    }
    else if(calculatedRisk < threshold){
      model.cell_status[includedCells[currCell]] = 0;
    }
    else{
      model.cell_status[includedCells[currCell]] = tieCellValue;
    }
       
    switch(model.cell_status[includedCells[currCell]]){
      case 1:
        stats.incorrectclass += totals[0][includedCells[currCell]];
        stats.correctclass += totals[1][includedCells[currCell]];
        stats.classhigh += totals[1][includedCells[currCell]];
        stats.misclasshigh += totals[0][includedCells[currCell]];
        break;
      case 0:
        stats.incorrectclass += totals[1][includedCells[currCell]];
        stats.correctclass += totals[0][includedCells[currCell]];
        stats.classlow += totals[0][includedCells[currCell]];
        stats.misclasslow += totals[1][includedCells[currCell]];
        break;
      case -1:
        break;
    };



  }

  stats.error = error_calculator->analysis_value(stats.classhigh, stats.misclasslow,
    stats.classlow, stats.misclasshigh);

  stats.balanced_error = fitness_calculator->analysis_value(stats.classhigh, stats.misclasslow,
    stats.classlow, stats.misclasshigh);
}

}

