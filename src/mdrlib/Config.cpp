//Config.cpp

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

#include "Config.h"

namespace mdr{


///
/// Constructor
///
Config::Config(){
  initialize();
}

///
/// Initializes variables
///
void Config::initialize(){
	balaccthresh = 100.0;
  combo_start = 1;
  combo_end = 1;
  crossval_interval = 5;
  modelbuild_interval = 10000;
  randomseed = 9245;
  notrecognized_response = 1;
  tiecell_value = 1;
  ptests = 0;
  ptests_done = 0;
  missing_value = -1;
  calc_threshold = true;
  bestpartition = false;
  show_rules = false;
  show_totals = false;
  random_shuffle = true;
  col_headers = false;
  power_model = true;
  showcombo_misclass = false;
  regress_test_on = false;
  mapfilename = "";
  checkpointfilename = "";
  checkpoint_interval = 0.0;
  input_format_string = "TEXT";
  fitness_type = BalancedAccuracy;
  biofilename = "";
  flush_models = false;
  perm_type = Combined;
}


int Config::model_size_start(){
  return combo_start;
}

///
/// Set starting model size
/// @param start starting combo
/// @throws MDRExcept when start not valid number
///
void Config::model_size_start(int start){
  if(start < 1)
    throw MDRExcept("Combo start must be a positive number");
  combo_start = start;
}


int Config::model_size_end(){
  return combo_end;
}

///
/// Sets end model size
/// @param end starting combo
/// @throws MDRExcept when end not valid number
///
void Config::model_size_end(int end){
  if(end < 1)
    throw MDRExcept("Combo end must be a positive number");
  combo_end = end;
}


int Config::num_crossval(){
  return crossval_interval;
}

///
/// Sets number of crossvalidations
/// @param cv number of crossvalidations
/// @throws MDRExcept when cv not valid number
///
void Config::num_crossval(int cv){
  if(cv < 1)
    throw MDRExcept("CROSSVALINTERVAL must be a positive number");
  crossval_interval = cv;
}

int Config::build_interval(){
  return modelbuild_interval;
}

///
/// Sets number of models to generate at a time
/// @param num number of models that are evaluated at a time in combination generator
/// @throws MDRExcept when num not valid number
///
void Config::build_interval(int num){
  if(num < 1)
    throw MDRExcept("MODELBUILDINTERVAL must be a positive number");
  modelbuild_interval = num;
}

int Config::random_seed(){
  return randomseed;
}

///
/// Sets value for random seed
/// @param seed
/// @throws MDRExcept when seed not valid number
///
void Config::random_seed(int seed){
  if(seed < 1)
    throw MDRExcept("RANDOMSEED must be a positive number");
  randomseed = seed;
}

int Config::models_to_keep(){
  return models_kept;
}

///
/// Sets value for number of models to keep
/// @param num
/// @throws MDRExcept when num not valid number
///
void Config::models_to_keep(int num){
  if(num < 1)
    throw MDRExcept("MODELSTOKEEP must be a positive number");
  models_kept = num;
}

void Config::set_balacc_thresh(float num){
	if(num < 0 || num > 1.0)
		throw MDRExcept("ACCTHRESHOLD must be between 0 and 1.0");
	// store as balanced error
	balaccthresh = 100.0 - num*100;
}

int Config::not_recognized_response(){
  return notrecognized_response;
}

///
///  When a pattern is encountered in testing that was not classified in
/// training, how should this indiviual treated?
/// -1: Unknown  0: Unaffected  1: Affected
/// @param value
/// @throws MDRExcept when value not valid number
///
void Config::not_recognized_response(int value){
  if(value < -1 || value > 1)
    throw MDRExcept("NOTRECOGNIZEDRESPONSE must be -1, 1, or 0");
  notrecognized_response = value;
}


int Config::tie_cell_value(){
  return tiecell_value;
}

///
/// When building a model, how should cells with equal affecteds and
/// unaffecteds be classified?
/// -1: Unknown  0: Unaffected  1: Affected
/// @param value
/// @throws MDRExcept when value is not -1, 0, or 1
///
void Config::tie_cell_value(int value){
  if(value < -1 || value > 1)
    throw MDRExcept("TIECELLVALUE must be -1, 1, or 0");
  tiecell_value = value;
}

int Config::num_ptests(){
  return ptests;
}

///
/// Set number of permutation tests
/// @param num
/// @throws MDRExcept when num is less than one
///
void Config::num_ptests(int num){
  if(num < 0)
    throw MDRExcept("PTESTS must be greater than or equal to 0");
  ptests = num;
}


int Config::num_ptests_completed(){
  return ptests_done;
}

///
/// Set number of permutation tests already completed.
/// @param value
/// @throws MDRExcept when value is less than one
///
void Config::num_ptests_completed(int value){
  if(value < 0)
    throw MDRExcept("PTESTDONE must be greater than 0");
  ptests_done = value;
}

int Config::missing_data_value(){
  return missing_value;
}

///
/// Set value for identifying missing data in data file
/// @param value
/// @throws MDRExcept when number is 0
///
void Config::missing_data_value(int value){
  if(value == 0)
    throw MDRExcept("MISSING cannot be 0");
  missing_value = value;
}

vector<string> Config::loci_always_included(){
  return always_included;
}

///
/// Set loci which will always be included in the models
/// @param included vector containing SNP names for loci to include
/// @throws MDRExcept when no loci in vector
///
void Config::loci_always_included(vector<string> included){
//   if(included.empty())
//     throw MDRExcept("ALWAYSINCLUDED must have at least one locus");
  always_included = included;
}

vector<string> Config::force_loci_included(){
  return force_loci;
}

///
/// Set loci which will be only model tested
/// @param included vector containing indexes for loci to include
/// @throws MDRExcept when no loci in vector
///
void Config::force_loci_included(vector<string> loci){
  force_loci = loci;
}



int Config::max_locus(){
  return max_locus_value;
}


///
/// Set maximum value for a genotype in data file
/// @param max
/// @throws MDRExcept when max less than 0
///
void Config::max_locus(int max){
  if(max < 1)
    MDRExcept("MAXLOCIVALUE must be greater than 0");
  max_locus_value = max;
}

///
/// Checks values within the object to see if any violate
/// the bounds of that parameter
/// @throws MDRExcept for any violation
///
void Config::check_parameters(){

  if(combo_end < combo_start)
    throw MDRExcept("COMBOEND must be greater than or equal to COMBOSTART");
  if(combo_start < 2 && regress_test_on)
    throw MDRExcept("COMBOSTART cannot be less than 2 when REGRESSTEST is ON");

}

}

