//Dataset.cpp

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

#include "Dataset.h"
#include <iostream>
#include <sstream>
#include <cstdlib>
using namespace std;

namespace mdr{

Dataset::Dataset(){
  initialize();
}


///
/// Initialize variables
///
void Dataset::initialize(){
  data = NULL;
  set_threshold = true;
  permodelthreshold = false;
  random_shuffle_on = true;
  contains_missing = false;
  crossval_interval = 1;
  max_genotype_value = 2;
  num_inds = 0;
  num_loci = 0;
  snp_names.push_back("status");
  positions.push_back(0);
  chromosomes.push_back("");
}

///
/// destructor
///
Dataset::~Dataset(){
  if(data != NULL){
    for(int i=0; i < num_inds; i++){
      delete [] data[i];
    }
    delete [] data;
  }
}

///
/// Initializes data array
/// @param nInds number of individuals in set
/// @param nLoci number of loci in set
/// @return
///
void Dataset::initialize_data(int nInds, int nLoci){
 
  num_loci = nLoci;
  num_inds = nInds;
 
  initialize_data();
}


///
/// Initializes data array
/// @param nInds number of individuals in set
/// @param nLoci number of loci in set
/// @return
///
void Dataset::initialize_data(){

  // Create the data array 
  // First index is reserved for the status of the individual
  // so need to add one to number of loci in set
  data = new char * [num_inds];
  for (int i=0; i < num_inds; i++){
    data[i] = new char[num_loci+1];
  }
   
}



///
/// Returns threshold for the loci combination.  When using same 
/// threshold for every model return the overall threshold in set.
/// @param lociCombination vector containing indexes for the 
/// loci that make up the model
/// @return threshold as float
/// 
float Dataset::determine_threshold(vector<unsigned int> & lociCombination){

  if(!permodelthreshold){
    return set_threshold;
  }

  int codedmissing = max_genotype_value + 1;
 
  vector<int> missing(2,0);
  vector<unsigned int>::iterator j;
  bool miss;
  
  for(int i=0; i<num_inds; i++){
    miss=false;
    for(j=lociCombination.begin(); j != lociCombination.end();
      j++){
      if(data[i][*j] == codedmissing){
        miss = true;
        break;
      }
    }
    if(!miss)
      missing[data[i][0]]++;
  }
  
  if(missing[0] == 0) // when no unaffected, threshold won't be used
    return 1;
  else
    return float(missing[1]) / missing[0];  
}



///
/// Randomize the affected status of individuals for permutation testing
/// @param  random number generator seed
/// @return 
///
void Dataset::permute_set(int seed){

  srand(seed);
  int affIndex;
  int numaffected=0;
  vector<int> indIndex(num_inds, 0);
  
  // zero out all status
  for(int i=0; i < num_inds; i++){
    if(data[i][0]){
      numaffected++;
    }
    data[i][0]=0;
    indIndex[i] = i;
  }

  int tempInt;

  // randomly select individuals to be affected
  for(int i=0; i<numaffected; i++){
    int numIndsLeft = num_inds-i;
    affIndex = rand() % numIndsLeft;
    data[indIndex[affIndex]][0]=1;
    tempInt = indIndex[affIndex];
    indIndex[affIndex] = indIndex[numIndsLeft-1];
    indIndex[numIndsLeft-1] = indIndex[affIndex];
  }
  
  split_data();
  
}


///
/// Split the data into training and testing sets for Cross-Validation
/// This function works by spliting the data into affected and unaffected sets, then added
/// these to training and testing sets in proportion
/// @param
/// @return
///
void Dataset::split_data(){ 

  // Split the Dataset into affecteds and unaffecteds
  vector<int> affecteds, unaffecteds, paired, unpaired;

  for(int i = 0; i < num_inds; i++){
    if(data[i][0] == 1)
      affecteds.push_back(i);
    else
      unaffecteds.push_back(i);
  }

  if(random_shuffle_on){ // If the RANDOMSHUFFLE was ON in the config file
    // Shuffle the Affecteds and Unaffecteds, randomizing the data
    shuffle_handles(affecteds);
    shuffle_handles(unaffecteds);
  }

  std::vector< std::vector<int> > bin;
  bin.resize(crossval_interval);

  if(crossval_interval == 1){
    training[0].clear();
    testing[0].clear();
    for(unsigned int i=0; i < affecteds.size(); i++){
      training[0].push_back(affecteds[i]);
    }
    for(unsigned int i=0; i < unaffecteds.size(); i++){
      training[0].push_back(unaffecteds[i]);
    }
  }
  else{
    if(random_shuffle_on){
      for(unsigned int i=0; i < affecteds.size(); i++){
        bin[i % crossval_interval].push_back(affecteds[i]);
      }
      for(unsigned int i=0; i < unaffecteds.size(); i++){
        bin[i % crossval_interval].push_back(unaffecteds[i]);
      }
    }
    else{
      bin_data(bin,affecteds);
      bin_data(bin,unaffecteds);
    }

    for(int testingcv=0; testingcv < crossval_interval; testingcv++){
      testing[testingcv].clear();
      training[testingcv].clear();
      for(int cv=0; cv < crossval_interval; cv++){
         if(cv == testingcv)
          testing[testingcv] = bin[cv];
        else{
          for(unsigned int i=0; i < bin[cv].size(); i++)
            training[testingcv].push_back(bin[cv][i]);
         }
      }
    }
   }
   
}

///
/// Bins data
///
void Dataset::bin_data(vector<vector<int> >& bin, vector<int>& inds){
  int binsize=int(inds.size()/crossval_interval);
  int remain=inds.size()%crossval_interval;
  int i=0;
  for(int currbin=0; currbin < crossval_interval; currbin++){
    for(int count=0; count<binsize; count++){
      bin[currbin].push_back(inds[i++]);
    }
    if(currbin < remain)
      bin[currbin].push_back(inds[i++]);
  }
}

///
/// Shuffles an array of handles by placing each item in a random position
/// @param handles - an STL vector to be randomly shuffled
/// @return
///
void Dataset::shuffle_handles(vector<int>& handles){

  int r, temp;

  srand(randseed);
  for(unsigned int i=0; i < handles.size(); i++){
    r = rand() % handles.size();
    temp = handles[i];
    handles[i] = handles[r];
    handles[r] = temp;
  }

}

///
/// Sets number of crossvalidation intervals and resizes the
/// training and testing vectors
/// @param cv
/// @return
///
void Dataset::set_crossval_interval(int cv){
  crossval_interval=cv;
  
  testing.resize(cv);
  training.resize(cv);
  split_data();
}

///
/// Fills SNP names with place holders only containing
/// reference numbers for the SNPs
/// 
void Dataset::reference_snp_names(){
  for(int i=1; i<=num_loci; i++){
    stringstream ss;
    ss << i;
    snp_names.push_back(ss.str());
    snp_index_map[ss.str()] = snp_names.size()-1;
    positions.push_back(i);
    chromosomes.push_back("1");
  }
}


///
/// Adds name of snp to vector and to map for reference later
/// @param name string
///
void Dataset::add_snp_name(string name){
  snp_names.push_back(name);
  snp_index_map[name] = snp_names.size()-1;
}


///
/// Do simple linear search to find 
/// the index of the SNP name in the snp_names vector
/// @param string SNP name
/// @returns index
/// @throws MDRExcept when snp_name is not found
///
int Dataset::locus_name_index(string snp_name){
  unsigned int num_snps = snp_names.size();
  for(unsigned int i=0; i<num_snps; i++){
    if(snp_names[i].compare(snp_name)==0){
      return i;
    }
  }
  throw MDRExcept("Unable to find SNP with name " + snp_name);
}


///
/// Initializes indexconverter that maps a multi-dimensional
/// array to a one-dimensional one
/// @param int combo_min Smallest combination size
/// @param int combo_max Largest combination size
/// 
void Dataset::initialize_converter(int combo_min, int combo_max){
  if(any_missing_data())
    indexConverter.set_genos_per_locus(get_max_locus_value()+2);
  else
    indexConverter.set_genos_per_locus(get_max_locus_value()+1);

  indexConverter.set_included_indexes(combo_min, combo_max,
    any_missing_data(), get_max_locus_value()+1);
}


}

