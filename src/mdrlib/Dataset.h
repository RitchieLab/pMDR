//Dataset.h

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

#ifndef __DATASET_H__
#define __DATASET_H__

#include <vector>
#include <string>
#include <map>
#include "MDRExcept.h"
#include "FlatIndex.h"

using namespace std;

namespace mdr{


/// Holds data for analysis and contains functionality for permutation
class Dataset{
  
  public:
    Dataset();
    
    ~Dataset();
    
    char** data;  
    
    vector<vector<int> > training, testing;  // The testing and training sets
    
    vector<string> snp_names;
    vector<string> chromosomes;
    vector<unsigned int> positions;
    
    map<string, unsigned int> snp_index_map;
  
    void add_snp_name(string name);
    
    void add_pos(unsigned int pos){positions.push_back(pos);}
    void add_chrom(string chr){chromosomes.push_back(chr);}
  
    float determine_threshold(vector<unsigned int> & lociCombination);
    
    void initialize_data(int nInds, int nLoci);
    
    void initialize_data();
    
    void permute_set(int seed);
    void split_data();
    
    void set_random_shuffle(bool shuffle_on){random_shuffle_on = shuffle_on;}
    void set_max_genotype_value(int max_value){max_genotype_value = max_value;}
    void set_num_inds(int nInds){num_inds=nInds;}
    void set_num_loci(int nLoci){num_loci=nLoci;}
    void set_crossval_interval(int cv);
    void set_permodelthreshold(bool permodel){permodelthreshold = permodel;}
    void set_entire_threshold(float thresh){set_threshold = thresh;}
    void set_contains_missing(bool missing){contains_missing = missing;}
    void set_missing_value(int miss_val){missing_value=miss_val;}
    void set_shuffle_seed(int seed){randseed=seed;}
    
    bool any_missing_data(){return contains_missing;}
    int get_num_loci(){return num_loci;}
    int get_num_inds(){return num_inds;}
    int get_max_locus_value(){return max_genotype_value;}
    int get_missing_value(){return missing_value;}
    int get_crossval_interval(){return crossval_interval;}
    int get_training_crossval_total(int cv){return training[cv].size();}
    int get_testing_crossval_total(int cv){return testing[cv].size();}
    float get_set_threshold(){return set_threshold;}
    
    int locus_name_index(string snp_name);
    
    void reference_snp_names();
    
    void initialize_converter(int combo_min, int combo_max);
    
    inline FlatIndex& converter(){return indexConverter;}
    
  private:

    FlatIndex indexConverter;
    
    void initialize();
    void shuffle_handles(vector<int>& handles);
    void bin_data(vector<vector<int> >& bin, vector<int>& inds);
    
    float set_threshold;
    bool permodelthreshold, random_shuffle_on, contains_missing;
    int max_genotype_value, num_inds, num_loci, crossval_interval, 
	    missing_value, randseed;
  
};

}

#endif

