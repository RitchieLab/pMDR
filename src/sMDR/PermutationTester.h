// PermutationTester.h

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


#ifndef __PERMUTATIONTESTER_H__
#define __PERMUTATIONTESTER_H__

#include "sAnalysis.h"
#include <LRInteract.h>

using namespace std;
using namespace mdr;


// sort floats in ascending order
class sort_ascending{
  public:
    bool operator()(float x, float y) const{
      return x > y;
    }
};


class PermutationTester{

  public:
    
    /// returns p value corresponding to the value passed
    float get_p_value(float value);
    
    /// runs permutations and creates tree from best results
    void run_permutations(unsigned int num_permutations, Dataset& set, Config& config, LogOutput& log_out);
    
    /// returns p value for the lr value
    float get_lr_p_value(float value);
    
  private:
  
    int total_perms;
    
    ModelTree permuted_models;
    
    vector<float> lr_results;
  
};


#endif

