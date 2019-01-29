//ParallelPerm.h

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


#ifndef __PARALLELPERM_H__
#define __PARALLELPERM_H__

#include "SlaveAnalysis.h"
#include "MasterAnalysis.h"
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

/// key is fitness and value is the position in the tree
typedef Utility::RBTree<float, unsigned int, floatLT > PermModelTree;
typedef Utility::RBTreeNode<float, unsigned int, floatLT > PermModelTreeNode;

/// Conducts permutation testing
class ParallelPerm{

  public:
  
    ParallelPerm(){};
  
    /// runs permutations for slave
    void run_slave_permutations(int num_ptests, Dataset& set, Config& config);
  
    /// runs permutations for master
    void run_master_permutations(int num_ptests, Dataset& set, Config& config, int nproc,
      LogOutput& log_out);
    
    /// returns p value corresponding to the value passed
    float get_p_value(float value, int rank);
    
    float get_lr_p_value(float value);
    
    void send_seed(int seed);
    
    int receive_seed();
  
  private:
  
    int total_perms;
    
    vector<PermModelTree> permuted_models;
    
    vector<float> lr_results;
    
    permutationTypes permtype;
    
};

#endif

