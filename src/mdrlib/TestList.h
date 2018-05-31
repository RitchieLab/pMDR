//TestList.h

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

#ifndef __TESTLIST_H__
#define __TESTLIST_H__

#include "TrainList.h"

using namespace std;

namespace mdr{


/// Holds results for the best models from the Test List
class TestList{

  public:

    ~TestList();

    vector<vector<ModelTree *> > result_trees;
    
    unsigned int num_crossvalidation_intervals(){return result_trees.size();}
    
    void Insert(float fitness, vector<int>& loci, int cv);
    
    void initialize_list(int cv, int max_loc_comb);
    
  private:
  

};



}

#endif

