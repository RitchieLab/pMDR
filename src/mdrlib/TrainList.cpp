//TrainList.cpp

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

#include "TrainList.h"

#include <iostream>
using namespace std;

namespace mdr{

///
/// Sets number of crossvalidatin intervals and creates
/// a tree for each crossvalidation interval and
/// size of model
/// @param cv
/// @param max_loc_comb Maximum size of model in number of loci
/// @param model_max Maximum number of models to keep
/// @return
///
void TrainList::initialize_list(int cv, int max_loc_comb, int model_max){

  vector<ModelTree*> temp_vector;
  temp_vector.resize(max_loc_comb+1);
  maximum_models = model_max;

  for(int curr_cv=0; curr_cv< cv; curr_cv++){
    result_trees.push_back(temp_vector);
    for(int order=0; order<=max_loc_comb; order++){
      result_trees[curr_cv][order] = new ModelTree(model_max);
    }
  }

}


///
/// Destructor must destroy trees
///
TrainList::~TrainList(){
  for(unsigned int curr_cv=0; curr_cv < result_trees.size(); curr_cv++){
    for(unsigned int order=0; order < result_trees[curr_cv].size(); order++){
      delete result_trees[curr_cv][order];
    }
  }

}

///
/// Clears trees
///
void TrainList::clear(){
  for(unsigned int curr_cv=0; curr_cv < result_trees.size(); curr_cv++){
    for(unsigned int order=0; order < result_trees[curr_cv].size(); order++){
      result_trees[curr_cv][order]->Clear();
    }
  }
}



///
/// Insert the model if better than worst model in tree
/// at indicated position, or if tree is not full yet
/// @param fitness
/// @param loci vector<int> of loci
/// @param cv Crossvalidation interval for model
/// @return
////
void TrainList::Insert(float fitness, Model& mod, int cv){
  result_trees[cv][mod.combination.size()]->Add(fitness, mod);
}


}

