//TrainList.h

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

#ifndef __TRAINLIST_H__
#define __TRAINLIST_H__

#include "rbtree.h"
#include <vector>
#include "Model.h"

using namespace std;

namespace mdr{

// can change ordering by changing this
struct floatLT {
        int operator()(const float l, const float r) const {
                if (l>r) return -1;
                if (l<r) return 1;
                return 0;
        }
};

typedef Utility::RBTree<float, Model, floatLT > ModelTree;
typedef Utility::RBTreeNode<float, Model, floatLT > ModelTreeNode;

/// Holds results for the best models in the train list
class TrainList{

  public:

    ~TrainList();

    vector<vector<ModelTree *> > result_trees;
    
    unsigned int num_crossvalidation_intervals(){return result_trees.size();}
    
    void Insert(float fitness, Model& mod, int cv);
    
    void initialize_list(int cv, int max_loc_comb, int model_max);
    
    int max_models_kept(){return maximum_models;}
    
    void clear();
    
  private:
    int maximum_models;

};



}

#endif

