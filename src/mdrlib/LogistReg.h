// LogistReg.h

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

#ifndef __LOGISTREG_H__
#define __LOGISTREG_H__

#include "FlatIndex.h"

namespace mdr{

///
/// Calculates logistic regression on a model.
/// Can run on multiple locus combinations.
/// Can calculate the likelihood ratio (1 degree of freedom for interaction models).
///

/// Runs logistic regression
class LogistReg{
  
  public:
  
    LogistReg();

    float run_lr(vector<vector<unsigned int> >& data, vector<unsigned int>& status, bool any_missing,
      int max_data_value, bool include_full_interaction=true);
    
  private:
    
    void initialize();
  
    void summarize_data(vector<vector<unsigned int> >& dataset, vector<unsigned int>& status,
      int max_data_value, bool any_missing, bool include_full_interaction);
    
    void initialize_summary(unsigned int currModelSize, bool include_full_interaction);
    
    void zero_summary(unsigned int array_size, unsigned int model_size, bool include_full_interaction=true);
    
    void initialize_interactions();
    
    unsigned int ix(int j,int k,int nCols);
 
    FlatIndex indexConverter;
 
    vector<vector<vector<unsigned int> > > interaction_lists;
    vector<vector< int> > includedIndexes;
    vector<vector<int> > summary_data;
    double PiD2;
    bool interaction_included;
    unsigned int modelSize, LociComboLimit, LociComboMin, maxIterations, defaultComboInterval, maxLocusValue,
      missingValue, numberGenosPerLocus;


};
}
#endif
