// MDR.h

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

#ifndef __MDR_H__
#define __MDR_H__

using namespace std;
#include "Model.h"
#include "Dataset.h"
#include "fitness.h"
#include "StatCreator.h"

namespace mdr{

class MDR{

  public:
  
    MDR();
    
    ~MDR();
    
    MDR(int crossval);
    
    MDR(int crossval, fitnessTypes type);
    
    void initialize(int crossval, fitnessTypes type);
    
    void test_model(vector<Model>& models, Dataset& set);
    
    void test_model(Model& curr_model, Dataset& set, int cv);
    
    void set_cv(int CrossVal);
    
    void set_tiecell_value(int value){tieCellValue = value;}
    
    void set_not_in_training(int value){notInTraining = value;}
    
    void full_model(Model& curr_model, Dataset& set);
    
    string fitness_name(){return fitness_calculator->get_analysis_name();}
    
    float display_fitness(float model_fitness){
      return fitness_calculator->convert_to_display(model_fitness);
    }
    
    Stat* get_calculator(){return fitness_calculator;}
    
    void test_model_complete(vector<Model>& models, Dataset& set);
    
  private:
  
    void fill_cv_totals(vector<vector<vector<int> > >& totals, Dataset& MasterData,
      vector<unsigned int>& loci);
  
    int set_cv_totals(vector<vector<vector<int> > >& totals, 
      vector<vector<int> >& cv_totals, int cv, int total_cv);
  
    void LinearMDR(vector<int>& current_set, Dataset& MasterData,
      vector<unsigned int>& loci, vector<vector<int> >& totals);
  
    void calculate_err(vector<vector<int> >& thresh_totals, 
        vector<vector<int> >& totals, float threshold, vector<int>& includedCells,
        model_stats& stats, Model& model);
        
    void check_included(Dataset& set);
    
    void set_calculators(fitnessTypes type);

    vector<vector<int> > cv_totals, zero_cv;
    vector<vector<vector<int> > > totals, zero_totals;
    vector<int> included_cells;
    vector<vector<int> > includedIndexes;

    int total_cv, tieCellValue, notInTraining;
    unsigned int last_model_size;
    
    StatCreator stat_creator;
    Stat * error_calculator;
    Stat * fitness_calculator;

};



}

#endif

