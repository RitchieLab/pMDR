//Analysis.h

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

#ifndef __ANALYSIS_H__
#define __ANALYSIS_H__

#include "TrainList.h"
#include "MDR.h"
#include "Config.h"
#include "LogOutput.h"
#include "ComboGenerator.h"
#include "OutputWriter.h"
#include "Checkpoint.h"

using namespace std;

namespace mdr{

class Analysis{

  public:

    Analysis(){generator=NULL;output_all=false;}

    virtual void set_combination_limits(int total_loci, int min_size, int max_size);

    virtual ~Analysis();

    vector<Model> get_best_models(Dataset& set, Config& config);

    void set_max_models(int max_training_models){max_train=max_training_models;}

    void set_parameters(Config& config, Dataset& set, bool log_all);

    void output_models(Dataset& set, Config& config, string base_out);

    void output_best(vector<Model> best_models, Dataset& set);

    Model get_single_best(vector<Model> best_models, Config& config);

    void output_single_best(Model& best_model, Dataset& set);

    void output_p_values(vector<Model>& models, Dataset& set, int p_tests);

    void output_lr_p_values(vector<Model>& models, Dataset& set, int p_tests);

    Stat* get_calculator(){return mdr.get_calculator();}

    void log_all_training(bool log_option){log_training = log_option;}
    bool log_all_training(){return log_training;}

    void set_generator_type(string biofilename, int idSize);
    
    void set_flushfn(string fn){flush_fn=fn;output_all=true;}

  protected:

    void test_single_model(vector<unsigned int>& loci, vector<Model>& models, Dataset& set,
      bool calc_prediction=false);

    void add_included_snps(vector<vector<unsigned int> >& combinations);

    double read_checkpoint(string checkpointfile);

    void write_checkpoint(float interval, string basefile, double completed_mods);
    
	void output_all_models(vector<Model>& models,std::vector<unsigned int>& combination);

    OutputWriter out_writer;
    Checkpoint check_pointer;
    MDR mdr;
    ComboGenerator* generator;
    vector<unsigned int> forced_loci;
    unsigned int total_cv;
    int max_train;
    unsigned int current_model_size;
    TrainList training_results;
    vector<unsigned int> always_included;
    bool log_training, output_all;

    float check_interval, error_threshold;
    string checkpointfile, basecheckpoint, flush_fn;

};

}

#endif

