//Config.h

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


#ifndef __CONFIG_H__
#define __CONFIG_H__

#include <string>
#include <vector>
#include "MDRExcept.h"
#include "fitness.h"

using namespace std;

namespace mdr{


class Config{

  public:

    Config();

    bool calc_threshold_each_model(){return calc_threshold;}
    void calc_threshold_each_model(bool ct){calc_threshold = ct;}

    bool display_partitions(){return bestpartition;}
    void display_partitions(bool part){bestpartition = part;}

    bool display_rules(){return show_rules;}
    void display_rules(bool show){show_rules = show;}

    bool display_totals(){return show_totals;}
    void display_totals(bool show){show_totals = show;}

    bool shuffle_data(){return random_shuffle;}
    void shuffle_data(bool shuffle){random_shuffle = shuffle;}

    bool data_with_headers(){return col_headers;}
    void data_with_headers(bool headers){col_headers = headers;}

    bool select_best_model(){return power_model;}
    void select_best_model(bool power){power_model = power;}

    bool show_all_misclass(){return showcombo_misclass;}
    void show_all_misclass(bool misclass){showcombo_misclass = misclass;}

    string mapfile(){return mapfilename;}
    void mapfile(string filename){mapfilename = filename;}

    string datafile(){return datafilename;}
    void datafile(string filename){datafilename = filename;}

    string check_point_file(){return checkpointfilename;}
    void check_point_file(string filename){checkpointfilename = filename;}

    float check_interval(){return checkpoint_interval;}
    void check_interval(float interval){checkpoint_interval = interval;}

    int model_size_start();
    void model_size_start(int start);

    int model_size_end();
    void model_size_end(int end);

    int num_crossval();
    void num_crossval(int cv);

    int build_interval();
    void build_interval(int num);

    int random_seed();
    void random_seed(int seed);

    int models_to_keep();
    void models_to_keep(int num);

    int not_recognized_response();
    void not_recognized_response(int value);

    int tie_cell_value();
    void tie_cell_value(int value);

    int num_ptests();
    void num_ptests(int num);

    int num_ptests_completed();
    void num_ptests_completed(int value);

    int missing_data_value();
    void missing_data_value(int value);

    int max_locus();
    void max_locus(int max);

    vector<string> loci_always_included();
    void loci_always_included(vector<string> included);

    vector<string> force_loci_included();
    void force_loci_included(vector<string> loci);

    void basename(string base){basefilename=base;}
    string basename(){return basefilename;}

    fitnessTypes fitness(){ return fitness_type;}
    void fitness(fitnessTypes fit){fitness_type = fit;}

    bool regress_test(){return regress_test_on;}
    void regress_test(bool do_test){regress_test_on = do_test;}

    string input_format(){return input_format_string;}
    void input_format(string i_format){input_format_string = i_format;}

    string get_biofilter_filename(){return biofilename;}
    void set_biofilter_filename(string name){biofilename = name;}

    void set_balacc_thresh(float num);
    float get_balacc_thresh(){return balaccthresh;}
    
    void setOutputAllModels(bool tf){flush_models=tf;}
    bool getOutputAllModels(){return flush_models;}

    void check_parameters();

  private:
    void initialize();

    string datafilename, mapfilename, checkpointfilename, basefilename, input_format_string,
      biofilename;

    bool calc_threshold, bestpartition, show_rules, show_totals, random_shuffle, col_headers,
      power_model, showcombo_misclass, regress_test_on, flush_models;

    int combo_start, combo_end, crossval_interval, modelbuild_interval, randomseed, models_kept,
      notrecognized_response, tiecell_value, ptests, ptests_done, missing_value, max_locus_value;

    float checkpoint_interval, balaccthresh;

    fitnessTypes fitness_type;

    vector<string> always_included, force_loci;
};

}

#endif

