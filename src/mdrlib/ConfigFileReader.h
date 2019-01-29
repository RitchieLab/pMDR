//ConfigFileReader.h

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

#ifndef __CONFIGFILEREADER_H__
#define __CONFIGFILEREADER_H__


using namespace std;

#include "Config.h"
#include "MDRExcept.h"
#include <map>
#include <fstream>

namespace mdr{


class ConfigFileReader{

  public:

    ConfigFileReader();

    Config read_config(string filename);

  private:

    void initialize();

    void add_keyword(string keyword, string value, Config& config);

    vector<unsigned int> convert_to_int_vect(string loci_string);

    vector<string> convert_to_string_vect(string snp_string);

    float convert_to_float(string value);

    bool on_or_off(string on_off);

    int convert_to_int(string value);

    string to_upper(string convert);

    string mapfilename, filename;

    enum configKeyWords{
      NoMatch,
      MaxLocusValue,
      CalcThreshold,
      ComboStart,
      ComboEnd,
      CrossValInterval,
      ModelBuildInterval,
      RandomSeed,
      Mapfilename,
      Filename,
      ModelsToKeep,
      BestPartition,
      AlwaysIncluded,
      ShowRules,
      ShowTotals,
      TrainBalWeight,
      RandomShuffle,
      NotRecognizedResponse,
      Colheaders,
      TieCellValue,
      ForceLoci,
      PtestNum,
      ShowComboMisclass,
      PowerModel,
      PTestsDone,
      Missing,
      CheckpointInterval,
      Checkpointfile,
      Fitness,
      RegressTest,
      InputFormat,
      BioFileName,
      BalAccThreshold,
      OutputAllModels,
      ModelFileName,
      PermutationMethod
    };

    map<string, configKeyWords> key_map;
    map<string, fitnessTypes> fitness_map;
    map<string, permutationTypes> permutation_map;

};


}

#endif

