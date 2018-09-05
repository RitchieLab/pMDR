//ConfigFileReader.cpp

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

#include "ConfigFileReader.h"
#include <sstream>

namespace mdr{


///
/// Constructor
///
ConfigFileReader::ConfigFileReader(){
  initialize();
}

///
/// Initalizes map for use in parsing configuration file
///
void ConfigFileReader::initialize(){
  key_map["FILENAME"] = Filename;
  key_map["ACCTHRESHOLD"] = BalAccThreshold;
  key_map["CALCTHRESHOLD"] = CalcThreshold;
  key_map["MAXLOCIVALUE"] = MaxLocusValue;
  key_map["COMBOSTART"] = ComboStart;
  key_map["COMBOEND"] = ComboEnd;
  key_map["CROSSVALINTERVAL"] = CrossValInterval;
  key_map["MODELBUILDINTERVAL"] = ModelBuildInterval;
  key_map["MODELSTOKEEP"] = ModelsToKeep;
  key_map["RANDOMSHUFFLE"] = RandomShuffle;
  key_map["RANDOMSEED"] = RandomSeed;
  key_map["NOTRECOGNIZEDRESPONSE"] = NotRecognizedResponse;
  key_map["TIECELLVALUE"] = TieCellValue;
  key_map["FORCELOCI"] = ForceLoci;
  key_map["SHOWCOMBOMISCLASS"]  = ShowComboMisclass;
  key_map["POWERMODEL"] = PowerModel;
  key_map["MISSING"] = Missing;
  key_map["PTESTS"] = PtestNum;
  key_map["PTEST-N"] = PtestNum;
  key_map["PTESTDONE"] = PTestsDone;
  key_map["ALWAYSINCLUDED"] = AlwaysIncluded;
  key_map["SHOWRULES"] = ShowRules;
  key_map["SHOWTOTALS"] = ShowTotals;
  key_map["MAPFILE"] = Mapfilename;
  key_map["MAPFILENAME"] = Mapfilename;
  key_map["TRAINBALWEIGHT"] = TrainBalWeight;
  key_map["COLHEADERS"] = Colheaders;
  key_map["SHOWPARTITION"] = BestPartition;
  key_map["BESTPARTITION"] = BestPartition;
  key_map["CHECKPOINTFILE"] = Checkpointfile;
  key_map["CHECKINTERVAL"] = CheckpointInterval;
  key_map["FITNESS"] = Fitness;
  key_map["REGRESSTEST"] = RegressTest;
  key_map["INPUTFORMAT"] = InputFormat;
  key_map["BIOFILTERNAME"] = BioFileName;
  key_map["OUTPUTALL"] = OutputAllModels;
  key_map["MODELFILE"] = ModelFileName;

  fitness_map["ERROR"] = Error;
  fitness_map["ODDSRATIO"] = OddsRatio;
  fitness_map["BALANCED"] = BalancedAccuracy;

}


///
/// Reads configuration file and returns Config object
/// @param configfilename string containing filename to read
/// @return Config
/// @throws MDRExcept on error
///
Config ConfigFileReader::read_config(string configfilename){

  ifstream infile;
  infile.open(configfilename.c_str());
  if(!infile.is_open())
    throw MDRExcept("Unable to open " + configfilename);

  Config config;
  config.models_to_keep(1);
  string keyword, value, line;

  while(!infile.eof()){
    getline(infile, line);
    stringstream ss(line);

    ss >> keyword;

    keyword = to_upper(keyword);
    if(line.find_first_of("ABCDEFGHIJKLMNOPQRSTUVWXYZ") == string::npos){
      continue;
    }

    if(keyword[0] == '/' || keyword[0] == '#')
      continue;

   // find start of rest of line
   value=line.substr(keyword.length(), line.length()-keyword.length());
   // find first part of value
   int pos = value.find_first_of("/\\-.ABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789abcdefghijklmnopqrstuvwxyz");
   if(pos > 0)
     value = value.substr(pos, value.length()-pos);

    // get value by taking the rest of the line
    if(key_map[keyword] != Filename && key_map[keyword] != Mapfilename
      && key_map[keyword] != AlwaysIncluded && key_map[keyword] != ForceLoci &&
      key_map[keyword] != Checkpointfile && key_map[keyword] != BioFileName && key_map[keyword] != ModelFileName)
      value = to_upper(value);
    add_keyword(keyword, value, config);
  }

  infile.close();

  // adjust for forced_loci presence
  if(config.force_loci_included().size() > 0){
    config.model_size_start(config.force_loci_included().size());
    config.model_size_end(config.force_loci_included().size());
  }


  return config;
}


///
/// Converts string to uppercase
/// @param convert string to convert to uppercase
/// @return converted string
///
string ConfigFileReader::to_upper(string convert){
  // convert keyword to uppercase
  unsigned int strpos, length = convert.length();
  for(strpos=0; strpos<length; strpos++){
    convert[strpos] = toupper(convert[strpos]);
  }
  return convert;
}


///
/// Sets value for the keyword passed in
/// @param keyword
/// @param value
/// @param config Config to add information to
///
void ConfigFileReader::add_keyword(string keyword, string value, Config& config){

  switch(key_map[keyword]){
    case NoMatch:
      throw MDRExcept(keyword + " is not a valid keyword in configuration file");
      break;
    case MaxLocusValue:
      config.max_locus(convert_to_int(value));
      break;
    case CalcThreshold:
      if(value.compare("SET")==0)
        config.calc_threshold_each_model(false);
      else if(value.compare("MODEL")==0)
        config.calc_threshold_each_model(true);
      else
        throw MDRExcept(keyword + " must be SET or MODEL.");
      break;
    case ComboStart:
      config.model_size_start(convert_to_int(value));
      break;
    case ComboEnd:
      config.model_size_end(convert_to_int(value));
      break;
    case CrossValInterval:
      config.num_crossval(convert_to_int(value));
      break;
    case ModelBuildInterval:
      config.build_interval(convert_to_int(value));
      break;
    case RandomSeed:
      config.random_seed(convert_to_int(value));
      break;
    case Mapfilename:
      config.mapfile(value);
      break;
    case Filename:
      config.datafile(value);
      break;
    case ModelsToKeep:
      config.models_to_keep(convert_to_int(value));
      break;
    case BestPartition:
      config.display_partitions(on_or_off(value));
      break;
    case AlwaysIncluded:
      config.loci_always_included(convert_to_string_vect(value));
      break;
    case ShowRules:
      config.display_rules(on_or_off(value));
      break;
    case ShowTotals:
      config.display_totals(on_or_off(value));
      break;
    case TrainBalWeight:

      break;
    case RandomShuffle:
      config.shuffle_data(on_or_off(value));
      break;
    case NotRecognizedResponse:
      config.not_recognized_response(convert_to_int(value));
      break;
    case Colheaders:
      config.data_with_headers(on_or_off(value));
      break;
    case TieCellValue:
      config.tie_cell_value(convert_to_int(value));
      break;
    case ForceLoci:
      config.force_loci_included(convert_to_string_vect(value));
      break;
    case PtestNum:
      config.num_ptests(convert_to_int(value));
      break;
    case ShowComboMisclass:
      config.show_all_misclass(on_or_off(value));
      break;
    case PowerModel:
      config.select_best_model(on_or_off(value));
      break;
    case PTestsDone:
      config.num_ptests_completed(convert_to_int(value));
      break;
    case Missing:
      config.missing_data_value(convert_to_int(value));
      break;
    case Checkpointfile:
      config.check_point_file(value);
      break;
    case CheckpointInterval:
      config.check_interval(convert_to_float(value));
      break;
    case RegressTest:
      config.regress_test(on_or_off(value));
      break;
    case InputFormat:
      config.input_format(to_upper(value));
      break;
    case BioFileName:
      config.set_biofilter_filename(value);
      break;
    case ModelFileName:
      config.setModelFileName(value);
      break;
    case OutputAllModels:
    	config.setOutputAllModels(on_or_off(value));
    	break;
    case Fitness:
      if(fitness_map.find(value) == fitness_map.end()){
        throw MDRExcept(value + " is not a valid value for FITNESS");
      }
      config.fitness(fitness_map[value]);
      break;
		case BalAccThreshold:
      config.set_balacc_thresh(convert_to_float(value));
      config.models_to_keep(100000);
      break;
  }

}


///
/// Splits string into tokens that will fit integer
/// vector that is returned
/// @param loci_string
/// @param vector of ints
///
vector<unsigned int> ConfigFileReader::convert_to_int_vect(string loci_string){
  stringstream ss(loci_string);
  unsigned int int_val;
  vector<unsigned int> loci;
  do{
    ss >> int_val;
    loci.push_back(int_val);
  }while(!ss.eof());

  if(loci.size() > 0){
    if(loci[loci.size()-1] == loci[loci.size()-2])
      loci.pop_back();
  }

  return loci;
}


///
/// Splits string into tokens and adds to
/// string vector that is returned
/// @param snp_string
/// @returns vector of strings
///
vector<string> ConfigFileReader::convert_to_string_vect(string snp_string){

  stringstream ss(snp_string);
  string snp_name;
  vector<string> names;

  do{
    ss >> snp_name;
    names.push_back(snp_name);
  }while(!ss.eof());

  if(names.size() > 1){
    if(names[names.size()-1].compare(names[names.size()-2])==0)
      names.pop_back();
  }

  return names;
}


///
/// Returns true when string matches 'ON' and false when not
/// @param on_off
/// @return bool
///
bool ConfigFileReader::on_or_off(string on_off){
  if(on_off.compare("ON")==0)
    return true;
  return false;
}

///
/// Converts string into float and returns
/// @param value string with number
/// @return new float
///
float ConfigFileReader::convert_to_float(string value){
  stringstream ss(value);
  float float_val;
  ss >> float_val;
  return float_val;
}


///
/// Converts string into integer and returns it
/// @param value string with number
/// @returns int
///
int ConfigFileReader::convert_to_int(string value){
  stringstream ss(value);
  int int_val;
  ss >> int_val;
  return int_val;
}


}

