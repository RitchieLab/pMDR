//BioComboGenerator.cpp
#include "BioComboGenerator.h"
#include "BioBinFileAccess.h"
#include "BioTextFileAccess.h"

namespace mdr{

///
/// Constructor initializes parameters
///
BioComboGenerator::BioComboGenerator(){
  initialize();
}

///
/// Destructor
///
BioComboGenerator::~BioComboGenerator(){
  FreeMemory();
}


///
/// Frees any memory 
///
void BioComboGenerator::FreeMemory(){
  if(_fileReader != NULL){
    delete _fileReader;
    _fileReader = NULL;
  }
}


///
/// Passes configuration parameters for this generator.  In this case all
/// it needs is the name of the biofilter file.
/// @param config map of strings
///
void BioComboGenerator::PassInitConfig(map<string, string> config){
  
  string biofilename, keyword;
  keyword = config.begin()->first;
  
  if(keyword.compare("BIOFILENAME") == 0){
    biofilename = config.begin()->second;
  }
  
  FreeMemory();
  
  // Need to check if file is binary and then create appropriate reader
  if(BioFileAccess::FileBinary(biofilename)){
    _fileReader = new BioBinFileAccess(biofilename);
  }
  else{
    _fileReader = new BioTextFileAccess(biofilename);
  }
  
  _total_models_in_file = _fileReader->GetRecordCount();
  
  _fileReader->Open();
  
}


///
/// Initializes variables for generator
/// @return none
///
void BioComboGenerator::initialize(){
  _current_model = 0;
  _fileReader = NULL;
  _current_index = 0;
//  idSize=11;
}


///
/// Creates initial state.  Resets current model to zero.
///
void BioComboGenerator::initialize_state(){
  _current_model = 0;
}


///
/// Returns number of models that would be returned by the AdvanceParameters
///
int BioComboGenerator::get_count_generated(){
  if(_current_index < _total_models_in_file)
    return ComboInterval;
  else
    return ComboInterval - (_current_index - _total_models_in_file);
}


///
/// Generates combinations in amount passed.  Doesn't change original interval.
/// @param new_interval Combinations to create
/// 
bool BioComboGenerator::GenerateCombinations(int new_interval){
  int old_interval = getComboInterval();
  SetComboInterval(new_interval);
  bool done = GenerateCombinations();
  SetComboInterval(old_interval);
  return done;
}

///
///  Reads combinations from the file.
/// @return Return true when all combinations completed
///
bool BioComboGenerator::GenerateCombinations(){
  
  bool done =false;
  ComboList.clear();
  
  vector<vector<string> > model_strings;
  
  _fileReader->GetModels(_current_index, model_strings, ComboInterval, idSize);
  
  _current_index += model_strings.size();

  if(_current_index >= _total_models_in_file)
    done = true;
  
  vector<vector<string> >::iterator outerIter;
  vector<string>::iterator innerIter;
  
  vector<unsigned int> model;
  map<string, unsigned int>::iterator snpMapIter;
  bool include;
  
  for(outerIter = model_strings.begin(); outerIter != model_strings.end(); outerIter++){
    model.clear();
    include = true;
    for(innerIter = outerIter->begin(); innerIter != outerIter->end(); innerIter++){
      snpMapIter = dataset->snp_index_map.find(*innerIter);
      if(snpMapIter == dataset->snp_index_map.end()){
        include = false;
        break;
      }
      else{
        model.push_back(snpMapIter->second);
      }
    }
    if(include){
      ComboList.push_back(model);
    }
    else{
      cout << "Missing Model";
      for(vector<string>::iterator missingIter=outerIter->begin(); missingIter != outerIter->end();
        missingIter++){
      cout << "\t" << *missingIter;
      }
      cout << endl;
    }
  }
  return done;
}

///
/// Advances index 
/// @param when all models finished return true
///
bool BioComboGenerator::AdvanceParameters(){
  _current_index += ComboInterval;
  // check to see if done
  if(_current_index >= _total_models_in_file)
    return true;
  else
    return false;
}

///
/// Returns total number of models in file
/// @return number of models
///
double BioComboGenerator::calc_combinations(){
    return double(_total_models_in_file);
}

///
/// Reads checkpoint information into generator
/// @param instream Input to read from
///
void BioComboGenerator::ReadCheckpoint(ifstream& instream){
  instream >> _current_index;
}
 
///
/// Writes checkpoint information from generator
/// @param outstream Stream for writing
///
void BioComboGenerator::WriteCheckpoint(ofstream& outstream){
  outstream << _current_index << endl;
}

}

