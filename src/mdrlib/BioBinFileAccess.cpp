//BioBinFileAccess.cpp

#include "BioBinFileAccess.h"
#include "Stringmanip.h"

using namespace std;

namespace mdr{

///
/// Alternative constructor
///
BioBinFileAccess::BioBinFileAccess(string Filename):BioFileAccess(Filename){
}


///
/// Get number of records in file from header information in the text file
/// @return number of models in file
///
unsigned int BioBinFileAccess::GetRecordCount(){
  
  if(!biofile.is_open()){
    OpenBinary();
  }
 
  // total number of models is second unsigned int in the file
  unsigned int totalModels;
  // skip first unsigned int
  biofile.seekg(sizeof(unsigned int));
  biofile.read((char*)&totalModels, sizeof(unsigned int));

  return totalModels;
}


///
/// Returns models -- specific number or to end of file
/// @param StartIndex
/// @param models
/// @param nModels Number of models to get
/// @return False when end of file reached
///
bool BioBinFileAccess::GetModels(int StartIndex, vector<vector<string> >& models, 
  int nModels, int idlength){
  
  // size of these are 11 for each model and then the implication index is 7 characters long
  // followed by a return character
  // so each line is 30 bytes long (each model)
  int record_size = 2 * sizeof(unsigned int) + sizeof(float);

  bool end_reached = false;

  // set position to next model (zero based)
  biofile.seekg(record_size*StartIndex + headeroffset);

  unsigned int rsone, rstwo;
  float implication;
  vector<string> model(2,"");
  
  for(int n=0; n<nModels; n++){
    if(biofile.eof()){
      end_reached = true;
      break;
    }
  
    biofile.read( (char *)&rsone, sizeof(unsigned int));
    biofile.read( (char *)&rstwo, sizeof(unsigned int));
    biofile.read( (char *)&implication, sizeof(float));
    model[0] = "rs"+ Stringmanip::itos(rsone);
    model[1] = "rs"+ Stringmanip::itos(rstwo);

    models.push_back(model);
  }
  
  return end_reached;
}

}



