//BioTextFileAccess.cpp

#include "BioTextFileAccess.h"
#include "Stringmanip.h"

using namespace std;

namespace mdr{

BioTextFileAccess::BioTextFileAccess(){
}

/// constructor
BioTextFileAccess::BioTextFileAccess(std::string Filename):BioFileAccess(Filename){
}


///
/// Get number of records in file from header information in the text file
/// @return number of models in file
///
unsigned int BioTextFileAccess::GetRecordCount(){
  
  if(biofile.is_open()){
    biofile.close();
  }

  biofile.open(filename.c_str());
  unsigned int totalModels;
  if(biofile.is_open()){
    biofile >> totalModels;
    biofile.close();
  }
 
  OpenBinary();

  return totalModels;
}


///
/// Opens file
/// @return true when opened successfully
///
bool BioTextFileAccess::Open() {
	ios_base::openmode OpenFlags = ios_base::in | ios_base::binary;
	   
	biofile.open(filename.c_str(), OpenFlags);
	
	OpenBinary();
	isopen = biofile.is_open();
	
 	if(isopen)
   	headeroffset = GetHeaderOffset();
	
	return ( isopen );
}

///
/// Returns header size for text file.  Just has to count the length
/// of the first line in the 
///
int BioTextFileAccess::GetHeaderOffset(){

  biofile.close();
  
  biofile.open(filename.c_str());
  // move to second unsigned int
  char headerline[256];
  biofile.getline(headerline,256);
  biofile.close();
  string hline(headerline);
  
  OpenBinary();
  
  return int(hline.length())+1;
}


///
/// Returns models -- specific number or to end of file
/// @param StartIndex
/// @param models
/// @param nModels Number of models to get
/// @return False when end of file reached
///
bool BioTextFileAccess::GetModels(int StartIndex, vector<vector<string> >& models, 
  int nModels, int idLength){
  
//  int rs_size=11, to_end_of_line = 8; // end of line is 7 for the implication and 1 for newline
  int to_end_of_line=1;
  
  // size of these are 11 for each model and then the implication index is 7 characters long
  // followed by a return character
  // so each line is 30 bytes long (each model)
  bool end_reached = false;
  // set position to next model (zero based)
  biofile.seekg((idLength*2+to_end_of_line) *StartIndex + headeroffset);
  int string_len=idLength+1;
  
  char firstRS[string_len], secondRS[string_len], endofline[to_end_of_line];
  vector<string> model(2,"");
  string r1, r2;
  
  for(int n=0; n<nModels; n++){
    biofile.read( (char *)&firstRS, idLength);
    firstRS[idLength]='\0';
    biofile.read( (char *)&secondRS, idLength);
    secondRS[idLength]='\0';

    biofile.read( (char *)&endofline, to_end_of_line);
    if(biofile.eof()){
      end_reached = true;
      break;
    }
    r1 = firstRS;
    r2 = secondRS;
    Stringmanip::trim(r1);
    model[0] = r1;

    Stringmanip::trim(r2);
    model[1] = r2;
//cout << "r1=" << r1 << " r2=" << r2 << endl;
    models.push_back(model);
  } 

  return end_reached;
}

}

