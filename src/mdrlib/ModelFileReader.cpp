//ModelFileReader.cpp

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

#include "ModelFileReader.h"
#include <iomanip>

using namespace std;

namespace mdr{

int ModelFileReader::convert_model_file(std::string modelfile, std::string newbiofilename){

  std::ifstream modStream(modelfile.c_str(), ios::in);
  if(!modStream.is_open()){
    throw MDRExcept("ERROR: Unable to open " + modelfile + "\n");
  }

  size_t maxSize=0;
  int nModels=0;
  string line, id1, id2;

  // determine size needed for ID length in new file
  do{
    getline(modStream, line);
    if(line.find_first_of("0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz") != string::npos){
     nModels++;
     stringstream ss(line); 
     ss >> id1 >> id2;

     if(id1.length() > maxSize)
       maxSize=id1.length();
     if(id2.length() > maxSize)
       maxSize=id2.length();
    }
  }while(!modStream.eof());
  
  modStream.close();

  // write new biofilter formatted list
  modStream.open(modelfile.c_str(), ios::in);
  ofstream outFile;
  outFile.open(newbiofilename.c_str());
  outFile << nModels << endl;
 
  do{
    getline(modStream, line);
    if(line.find_first_of("0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz") != string::npos){
     stringstream ss(line);
     ss >> id1 >> id2;
     outFile << setw(maxSize) << id1 << setw(maxSize) << id2 << endl;
    }
  }while(!modStream.eof());

  modStream.close();
  outFile.close();

  return maxSize;
}


}

