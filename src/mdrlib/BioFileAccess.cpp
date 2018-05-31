//BioFileAccess.cpp

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

#include "BioFileAccess.h"
#include "MDRExcept.h"
using namespace std;

namespace mdr{

BioFileAccess::BioFileAccess(){
}

BioFileAccess::BioFileAccess(string Filename ){
  filename = Filename;
}

void BioFileAccess::SetFileName(string fname){
  filename = fname;
}

///
/// check whether file is binary and sets offset for determining position
/// @param Filename
/// @return
///
bool BioFileAccess::FileBinary(string Filename){

  // if it is binary the first unsigned int will equal 0
  // otherwise it will be some other number that forms part of
  // the count in the text file
  ifstream reader;
  bool isbinary=false;
  unsigned int binaryCheck;
  reader.open(Filename.c_str(), ios::binary);
  if(!reader.is_open()){
    throw MDRExcept("Unable to open bio filter file " + Filename);
  }
  reader.read((char*)&binaryCheck, sizeof(unsigned int));
  
  if(binaryCheck==0){
    isbinary = true;
  }
  
  reader.close();
  return isbinary;
}


///
/// Returns header size for binary file.  It has a single 
/// unsigned int that indicates it is binary and then the second 
/// unsigned int that contains the number of models.
///
int BioFileAccess::GetHeaderOffset(){
  // move to second unsigned int
  return 2 * sizeof(unsigned int);
}


///
/// Opens file
/// @return true when opened successfully
///
bool BioFileAccess::Open() {
        if(biofile.is_open()){
          biofile.close();
        }
	ios_base::openmode OpenFlags =ios_base::in | ios_base::binary;
	
	OpenBinary();
	isopen = biofile.is_open();
	
 	if(isopen)
   	headeroffset = GetHeaderOffset();
	
	return ( isopen );
}


///
/// opens file in binary mode
///
void BioFileAccess::OpenBinary(){
	ios_base::openmode OpenFlags =ios_base::in | ios_base::binary;
	   
	biofile.open( filename.c_str(), OpenFlags);
}


///
/// Destructor
///
BioFileAccess::~BioFileAccess() {
  if (isopen)
	  Close();
}

///
/// Closes file
///
void BioFileAccess::Close() {
	biofile.close();
}

}


