//BioFileAccess.h

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

#ifndef __BIOFILEACCESS_H__
#define __BIOFILEACCESS_H__

#include <fstream>
#include <vector>
#include <iostream>

namespace mdr{

class BioFileAccess{

  public:
  
    /// constructor
    BioFileAccess();
  
    /// constructor
    BioFileAccess(std::string Filename);
    
    /// destructor
    virtual ~BioFileAccess();
    
    /// sets file name
    void SetFileName(std::string fname);
   
  	/// opens file
  	virtual bool Open();
  	
  	/// reads an index
	  virtual bool GetModels(int StartIndex, std::vector<std::vector<std::string> >& models, int nModels,int idlength=11)=0;
  	
  	/// Gets header length for use as an offset
  	int GetHeaderOffset();
  	
  	/// close file
  	void Close();
  	
  	/// get number of records in file
	  virtual unsigned int GetRecordCount()=0;

    /// check whether file is binary
    static bool FileBinary(std::string Filename);

  protected:
  
    /// opens file in binary mode
    void OpenBinary();
  
    std::ifstream biofile;
  	bool isopen;
	  std::string filename;
	  int headeroffset;

};

}

#endif
