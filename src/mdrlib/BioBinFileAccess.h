//BioBinFileAccess.h

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

#ifndef __BioBinFileAccess_H__
#define __BioBinFileAccess_H__

#include "BioFileAccess.h"

namespace mdr{

class BioBinFileAccess: public BioFileAccess{

  public:
  
    /// constructor
    BioBinFileAccess(){}
  
    /// constructor
    BioBinFileAccess(std::string Filename);

  	/// reads the indicated number of models starting at the index given
	  bool GetModels(int StartIndex, std::vector<std::vector<std::string> >& models, int nModels,int idlength=11);
  	 	
  	/// get number of records in file
	  unsigned int GetRecordCount();


};

}

#endif
