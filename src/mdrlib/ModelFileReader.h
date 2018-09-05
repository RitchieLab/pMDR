//ModelFileReader.h

//     Copyright 2018 Marylyn Ritchie

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

#ifndef __MODELFILEREADER_H__
#define __MODELFILEREADER_H__

#include "MDRExcept.h"
#include <iostream>
#include <sstream>
#include <fstream>


namespace mdr{

class ModelFileReader{

  public:

    /// returns size of strings used in file
    static int convert_model_file(std::string modelfile, std::string newbiofilename);
  
  private:
  
};


}

#endif

