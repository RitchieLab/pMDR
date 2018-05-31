//Checkpoint.h

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

#ifndef __CHECKPOINT_H__
#define __CHECKPOINT_H__

#include <fstream>
#include <iostream>
#include "ComboGenerator.h"
#include "TrainList.h"
#include "MDRExcept.h"

namespace mdr{

class Checkpoint{

  public:
    Checkpoint(){}
    Checkpoint(string checkfilename){filename=checkfilename;}

    void write_checkpoint(TrainList& models, ComboGenerator* generator, string filename,
      float interval, double mods_complete);
    
    double read_checkpoint(TrainList& models, ComboGenerator* generator, string filename);

    double read_checkpoint(TrainList& models, ComboGenerator* generator);

    string checkpoint_name(){return filename;}
    
    void checkpoint_name(string file){filename = file;}

  private:
    
    ofstream& open_writing();
    
    ifstream& open_reading();
    
    void close_file();
    
    string filename;

};

}






#endif
