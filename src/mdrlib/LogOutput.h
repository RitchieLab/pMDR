//LogOutput.h

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

#ifndef __LOG_H__
#define __LOG_H__

#include <fstream>
#include <iostream>
#include "Model.h"
#include "Dataset.h"
#include "Stat.h"

namespace mdr{

class LogOutput{

  public:
    LogOutput();
    LogOutput(string filename);
    ~LogOutput();

    void open_log(string filename);
    void close_log();

    void add_line(string logline);
    void add_time();
    void report_completed(double completed_mods, double total_mods);
    void add_perm_model(Model& model, int numPerm, Dataset& set, Stat* calculator, 
      bool include_lr=false);
    
    void total_models(double total_mods);
    void add_message(string message);

  private:
    
    void append_log();
  
    ofstream logstream;
    string logfilename;

};

}
#endif

