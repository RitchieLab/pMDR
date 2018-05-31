// LogOutput.cpp


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

#include "LogOutput.h"
#include <time.h>
#include <sstream>
#include <iomanip>


namespace mdr{

LogOutput::LogOutput(){
}


LogOutput::LogOutput(string filename){
  open_log(filename);
}

LogOutput::~LogOutput(){
  close_log();
}


void LogOutput::open_log(string filename){
  logfilename = filename;
  logstream.open(logfilename.c_str(), ofstream::out);
}


void LogOutput::append_log(){
  if(logstream.is_open()){
    logstream.close();
  }
  logstream.open(logfilename.c_str(), ofstream::app);
}


void LogOutput::close_log(){
  logstream.close();
}

void LogOutput::add_time(){ 
  time_t rawtime;
  struct tm * timeinfo;
  time ( &rawtime );
  timeinfo = localtime ( &rawtime );
  logstream << asctime (timeinfo);
}

void LogOutput::add_line(string logline){
  logstream << logline << endl;
}


void LogOutput::report_completed(double completed_mods, double total_mods){
  double comp;

  append_log();

  if(completed_mods < total_mods){
    comp = completed_mods / total_mods; 
  }
  else{
    comp = 1.0;
  }
  logstream.precision(2);
  logstream << "Completed " << fixed << comp * 100 << "% at ";
  add_time();
  logstream << endl;
  
  close_log();
  
}


void LogOutput::add_perm_model(Model& model, int permNum, Dataset& set,
  Stat* calculator, bool include_lr){
  
  stringstream comb;

  append_log();

  logstream << left << setw(5)<< permNum;
  
  comb << "[";
  unsigned int loc;
  for(loc=0; loc < model.combination.size()-1; loc++)
    comb << set.snp_names[model.combination[loc]] << "_";
  comb << set.snp_names[model.combination[loc]]  << "]";
  
  logstream << left << setw(14) << comb.str();
  
  comb.str("");
  comb.flags(ios::showpoint);
  // output training err
  comb << "TrainErr: " << setprecision(4) << model.training.error;
  
  logstream << setw(17) << comb.str();
  
  // output avg bal acc
  comb.str("");
  comb << "Avg " << calculator->get_analysis_name() << ": " << calculator->convert_to_display(model.get_balpredavg())
    << " ";
  logstream << setw(25) << comb.str();
  
  // output avg prediction err
  comb.str("");
  comb << "Avg Pred Err: " << model.get_predictavg();
  logstream << setw(20) << comb.str();
  
  // output lr reduced model result if needed
  if(include_lr){
    comb.str("");
    comb << "LR stat: " << model.get_interact_llr();
    logstream << setw(20) << comb.str();
  }
  
  logstream << endl;
  
  close_log();
  
}

void LogOutput::total_models(double total_mods){
  logstream << "Number of models: "<< total_mods << "\n";
}


}

