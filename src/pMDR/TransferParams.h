//TransferParams.h

//     Copyright 2008 Marylyn Ritchie

//     This file is part of pMDR.
// 
//     pMDR is free software: you can redistribute it and/or modify
//     it under the terms of the GNU General Public License as published by
//     the Free Software Foundation, either version 3 of the License, or
//     (at your option) any later version.
// 
//     pMDR is distributed in the hope that it will be useful,
//     but WITHOUT ANY WARRANTY; without even the implied warranty of
//     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//     GNU General Public License for more details.
// 
//     You should have received a copy of the GNU General Public License
//     along with pMDR.  If not, see <http://www.gnu.org/licenses/>.


#ifndef __TRANSFERPARAMS_H__
#define __TRANSFERPARAMS_H__

#include "mpi.h"
#include <Dataset.h>
#include <Config.h>

using namespace std;
using namespace mdr;

/// transfers params and data between master and slaves
class TransferParams{

  public:
  
    TransferParams();
  
    /// send parameters to slaves
    void send_params(Config& config, Dataset& set);
    
    /// receives parameters from master
    vector<unsigned int> receive_params(Config& config, Dataset& set);

    /// send data
    void send_data(Dataset& set);
    
    /// receive data
    void receive_data(Dataset& set);
  
  private:
    int colsplit;

};

#endif

