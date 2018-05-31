//pKnuthComboGenerator.h

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

#ifndef __PKNUTHCOMBOGENERATOR_H__
#define __PKNUTHCOMBOGENERATOR_H__

#include "pComboTransfer.h"
#include <KnuthComboGenerator.h>

namespace mdr{

class pKnuthComboGenerator: public pComboTransfer, public KnuthComboGenerator{

  public:
  
    /// For sending combo parameters
    void SendParameters(int node, int srParameters);
    
    /// Send stop
    void SendStop(int node, int srParameters);

    /// Sends signal to slaves to send models for checkpoint
    void SendModelSignal(int node, int srParameters);
    
    /// For receiving combo parameters
    int ReceiveParameters(int srParameters);

};

}

#endif
