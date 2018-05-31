//OddsRatioStat.h

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

#ifndef __ODDSRATIOSTAT_H__
#define __ODDSRATIOSTAT_H__

using namespace std;

namespace mdr{
///
/// Abstract base class that calculates stat to be used in ranking models
///
class OddsRatioStat: public Stat{

  public:
  
    OddsRatioStat(){set_analysis_name("Odds Ratio");}
  
    /// invert result as MDR assumes smaller values are better
    virtual float analysis_value(unsigned int true_positives, unsigned int false_positives, unsigned int true_negatives, 
      unsigned int false_negatives){
        return 1/((float(true_positives)+0.5)*(float(true_negatives)+0.5)/
          ((float(false_positives)+0.5) * (float(false_negatives)+0.5)));
      }
    
    virtual float convert_to_display(float calculated_value){ return 1/calculated_value;}
    
  private:

};



}

#endif
