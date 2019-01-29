//fitness.h

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

#ifndef __FITNESS_H__
#define __FITNESS_H__

namespace mdr{

/// Lists types of fitness available 
enum fitnessTypes{
  BalancedAccuracy,
  Error,
  OddsRatio
};

/// Lists types of permutation methods available
enum permutationTypes{
  TopRank,
  Ranked,
  Combined
};

}

#endif
