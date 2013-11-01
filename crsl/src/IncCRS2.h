//============================================================================
//   Copyright (c) 2013 Southwest Jiaotong University.
//        All rights reserved.
//
//   Licensed under the Apache License, Version 2.0 (the "License");
//   you may not use this file except in compliance with the License.
//   You may obtain a copy of the License at
//
//     http://www.apache.org/licenses/LICENSE-2.0
//
//   Unless required by applicable law or agreed to in writing, software
//   distributed under the License is distributed on an "AS IS" BASIS,
//   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
//   See the License for the specific language governing permissions and
//   limitations under the License.
//
//   For more about this software visit:
//
//     http://cs.gsu.edu/~jbzhang
//============================================================================
// Name        : IncCRS2.cpp
// Created on  : Aug. 26, 2012 at Home (Atlanta, GA 30338, USA)
// Modified on : Apr. 17, 2013 at SWJTU
// Description : Calculate Approximations based on IncCRS
//============================================================================

#ifndef __INCCRS2_H_
#define __INCCRS2_H_

#include "CRS2.h"

namespace CRS2{

extern uint _n, _d;
extern uint _n_add, _d_add;
extern uint _n_del, _d_del;

/************************add objects********************************/
void AddInData(const char * fname, const vector<int> data_type, vector<CompositeData> & data);
void AddUpdateMatrix(const vector<CompositeData> data, VVI & RM, VVI & DM, VS & deci_val);
/// Update diagonal matrix Lambda when adding objects
void AddUpdateLambda(const VVI RM, VI & LM);

void AddUpdateOmega(const VVI & RM, const VVI & DM, VVI & OmegaM);

/************************del objects********************************/
void DelUpdateOmega(const VVI & RM, const VVI & DM, VVI & OmegaM, const uint del_Num);

void DelUpdateLambda(const VVI RM, VI & LM, const uint del_Num);

void DelUpdateMatrix(vector<CompositeData> & data, VVI & RM, VVI & DM, VS & deci_val, const uint del_Num);

void DelInData(vector<CompositeData> & data, const uint del_Num); ///del_Num: 0 ~ data.size()
}
#endif // __INCCRS2_H_ 
