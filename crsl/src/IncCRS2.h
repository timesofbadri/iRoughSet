//============================================================================
// Name        : IncCompositeRoughSets.h
// Author      : Junbo Zhang
// Version     : 0.6
// Created on  : Aug. 26, 2012 at Home (Atlanta, GA 30338, USA)
// Modified on : Apr. 17, 2013 at SWJTU
// Copyright   : Your copyright notice
// Description : Calculate Approximations based on IncCRS
//============================================================================

#ifndef INCCRS2_H_
#define INCCRS2_H_

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
#endif /* INCCRS2_H_ */
