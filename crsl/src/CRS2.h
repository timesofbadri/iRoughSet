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
// Name        : CRS2.cpp
// Created on  : Aug. 26, 2012 at Home (Atlanta, GA 30338, USA)
// Modified on : Apr. 17, 2013 at SWJTU
// Description : Calculate Approximations based on CRS
//============================================================================

#ifndef __CRS2_H_
#define __CRS2_H_

//#pragma warning(disable: 4786)  // Disable warning for too long decorated name

/************************Standard C++ header files********************************/
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <deque>
#include <list>
#include <set>
#include <map>
#include <algorithm>
#include <functional>
#include <numeric>
#include <cmath>
#include <memory.h>
#include <time.h>
#include <stdlib.h>
#include <cstring>
#include <ctime>
/************************Standard C++ namespace**********************************/
using namespace std;

namespace CRS2{

/************************User-defined class**************************************/
class SetValuedData
{
public:
	vector<string> val;
	//string raw_val;
	SetValuedData(){};
	SetValuedData(string raw_val)
	{
		val.clear();
		unsigned int i=0, s=0;
		for(i=0;i<raw_val.size();++i)
		{
			if(raw_val[i] == ',')
			{
				val.push_back( raw_val.substr(s, i-s) );
				s = i+1;
			}
		}
		if(s<raw_val.size()) val.push_back(raw_val.substr(s));
	}
	bool operator == (const SetValuedData & C) const;
};

class MissingData
{
public:
	string val;
	bool operator == (const MissingData & C) const;
};

class CompositeData
{
public:
	CompositeData() {};
	vector<string> c; // Categorical data : 1
	vector<MissingData> m; // Missing data : 2
	vector<SetValuedData> sv; // Set-valued data : 3
	vector<double> n; // Numerical data : 4
	string d; // Decision : 0
	bool operator == (const CompositeData &c) const;
	void clear()
	{
		c.clear();
		m.clear();
		sv.clear();
		n.clear();
	}
};

/************************typedef**************************************/
typedef vector<int> VI;
typedef vector<double> VD;
typedef vector<VI > VVI;
typedef vector<string> VS;
typedef map<string, VI > STR2VI_MAP;
typedef vector<VD > VVD;


/************************Main functions**************************************/
void InDataType(const char * fname, vector<int> & data_type);

void InData(const char * fname, const vector<int> data_type, vector<CompositeData> & data);

///Calculate relation matrix and decision matrix
///RM: n * n
///DM: d * n
///deci_val: d * 1
//void CalcMatrix(const vector<CompositeData> data, VVI & RM, VVI & DM, STR2VI_MAP & DMap, VS & deci_val);
void CalcMatrix(const vector<CompositeData> data, VVI & RM, VVI & DM, VS & deci_val);

/// Calc diagonal matrix Lambda
void CalcLambda(const VVI RM, VI & LM);
/// Calc middle matrix Omega
/// Omega: n * d
void CalcOmegaMatrix(const VVI & RM, const VVI & DM, VVI & OmegaM);

/// Calc matrix H
//void CalcH(const VVI RM, const VVI DM, const VI LM, VVD & H);
void CalcH(const VVI OmegaM, const VI LM, VVD & H);

void OutputAx(const VVD H, const VS deci_val, const double alpha, const double beta, const char * fname);

void OutputMatrix(const VVI matrix);
void OutputMatrix(const VVD matrix);
}
#endif // __CRS2_H_ 
