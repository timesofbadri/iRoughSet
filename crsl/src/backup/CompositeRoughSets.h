//============================================================================
// Name        : CompositeRoughSets.h
// Author      : Junbo Zhang
// Version     : 0.6
// Created on  : Aug. 26, 2012 at Home (Atlanta, GA 30338, USA)
// Modified on : Apr. 17, 2013 at SWJTU
// Copyright   : Your copyright notice
// Description : Calculate Approximations based on CRS
//============================================================================

#ifndef COMPOSITEROUGHSETS_H_
#define COMPOSITEROUGHSETS_H_

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
//void CalcMatrix(const vector<CompositeData> data, VVI & RM, VVI & DM, STR2VI_MAP & DMap, VS & deci_val);
void CalcMatrix(const vector<CompositeData> data, VVI & RM, STR2VI_MAP & DMap, VS & deci_val);

/// Calc diagonal matrix Lambda
void CalcLambda(const VVI RM, VI & LM);

/// Calc matrix H
//void CalcH(const VVI RM, const VVI DM, const VI LM, VVD & H);
void CalcH(const VVI RM, STR2VI_MAP DMap, const VI LM, VVD & H);

void OutputAx(const VVD H, const VS deci_val, const double alpha, const double beta, const char * fname);

void OutputMatrix(const VVI matrix);
void OutputMatrix(const VVD matrix);

#endif /* COMPOSITEROUGHSETS_H_ */
