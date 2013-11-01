//============================================================================
// Name        : CompositeRoughSets.cpp
// Author      : Junbo Zhang
// Version     : 0.6
// Created on  : Aug. 26, 2012 at Home (Atlanta, GA 30338, USA)
// Modified on : Apr. 17, 2013 at SWJTU
// Copyright   : Your copyright notice
// Description : Calculate Approximations based on CRS
//============================================================================

#include "CompositeRoughSets.h"

bool SetValuedData::operator == (const SetValuedData & C) const
{
	unsigned int i,j;
	for(i=0;i<(*this).val.size();++i)
	{
		for(j=0;j<C.val.size();++j)
		{
			if((*this).val[i] == C.val[j]) return true;
		}
	}
	return false;
}

bool MissingData::operator ==(const MissingData & C) const
{
	if((*this).val == "?") return true;
	if(C.val == "?") return false;
	if((*this).val == "*" || C.val == "*") return true;
	if((*this).val == C.val) return true;
	return false;
}

///Priority level: C > M > S > N (1 > 2 > 3 > 4)
bool CompositeData::operator == (const CompositeData & C) const
{
	unsigned int i;
	/// Categorical data
	for(i=0;i<(*this).c.size();++i)
	{
		if((*this).c[i] != C.c[i]) return false;
	}
	/// Missing data
	for(i=0;i<(*this).m.size();++i)
	{
		if( !((*this).m[i] == C.m[i]) ) return false;
	}
	/// Set-valued data
	for(i=0;i<(*this).sv.size();++i)
	{
		if( !((*this).sv[i] == C.sv[i]) ) return false;
	}
	///Numerical data, 2-norm, delta
	double dis = 0;
	for(i=0;i<(*this).n.size();++i)
	{
		dis += ((*this).n[i] - C.n[i]) * ((*this).n[i] - C.n[i]);
	}
	dis = sqrt(dis);
	//if(dis > gDELTA) return false;
	if(dis > 0.15) return false;
	return true;
}

void InData(const char * fname, const vector<int> data_type, vector<CompositeData> & data)
{
	ifstream fin(fname);
	// file opened?
	if (! fin)
	{
	   // NO, abort program
	   cerr << "can't open input file \"" << fname << "\""
	        << endl;
	   exit(EXIT_FAILURE);
	}

	data.clear();
	unsigned int i;
	string tc, td, t_str;
	MissingData tm;
	//SetValuedData tsv;
	double tn;

	while(!fin.eof())
	{
		CompositeData cTmp;
		cTmp.clear();
		for(i=0;i<data_type.size();++i)
		{
            switch(data_type[i])
			{
			case 1:
				fin >> tc;
				cTmp.c.push_back(tc);
				break;
			case 2:
				fin >> tm.val;
				cTmp.m.push_back(tm);
				break;
			case 3:
				fin >> t_str;
				cTmp.sv.push_back(SetValuedData(t_str));
				break;
			case 4:
				fin >> tn;
				cTmp.n.push_back(tn);
				break;
			case 0:
				fin >> td;
				cTmp.d = td;
				break;
			default:
				break;
			}
		}
        if(fin.eof()) break;
		data.push_back(cTmp);
	}
}

///Read data type
void InDataType(const char * fname, vector<int> & data_type)
{
	ifstream fin(fname);
	data_type.clear();
	string id, type;
	while(fin>>id>>type)
	{
		switch(type.at(0))
		{
			case 'C': data_type.push_back(1); break;
			case 'M': data_type.push_back(2); break;
			case 'S': data_type.push_back(3); break;
			case 'N': data_type.push_back(4); break;
			case 'D': data_type.push_back(0); break;
			default : data_type.push_back(1); break;
		}
	}
	fin.close();
}

///Calculate relation matrix and decision matrix
/*void CalcMatrix(const vector<CompositeData> data, VVI & RM, VVI & DM, STR2VI_MAP & DMap, VS & deci_val)
{
	unsigned int i, j;
	//Calc Relation Matrix
	RM.clear();
	VI vbTmp(data.size(),0);
	RM.insert(RM.begin(),data.size(),vbTmp);
	for(i=0; i<data.size(); ++i)
	{
		RM[i][i] = 1;
		for(j=0; j<data.size(); ++j)
		{
			if(data[i] == data[j]) RM[i][j] = 1;
		}
	}

	//Calc Relation Matrix
	DMap.clear();
	STR2VI_MAP::iterator it;
	for(i=0; i<data.size(); ++i)
	{
		it = DMap.find(data[i].d);
		if(it == DMap.end() )
		{
			VI tmp;
			tmp.push_back(i);
			DMap.insert(STR2VI_MAP::value_type(data[i].d, tmp) );
		}
		else
		{
			((*it).second).push_back(i);
		}
	}

	DM.clear();
	DM.insert(DM.begin(),DMap.size(),vbTmp);
	deci_val.clear();
	for(it = DMap.begin(), i=0; it!=DMap.end(); ++it, ++i)
	{
		deci_val.push_back((*it).first);
		for(j=0;j<(*it).second.size();++j)
		{
			DM[i][(*it).second[j] ] = 1;
		}
	}
}*/

void CalcMatrix(const vector<CompositeData> data, VVI & RM, STR2VI_MAP & DMap, VS & deci_val)
{
	unsigned int i, j;
	//Calc Relation Matrix
	RM.clear();
	VI vbTmp(data.size(),0);
	RM.insert(RM.begin(),data.size(),vbTmp);
	for(i=0; i<data.size(); ++i)
	{
		for(j=0; j<data.size(); ++j)
		{
			if(data[i] == data[j]) RM[i][j] = 1;
		}
		RM[i][i] = 1;
	}

	//Calc Relation Matrix
	DMap.clear();
	STR2VI_MAP::iterator it;
	for(i=0; i<data.size(); ++i)
	{
		it = DMap.find(data[i].d);
		if(it == DMap.end() )
		{
			VI tmp;
			tmp.push_back(i);
			DMap.insert(STR2VI_MAP::value_type(data[i].d, tmp) );
		}
		else
		{
			((*it).second).push_back(i);
		}
	}

	for(it = DMap.begin(); it!=DMap.end(); ++it)
	{
		deci_val.push_back((*it).first);
	}
}


/// Calc diagonal matrix Lambda
void CalcLambda(const VVI RM, VI & LM)
{
	LM.clear();
	unsigned int i,j;
	int sum = 0;
	for(i=0;i<RM.size();++i)
	{
		sum = 0;
		for(j=0;j< RM[i].size();++j)
		{
			sum += RM[i][j];
		}
		LM.push_back(sum);
	}
}

/// Calc matrix H
/*void CalcH(const VVI RM, const VVI DM, const VI LM, VVD & H) // old method
{
	if(RM.size() == 0 || DM.size() == 0 || RM.size() != DM[0].size())
	{
		cerr << "Size of relation matrix and decison matrix is 0 or not equal !" << endl;
		exit(-1);
	}
	H.clear();
	VD vdTmp(DM[0].size(),0);
	H.insert(H.begin(),DM.size(),vdTmp);

	unsigned int i,j,k;
	int sum = 0;
	for(i=0;i<H.size();++i)
	{
		for(j=0;j<RM.size();++j)
		{
			sum = 0;
			for(k=0;k<H[i].size();++k)
			{
				sum += (RM[j][k] & DM[i][k]);
			}
			H[i][j] = 1.0 * sum / LM [j];
		}
	}
}
*/
void CalcH(const VVI RM, STR2VI_MAP DMap, const VI LM, VVD & H)
{
	if(RM.size() == 0 || DMap.size() == 0)
	{
		cerr << "Size of relation matrix and decison matrix is 0 !" << endl;
		exit(-1);
	}
	H.clear();
	VD vdTmp(RM.size(),0);
	H.insert(H.begin(),DMap.size(),vdTmp);

	unsigned int i,j,k;
	STR2VI_MAP::iterator it;
	int sum = 0;

	for(it=DMap.begin(), i=0; it!= DMap.end(); ++it, ++i)
	{
		for(j=0;j<RM.size();++j)
		{
			sum = 0;
			for(k=0;k<(*it).second.size();++k)
			{
				sum += RM[j][(*it).second[k]];
			}
			H[i][j] = 1.0 * sum / LM [j];
		}
	}
}

///Ouput Approximations
void OutputAx(const VVD H, const VS deci_val, const double alpha, const double beta, const char * fname)
{
	if(H.size() == 0)
	{
		cerr << "Size of H is 0!" << endl;
		exit(-1);
	}
	ofstream fout(fname);
	unsigned int i,j;
	int sum = 0;
	for(i=0;i<H.size();++i)
	{
		fout << "lower approximation of the decision\t" << deci_val[i] << endl;
		for(j=0;j<H[i].size();++j)
		{
			if(H[i][j] >= alpha)
			{
				fout << j << " ";
				++ sum;
			}
		}
		fout << endl;
	}

	for(i=0;i<H.size();++i)
	{
		fout << "upper approximation of the decision\t" << deci_val[i] << endl;
		for(j=0;j<H[i].size();++j)
		{
			if(H[i][j] > beta)
			{
				fout << j << " ";
			}
		}
		fout << endl;
	}

	fout << endl << "number of positive region:\t" << sum << endl;
}


void OutputMatrix(const VVI matrix)
{
	unsigned int i,j;
	for(i=0;i<matrix.size();++i)
	{
		for(j=0;j<matrix[i].size();++j)
		{
			cout << matrix[i][j] <<" ";
		}
		cout << endl;
	}
}

void OutputMatrix(const VVD matrix)
{
	unsigned int i,j;
	for(i=0;i<matrix.size();++i)
	{
		for(j=0;j<matrix[i].size();++j)
		{
			cout << matrix[i][j] <<" ";
		}
		cout << endl;
	}
}
