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

#include "IncCRS2.h"

namespace CRS2 {
uint _n, _d;
uint _n_add, _d_add;
uint _n_del, _d_del;
/************************add objects********************************/
void AddInData(const char * fname, const vector<int> data_type,
		vector<CompositeData> & data) {
	ifstream fin(fname);
	if (!fin) {
		// NO, abort program
		cerr << "can't open input file \"" << fname << "\"" << endl;
		exit(EXIT_FAILURE);
	}

	unsigned int i;
	string tc, td, t_str;
	MissingData tm;
	double tn;

	while (!fin.eof()) {
		CompositeData cTmp;
		cTmp.clear();
		for (i = 0; i < data_type.size(); ++i) {
			switch (data_type[i]) {
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
		if (fin.eof())
			break;
		data.push_back(cTmp);
	}
}

void AddUpdateMatrix(const vector<CompositeData> data, VVI & RM, VVI & DM,
		VS & deci_val) {
	unsigned int i, j, n, k;
	n = RM.size();
	_n = RM.size();
	for (i = 0; i < n; ++i) {
		for (j = n; j < data.size(); ++j) {
			if (data[i] == data[j])
				RM[i].push_back(1);
			else
				RM[i].push_back(0);
		}
	}

	for (j = n; j < data.size(); ++j) {
		VI viTmp;
		viTmp.clear();
		for (i = 0; i < data.size(); ++i) {
			if (data[j] == data[i])
				viTmp.push_back(1);
			else
				viTmp.push_back(0);
		}
		//viTmp[j][j] = 1;
		RM.push_back(viTmp);
	}
	_n_add = RM.size();
	///update decision map when adding object

	_d = deci_val.size();
	//uint kk = 0;
	//cout << "DD:\t" << _d <<"\t" << data.size() - n << endl;
	for (i = n; i < data.size(); ++i) {
		for (j = 0; j < deci_val.size(); ++j) {
			if (data[i].d == deci_val[j]) {
				//cout << "id:\t" << data[i].d << endl;
				break;
			}
		}

		if (j != deci_val.size()) /// have not formed new decision
				{
			VI viTmp(DM[0].size(), 0);
			viTmp[j] = 1;
			DM.push_back(viTmp);
			viTmp.clear();
			//cout << kk++ <<"\t"<<DM[0].size() << "\t" <<j << endl;
		} else {
			for (k = 0; k < DM.size(); ++k) {
				DM[k].push_back(0);
			}
			VI viTmp(DM[0].size(), 0);
			viTmp[j] = 1;
			DM.push_back(viTmp);
			deci_val.push_back(data[i].d);
			viTmp.clear();
			//cout << kk++ <<"\t --- "<<DM[0].size() << endl;
		}

	}

	_d_add = deci_val.size();

	/*STR2VI_MAP::iterator it;
	 for(i=n; i<data.size(); ++i)
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
	 deci_val.clear();
	 for(it = DMap.begin(); it!=DMap.end(); ++it)
	 {
	 deci_val.push_back((*it).first);
	 }*/
}

/// Update diagonal matrix Lambda when adding objects
void AddUpdateLambda(const VVI RM, VI & LM) {
	unsigned int i, j, n;
	n = LM.size();

	for (i = 0; i < n; ++i) {
		for (j = n; j < RM.size(); ++j) {
			LM[i] += RM[i][j];
		}
	}

	int sum = 0;
	for (i = n; i < RM.size(); ++i) {
		sum = 0;
		for (j = 0; j < RM[i].size(); ++j) {
			sum += RM[i][j];
		}
		LM.push_back(sum);
	}
}

void AddUpdateOmega(const VVI & RM, const VVI & DM, VVI & OmegaM) {
	//cout << OmegaM[i.size() << "\t" << _d_add << endl;
	if (RM.size() != _n_add || DM.size() != _n_add || OmegaM.size() != _n) {
		cerr << "Size is error in AddUpdateOmega!" << endl;
		exit(-1);
	}
	uint i, j, k;
	if (_d_add > _d) {
		for (i = 0; i < OmegaM.size(); ++i) {
			OmegaM[i].insert(OmegaM[i].end(), _d_add - _d);
		}
	}
	VI viTmp(_d_add, 0);
	for (i = _n; i < _n_add; ++i) {
		OmegaM.push_back(viTmp);
	}

	///left-top: (0,0) ~ (_n - 1, _d -1)
	/// Omega + P * Q'
	for (i = 0; i < _n; ++i)
		for (j = 0; j < _d; ++j) {
			for (k = _n; k < _n_add; ++k) {
				OmegaM[i][j] += RM[i][k] * DM[k][j];
			}
		}
	///right-top: (0,_d) ~ (_n - 1, _d - 1)
	///RM*P' + P*R'
	for (i = 0; i < _n; ++i) {
		for (j = _d; j < _d_add; ++j) {
			for (k = 0; k < _n_add; ++k) {
				OmegaM[i][j] += RM[i][k] * DM[k][j];
			}
		}
	}
	///left-bottom: (_n, 0) ~ (_n_add - 1, _d - 1)
	///Q*DM + R* Q'
	for (i = _n; i < _n_add; ++i)
		for (j = 0; j < _d; ++j) {
			for (k = 0; k < _n_add; ++k) {
				OmegaM[i][j] += RM[i][k] * DM[k][j];
			}
		}
	///right-bottom: (_n,_d) ~ (_n_add - 1, _d_add - 1)
	for (i = _n; i < _n_add; ++i)
		for (j = _d; j < _d_add; ++j) {
			for (k = 0; k < _n_add; ++k) {
				OmegaM[i][j] += RM[i][k] * DM[k][j];
			}
		}
}

/************************del objects********************************/

void DelUpdateOmega(const VVI & RM, const VVI & DM, VVI & OmegaM,
		const uint del_Num) {
	if (del_Num > OmegaM.size() || del_Num < 0) {
		// NO, abort program
		cerr << "del_Num is out of size in [DelUpdateOmega]! " << endl;
		exit(EXIT_FAILURE);
	}

	uint n = del_Num;
	uint i, j;
	OmegaM.erase(OmegaM.end() - n, OmegaM.end());
	uint DelNum = 0;
	for (i = 0; i < OmegaM.size(); ++i) {
		for (j = 0; j < OmegaM[i].size(); ++j) {
			if (OmegaM[i][j] != 0)
				break;
		}
		if (j == OmegaM[i].size()) {
			++DelNum;
		}
	}
}

void DelUpdateLambda(const VVI RM, VI & LM, const  uint del_Num) {
	if (del_Num > LM.size() || del_Num < 0) {
		// NO, abort program
		cerr << "del_Num is out of size in [DelUpdateLambda]! " << endl;
		exit(EXIT_FAILURE);
	}
	uint n = del_Num;
	uint i, j;
	LM.erase(LM.end() - n, LM.end());
	for (i = 0; i < LM.size(); ++i) {
		for (j = n; j < RM.size(); ++j) {
			LM[i] -= RM[i][j];
		}
	}
}

void DelUpdateMatrix(vector<CompositeData> & data, VVI & RM, VVI & DM,
		VS & deci_val, const uint del_Num) {
	if (del_Num > RM.size() || del_Num < 0) {
		// NO, abort program
		cerr << "del_Num is out of size in [DelUpdateMatrix]! " << endl;
		exit(EXIT_FAILURE);
	}
	unsigned int n = del_Num;
	//unsigned int del_num = data.size() - n;
	unsigned int i, j;

	/// delete from Relation Matrix & Decision Matrix
	RM.erase(RM.end() - n, RM.end());
	DM.erase(DM.end() - n, DM.end());
	for (i = 0; i < RM.size(); ++i) {
		RM[i].erase(RM[i].end() - n, RM[i].end());
		//for (j = 0; j < n; ++j)
		//	RM[i].pop_back();
	}
	/// delete from Decision Matrix
	uint DelNum = 0;
	for (i = 0; i < DM.size(); ++i) {
		for (j = 0; j < DM[i].size(); ++j) {
			if (DM[i][j] != 0)
				break;
		}
		if (j == DM[i].size()) {
			DM.erase(DM.begin() + i - DelNum);
			deci_val.erase(deci_val.begin() + i - DelNum);
			++DelNum;
		}
	}
}

void DelInData(vector<CompositeData> & data, uint del_Num) ///del_Num: 0 ~ data.size()
		{
	if (del_Num > data.size() || del_Num < 0) {
		// NO, abort program
		cerr << "del_Num is out of size in [DelInData]! " << endl;
		exit(EXIT_FAILURE);
	}
	/// delete from data
	data.erase(data.end() - del_Num, data.end());
}
}

