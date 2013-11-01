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
// Name        : crsl.cpp
// Created on  : Aug. 26, 2012 at Home (Atlanta, GA 30338, USA)
// Modified on : Apr. 17, 2013 at SWJTU
// Description : Calculate Approximations based on IncCRS
//============================================================================

#include "parse_args.h"
#include "CRS2.h"
#include "IncCRS2.h"

#include <time.h>

double time_diff(timespec start, timespec end) {
	timespec temp;
	if ((end.tv_nsec - start.tv_nsec) < 0) {
		temp.tv_sec = end.tv_sec - start.tv_sec - 1;
		temp.tv_nsec = 1000000000 + end.tv_nsec - start.tv_nsec;
	} else {
		temp.tv_sec = end.tv_sec - start.tv_sec;
		temp.tv_nsec = end.tv_nsec - start.tv_nsec;
	}
	return 1.0 * temp.tv_sec + 1.0 * temp.tv_nsec * (1e-9);
}

parameter m_parameter;
double gDELTA = 0.15;
vector<CRS2::CompositeData> gDATA;
CRS2::VI gDATA_TYPE;
CRS2::VVI gRM; // gRelationMatrix;
//STR2VI_MAP gDMap; // gDecisionMap
CRS2::VVI gDM; // gDecisionMatrix
CRS2::VVI gOmegaM;
CRS2::VS gDeci_val; // decision value
CRS2::VI gLM; // diagonal matrix Lambda
CRS2::VVD gH; // result matrix H, which is also accuracy matrix
//VVD gCov;			// coverage matrix, which is used to induce rule

void TestPrint();
//void time_elapsed();
timespec time1, time2;

void static_work();
void add_work();
void del_work();

int main(int argc, char ** argv) {
	parse_command_line(argc, argv, m_parameter);
	switch (m_parameter.model[0]) {
	case 's':
		static_work();
		break;
	case 'a':
		add_work();
		break;
	case 'd':
		del_work();
		break;
	default:
		exit_with_help();
		break;
	}
	return 0;
}

///static, matrix method used
void static_work() {
	//long t1 = clock();
	clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &time1);
	CRS2::InDataType(m_parameter.data_type, gDATA_TYPE);
	CRS2::InData(m_parameter.data, gDATA_TYPE, gDATA);
	CRS2::CalcMatrix(gDATA, gRM, gDM, gDeci_val);
	CRS2::CalcLambda(gRM, gLM);
	CRS2::CalcOmegaMatrix(gRM, gDM, gOmegaM);
	CRS2::CalcH(gOmegaM, gLM, gH);
	CRS2::OutputAx(gH, gDeci_val, m_parameter.alpha, m_parameter.beta,
			m_parameter.output);

	//long t2 = clock();
	//printf("elapsed time:\t%f\tseconds [static]\n",
	//		(t2 - t1) * 1.0 / CLOCKS_PER_SEC);
	clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &time2);
	printf("Elapsed time (in second) [Static]:\t%lf\n",
			time_diff(time1, time2));
}

///dynamic, incremental method used
void add_work() {
	CRS2::InDataType(m_parameter.data_type, gDATA_TYPE);
	CRS2::InData(m_parameter.data, gDATA_TYPE, gDATA);
	CRS2::CalcMatrix(gDATA, gRM, gDM, gDeci_val);
	CRS2::CalcLambda(gRM, gLM);
	CRS2::CalcOmegaMatrix(gRM, gDM, gOmegaM);
	//CalcH(gRM, gDMap, gLM, gH);
	//cout << 1 << endl;
	//time start
	//long t1 = clock();
	clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &time1);
	CRS2::AddInData(m_parameter.data_add, gDATA_TYPE, gDATA);
	CRS2::AddUpdateMatrix(gDATA, gRM, gDM, gDeci_val);
	//cout << gDeci_val.size() << endl;
	gDATA.clear();
	CRS2::AddUpdateLambda(gRM, gLM);
	CRS2::AddUpdateOmega(gRM, gDM, gOmegaM);
	gRM.clear();
	gDM.clear();
	CRS2::CalcH(gOmegaM, gLM, gH);
	CRS2::OutputAx(gH, gDeci_val, m_parameter.alpha, m_parameter.beta,
			m_parameter.output);
	//long t2 = clock();
	//printf("elapsed time:\t%f\tseconds [add]\n",
	//		(t2 - t1) * 1.0 / CLOCKS_PER_SEC);
	clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &time2);
	printf("Elapsed time (in second) [add]:\t%lf\n", time_diff(time1, time2));
}

void del_work() {
	CRS2::InDataType(m_parameter.data_type, gDATA_TYPE);
	CRS2::InData(m_parameter.data, gDATA_TYPE, gDATA);
	CRS2::CalcMatrix(gDATA, gRM, gDM, gDeci_val);
	CRS2::CalcLambda(gRM, gLM);
	CRS2::CalcOmegaMatrix(gRM, gDM, gOmegaM);
	//CalcH(gRM, gDMap, gLM, gH);

	//time start
	clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &time1);
	//long t1 = clock();
	uint delNum = m_parameter.del_ratio * gRM.size() / 100;
	CRS2::DelUpdateOmega(gRM, gDM, gOmegaM, delNum);
	CRS2::DelUpdateLambda(gRM, gLM, delNum);
	CRS2::DelUpdateMatrix(gDATA, gRM, gDM, gDeci_val,delNum);
	CRS2::DelInData(gDATA, delNum);
	CRS2::CalcH(gOmegaM, gLM, gH);
	CRS2::OutputAx(gH, gDeci_val, m_parameter.alpha, m_parameter.beta,
			m_parameter.output);
	//long t2 = clock();
	//printf("elapsed time:\t%f\tseconds [del]\n",
	//		(t2 - t1) * 1.0 / CLOCKS_PER_SEC);
	clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &time2);
	printf("Elapsed time (in second) [del]:\t%lf\n", time_diff(time1, time2));
}

/*
 void time_elapsed() {
 long t1 = clock();
 printf("elapsed time:\t%f seconds\n", t1 * 1.0 / CLOCKS_PER_SEC);
 }*/

void TestPrint() {
	cout << "RM" << endl;
	CRS2::OutputMatrix(gRM);
	cout << "H" << endl;
	CRS2::OutputMatrix(gH);
	cout << "Set-valuedData" << endl;
	for (uint i = 0; i < gDATA.size(); ++i) {
		for (uint j = 0; j < gDATA[i].sv.size(); ++j) {
			for (uint k = 0; k < gDATA[i].sv[j].val.size(); ++k)
				cout << gDATA[i].sv[j].val[k] << ",";
			cout << "\t";
		}
		cout << endl;
	}
}


// vim: ts=4 sw=4 sts=4 smarttab smartindent expandtab


