//============================================================================
// Name        : parse_args.h
// Author      : Junbo Zhang
// Version     : 0.6
// Created on  : Aug. 26, 2012 at Home (Atlanta, GA 30338, USA)
// Modified on : Apr. 17, 2013 at SWJTU
// Copyright   : Your copyright notice
// Description : Parse args
//============================================================================

#ifndef PARSE_ARGS_H_
#define PARSE_ARGS_H_

#include <iostream>
#include <stdlib.h>
#include <string.h>
//#include <time.h>
//#include <sys/time.h>

#ifndef MAX_PATH
#define MAX_PATH 1024
#endif //MAX_PATH

typedef unsigned int uint;

struct parameter //parameter config
{
	char   data_type[MAX_PATH]; // the path of input file name about data type
	char   data[MAX_PATH];      // the path of input file name
	char   data_add[MAX_PATH];  // the path of added input file name
	int    del_ratio;           // delete ratio: 0 ~ 100
	char   output[MAX_PATH];
	double alpha, beta, delta;  // x, y, z
	char   model[MAX_PATH];     //
};

void exit_with_help();

/** parse_command_line()
 *  Parse the user arguments
 */
void parse_command_line(const int argc, char **argv, parameter & m_parameter);


#endif /* PARSE_ARGS_H_ */
