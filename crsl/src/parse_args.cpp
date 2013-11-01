//============================================================================
// Name        : parse_args.cpp
// Author      : Junbo Zhang
// Version     : 0.6
// Created on  : Aug. 26, 2012 at Home (Atlanta, GA 30338, USA)
// Modified on : Apr. 17, 2013 at SWJTU
// Copyright   : Your copyright notice
// Description : Parse args
//============================================================================

#include "parse_args.h"

void exit_with_help()
{
	std::cerr << "Usage:\n"
			<< "-m model {static, add, del} \n"
			<< "-t data type file\n"
			<< "-i input data file\n"
			<< "-a added data file\n"
			<< "-d delete radio [0,100]\n"
			<< "-o output file\n"
			<< "-x alpha (beta, 1]\n"
			<< "-y beta [0,alpha) \n"
			<< "-z delta (0,1] \n"
			<< std::endl;
    exit(-1);
}

/** parse_command_line()
 *  Parse the user arguments
 */
void parse_command_line(const int argc, char **argv, parameter & m_parameter)
{
	m_parameter.alpha = 1;
	m_parameter.beta = 0;
	m_parameter.del_ratio = 10;
	m_parameter.delta = 0.15;
	int i=1;
	if (i>=argc)
		exit_with_help();
	for (i=1;i<argc;++i)
	{
		if(argv[i][0] != '-') break;
		if(++i>=argc)
			exit_with_help();
		switch(argv[i-1][1])
		{
		case 'm':
			strcpy(m_parameter.model, argv[i]);
		case 't':
			strcpy(m_parameter.data_type, argv[i]);
			break;
		case 'i':
			strcpy(m_parameter.data, argv[i]);
			break;
		case 'a':
			strcpy(m_parameter.data_add, argv[i]);
			break;
		case 'd':
			m_parameter.del_ratio = atoi(argv[i]);
			break;
		case 'o':
			strcpy(m_parameter.output, argv[i]);
			break;
		case 'x':
			m_parameter.alpha = atof(argv[i]);
			break;
		case 'y':
			m_parameter.beta = atof(argv[i]);
			break;
		case 'z':
			m_parameter.delta = atof(argv[i]);
			break;
		default:
			std::cerr << "Unknown option: -"<<argv[i-1][1]<< std::endl;
			//fprintf(stderr,"Unknown option: -%c\n", argv[i-1][1]);
			exit_with_help();
			break;
		}
	}
}
