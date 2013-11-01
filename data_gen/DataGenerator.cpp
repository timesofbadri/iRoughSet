/** \file DataGenerator.cpp
    \brief  CRS
* @author 张钧波
* @version 1.0
* @date 2012.08.30
* @language English & Chinese
*/

#pragma once

/************************标准C++必要的头文件************************************/
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
//#include <windows.h>
/************************使用标准C++命名空间**********************************/
using namespace std;

#define M_MAX 2147483647

struct m_parameter //parameter config IS=(U, A, V, f)
{
	 unsigned int    U;		  //对象个数
	// unsigned int    A;       //属性个数
	 unsigned int    V;		  //最大值
	// double	Ratio;            //Ratio 集值比例
	 unsigned int	 A_SV;    //集值属性
	 unsigned int	 A_C;     //符号型属性
	 unsigned int	 A_N;     //数值属性
	 unsigned int	 A_M;     //missing属性
     unsigned int    D;       //决策个数，0表示不需要决策
     char Name[81];
};

void exit_with_help()
{
    printf(
		"Usage: \n"
		"-u Size of universe >0 \n"
    	"-c Number of categorical attribute >0 \n"
    		"-n Number of numerical attribute >0 \n"
    		"-m Number of missing attribute >0 \n"
    		"-s Number of set-valued attribute >0 \n"
		"-v Value region >0\n"
		"-d Size of decision  >0 \n"
                "-o the output data file [default: test1]\n"
		);
	exit(1);
}

void parse_command_line(int argc, char **argv, m_parameter & mp)
{
	///default
	mp.U = 24;
	mp.V = 3;
	mp.A_C = 0;
	mp.A_N = 0;
	mp.A_SV = 0;
	mp.A_M = 0;
	mp.D = 0;
        strcpy(mp.Name, "test1");
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
		case 'u':
			mp.U = atoi(argv[i]);
			if(mp.U <= 0) {
				printf("Parameter error: -%c\n", argv[i-1][1]);
				exit_with_help();
			}
			break;
		case 'v':
			mp.V = atoi(argv[i]);
			if(mp.V <= 0) {
				printf("Parameter error: -%c\n", argv[i-1][1]);
				exit_with_help();
			}
			break;
		case 'c':
			mp.A_C = atoi(argv[i]);
			break;
		case 's':
			mp.A_SV = atoi(argv[i]);
			break;
		case 'm':
			mp.A_M = atoi(argv[i]);
			break;
		case 'n':
			mp.A_N = atoi(argv[i]);
			break;
		case 'd':
			mp.D = atoi(argv[i]);
			if(mp.D < 0) {
				printf("Parameter error: -%c\n", argv[i-1][1]);
				exit_with_help();
			}
			break;
                case 'o':
                        strcpy(mp.Name, argv[i]);
                        break;
		default:
			printf("Unknown option: -%c\n", argv[i-1][1]);
			exit_with_help();
		}
	}
}

void GenerateCategorical(unsigned int V)
{
	printf("%d", rand()%V);
}

void GenerateNumerical(unsigned int V_max)
{
	printf("%lf", (double)(rand()%V_max) / V_max);
}

void GenerateMissing(unsigned int V)
{
	int flag = rand()%V;
	if(flag%3 == 0)
		printf("*");
	else if(flag%3 ==1)
		printf("?");
	else printf("%d", rand()%V);
}

void GenerateSet_valued(unsigned int V)
{
	unsigned int i;
	int v=rand()%V+1;
	set<unsigned int> G;
	for (i=0;i<v;++i)
	{
		G.insert(rand()%V);
	}
	set<unsigned int>::iterator sit=G.begin();
	printf("%d",*sit);
	for (++sit;sit!=G.end();++sit)
	{
		printf(",%d",*sit);
	}
}

int main(int argc, char **argv)
{
	m_parameter mp;
	parse_command_line(argc,argv,mp);
	srand((unsigned)time(NULL));
        
        FILE *fData=0;
        char fD[81];
        strcpy(fD, mp.Name);
        strcat(fD, ".data");
        if(0 == ( fData = freopen(fD, "w" ,stdout)))
        {
            printf("Cannot open file.\n");
            return 0;
        }
	unsigned int i,j;
	for (i=0;i<mp.U;++i)
	{
		if (mp.D != 0)
		{
			GenerateCategorical(mp.D);
		}		
		for(j=0;j<mp.A_C;++j)
		{
			printf("\t");GenerateCategorical(mp.V);
		}
		for(j=0;j<mp.A_M;++j)
		{
			printf("\t");GenerateMissing(mp.V);
		}
		for(j=0;j<mp.A_SV;++j)
		{
			printf("\t");GenerateSet_valued(mp.V);
		}
		for(j=0;j<mp.A_N;++j)
		{
			printf("\t");GenerateNumerical(1000);
		}

		printf("\n");
	}
        fclose(fData);

        FILE *fName=0;
        char fN[81];
        strcpy(fN, mp.Name);
        strcat(fN, ".names");
        if(0 == ( fName = freopen(fN, "w" ,stdout)))
        {
            printf("Cannot open file.\n");
            return 0;
        }
	if (mp.D != 0)
            printf("class: D\n");
	for(j=0;j<mp.A_C;++j)
		printf("a: C\n");
	for(j=0;j<mp.A_M;++j)
		printf("a: M\n");
	for(j=0;j<mp.A_SV;++j)
		printf("a: S\n");
	for(j=0;j<mp.A_N;++j)
		printf("a: N\n");
        fclose(fName);

	return 0;
}
