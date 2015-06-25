#ifndef _ordering_io_hpp_included_
#define _ordering_io_hpp_included_

// c++ header files
#include <fstream>
#include <iostream>
#include <string>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <vector>
#include <ctime>
// #include <cpgplot.h>
#include <map>

// my header files
#include "ordering-tools.hpp"

// namespaces
using namespace std;

vector<int> OrderVector(int l1,int l2,int w,int n);

vector<int> GetNodeOrder(const vector<int>& kernelOrder, const vector<int>& klines,
			  			 int netSize, int option);

vector< vector<int> > GetKernels(vector<double>& similarityMatrix,
				  				 vector<int>& originalOrder,
				  				 vector <int>& klines, int& n,
				  				 vector<int>& translationTable );

void ReadMatrix(const int format, const char* fileName,
		 vector<int>& transTable , vector<double>& sim );

void ReadTranslationTableMatrix(const char* fileName, vector<int>& transTable);

void ReadData(const int format,
			  vector<int>& transTable,
			  vector<double>& sim);

void ReadData(double& temperature,double& factor,unsigned int& iter,
	      	  const unsigned int format,
	      	  vector<int>& translationTable,
	      	  vector<double>& similarityMatrix,
	      	  double& minenergy);

void PrintMatrix (const char* fileName,
				  const vector<double>& sim,
				  const vector<int>& kernelOrder,
		  		  const vector<int>& klines,
		  		  int n,
		  		  int mode);

void PrintTranstable(const char* fileName,
					 const vector<int>& translationTable,
					 const vector<int>& order,
					 const vector<vector<int> >& kernels);

void PrintValue(const char* fileName, const double value);

#endif
