#ifndef _ordering_annealing_hpp_included_
#define _ordering_annealing_hpp_included_

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
#include "ordering-io.hpp"

// namespaces
using namespace std;


//********************************************///

int ChangeSegmentOrder(vector<int>& order, int line1, int line2, int width);

int MakeChangeVectors(const vector<int>& order, vector<int>& change, 
		      vector<int>& nochange, int line1, int line2, 
		      int width, int n, const vector<int>& klines);

double ComputeEnergy(const vector<double>& matrix, int netSize, const vector<int>& kernelOrder,
		         	 const vector<int>& klines);

double EnergyChangeReOrd(const vector<double>& matrix,int l1,int l2,int w,int n, 
			 const vector<int>& kernelOrder,const vector<int>& klines);

double AnnealStep(double temperature, const vector<double>& matrix, int netSize,
		  		  vector <int>& order,int nblocks, const vector<int>& klines);

double AnnealStepFree(double temperature, const vector<double>& matrix, int netSize,
		      	      vector <int>& order,int nblocks, const vector<int>& klines);

void AnnealIter(double temperature,
				const vector<double>& matrix,
                const vector<int>& translationTable,
                const vector<vector<int> >& kernels,
                int nsteps,
                int n,
                vector <int>& order,
                double& energy,
                double& minenergy,
                int nblocks,
                const vector<int>& kline);
  
void InitialTemp(double& temperature,
		 		 const vector<double>& matrix,
                 const vector<int>& translationTable,
                 const vector<vector<int> >& kernels,
                 int nsteps,
                 int n,
                 vector <int>& order,
                 int nblocks,
                 const vector<int>& kline);

double ExhaustiveStep(vector<double>& matrix,
		      		  int n,
		      		  vector<int>& order,
		      		  const vector<int>& klines,
		      		  int row,
		      		  int nblocks);
  
void ExhaustiveIter(vector<double>& matrix,
				    int n,
				    vector<int>& order,
				    double& energy,
				    const vector<int>& kline,
				    int nblocks);

#endif
