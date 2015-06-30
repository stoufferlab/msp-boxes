
//c++ header files
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

// my header files
#include "ordering-annealing.hpp"
#include "ordering-io.hpp"
#include "ordering-tools.hpp"

// namespaces
using namespace std;

/****************************

This program reorders the elements of a matrix making first a blocking
into kernels

********************************************/

int main(int argc, char** argv)  {

  unsigned int format, nblocks, restart;
  char *inFile, *translationFile;

  ORDERINGTOOLS::OptionParserExhaustive(argc, argv, inFile, translationFile, format, nblocks, restart);

  vector<double>  similarityMatrix;
  vector<int> translationTable;
  vector<int> kernelOrder,klines;

  double energy;

  if(restart){
    ReadData(format, translationTable, similarityMatrix);
  }else{
    ReadMatrix(format, inFile, translationTable, similarityMatrix);
  }
   
  vector<vector<int> > kernels;
  unsigned int netSize = translationTable.size();
  int nkernels=netSize;
  kernels = GetKernels(similarityMatrix, kernelOrder, klines, nkernels, translationTable);
  
  cout<<"Reduced matrix size "<<nkernels<<endl;
  cout<<"Reduced matrix size "<<kernelOrder.size()<<endl;
  cout<<"Reduced matrix size "<<klines.size()<<endl;

  energy = ComputeEnergy(similarityMatrix, netSize, kernelOrder, klines);
  cout.precision(16);
  cout<<"e: "<<energy<<endl;

  ExhaustiveIter(similarityMatrix, netSize, kernelOrder, energy, klines, nblocks);
  cout<<"After exhaustive "<<energy<<endl;
   
  PrintMatrix( "coclas-exhaustive.dat", similarityMatrix, kernelOrder, klines, netSize, 0);
  PrintTranstable( "transtable-exhaustive.dat", translationTable, kernelOrder, kernels);
  return 0;
   
}


