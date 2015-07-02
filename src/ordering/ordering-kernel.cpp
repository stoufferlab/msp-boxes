
// c++ header files
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <ctime>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

// my header files
#include "ordering-tools.hpp"
#include "ordering-io.hpp"
#include "ordering-annealing.hpp"

// namespaces
using namespace std;

// constants
#define EPS 1e-8

/***************************88

This program reorders the elements of a matrix making first a blocking
into kernels (so it is better that they are hierarchically clustered
first. The program will not perform hierarhchical clustering itself.
Then, the algorith repositions bloks selected at random, the width of
the block is selected so that at low temperature is of size one.

********************************************///


////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////

int main(int argc, char** argv)  {

  time_t start,now;
  time(&start);
  
  unsigned int randomization, format, restart;
  double factor, temperature;
  char *inFile;
  unsigned long seed;

  ORDERINGTOOLS::OptionParser(argc, argv, inFile, temperature, factor, randomization, format, seed, restart);

  //initialize random number generator
  if(seed == 0)
    srand(time(NULL));
  else
    srand(seed);
     
  vector<double> similarityMatrix;
  vector<int> translationTable;
  vector<int> kernelOrder,klines;

  //int i,j,i1,j1;
  unsigned int iter=0;
  double energy,energyold,minenergy=-1;
  //double x;
  int nblocks = 1;

  ofstream fout;

  if(restart){
    ReadData(temperature, factor, iter, format, translationTable, similarityMatrix, minenergy);
    ++iter;
    temperature *= factor;
    fout.open("outt.dat",ifstream::app);
  }else{
    ReadMatrix(format, inFile, translationTable, similarityMatrix);
    PrintValue("temperature-factor.dat",factor);
    fout.open("outt.dat");
  }

  vector<vector<int> > kernels;
  unsigned int netSize = translationTable.size();
  int nkernels=netSize;
  kernels = GetKernels( similarityMatrix, kernelOrder, klines, nkernels, translationTable);

  cout<<"Reduced matrix size "<<nkernels<<endl;
  cout<<"Reduced matrix size "<<kernelOrder.size()<<endl;
  cout<<"Reduced matrix size "<<klines.size()<<endl;
 
  if(randomization){
    kernelOrder = RandomizeOrder(nkernels, nblocks);
    cout<<"Randomized\n";
  }
    
  int count=0;
  int Nsteps = (int)(0.01*nkernels*nkernels/nblocks/nblocks);
  if (Nsteps < 10) Nsteps = 10;//at least 10 steps in annealing process
  //cout << "Nteps: "<<Nsteps<<endl;

  energyold = energy = ComputeEnergy( similarityMatrix, netSize, kernelOrder, klines);
  if(minenergy == -1)
    minenergy = energy;

  cout<<"e: "<<energy<<endl;

  // if the temperature has not been specified in some manner, estimate a good starting point
  if(temperature == -1){
    InitialTemp(temperature, similarityMatrix, translationTable, kernels, 50, netSize, kernelOrder, nblocks, klines);
  }

  if (nkernels>2){

    while(count<20){

      AnnealIter( temperature, similarityMatrix, translationTable, kernels, Nsteps,
		  netSize, kernelOrder, energy, minenergy, nblocks, klines);
     
      fout.precision(16);
      fout<<iter<<" "<<energy<<" "<<temperature<<" "<<minenergy<<endl;
      fout.flush();

      if( fabs( energyold - energy ) < EPS ){
	count ++;
      }
      else{
	count = 0;
      }
     
      energyold = energy;

      //PrintMatrix( "coclas-ordert.dat", similarityMatrix, kernelOrder, klines, netSize, 0);
      PrintTranstable( "transtable-ordert.dat", translationTable, kernelOrder, kernels);
      PrintValue( "temperature.dat", temperature);
      PrintValue( "iteration.dat", iter);

      time(&now);

      //if(difftime(now,start) > 5)
	    //  exit(1);

      iter ++;
      temperature *= factor;

    }
  }

  cout<<"e: "<<energy<<endl;
   
//print matrix in format for level.py

//     PrintMatrix(fileName,sim,order,klines,n,1);
  PrintMatrix( "coclas-final.dat", similarityMatrix, kernelOrder, klines, netSize, 0);
  PrintTranstable( "transtable-final.dat", translationTable, kernelOrder, kernels);
  return 0;
 
}


