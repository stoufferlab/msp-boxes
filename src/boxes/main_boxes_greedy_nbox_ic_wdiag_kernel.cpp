#include <fstream>
#include <iostream>
#include <string>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <vector>
#include <algorithm>
using namespace std;

#include "greedy.hpp"
//In this version, the breaking point can only be between kernels

int main(int argc, char** argv)  {

   if(argv[1]==NULL){
     printf("Please write in the command line:\n");
     printf("1--the size of the network\n");
     printf("2--the datafile that contains the similarity matrix\n");
     printf("3-- the format of the data file\n");
     printf("format '0':  if the data is stored i,j,x\n");
     printf("format '1':  if the data is stored x11,..,..,..,..,x1n\n");
     printf("                          xn1,..,..,..,..,xnn\n");
     printf("4--the minimum node label\n");
     exit(0);
    } 

  int netSize = atoi(argv[1]);
  char *inFile = argv[2];
  int type = atoi(argv[3]),ii[6];
  int nodeMin = atoi(argv[4]);

  vector <double> similarityMatrix(netSize*netSize,0.);
  double average=0.,variance=0.;

  ReadMatrix( type, inFile, &similarityMatrix, nodeMin,	
	      netSize, &average, &variance);
  
  vector <int> kernelRows;
  kernelRows = GetKernelRows (similarityMatrix, netSize);
  int reducedNetSize = kernelRows.size();//number of different kernels

//   cout<< "number of kernels "<<reducedNetSize<<endl;
  /************************************************************/

  vector<vector<int> > partition, partMin, partOne;//nodeini and nodelast of each box
  vector <vector<vector<int> > > partHistory;
  vector <double> lSqHistory;
  vector<double > leastSquares(netSize,0.);
  double leastSq, leastSqMin, leastSqOne;

  vector<int> label(1,1);//store whether this partition was split 
  //or not in the previous iteration
  
  //partitions are stored by the kernel number, 
  //in kernel Rows we have the info about the "reduce size of the matrix 
  //to explore"
  
  InitializePartition(1, reducedNetSize, &partition);
  leastSqOne = leastSqMin = leastSq = (variance-average*average);
  lSqHistory.push_back(leastSqOne);
  partMin = partOne = partition;
  partHistory.push_back(partOne);
  
//   cout<<"Finding best partition... \n";

  FindOptimalPartition( &partition, &label, similarityMatrix, netSize,
			&leastSq, &partHistory, &lSqHistory, kernelRows);
  partMin = partition;
  leastSqMin = leastSq ;
  
  FindOptimalPartition(&partition, &label, similarityMatrix, netSize, 
		       &leastSq, &partHistory, &lSqHistory, partition.size() + 10,
		       kernelRows);
  

  while (leastSq < leastSqMin){
//     cout<<"minimizing\n";
    leastSqMin = leastSq;
    partMin = partition;
    FindOptimalPartition(&partition, &label, similarityMatrix, netSize, 
			 &leastSq, &partHistory, &lSqHistory, partition.size() + 10,
			 kernelRows);
  
  }

  partMin = GetBestPartitionBIC( &leastSqMin, partHistory, lSqHistory,
  				 (netSize*netSize + netSize)/2);

  int i; //printing output
  int iniNode, endNode;
  for(i=0;i<partMin.size();i++){
    iniNode = kernelRows[ partMin[i][0]];
    endNode = netSize - 1;
    if (i != partMin.size() -1 )
      endNode = kernelRows[ partMin[i + 1][0]] - 1;
    cout<<iniNode+ nodeMin<<"--"<<endNode + nodeMin<<" ";
  }
  cout<<endl;
  
  return 0;
  

}











