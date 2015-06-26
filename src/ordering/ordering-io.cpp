#ifndef _ordering_io_cpp_included_
#define _ordering_io_cpp_included_

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
#include "ordering-io.hpp"
#include "ordering-tools.hpp"
#include "ordering-templated.hpp"

// namespaces
using namespace std;


//********************************************///

vector<int> OrderVector(int l1,int l2,int w,int n){

  int i;
  vector<int> ordvec(n);

  for(i=0;i<n;i++)
    ordvec[i] = i;

  for(i=0;i<w;i++){

    Swap( ordvec[ l1 + i ], ordvec[ l2 + i ]);

  }

  return ordvec;

}


/****************************************************************/

 vector<int> GetNodeOrder(const vector<int>& kernelOrder,
                          const vector<int>& klines,
			   int netSize, int option){
   //kernelOrder is kernel order
   //netSize is total number of nodes <= KernelOrder.size()
   //option:
   //0: means  nodeOrder[originalposition]=currentposition
   //1: means nodeOrder[currentposition]=originalposition

   int i,j,k,jmax;
   vector<int> nodeOrder(netSize);

   switch(option){

   case 0:

     k=0;
     for(i=0;i<kernelOrder.size();i++){
       
       jmax=netSize-klines[kernelOrder[i]];
       if(kernelOrder[i]!=kernelOrder.size()-1) 
	 jmax=klines[kernelOrder[i]+1]-klines[kernelOrder[i]];
       
       for(j=0;j<jmax;j++){
	 nodeOrder[klines[kernelOrder[i]]+j]=k;
 	 k++;
       }
       
     }

     break;
   case 1:

     k=0;

     for(i=0;i<kernelOrder.size();i++){
       jmax=netSize-klines[kernelOrder[i]];

       if(kernelOrder[i]!=kernelOrder.size()-1) 
	 jmax=klines[kernelOrder[i]+1]-klines[kernelOrder[i]];
       for(j=0;j<jmax;j++){
	 nodeOrder[k]=klines[kernelOrder[i]]+j;

	 k++;
       }
       
     }
   }

   return nodeOrder;

 }

/*********************************************************************/

vector< vector<int> > GetKernels( vector<double>& similarityMatrix,
				  vector<int>& originalOrder,
				  vector <int>& klines, int& n,
				  vector<int>& translationTable ){
  int nn=n,i,j;
  vector<int> kernelsize,knodelist;
  vector<vector<int> > kernelList;//list of nodes inside each kernel
  vector<int> kernelOrder;//order of the kernels
  vector<int> kkline;//list of lines at which each kernel starts
  vector<int> assignedNodes;//previously assigned nodes
  vector<double> similarityMatrixNew (nn*nn, 0.); //we create new matrix reordered by kernels
  vector<int> translationTableNew(nn);//and update translation table to new order

  int kernelCount, nodeCount;
  nodeCount =0;
  kernelCount =0;
  int k;
  for (i=0;i<nn; i++){
//     cout<<i<<" "<<nn<<endl;
    if (Find(i, assignedNodes) ==-1){
      knodelist.push_back(nodeCount);
      kkline.push_back(nodeCount);
      kernelOrder.push_back(kernelCount);
      assignedNodes.push_back(i);
      translationTableNew[nodeCount]=  translationTable[i];	
      nodeCount ++;
      for( j= i + 1; j < nn; j++){

	if(similarityMatrix[ i + j*nn ] == 1.){
	  knodelist.push_back(nodeCount);
	  assignedNodes.push_back(j);
	  translationTableNew[nodeCount] = translationTable[j];	
	  nodeCount ++;
	  
	}
      }
      kernelCount++;
      kernelList.push_back( knodelist );
      kernelsize.push_back( knodelist.size() );
      knodelist.clear();

    }
    
  }

  cout<<"Count "<<nodeCount <<endl;
  cout<<assignedNodes.size()<<endl;
  //Construct new matrix
  for (i=0;i<assignedNodes.size();i++){
    for (j=i;j<assignedNodes.size();j++){
      similarityMatrixNew[i + nn*j]= similarityMatrixNew[j + nn*i] 
	= similarityMatrix[assignedNodes[i] + nn*assignedNodes[j]];
    }
  }


  n = kernelList.size();
  originalOrder = kernelOrder;
  klines = kkline;
  similarityMatrix = similarityMatrixNew;
  translationTable = translationTableNew;

  cout<<"size "<<klines.size()<<endl;
  ofstream fout( "kernel-size.dat" );
  ofstream fout1( "kernel-list.dat" );
  
  int count = 0;
  for(i = 0; i < kernelsize.size(); i++){
    fout<<kernelsize[i]<<endl;
    count += kernelList[i].size();
    for( j=0; j<kernelList[i].size(); j++)
      fout1<<kernelList[i][j]<<" ";
    fout1<<endl;
  }
  
  fout.close();
  fout1.close();


  cout<<"Total n of nodes "<<count<<endl;

  return kernelList;

}

/**************************************************************/

void ReadMatrix(const int format, const char *fileName,
		vector<int>& transTable , vector<double>& sim ) {

  ifstream gin;
  gin.open(fileName);
  double x;
  int n = transTable.size();
  int i, j, count;
  vector <int> transtableTemp;

  cout<<"Reading matrix....";

  switch( format ) { 
    
  case 0 :
    int ii [2];
    
    // determine how many species there are
    count = 0;
    while(gin>>ii[0]>>ii[1]>>x){
      
      i = Find(ii[0], transtableTemp);
      if( i == -1){
	i = count;
	transtableTemp.push_back(ii[0]);
	count ++;
      }
      j = Find(ii[1], transtableTemp);
      if( j == -1){
	j = count;
	transtableTemp.push_back(ii[1]);
	count ++;
      }
    }
    gin.close();

    n = transtableTemp.size();

    sim.clear();
    for(i=0;i<n;i++)
      for(j=0;j<n;j++)
	sim.push_back(0);
    
    gin.open(fileName);
    while(gin>>ii[0]>>ii[1]>>x){
      i = Find(ii[0], transtableTemp);
      j = Find(ii[1], transtableTemp);

      sim[ i + n * j ] = x;
    }

    break;

  case 1:
    char c;
    count = 0;

    // determine how many species there are
    while(gin.get(c)){
      if(c == '\n')
	++count;
    }
    gin.close();
    
    n = count;

    transtableTemp.clear();
    sim.clear();
    for(i=0;i<n;i++){
      transtableTemp.push_back(i);
      for(j=0;j<n;j++)
	sim.push_back(0);
    }

    gin.open(fileName);
    for( i=0; i<n; i++ ){
      
      for( j=0; j<n; j++){
	
	gin >> x;
	 sim[ i + n * j ] = x;
	 
      }
    }

    break;
  } 
  
  transTable = transtableTemp;
  gin.close();
  cout<<"Read\n";
     
}
    
/**************************************************************/

void ReadTranslationTable(const char *fileName, vector<int>& transTable) {

  unsigned int i,j;
  transTable.clear();

  ifstream gin;
  gin.open(fileName);

  cout<<"Reading translation table....";

  while(gin>>i>>j){
    transTable.push_back(j);
  }

  gin.close();

  cout<<"Read\n";
     
}
    
/**************************************************************/

void ReadData(const int format, vector<int>& translationTable , vector<double>& similarityMatrix){

  ifstream inFile;

  if(!fileExists("coclas-final.dat")){
    cout << "Attempting to read in coclassification data but the file coclas-final.dat does not exist.\n";
    cout << "Exiting program.\n";
    exit(1);
  }else{
    ReadMatrix(format, "coclas-final.dat", translationTable, similarityMatrix);
  }

  if(!fileExists("transtable-final.dat")){
    cout << "Attempting to read in the translation table but the file transtable-final.dat does not exist.\n";
    cout << "Exiting program.\n";
    exit(1);
  }else{
    ReadTranslationTable("transtable-final.dat", translationTable);
  }

}

void ReadData(double& temperature,double& factor,unsigned int& iter,
	      const unsigned int format, vector<int>& translationTable,vector<double>& similarityMatrix,
	      double& minenergy){
  ifstream inFile;

  if(fileExists("coclas-final.dat") && fileExists("transtable-final.dat")){
    cout << "This system has already converged on a final solution.\n";
    cout << "Exiting program.\n";
    exit(1);
  }

  if(!fileExists("temperature.dat")){
    cout << "Attempting to read in the previous temperature but the file temperature.dat does not exist.\n";
    cout << "Exiting program.\n";
    exit(1);
  }else{
    inFile.open("temperature.dat");
    inFile >> temperature;
    inFile.close();
  }

  if(!fileExists("temperature-factor.dat")){
    cout << "Attempting to read in the temperature factor but the file temperature-factor.dat does not exist.\n";
    cout << "Exiting program.\n";
    exit(1);
  }else{
    inFile.open("temperature-factor.dat");
    inFile >> factor;
    inFile.close();
  }

  if(!fileExists("coclas-ordert.dat")){
    cout << "Attempting to read in the coclassification data but the file coclas-ordert.dat does not exist.\n";
    cout << "Exiting program.\n";
    exit(1);
  }else{
    ReadMatrix(format, "coclas-ordert.dat", translationTable, similarityMatrix);
  }

  if(!fileExists("transtable-minimum.dat")){
    cout << "Attempting to read in the translation table but the file transtable-minimum.dat does not exist.\n";
    cout << "Exiting program.\n";
    exit(1);
  }else{
    ReadTranslationTable("transtable-minimum.dat", translationTable);
  }

  if(!fileExists("energy-minimum.dat")){
    cout << "Attempting to read in the minimum energy but the file energy-minimum.dat does not exist.\n";
    cout << "Exiting program.\n";
    exit(1);
  }else{
    inFile.open("energy-minimum.dat");
    inFile >> minenergy;
    inFile.close();
  }

  if(!fileExists("iteration.dat")){
    cout << "Attempting to read in the iteration but the file iteration.dat does not exist.\n";
    cout << "Exiting program.\n";
    exit(1);
  }else{
    inFile.open("iteration.dat");
    inFile >> iter;
    inFile.close();
  }

}


/**************************************************************/

void PrintMatrix (const char* fileName, const vector<double>& sim, const vector<int>& kernelOrder,
		              const vector<int>& klines, int n, int mode){
  //mode
  //0: prints in matrix format
  //1: prints in column format (position_node1 position_node2 sim[node1][node2])

  int i,j,i1,j1;

  vector<int> nodeOrder = GetNodeOrder( kernelOrder, klines, n, 1);

//   cout<<"printing"<<endl;

  ofstream fo(fileName); 

  switch(mode){

  case 0://all in a row

    for(i=0;i<n;i++){
      i1=nodeOrder[i];
      for(j=0;j<n;j++){
	j1=nodeOrder[j];
	fo<<sim[i1+j1*n]<<" ";
      }
      fo<<endl;
    }
    break;

  case 1://coclas style
    
    for(i=0;i<n;i++){
      i1=nodeOrder[i];
     
      for(j=0;j<n;j++){
	j1=nodeOrder[j];
	fo<<i <<" "<<j<<" "<<sim[i1+j1*n]<<endl;
      }
   
    }
    break;

  }

  fo.close();

}

void PrintTranstable(const char* fileName, const vector<int>& translationTable, const vector<int>& order, const vector<vector<int> >& kernels){
  unsigned int i,j,k,l;

  k=0;
  ofstream outFile(fileName);
  
  for(i=0;i<order.size();i++){
    l=order[i];
     
    for(j=0;j<kernels[l].size();j++){
      outFile<<k<<" "<<translationTable[kernels[l][j]]<<endl;
      k++;
    }
  }   

  outFile.close();
}

void PrintValue(const char* fileName, const double value){
  ofstream outFile(fileName);
  outFile.precision(16);
  outFile << value << endl;
  outFile.close();
}

#endif
