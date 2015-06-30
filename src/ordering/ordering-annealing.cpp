#ifndef _ordering_annealing_cpp_included_
#define _ordering_annealing_cpp_included_

// c++ header files
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <ctime>
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

/*****************************

This program reorders the elements of a matrix making first a blocking
into kernels (so it is better that they are hierarchically clustered
first. The program will not perform hierarhchical clustering itself.
Then, the algorith repositions bloks selected at random, the width of
the block is selected so that at low temperature is of size one.

********************************************///

int ChangeSegmentOrder(vector<int>& order, int line1, int line2, int width){

  unsigned int i, j;
  vector<int> pos(width);
  
  for(i=0;i<width;i++){
    pos[i] = order[line1+i];
  }

  if(line1 < line2){
    j=0;
    for( i=line1+width; i<line2+width; i++){
      order[ line1 + j ] = order[i];
      j++;
    }
    for(i=0;i<width;i++)
      order[ line2 + i ] = pos[i];
  }else
    if( line1 > line2 ){
      for( i=line1-1; i>=line2; i--){
        order[ i + width ] = order[i];
      }

      for( i=0; i<width; i++)
        order[ line2 + i ] = pos[i];
    }

  return 1;
}      

/****************************************************************/

int MakeChangeVectors(const vector <int>& order,
                      vector<bool>& change,
                      int line1,
                      int line2, 
		                  int width,
                      int n,
                      const vector<int>& klines){

  int i,in,j,jmax, nk = klines.size();
  vector<int> pos;
  
  for(i=0; i<width; i++){
    in = order[line1+i];
    jmax = n - klines[in];
    
    if(in != nk-1)
      jmax = klines[in+1] - klines[in];

    for(j=0;j<jmax;j++)
      change[klines[in] + j] = true;
  }
    
  if(line1<line2){
    for(i=line1+width; i<line2+width; i++){
      in = order[i];
      jmax = n-klines[in];
      
      if( in != nk-1 )
        jmax = klines[in+1] - klines[in];

      for( j=0; j<jmax; j++)
        change[klines[in] + j] = true;
    }
  }else
    if(line1 > line2){
      for( i=line2; i<line1; i++){
        in = order[i];
        jmax = n - klines[in];

        if( in != nk-1 )
          jmax = klines[in+1] - klines[in];
      
        for(j=0; j<jmax; j++)
          change[klines[in] + j] = true;
      }
    }

  /*for(i=0; i<n; i++){
    if(find(change.begin(),change.end(),i) == change.end()){
      nochange.push_back(i);
    }
  }*/

  return 1;
}      

/*************************************************************/

double ComputeEnergy(const vector<double>& mat, int netSize, const vector<int>& kernelOrder,
		                 const vector<int>& klines){
  //kernelOrder: vector with the order of kernels
  //klines: line in original matrix where kernel starts
  int i,j;
  double energy = 0.;
  vector<int> nodeOrder = GetNodeOrder( kernelOrder, klines, netSize, 1);
  
  for(i=0;i<netSize-1;i++){
    for(j=i+1;j<netSize;j++){
      energy += double(abs(i-j)) * mat[ nodeOrder[i] + netSize* nodeOrder[j] ];
      // cout<<i<< " "<<j<<" "<< nodeOrder[i] << " "<< netSize<<endl;
      // cout<<i<< " "<<j<<" "<<mat[ nodeOrder[i] + netSize* nodeOrder[j] ]<<endl;
    }
  }

  energy = 2*energy/double(netSize);

  return energy;

}

/*****************************************************************/
/*****************************************************************/

double EnergyChangeReOrd(const vector<double>& mat,
                         int l1,
                         int l2,
                         int w,
                         int n, 
			                   const vector<int>& kernelOrder,
                         const vector<int>& klines){

  double deltaen = 0;
  unsigned int i,j;
  
  // figure out what the new order would be
  vector<int> newKernelOrder(kernelOrder);
  ChangeSegmentOrder(newKernelOrder,l1,l2,w);

  if(newKernelOrder != kernelOrder){
    //make a vector that is true for the rows that change position and false for the rows that don't
    vector<bool> lchange(n,false);
    MakeChangeVectors(kernelOrder, lchange, l1, l2, w, n, klines);

    //given a matrix line, it tells me the position of the nodes in the
    //"virtual" ordered matrix
    vector<int> oldorder = GetNodeOrder(kernelOrder,klines,n,0);
    vector<int> neworder = GetNodeOrder(newKernelOrder,klines,n,0);
  
    // compute the change in energy between the changed lines and fixed lines
    for(i=0;i<lchange.size();i++){
      for(j=i+1;j<lchange.size();j++){
        if(lchange[i] || lchange[j]){
          deltaen+=mat[i+j*n]*(double)(abs(neworder[i]-neworder[j])-abs(oldorder[i]-oldorder[j]));
        }
      }
    }
  
    // add the appropriate scaling
    deltaen = 2*deltaen/double(n);
  }
  
  return deltaen;
}

/**************************************************************/

double AnnealStep(double temperature,
                  const vector<double>& matrix,
                  int netSize,
		              vector <int>& order,
                  int nblocks,
                  const vector<int>& klines)
  //Only reordering matrix rordering matrix
{
  
  int line1,line2,width, nkernels=klines.size();
  double deltaE;

//    double energy,energyold;
//    energyold=energy=ComputeEnergy(*mat, n, *order, klines);

   //1 pick two random numbers (starting line that we are swapping)
  
  line1 = line2 = width = 0;
  
  //"Computing line and width\n";

  width = (int)( fabs(gauss(sqrt(temperature)*.05*nkernels) ) ) + nblocks;
//   width = nblocks;
  width -= (width%nblocks);
  

  line1 = width/2 + (int)( (nkernels-width)*ran1f() );
  line1 -= width/2;
  line1 -= (line1%nblocks);
  
  if ( line1 < 0 ) line1 = 0;
  

  while( width > nkernels-line1 ){
    
    width = (int)( fabs( gauss( sqrt(temperature)*.05*nkernels ) ) ) + nblocks;
    width -= (width%nblocks);
      
  }

  line2=(int)(gauss(sqrt(temperature)*.1));
  
  
  if( line2 < 0 ){

    line2 -= line1;
    if( line2 < 0 ) line2=0;
    
  }
  else{
    
    line2 += line1 + width;
    if ( line2 > nkernels-width ) line2 = nkernels-width;
    
  }
  
  line2 -= (line2%nblocks);


  deltaE = EnergyChangeReOrd(matrix, line1, line2, width, netSize, order, klines);
 
  double prob;

  prob = exp( -deltaE/temperature );

  if( ranf() < prob ){
    
    ChangeSegmentOrder( order, line1, line2, width);
    return deltaE;
  }	 
  
  return 0;
  
  
}

// Make an annealing step of reordering the matrix but without
// worrying about acceptance temperature
//
// This is in order to get a good starting temperature
//
double AnnealStepFree(double temperature,
                      const vector<double>& matrix,
                      int netSize,
		                  vector<int>& order,int nblocks,
                      const vector<int>& klines)
  //Only reordering matrix rordering matrix
{
  
  int line1,line2,width, nkernels=klines.size();
  double deltaE;

//    double energy,energyold;
//    energyold=energy=ComputeEnergy(*mat, n, *order, klines);

   //1 pick two random numbers (starting line that we are swapping)
  
  line1 = line2 = width = 0;
  
  //"Computing line and width\n";

  width = (int)( fabs(gauss(sqrt(temperature)*.05*nkernels) ) ) + nblocks;
//   width = nblocks;
  width -= (width%nblocks);
  

  line1 = width/2 + (int)( (nkernels-width)*ran1f() );
  line1 -= width/2;
  line1 -= (line1%nblocks);
  
  if ( line1 < 0 ) line1 = 0;
  

  while( width > nkernels-line1 ){
    
    width = (int)( fabs( gauss( sqrt(temperature)*.05*nkernels ) ) ) + nblocks;
    width -= (width%nblocks);
      
  }

  line2=(int)(gauss(sqrt(temperature)*.1));
  
  
  if( line2 < 0 ){

    line2 -= line1;
    if( line2 < 0 ) line2=0;
    
  }
  else{
    
    line2 += line1 + width;
    if ( line2 > nkernels-width ) line2 = nkernels-width;
    
  }
  
  line2 -= (line2%nblocks);


  deltaE = EnergyChangeReOrd(matrix, line1, line2, width, netSize, order, klines);
 
  ChangeSegmentOrder( order, line1, line2, width);
  return deltaE;
  
}

/**************************************************************/

void AnnealIter(double temperature,
                const vector<double>& matrix,
                const vector<int>& translationTable,
                const vector<vector<int> >& kernels,
		            int nsteps,int n,
                vector<int>& order,double& energy,double& minenergy,int nblocks, 
		            const vector<int>& kline){
  
  int step;
  double deltaen;

  for( step=0; step<nsteps; step++ ){
    deltaen = AnnealStep( temperature, matrix, n, order, nblocks, kline);
    energy += deltaen;

    if(energy < minenergy ){
      minenergy = energy;
      //PrintMatrix( "coclas-minimum.dat", matrix, order, kline, n, 0);
      PrintTranstable( "transtable-minimum.dat", translationTable, order, kernels);
      PrintValue( "energy-minimum.dat", minenergy);
    }
  }
   	     
}

/**************************************************************/

void InitialTemp(double& temperature,
            		 const vector<double>& matrix,
                 int nsteps,
                 int n,
                 const vector<int>& order,
                 int nblocks,
                 const vector<int>& kline){
  
  unsigned int i,step;
  double deltaen;
  double initialProbability = 0.95;
  vector<int> neworder(order);

  temperature = ComputeEnergy(matrix, n, neworder, kline);
  for(i=0;i<20;i++){
    deltaen = 0;
    for( step=0; step<nsteps; step++ ){
        deltaen += fabs(AnnealStepFree(temperature, matrix, n, neworder, nblocks, kline));
    }

    deltaen /= double(nsteps);
    temperature = deltaen/(-log(initialProbability));
  }

}

/*****************************************************************/
/*****************************************************************/

double ExhaustiveStep(vector<double>& matrix,int n,
            		      vector<int>& order,
                      const vector<int>& klines,
                      int row, int nblocks)
{
  
  int line1, line2, width, nk=klines.size();
  double deltaE,maxdeltaE=0;
  int line=-1,wl=0;

  line1=row;

  //cout<<"Iterating \n"<<endl;
  for (width = nblocks; width>0 ; width--){
    if(line1+width < nk){
      for (line2=0;line2<= nk-width;line2++){
        if(line2<line1 || line2 >= line1+width){
          deltaE = EnergyChangeReOrd(matrix, line1, line2, width, n, order, klines);

          if(deltaE<maxdeltaE){
            maxdeltaE = deltaE;
            line = line2;
            wl = width;
          }
        }
      }
    }
  }
      
  if(line !=-1){
    ChangeSegmentOrder(order,line1,line,wl);
    return maxdeltaE;
  }
  else{
    return 0;
  }
}


void ExhaustiveIter(vector<double>& matrix,int n,
		                vector<int>& order,
                    double& energy,
                    const vector<int>& kline,int nblocks){

  int step,nsteps=kline.size();
  double deltaen;
  int line;

  for(step=0;step<nsteps;step++){
    ///we try to change the kernel originally at position 0, line is ist current position
    line = find(order.begin(),order.end(),step) - order.begin();
    deltaen = ExhaustiveStep( matrix, n, order, kline, line, nblocks);
    energy+=deltaen;
  }
}

#endif
