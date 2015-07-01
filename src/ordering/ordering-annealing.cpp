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
#include <iterator>
#include <string>
#include <utility>
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

*****************************/

int ChangeSegmentOrder(vector<int>& order, int line1, int line2, int width){

  if( line1 < line2 )
    rotate(order.begin()+line1,order.begin()+line1+width,order.begin()+line2+width);
  else
    rotate(order.begin()+line2,order.begin()+line1,order.begin()+line1+width);

  return 1;
}      

/****************************************************************/

int MakeChangeVectors(const vector <int>& order,
                      vector<int>& change,
                      vector<int>& nochange,
                      int line1,
                      int line2,
                      int width,
                      int n,
                      const vector<int>& klines){

  int i,in,in2,j,jmax, nk = klines.size();
  vector<bool> changed(n,false);
  
  change.clear();
  change.reserve(n);

  for(i=line1; i<line1+width; i++){
    in = order[i];
    
    if(in != nk-1)
      jmax = klines[in+1] - klines[in];
    else
      jmax = n - klines[in];  

    for(j=0;j<jmax;j++){
      change.push_back( klines[in] + j );
      changed[klines[in]+j] = true;
    }

  }
    
  if(line1<line2){

    for(i=line1+width; i<line2+width; i++){
      in = order[i];
      
      if(in != nk-1)
        jmax = klines[in+1] - klines[in];
      else
        jmax = n-klines[in];

      for(j=0; j<jmax; j++){
        change.push_back( klines[in] + j );
        changed[klines[in]+j] = true;
      }

    }

  }

  if( line1 > line2 ){

    for(i=line2; i<line1; i++){
      in = order[i];
      
      if(in != nk-1)
        jmax = klines[in+1] - klines[in];
      else
        jmax = n - klines[in];
      
      for(j=0; j<jmax; j++){
        change.push_back( klines[in] + j );
        changed[klines[in]+j] = true;
      }

    }

  }

  nochange.clear();
  nochange.reserve(n-change.size());
  for(i=0;i<n;i++){
    if(!changed[i]){
      nochange.push_back(i);
    }
  }

  return 1;

}      

/*************************************************************/

double ComputeEnergy(const vector<double>& mat, int netSize, const vector<int>& kernelOrder,
		                 const vector<int>& klines){
  //kernelOrder: vector with the order of kernels
  //klines: line in original matrix where kernel starts
  int i,j;
  double energy = 0.;
  vector<int> nodeOrder = GetNodeOrder(kernelOrder, klines, netSize, 1);
  
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

double EnergyChangeReOrd(const vector<double>& mat,int l1,int l2,int w,int n, 
			                   const vector<int>& kernelOrder,
                         const vector<int>& klines,
                         vector<int>& newKernelOrder){

  double deltaen = 0;
  unsigned int i,j;
  int line1,line2;
  int do1, do2, oo1, no1;

  // determine what the new order will be
  ChangeSegmentOrder(newKernelOrder,l1,l2,w);
  
    // only calculate a change in energy when some lines actually changed
  if(newKernelOrder != kernelOrder){

    // determine what the node by node order is based on their kernels
    vector<int> oldorder,neworder;
    oldorder=GetNodeOrder(kernelOrder,klines,n,0);
    neworder=GetNodeOrder(newKernelOrder,klines,n,0);

    //make two vectors, one with the rows that
    //change position and one with the rows
    //that keep the same position
    vector<int> lchange,lfix;
    MakeChangeVectors(kernelOrder, lchange, lfix, l1, l2, w, n, klines);

    // calculate the change in energy only for pairs that are altered by change    
    for(i=0;i<lchange.size();i++){
      line1 = lchange[i];
      oo1 = oldorder[line1];
      no1 = neworder[line1];
      for(j=0;j<lfix.size();j++){
        line2=lfix[j];
        do1 = abs(oo1-oldorder[line2]);
        do2 = abs(no1-neworder[line2]);
        deltaen+=mat[line1+line2*n]*(do2-do1);
      }

      for(j=i+1;j<lchange.size();j++){
        line2=lchange[j];
        do1 = abs(oo1-oldorder[line2]);
        do2 = abs(no1-neworder[line2]);
        deltaen+=mat[line1+line2*n]*(do2-do1);
      }
    }

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
                  const vector<int>& klines,
                  bool freeStep)
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

  vector<int> neworder(order);
  deltaE = EnergyChangeReOrd(matrix, line1, line2, width, netSize, order, klines, neworder);
 
  if(freeStep || ranf() < exp( -deltaE/temperature )){
    order = neworder;
    return deltaE;
  }

  return 0;
}

/**************************************************************/

void AnnealIter(double temperature,
                const vector<double>& matrix,
                const vector<int>& translationTable,
                const vector<vector<int> >& kernels,
		            int nsteps,int n,
                vector<int>& order,double& energy,double& minenergy,int nblocks, 
		            const vector<int>& kline){
  
  int step,i1,ii1;
  double deltaen;

  for( step=0; step<nsteps; step++ ){
    deltaen = AnnealStep(temperature, matrix, n, order, nblocks, kline, false);
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
                 const vector<int>& translationTable,
		             const vector<vector<int> >& kernels,
                 int nsteps,int n,
                 vector<int>& order,
                 int nblocks,
                 const vector<int>& kline){
  
  int step;
  double deltaen;
  double initialProbability = 0.95;
  int i;
  vector<int> itorder(order);

  temperature = ComputeEnergy(matrix, n, itorder, kline);
  for(i=0;i<20;i++){
    deltaen = 0;
    for( step=0; step<nsteps; step++ ){
        deltaen += fabs(AnnealStep(temperature, matrix, n, itorder, nblocks, kline, true));
    }

    deltaen /= double(nsteps);
    temperature = deltaen/(-log(initialProbability));
  }

  order = itorder;
}

/*****************************************************************/
/*****************************************************************/

double ExhaustiveStep(const vector<double>& matrix,int n,
            		      vector<int>& order,
                      const vector<int>& klines,
                      int row, int nblocks)
{
  
  int line1, line2, width, nk=klines.size();
  double deltaE,maxdeltaE=0;
  int line=-1,wl=0;
  vector<int> maxorder, neworder(order);

  line1=row;

  //cout<<"Iterating \n"<<endl;
  for (width = nblocks; width>0 ; width--){

    if( line1+width < nk ){

      for (line2=0;line2<= nk-width;line2++){

        if(line2<line1 || line2 >= line1+width){

          deltaE = EnergyChangeReOrd(matrix, line1, line2, width, n, order, klines, neworder);

          if(deltaE<maxdeltaE){
            maxdeltaE = deltaE;
            line = line2;
            wl = width;
            maxorder = neworder;
          }
        }
      }
    }
  }
      
  if(line !=-1){
    order = maxorder;
    return maxdeltaE;
  }
  else{
    return 0;
  }
}


void ExhaustiveIter(const vector<double>& matrix,
                    int n,
		                vector<int>& order,
                    double& energy,
                    const vector<int>& kline,
                    int nblocks){

  int step,nsteps=kline.size();
  double deltaen;
  int line;

  for(step=0;step<nsteps;step++){
    ///we try to change the kernel originally at position 0, line is ist current position
    line = find(order.begin(),order.end(),step) - order.begin();
    deltaen = ExhaustiveStep(matrix, n, order, kline, line, nblocks);
    energy+=deltaen;
  }
}



#endif
