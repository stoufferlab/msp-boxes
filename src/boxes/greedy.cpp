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

#define PI (4*atan(1.))

/********************************************************************************/
/*random number generators and miscellaneous fucntions***************************/
/********************************************************************************/

//uniform number between [0,1]
double ranf(){

  return (double)(rand())/(RAND_MAX);

}

//uniform number between [0,1)
double ran1f(){

  return (double)(rand())/(RAND_MAX+1.);

}

//generates Gaussian distributed variable
double gauss(double variance){
  

  double x1=ranf();
  double x2=ranf();

  return variance*sqrt( - 2.*log(x1) )* cos( 2.* PI* x2 );

}


int ipow(int a, int b){

  if( b ==0)
    return 1;
  else
    return ipow(a, b-1);

}
int FindNumber(int n,vector<int> nlist){
  
  int i;
 
  for(i=0;i<nlist.size();i++)
    if(nlist[i]==n)return i;
 
  return -1;
}



/********************************************************************************/
/*Information Criterion**********************************************************/
/********************************************************************************/

//return Akaike's criterion
double GetAIC ( double L, int n, int kp){
  //L: sum of residual squares
  //n = number of data points
  //k = number of parameters
  int k = 2*kp + 1;
  return n * log(L) + 2 *  k  + (double)( 2 * k * ( k + 1 ) ) / (double)(n - k - 1);

}  

//returns Bayesian Information Criterion
double GetBIC ( double L, int n, int kp){
  //L: sum of residual squares
  //n = number of data points
  //k = number of parameters
  int k = 2 * kp + 1;
  return n * log(L) + k * log( (double) n );
	
}  

/* Gets the best partition for partitionHistory according to Aiakike's
   information criterion*/
vector<vector<int> > GetBestPartitionAIC(double *leastSq,vector < vector< vector <int> > > partitionHist, 
					 vector< double > lSqHist, int n)
{

  int i;
  double min=10000.,max=0, aic, leastSquares;

  for (i=0;i<lSqHist.size();i++){
    aic = GetAIC ( lSqHist[i], n , partitionHist[i].size()); 

    if(aic < min) min = aic;
    if(aic >max) max = aic;;
//        cout<<i+1<<" "<<chist[i]<<" "<<min<<" "<<max<<endl;
    cout <<"AIC "<< partitionHist[i].size()<< " "<< aic<<endl;
  }

  for (i=0;i<lSqHist.size();i++){
    aic = GetAIC ( lSqHist[i], n , partitionHist[i].size()); 
    leastSquares = (aic - min);
//        cout<<"chirel "<<chi<<endl;
    if (leastSquares == 0.){
      *leastSq = aic;
      
      return partitionHist[i];
    }
  }
  
  return partitionHist[partitionHist.size()-1];
}

/* Gets the best partition for partitionHistory according to the Bayesian
   information criterion*/
vector<vector<int> > GetBestPartitionBIC(double *leastSq,vector < vector< vector <int> > > partitionHist, 
					 vector< double > lSqHist, int n)
{

  int i, part = -1;
  double min = 10000., max = 0., bic, leastSquares;

   for (i=0;i<lSqHist.size();i++){
    bic = GetBIC ( lSqHist[i], n , partitionHist[i].size()); 

    if(bic < min) {
      min = bic;
      part = i;
    }
    if(bic >max) max = bic;;

  }

  if (part != -1)
    return partitionHist[part];

  return partitionHist[partitionHist.size()-1];
}
 
  


/********************************************************************************/
/*Functions for Data Processing**************************************************/
/********************************************************************************/

//Reads input matrix from filename using two possible formats
//0: column format x y mat_element
//1: matrix format mat_el_11 mat_el_12 .... mat_el_1netSize
//ad computes the average and variance of the upper triangle including the 
//diagonal
//OK
void ReadMatrix( int format, char *fileName,
	vector <double> *sim, int nodeMin,	
	int netSize, double *average, double *variance) {
  
  
  ifstream gin(fileName);
  double x;
  
//   cout<<"Reading....\n";

  switch( format ) { 
    
  case 0 : 
    int ii [6];
    
    while(gin>>ii[0]>>ii[1]>>x){		
      	(*sim)[ ii[0] - nodeMin + netSize * (ii[1] - nodeMin) ] = x;
      	(*sim)[ ii[1] - nodeMin + netSize * (ii[0] - nodeMin) ] = x;
	(*average) += x* (1+bool(ii[0]==ii[1]));		
	(*variance) += x*x * (1+bool(ii[0]==ii[1]));		
    }
    break;
  case 1:

    int i,j;

     for( i=0; i<netSize; i++)
       for( j=0; j<netSize; j++){
	 gin >> x;
	  (*sim)[i+netSize*j]=x;
	  (*sim)[j+netSize*i]=x;
	  (*average)+=x* (1+bool(i==j));
	  (*variance)+=x*x* (1+bool(i==j));

       }
     break;
 } 
  
 gin.close();

 (*average) /= (double)( netSize*netSize + netSize ) ;
 (*variance) /= (double)( netSize*netSize + netSize ); 
//   cout<<"Read\n";
     
}

/* Gets starting row of each kernel in the matrix. Because breaking points
   cannot be found inside kernels, it speeds up the process of finding
   the boxes.*/
vector< int > GetKernelRows( vector<double> similarityMatrix, int netSize){

  vector<int> kernelRows;
  int row, col, nrow;
//   cout<<"Netsize "<<netSize<<endl;
  nrow = row = 0;
  kernelRows.push_back( row );
  while (row < netSize -1){
//     cout<<row<<endl;
    for( col=row+1; col<netSize; col++){
      row = col;
      if( similarityMatrix[ nrow + col*netSize ] != 1.){
	kernelRows.push_back( col );
	nrow = col;
	break;
      }
    }
  }
//   cout<<"Done\n";
  return kernelRows;
}

 
/********************************************************************************/
/*Functions for computing Least Sqaures to box model*/
/********************************************************************************/

double ComputeLeastSquares(vector<double> similarityMatrix,
			   vector< vector<int> > partition,
			   int netSize, vector<int> kernelRows){

  int i,j,k;
  int row, col;
  int numberBoxes=partition.size();
  double leastSquares=0.;
  double boxAverage=0., boxAverageSq=0.;
  double AverageTotal=0., AverageSqTotal=0.;
  double outsideAverage = 0., outsideAverageSq = 0.;
  int inCount=0, outCount=0,totalCount=0;
 
  int nodesInKerneli;
  int nodesInKernelk;
//   cout<<"pp: "<<partition[0][0]<<" "<<partition[0][1]<<endl;
//   cout<<"pp: "<<partition[1][0]<<" "<<partition[1][1]<<endl;

  for(j=0;j<numberBoxes;j++){
    inCount=0;
    boxAverage=0;
    boxAverageSq=0;
    for(i=partition[j][0];i<=partition[j][1];i++){
      row = kernelRows[i];
      nodesInKerneli = netSize - row;
      if (i < kernelRows.size()-1)
	nodesInKerneli = kernelRows[i+1] - row;
      //cout<<"kernel "<<i<<" "<<nodesInKerneli<<endl;

      inCount += nodesInKerneli * (nodesInKerneli + 1)/2;
      boxAverage += similarityMatrix[row+row*netSize] * nodesInKerneli*(nodesInKerneli + 1)/2;
      boxAverageSq += similarityMatrix[row+row*netSize]*similarityMatrix[row+row*netSize] * nodesInKerneli*(nodesInKerneli + 1)/2;

      for(k=i+1;k<kernelRows.size();k++){
	col = kernelRows[k];
	nodesInKernelk = netSize - col;
	if (k < kernelRows.size()-1)
	  nodesInKernelk = kernelRows[k+1] - col;
	//cout<<"kernel "<<k<<" "<<nodesInKernelk<<endl;

	if(k<=partition[j][1]){
	  inCount += nodesInKerneli * nodesInKernelk;
	  boxAverage += similarityMatrix[row+col*netSize] * nodesInKerneli*nodesInKernelk;
	  boxAverageSq += similarityMatrix[row+col*netSize]*similarityMatrix[row+col*netSize] * nodesInKerneli*nodesInKernelk;
	}
	else{
	  outCount += nodesInKerneli * nodesInKernelk;
	  outsideAverage += similarityMatrix[row+col*netSize] * nodesInKerneli*nodesInKernelk;
	  outsideAverageSq += similarityMatrix[row+col*netSize]*similarityMatrix[row+col*netSize] * nodesInKerneli*nodesInKernelk;
	}
      }
    }
    
    totalCount += inCount;
    AverageSqTotal += boxAverageSq;
    if(inCount>0)
      AverageTotal +=  boxAverage * boxAverage/double(inCount);
    //cout<<boxAverage<<" "<< inCount<<endl;

  }

  //cout<<outsideAverage<<" "<<outsideAverageSq<<" "<< outCount<<endl;

  //cout<<"boxo "<<avto<<" "<<counto<<endl;

  totalCount += outCount;
  AverageSqTotal += outsideAverageSq;
  if(outCount>0)
    AverageTotal += outsideAverage * outsideAverage/double(outCount);

  leastSquares = AverageSqTotal - AverageTotal;

  return leastSquares/(double)(totalCount);

}

double ComputeLeastSquares(vector<double> similarityMatrix, vector<vector <int> > partition,
			   vector<double> *mean,vector<int> *nCount,int netSize,
			   vector<int> kernelRows){

  int i,j,k;
  int numberBoxes=partition.size();
  double leastSquares=0;
  double boxAverage=0., boxAverageSq=0;
  double outAverage=0., outAverageSq=0;
  double AverageTotal=0., AverageSqTotal=0.;
  int inCount=0, outCount=0, totalCount = 0;

  int nodesInKerneli;
  int nodesInKernelk;
  int row, col;
  vector<double> tempMean;
  vector<int> tempCount;
  mean->clear();
  nCount->clear();

  for(j=0;j<numberBoxes;j++){
    inCount=0;
    boxAverage=0;
    boxAverageSq=0;
//     cout<<partition[j][0]<<" "<<partition[j][1]<<endl;
    
    for(i=partition[j][0];i<=partition[j][1];i++){
      row = kernelRows[i];
      nodesInKerneli = netSize - row;
      if (i < kernelRows.size()-1)
	nodesInKerneli = kernelRows[i+1] - row;
      //consider contribution of own kernel
      inCount += nodesInKerneli * (nodesInKerneli + 1)/2;
//       cout<<"kernel "<<i<<" "<<nodesInKerneli<<" "<<inCount<<endl;
      boxAverage += similarityMatrix[row+row*netSize]
	*nodesInKerneli*(nodesInKerneli + 1)/2;
      boxAverageSq += 
	similarityMatrix[row+row*netSize]*similarityMatrix[row+row*netSize] 
	*nodesInKerneli *(nodesInKerneli + 1)/2;
      
      for(k=i+1;k<kernelRows.size();k++){
	col = kernelRows[k];
	nodesInKernelk = netSize - kernelRows[k];
	if (k < kernelRows.size()-1)
	  nodesInKernelk = kernelRows[k+1]-kernelRows[k];
	if(k<=partition[j][1]){
	  inCount += nodesInKerneli* nodesInKernelk;
	  boxAverage += similarityMatrix[row+col*netSize]
	    *nodesInKerneli*nodesInKernelk;
	  boxAverageSq += 
	    similarityMatrix[row+col*netSize]*similarityMatrix[row+col*netSize] 
	    *nodesInKerneli *nodesInKernelk;
// 	  cout<<"kernel "<<k<<" "<<nodesInKernelk<<" "<<inCount<<endl;
	}
	else{
	  outCount += nodesInKerneli* nodesInKernelk;
	  outAverage += similarityMatrix[row+col*netSize] 
	    *nodesInKerneli * nodesInKernelk;
	  outAverageSq +=
	    similarityMatrix[row+col*netSize]*similarityMatrix[row+col*netSize]
	    *nodesInKerneli * nodesInKernelk;
// 	  cout<<similarityMatrix[row+col*netSize]<< " "
// 	      <<nodesInKerneli * nodesInKernelk<<endl;
	}
      }
    }

    tempMean.push_back(boxAverage);  
//     cout<<"box av "<<boxAverage<<" "<<double(inCount)<<endl;
    tempCount.push_back(inCount);  

 
    totalCount += inCount;
    AverageSqTotal += boxAverageSq;
    if(inCount>0)
      AverageTotal +=  boxAverage* boxAverage/double(inCount);

//     cout<<av2tt<<" "<<avtt<<" "<<tcount<<endl;
    
    
//     cout<<fluct<<endl;
  }

  tempMean.push_back( outAverage);  
  tempCount.push_back(outCount);  

  if(outCount >0){
    AverageSqTotal += outAverageSq ;
    AverageTotal += outAverage*outAverage/double(outCount);
  }

  *mean = tempMean;
  *nCount = tempCount;
   
  leastSquares = AverageSqTotal - AverageTotal;

  totalCount += outCount;
//   cout<<"out Av "<<outAverage<<" "<<double(outCount)<<endl;
//   cout<<"total "<<totalCount<<" "<<AverageSqTotal<<" "<<AverageTotal<<endl;
  

//  cout<<leastSquares<<" "<<n<<endl;
//    return leastSquares;
//   cout<<"leastSquares "<<leastSquares<<"---  "<<av2tt<<" "<<avtt<<endl;

//    cout<<leastSquares<<" "<<tcount<< " "<<(n*n-n)/2<<endl;
  return leastSquares/double(totalCount);
//    return chi/double(n*(n+1)/2);

}

/********************************************************************************/
/*Generates initial partition of nbox boxes**************************************/
/********************************************************************************/
//OK
void InitializePartition(int nbox, int netSize, vector<vector<int> > *partition){
  //netSize is the  total number of places where a breaking point can be placed +1

  vector<int> tempPart;//store temporaty partition
  int i, node;
  partition->clear();
  node=0;
  tempPart.push_back(node); //initial box starts at node 0

  
  for(i=0;i<nbox-1;i++){

    node+=int((netSize-1-(nbox-i-1)-node+1)*ran1f()); //following node

    tempPart.push_back(node);
    (*partition).push_back(tempPart);
    tempPart.clear();
    node++;
    tempPart.push_back(node);//initial node of second box
  }

  tempPart.push_back(netSize-1); //last box ends at n-1
  (*partition).push_back(tempPart);


}    
  
/*Finds optimal partition using greedy algorithm The algorithm will
try to find iteratively the best partition with two boxes by
evaluation of leastSq(2) for each breaking point bPoint, and selecting
the bPoint for which leastSq is minimal.  Then, the algorithm will add
another breaking point, if it finds a partition for which
leastSq(3)<leastSq(2). It will iterate the process until it cannot
find a better partition.  The function also stores the history of the
best partition found for each number of boxes
*/
//OK
void FindOptimalPartition(vector<vector<int> > *partition,vector<int>*label,
			  vector<double> simMatrix, int netSize,
			  double *leastSqOld,vector < vector< vector <int> > > *partitionHistory,
			  vector< double > *leastSqHistory, vector <int> kernelRows){

  int i,count=0;//count will be set to one when all the partition did not change in the iteration before
  int nbox = 0;
//   cout<<"1 box. LeastSquares = "<<*leastSqOld<<endl;

  while (count==0){
    
    SplitPartition(partition, label, simMatrix , netSize, 
		   leastSqOld, kernelRows);
    count=1;

    for(i=0;i<partition->size();i++)
      if((*label)[i]==1){
	count=0;
 	break;
     }
    nbox ++;
//     cout<< nbox +1<< " boxes. LeastSquares = "<<*leastSqOld<<endl;
    
    leastSqHistory->push_back(*leastSqOld);
    partitionHistory->push_back(*partition);
  }

}

/*Finds optimal partitions using greedy algorithm. The algorithm will
try to find the best partitions for nboxini +1, nboxini+ nbox. nboxini
is the number of boxes in the partition given as initial partition for
the search. The starting point for the best partition of nboxes, i
sthe best partition found for nboxes-1. Note that it is not required
that leastSq(noboxes)<leastSq(nboxes-1). The function also stores the
history of the best partition found for each number of boxes
*/

void FindOptimalPartition(vector<vector<int> > *partition,vector<int>*label,
			  vector<double> simMatrix,int netSize, double *leastSqOld,
			  vector < vector< vector <int> > > *partitionHist, 
			  vector< double > *leastSqHist, int nBox, vector<int> kernelRows)
{

  int i;
  
  int min = kernelRows.size()-1;
  if(nBox < min) min = nBox;
//   cout<<min<<endl;

  while (partition->size()<min){

    SplitPartitionAlways( partition, label, simMatrix, netSize, 
			  leastSqOld, kernelRows);

    leastSqHist->push_back(*leastSqOld);
    partitionHist->push_back(*partition);

  }


}

/*********************************************************************/
/***Splitting box functions*******************************************/
/*********************************************************************/
/*Finds the best split.  The labels vector keeps track of which is the
  box that was split (1). If all the labels are set to zero, then no
  boxes were split. No boxes are split when lestSq cannot be minimized
  by the split of any individual box into two.
*/
void SplitPartition(vector<vector<int> > *partition, vector<int>*label,
		    vector<double> simMatrix, int netSize, double *leastSqOld,
		    vector<int> kernelRows){

  int i,j;
  vector<vector<int> >partitionNew;//stores initial row of partition and final node of partition
  vector<int>labelNew;


  vector<int> result, dummy;///two elements first 0 no solit 1 yes split|| second -1 or the node boundary

  vector< vector <int> > splitResult;
  double leastSq = *leastSqOld, leastSqMin;
  int nSplit = -1;//breaking point

  for(i=0;i<partition->size();i++){
    
    result = GetBestSplit(*partition, i, simMatrix, netSize,
			  *leastSqOld, &leastSqMin, kernelRows);

    splitResult.push_back(result);

//      cout<<i<<" "<<chimin<<endl;
    if(leastSqMin < leastSq - 1e-8){
      nSplit = i;
      leastSq = leastSqMin;
    }

//      cout<<i<<" "<<partition->size()<<" "<<result[1]<<" "
//  	<<result[0]<<" "<<chimin<<endl;

  }

  for(i=0;i<partition->size();i++){
    
    if(i == nSplit){
      //     if(splitResult[i][0]==1){
      
      dummy.clear();
      dummy.push_back((*partition)[i][0]);
      dummy.push_back(splitResult[i][1]);
      labelNew.push_back(1);
      partitionNew.push_back(dummy);
      
      dummy.clear();
      dummy.push_back(splitResult[i][1]+1);
      dummy.push_back((*partition)[i][1]);
      partitionNew.push_back(dummy);
      labelNew.push_back(1);
      
    }
    else{
      partitionNew.push_back((*partition)[i]);
      labelNew.push_back(0);
      
    }
    
  }
  
  *partition = partitionNew;
  *label = labelNew;
  
  if (nSplit!=-1){
    double leastSquares = ComputeLeastSquares(simMatrix, partitionNew,
					      netSize, kernelRows);
    //   cout<<"chi "<<chi<<endl;
    
    *leastSqOld = leastSq;  
  }
  
}

/*Finds the best split for a certain number of boxes. The labels
  vector keeps track of which is the box that was split (1).
*/void SplitPartitionAlways
(vector<vector<int> > *partition, vector<int>*label,
 vector<double> simMatrix, int netSize, double *leastSqOld,
 vector<int> kernelRows){

  int i,j;
  vector<vector<int> >partitionNew;
  vector<int>labelNew;


  vector<int> result,dummy;//two elements:
                          // first element:  0 no split 1 yes split
                          // second element: -1 or the node boundary

  vector< vector <int> > splitResult;
  double leastSq=10000., leastSqMin;
  int nSplit=-1;//keeps track of the box that is split into two.

  for(i=0;i<partition->size();i++){

    result = GetAnyBestSplit( *partition, i, simMatrix, netSize,
			      *leastSqOld, &leastSqMin, kernelRows);

    splitResult.push_back(result);

//     cout<<i<<" "<<chimin<<endl;
    if((*partition)[i][1]!= (*partition)[i][0])
      if( leastSqMin < leastSq - 1e-8){
	nSplit = i;
	leastSq = leastSqMin;
      }

//     cout<<i<<" "<<partition->size()<<" "<<result[1]<<" "
// 	<<result[0]<<" "<<chimin<<endl;

  }

  for(i = 0; i<partition->size(); i++){
    
    if(i == nSplit){
      //     if(splitResult[i][0]==1){
      
      dummy.clear();
      dummy.push_back((*partition)[i][0]);
      dummy.push_back(splitResult[i][1]);
      labelNew.push_back(1);
      partitionNew.push_back(dummy);
      
      dummy.clear();
      dummy.push_back(splitResult[i][1]+1);
      dummy.push_back((*partition)[i][1]);
      partitionNew.push_back(dummy);
      labelNew.push_back(1);
      
    }
    else{
      partitionNew.push_back((*partition)[i]);
      labelNew.push_back(0);
      
    }
    
  }
  
 
  
  *partition=partitionNew;
  *label=labelNew;
  
  if (nSplit != -1){
     double leastSquares = ComputeLeastSquares(simMatrix, partitionNew,
					       netSize, kernelRows);
   
    
     *leastSqOld = leastSq;
  }
  
}
  

/*********************************************************************/
/*********************************************************************/

/*********************************************************************/
/**functions that find best split*************************************/
/*********************************************************************/
/*Function that finds the best split inside a box and checks wether
  the least Squares in minimized or not. If it is not, there is no
  split returned (that is, all the labels are set to zero. Breaking
  points are considered to be between kernels.*/
vector<int >  GetBestSplit(vector<vector<int> > partition, int nPart,
			   vector<double> simMatrix, int netSize, 
			   double leastSqOld, double *leastSqMin,
			   vector<int> kernelRows){


  vector<int> kernelChange, kernelStay1, kernelStay2; //nodes that
						      //change
						      //position
                                                   //nodes that do not
                                                   //change (top/bottom)
  vector<vector<int> > partitionTemp, partitionOld;
  //we split box npart into two, putting the first kernel on its own
  //box)
  
  partitionOld = partitionTemp = 
    AddPartition(partition, nPart, partition[nPart][0]);
  vector<double> meanChange;
  vector<double> mean;
  vector<int> count;

  double leastSq = ComputeLeastSquares
    (simMatrix, partitionTemp, &mean, &count, netSize, kernelRows);
//   cout<<"leastSq: "<<leastSq<<endl;

  //mean is avector that contains the means inside each one of the
  //boxes count is a vector that contains the elements inside each one
  //of the boxes and of the 'outside box'

  int nodeMin = partition[nPart][0];
  int j;
  int nBox1;
  int nBox2;
  int width;
  double leastSqTemp, deltaLeastSq=0;
  //initialization of leastSquares value
  *leastSqMin = leastSq;   
  
   
  //loop that moves breaking point from j to to end of partition -1
  for(j=partition[nPart][0]+1;j<partition[nPart][1];j++){
      
    kernelStay1.clear();
    kernelStay2.clear();
    kernelChange.clear();

    kernelStay1.push_back(j+1);
    kernelStay1.push_back(partition[nPart][1]);
    kernelStay2.push_back(partition[nPart][0]);
    kernelStay2.push_back(j-1);
    kernelChange.push_back(j);
    kernelChange.push_back(j);
    
    ComputeChanges(simMatrix, kernelChange, kernelStay1, kernelStay2, 
		   netSize, &meanChange, kernelRows);

    //count lateral size of boxes before moving breaking point
    nBox1 = GetCount(partition[nPart][0], j-1, kernelRows, netSize);
    nBox2 = GetCount(j, partition[nPart][1], kernelRows, netSize);
    
    width = kernelRows[j+1] - kernelRows[j];//number of nodes in
					    //kernel that changes form
					    //one box to the other
//     cout<<"nnodes 1 "<<nBox1<<endl;
//     cout<<"nnodes 2 "<<nBox2<<endl;
//     cout<<"width "<<width<<endl;

//     cout<<nBox1<<" "<<nBox2<<endl;
    deltaLeastSq = ComputeDeltaLeastSq(netSize, width, nPart+1, -1, meanChange,
			      mean, count, nBox2, nBox1 );
    
    UpdatePartition( width, nPart+1, -1, &partitionOld, &mean,
		     meanChange, &count, kernelRows, netSize);
//     leastSqTemp = ComputeLeastSquares(simMatrix, partitionOld, netSize, kernelRows);

    leastSq += deltaLeastSq;
//     cout<< "Comparison "<<leastSq <<" "<< leastSqTemp<<endl;
    

    if(leastSq < *leastSqMin - 1e-8){
      *leastSqMin = leastSq;
      nodeMin=j;
//       cout<<chimin<<" "<<nodemin<<endl;
    }
    
    
  }
  
  vector<int> result;


  if (*leastSqMin < leastSqOld - 1e-8){

    result.push_back(1);
    result.push_back(nodeMin);
    result.push_back(nodeMin);
    
  }
  else{

    result.push_back(0);
    result.push_back(-1);
    *leastSqMin=10000;
  }
 
  return result;


}

/*Function that finds the best split inside a box, regaaardless of
  whether the least Squares in minimized with respect the partition
  with one box less or not. Breaking points are considered to be
  between kernels.*/
vector<int >  GetAnyBestSplit(vector<vector<int> > partition,int nPart,
			      vector<double> simMatrix, int netSize, 
			      double leastSqOld, double *leastSqMin,
			      vector<int> kernelRows){


  vector<int> kernelChange, kernelStay1, kernelStay2; //nodes that
						      //change
						      //position
                                                   //nodes that do not
                                                   //change (top/bottom)
  vector<vector<int> > partitionTemp, partitionOld;
  //we split box npart into two, putting the first kernel on its own
  //box)
  
  partitionOld = partitionTemp = 
    AddPartition(partition, nPart, partition[nPart][0]);
  vector<double> meanChange;
  vector<double> mean;
  vector<int> count;

  double leastSq = ComputeLeastSquares
    (simMatrix, partitionTemp, &mean, &count, netSize, kernelRows);
  //mean is avector that contains the means inside each one of the
  //boxes count is a vector that contains the elements inside each one
  //of the boxes and of the 'outside box'

  int nodeMin = partition[nPart][0];
  int j;
  int nBox1;
  int nBox2;
  int width;
  double leastSqTemp, deltaLeastSq=0;
  //initialization of leastSquares value
  *leastSqMin = leastSq;   


 for(j=partition[nPart][0]+1;j<partition[nPart][1];j++){
      
    kernelStay1.clear();
    kernelStay2.clear();
    kernelChange.clear();

    kernelStay1.push_back(j+1);
    kernelStay1.push_back(partition[nPart][1]);
    kernelStay2.push_back(partition[nPart][0]);
    kernelStay2.push_back(j-1);
    kernelChange.push_back(j);
    kernelChange.push_back(j);
    
    ComputeChanges(simMatrix, kernelChange, kernelStay1, kernelStay2, 
		   netSize, &meanChange, kernelRows);

    //count lateral size of boxes before moving breaking point
    nBox1 = GetCount(partition[nPart][0], j-1, kernelRows, netSize);
    nBox2 = GetCount(j, partition[nPart][1], kernelRows, netSize);
    width = kernelRows[j+1] - kernelRows[j];//number of nodes in
					    //kernel that changes form
					    //one box to the other

    deltaLeastSq = ComputeDeltaLeastSq(netSize, width, nPart+1, -1, meanChange,
			      mean, count, nBox2, nBox1 );
    
    UpdatePartition( width, nPart+1, -1, &partitionOld, &mean,
		     meanChange, &count, kernelRows, netSize);
//     leastSqT=ComputeLeastSquares(simMatrix, partitionOld, netSize, kernelRows);

    leastSq += deltaLeastSq;
//      cout<< chi <<" "<<ComputeLeastSquares(sim,partitionold,n)<<endl;
    

    if(leastSq < *leastSqMin - 1e-8){
      *leastSqMin = leastSq;
      nodeMin=j;
//       cout<<chimin<<" "<<nodemin<<endl;
    }
    
    
  }
  
  vector<int> result;
  
  result.push_back(1);
  result.push_back(nodeMin);
  result.push_back(nodeMin);

  return result;


}


/*Adds a new breaking point to the box partition[nPart], putting the
   kernel 'kernel' into a new box*/
vector<vector<int> > AddPartition
(vector<vector<int> > partition, int nPart, int kernel){

  vector<int> dummy;
  int j;
  vector<vector<int> > partitionnew;

  for(j=0;j<partition.size();j++){
    if(j==nPart){
      dummy.clear();
      dummy.push_back(partition[j][0]);
      dummy.push_back(kernel);
      partitionnew.push_back(dummy);
      dummy.clear();
      dummy.push_back(kernel+1);
      dummy.push_back(partition[j][1]);
      partitionnew.push_back(dummy);
    }
    else{
      partitionnew.push_back(partition[j]);
    }
  }
  return partitionnew;

}

/*Function that computes the changes in mean, when changing the
  breaking point*/

void ComputeChanges
(vector<double> simMat,vector<int> kernelChange,
 vector<int> kernelStay1, vector<int> kernelStay2,
 int netSize, vector<double> *mean, vector<int>kernelRows){
  
  double average=0., deltaLeastSq=0., avv=0.;
  int i,j;
  int col, row;
  int nodesInKerneli, nodesInKernelj;

  mean->clear();

//     cout<<"ch "<<nodechange[0]<<" "<<nodechange[1]<<endl;
//     cout<<"ch "<<nodestay1[0]<<" "<<nodestay1[1]<<endl;
//     cout<<"ch "<<nodestay2[0]<<" "<<nodestay2[1]<<endl;

  /*For each kernel we have to compute the width of some kernel*/
//   if(kernelChange[1] <kernelRows.size()-1)
//     cout << "change "<<kernelRows[kernelChange[0]] << " "<< kernelRows[kernelChange[1]+1]-1<<endl;
//   else
//     cout << "change "<<kernelRows[kernelChange[0]] << " "<< netSize -1<<endl;
//   if(kernelStay1[1] < kernelRows.size()-1)
//     cout << "stay1 "<<kernelRows[kernelStay1[0]] << " "<< kernelRows[kernelStay1[1]+1] - 1<<endl;
//   else
//     cout << "stay1 "<<kernelRows[kernelStay1[0]] << " "<< netSize -1<<endl;
//   if(kernelStay2[1] < kernelRows.size()-1)
//     cout << "stay2 "<<kernelRows[kernelStay2[0]] << " "<< kernelRows[kernelStay2[1]+1] - 1<<endl;
//   else
//     cout << "stay2 "<<kernelRows[kernelStay2[0]] << " "<< netSize -1<< endl;

  for(i=kernelChange[0];i<=kernelChange[1];i++){
    row = kernelRows[i];
    nodesInKerneli = netSize - row;
    if (i < kernelRows.size()-1)
      nodesInKerneli = kernelRows[i+1] - row;
    average += simMat[row + netSize*row] 
      * nodesInKerneli* (nodesInKerneli + 1)/2 ;
    avv += simMat[row + netSize*row]
      * nodesInKerneli* (nodesInKerneli +1)/2 ;
    for(j=i+1;j<=kernelChange[1];j++){
      col = kernelRows[j];
      nodesInKernelj = netSize - kernelRows[j];
      if (j < kernelRows.size()-1)
	nodesInKernelj = kernelRows[j+1]-kernelRows[j];
	
      average += simMat[row + netSize*col] 
	* nodesInKerneli* nodesInKernelj ;
      avv += simMat[row + netSize*col]
	* nodesInKerneli* nodesInKernelj ;
    }

    for(j=kernelStay1[0];j<=kernelStay1[1];j++){
      col = kernelRows[j];
      nodesInKernelj = netSize - kernelRows[j];
      if (j < kernelRows.size()-1)
	nodesInKernelj = kernelRows[j+1]-kernelRows[j];
      
      average += simMat[row + netSize*col] 
	* nodesInKerneli* nodesInKernelj;
//        cout<<"ch: "<<i<<" "<<j<<" "<<av<<" "<<avv<<" "<<count1<<
//   	" "<<simMat[i+j*n]<<" "<<simMat[j+i*n]<<endl;
    }
    

    for(j=kernelStay2[0];j<=kernelStay2[1];j++){
      col = kernelRows[j];
      nodesInKernelj = netSize - kernelRows[j];
      if (j < kernelRows.size()-1)
	nodesInKernelj = kernelRows[j+1]-kernelRows[j];
      
      avv += simMat[row + netSize*col] 
	* nodesInKerneli* nodesInKernelj;

    }
  }
//   cout <<"changes "<< -average<<" "<<avv<<endl;
  mean->push_back(-average);
  mean->push_back(avv);
  mean->push_back(average-avv);

 

}


/*Function that returns the lateral size of a box with 
kernels form kerenlIni to kernelEnd*/

int GetCount
(int kernelIni, int kernelEnd, vector<int> kernelRows, int netSize){

  int nodeCount = 0;
  int kernel;
  int kernelSize;
  
  for (kernel = kernelIni; kernel<=kernelEnd;kernel++){
    kernelSize = netSize - kernelRows[kernel];
    if(kernel < kernelRows.size() -1)
      kernelSize = kernelRows[kernel +1 ] - kernelRows[kernel];
    nodeCount += kernelSize;
      
      }
  return nodeCount;

}
  
/*Computes the change in the least Squares introduced by changing one
  kernel form one box to the neoiughboring one*/
double ComputeDeltaLeastSq
(int netSize, int width, int box, int neig, vector<double> meanChange, 
 vector<double> mean, vector<int> count,int boxSize1,int boxSize2 ){


  double deltaEn = 0.;
  //compute the sizes of the elements after moving one kernel form one box to another
  int nbox1p = count[box] - width*(2*boxSize1 - width + 1)/2;
  int nbox2p = count[box+neig] + width*(2*boxSize2 + width + 1)/2;
  int nboxop = count[count.size()-1] - (width*(boxSize2 - boxSize1 + width));
//   cout<< nbox1p << " " <<nbox2p<< " "<<nboxop<<endl;
//   cout<< count[box] << " " <<count[box+neig]<< " "<<count[count.size()-1]<<endl;
  
  double deltaMeanBox1 = 0.;
  if (nbox1p >0 )
    deltaMeanBox1 += (mean[box]+meanChange[0])*(mean[box]+meanChange[0])/double(nbox1p);
//   cout <<"deltaMeanBox1 "<<deltaMeanBox1 <<endl;
  if (count[box]>0 )
    deltaMeanBox1 += -mean[box]*mean[box]/double(count[box]);
//   cout <<"deltaMeanBox1 "<<deltaMeanBox1 <<endl;

  double deltaMeanBox2= 0.;

  if( nbox2p >0)
    deltaMeanBox2 += (mean[box+neig]+meanChange[1])*(mean[box+neig]+meanChange[1])
      /double(nbox2p);
//   cout <<"deltaMeanBox2 "<<deltaMeanBox2 <<endl;
  if (count[box+neig]>0)
    deltaMeanBox2 += -mean[box+neig]*mean[box+neig]/double(count[box+neig]);
//   cout <<"deltaMeanBox2 "<<deltaMeanBox2 <<endl;


  double deltaMeanOut= 0.;
  if(nboxop >0 )
    deltaMeanOut += (mean[mean.size()-1]+meanChange[2])*(mean[mean.size()-1]+meanChange[2])
      /double(nboxop);
//   cout <<"deltaMeanOut "<<deltaMeanOut <<endl;
  if(count[count.size()-1] > 0)
    deltaMeanOut +=  -mean[mean.size()-1]*mean[mean.size()-1]/double(count[count.size()-1]);
//   cout <<"deltaMeanOut "<<deltaMeanOut <<endl;

  double deltaMean1=-deltaMeanBox1;
  double deltaMean2=-deltaMeanBox2;
  double deltaOut=-deltaMeanOut;
//   cout <<deltaMean1<< " "<<deltaMean2<<" "<<deltaOut<<endl;
  deltaEn = deltaMean1 + deltaMean2 + deltaOut;

//   cout<< deltaEn<<" "<<(double)((netSize*netSize+netSize)/2)<<endl;

  return deltaEn/(double)((netSize*netSize+netSize)/2);


} 
  
/*Updates partition: meaning it keeps track of how the average and
  element count inside each box change when the kernel at the boundary
  is moved from one box to the neighboring one*/

void    UpdatePartition
(int width, int box, int neig, vector<vector<int> > *partition,
 vector<double> *mean, vector<double> meanChange, 
 vector<int> *count, vector<int> kernelRows, int netSize){

  
  int nBox1 = GetCount((*partition)[box][0],(*partition)[box][1], kernelRows, netSize);
  int nBox2 = GetCount((*partition)[box + neig][0], (*partition)[box+neig ][1], 
		       kernelRows, netSize);

  //compute new number of elements inside the boxes
  int nbox1p = (*count)[box] - width*(2*nBox1 - width + 1)/2;
  int nbox2p = (*count)[box+neig] + width*(2*nBox2 + width + 1)/2;
  int nboxop = (*count)[count->size()-1] - (width*(nBox2 - nBox1 + width));
  
  int kernelWidth= 1; //number of kernels moved

  switch(neig){

  case 1:
    
    (*partition)[box][1]=(*partition)[box][1]-kernelWidth;
    (*partition)[box+neig][0]=(*partition)[box][1]+1;
    break;

  case -1:

    (*partition)[box+neig][1]=(*partition)[box][0]+kernelWidth-1;
    (*partition)[box][0]=(*partition)[box][0]+kernelWidth;
    break;

  }
 
  (*mean)[box]+=meanChange[0];
  (*mean)[box+neig]+=meanChange[1];
  (*mean)[mean->size()-1]+=meanChange[2];

  (*count)[box]=nbox1p;
  (*count)[box+neig]=nbox2p;
  (*count)[count->size()-1]=nboxop;


}
