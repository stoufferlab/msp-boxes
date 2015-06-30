#ifndef _ordering_tools_cpp_included_
#define _ordering_tools_cpp_included_

// c++ header files
#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iostream>
#include <vector>

// my header files
#include "ordering-tools.hpp"

// namespaces
using namespace std;

// constants
double PI=4*atan(1.);

/////////////////////////////////////////////////////////

void ORDERINGTOOLS::StandardHelp(){
  cout << "Usage: executable [OPTIONS] -f SOURCE\n";
  cout << "     : executable -R\n\n";

  cout << "Mandatory arguments to long options are mandatory for short options too.\n";
  cout << "  -d, --data-format            the format of the data\n";
  cout << "                                 0: the data is stored i1 j1 x\n";
  cout << "                                 1: the data is stored x11 x12 ... x1N\n";
  cout << "                                                       xN1 xN2 ... xNN\n";
  cout << "                                 (default: 1)\n";
  cout << "  -f, --data-file              the data file either in matrix form or three\n";
  cout << "                                 column form (see -d option for details)\n";
  cout << "                                 WARNING: One of -f or -R must be used\n";
  cout << "  -r, --randomized-start       whether or not to randomize the network\n";
  cout << "                                 initially\n";
  cout << "                                 0: no\n";
  cout << "                                 1: yes\n";
  cout << "                                 (default: 0)\n";
  cout << "  -R, --restart                whether to restart from earlier simulation\n";
  cout << "                                 WARNING: All other options will be ignored\n";
  cout << "                                 WARNING: One of -f or -R must be used\n";
  cout << "  -s, --random-seed            the random seed (default: seeded by current\n";
  cout << "                                 time)\n";
  cout << "  -t, --temperature-factor     the factor to multiply the temperature by\n";
  cout << "                                 (default: 0.95)\n";
  cout << "  -T, --initial-temperature    the initial temperature in simulated annealing\n";
  cout << "                                 (default: calculated such that 95% of moves\n";
  cout << "                                 initially accepted)\n";
  cout << endl;
  cout << "      --help                   display this help and exit\n";
  exit(0);
}

void ORDERINGTOOLS::OptionParser(int argc,
				 char *argv[],
			 	 char* &filename,
				 double &temperature,
				 double &factor,
				 unsigned int &randomization,
				 unsigned int &format,
				 unsigned long &randomSeed,
				 unsigned int &restart){
  
  filename = " ";
  temperature = -1;
  factor = 0.95;
  randomization = 0;
  format = 1;
  randomSeed = 0;
  restart = 0;

  for(int i=1;i<argc;++i){
    if(strcmp(argv[i],"--help") == 0){
      StandardHelp();
    }else
      if(strcmp(argv[i],"-R") == 0 || strcmp(argv[i],"--restart") == 0){
	restart = 1;
      }else
	if(strcmp(argv[i],"-f") == 0 || strcmp(argv[i],"--data-file") == 0){
	  filename = argv[++i];
	}else
	  if(strcmp(argv[i],"-d") == 0 || strcmp(argv[i],"--data-format") == 0){
	    format = int(strtol(argv[++i],NULL,10));
	  }else
	    if(strcmp(argv[i],"-r") == 0 || strcmp(argv[i],"--randomized-start") == 0){
	      randomization = 1;
	    }else
	      if(strcmp(argv[i],"-s") == 0 || strcmp(argv[i],"--random-seed") == 0){
		randomSeed = strtol(argv[++i],NULL,10);
	      }else
		if(strcmp(argv[i],"-t") == 0 || strcmp(argv[i],"--temperature-factor") == 0){
		  factor = strtod(argv[++i],NULL);
		}else
		  if(strcmp(argv[i],"-T") == 0 || strcmp(argv[i],"--initial-temperature") == 0){
		    temperature = strtod(argv[++i],NULL);
		  }
  }

  if(restart){
    temperature = -1;
    randomization = 0;
  }else{
    if(!restart && filename == " "){
      cout << "You must specify an input filename (with -f) or to restart (with -R)\n";
      cout << "Try --help for more information.\n";
      exit(1);
    }
  }

  if(format != 0 && format != 1){
    cout << "You have submitted an invalid option for the file format\n";
    cout << "Try --help for more information.\n";
    exit(1);
  }

}

// Help for the exhaustive ordering code
void ORDERINGTOOLS::ExhaustiveHelp(){
  cout << "Usage: executable [OPTIONS] -f SOURCE\n";
  cout << "     : executable -R\n\n";

  cout << "Mandatory arguments to long options are mandatory for short options too.\n";
  cout << "  -d, --data-format            the format of the data\n";
  cout << "                                 0: the data is stored i1 j1 x\n";
  cout << "                                 1: the data is stored x11 x12 ... x1N\n";
  cout << "                                                       xN1 xN2 ... xNN\n";
  cout << "                                 (default: 1)\n";
  cout << "  -f, --data-file              the data file either in matrix form or three\n";
  cout << "                                 column form (see -d option for details)\n";
  cout << "                                 WARNING: One of -f or -R must be used\n";
  cout << "  -n, --block-size             the size of the blocks to swap for the\n";
  cout << "                                 exhaustive optimatization (default: 1)\n";
  cout << "  -R, --restart                whether to restart from earlier simulation\n";
  cout << "                                 WARNING: All other options will be ignored\n";
  cout << "                                 WARNING: One of -f or -R must be used\n";
  cout << "  -t, --translation-table      the translation table for the matrix order\n";
  cout << endl;
  cout << "      --help                   display this help and exit\n";
  exit(0);
}

void ORDERINGTOOLS::OptionParserExhaustive(int argc,
					   char *argv[],
					   char* &filename,
					   char* &translationFilename,
					   unsigned int &format,
					   unsigned int &nblocks,
					   unsigned int &restart){

  filename = " ";
  translationFilename = " ";
  format = 1;
  nblocks = 1;
  restart = 0;

  for(int i=1;i<argc;++i){
    if(strcmp(argv[i],"--help") == 0){
      ExhaustiveHelp();
    }else
      if(strcmp(argv[i],"-R") == 0 || strcmp(argv[i],"--restart") == 0){
	restart = 1;
      }else
	if(strcmp(argv[i],"-f") == 0 || strcmp(argv[i],"--data-file") == 0){
	  filename = argv[++i];
	}else
	  if(strcmp(argv[i],"-t") == 0 || strcmp(argv[i],"--translation-table") == 0){
	  translationFilename = argv[++i];
	  }else
	    if(strcmp(argv[i],"-d") == 0 || strcmp(argv[i],"--data-format") == 0){
	      format = int(strtol(argv[++i],NULL,10));
	    }else
	      if(strcmp(argv[i],"-n") == 0 || strcmp(argv[i],"--block-size") == 0){
		nblocks = int(strtol(argv[++i],NULL,10));
	      }
  }

  if(restart){
    if(filename != " " || translationFilename != " "){
      if(filename != " ")
	cout << "You cannot specify restart (-R) and an input filename (-f)\n";
      if(translationFilename != " ")
	cout << "You cannot specify restart (-R) and an input translation filename (-t)\n";
      cout << "Try --help for more information.\n";
      exit(1);
  }else{
      if(!restart && filename == " "){
	cout << "You must specify an input filename (with -f) or to restart (with -R)\n";
	cout << "Try --help for more information.\n";
	exit(1);
      }
    }
  }

  if(format != 0 && format != 1){
    cout << "you have submitted an invalid option for the file format\n";
    cout << "Try --help for more information.\n";
    exit(1);
  }

}


//uniform number (0,1]
double ranf(){

  return (double)(rand())/(RAND_MAX);

}

//uniform number [0,1)
double ran1f(){

  return (double)(rand())/(RAND_MAX+1.);

}

//generate a gaussian random variable with mean 0 and specified variance
double gauss(double variance){
  

  double x1=ranf();
  double x2=ranf();

  return variance*sqrt( - 2.*log(x1) )* cos( 2.* PI* x2 );

}

// return a randomized ordering of kernels
vector <int> RandomizeOrder(int n, int nblocks){

  vector<int> order, ord(n);
  int i, j, k=0;

  cout<<"Randomizing order\n";

  while( order.size() < n ){

    i = (int)(n*ran1f());
    i -= (i%nblocks);

    while(find(order.begin(),order.end(),i)!=order.end()){
      i = (int)(n*ran1f());
      i -= (i%nblocks);
    }    

    for( j=0; j<nblocks; j++){
      order.push_back( i + j );
      ord[ i + j ] = k;
      k ++;
//       cout<<"a "<<i+j<<" "<<order[i+j]<<endl;
    }
   
  }
  return ord;
}

bool fileExists(const char* fileName)
{
  fstream fin;
  fin.open(fileName);
  if( fin.is_open() )
  {
    fin.close();
    return true;
  }
  fin.close();
  return false;
}

#endif
