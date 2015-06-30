#ifndef _ordering_tools_hpp_included_
#define _ordering_tools_hpp_included_

// c++ header files
#include <vector>

// my header files

// namespaces
using namespace std;

// constants
// #define seed  243777//19371

namespace ORDERINGTOOLS {

  void StandardHelp();

  void OptionParser(int argc,
		    char *argv[],
		    char* &filename,
		    double &temperature,
		    double &factor,
		    unsigned int &randomization,
		    unsigned int &format,
		    unsigned long &randomSeed,
		    unsigned int &restart);

  void ExhaustiveHelp();

  void OptionParserExhaustive(int argc,
			      char *argv[],
			      char* &filename,
			      char* &translationFilename,
			      unsigned int &format,
			      unsigned int &nblocks,
			      unsigned int &restart);

}

double ranf();

double ran1f();

double gauss(double variance);

vector <int> RandomizeOrder(int n, int nblocks);
bool fileExists(const char* fileName);

#endif
