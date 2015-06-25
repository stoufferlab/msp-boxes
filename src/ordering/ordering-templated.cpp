#ifndef _ordering_templated_cpp_included_
#define _ordering_templated_cpp_included_

// c++ header files

// my header files

// namespaces
using namespace std;

// swap the values of two variables
template <class T>
void Swap(T& a,T& b)
{
  T q=b;
  b=a;
  a=q;
}

// search for an object in a vector
template <class T>
int Find(T n, vector<T> nlist){

  int i;

  for(i=0;i<nlist.size();i++)
    if(nlist[i]==n)return i;

  return -1;
}

#endif
