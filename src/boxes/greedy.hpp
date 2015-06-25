
double ranf();
double ran1f();
double gauss(double variance);
int ipow(int a, int b);
double GetAIC ( double L, int n, int kp);
double GetBIC ( double L, int n, int kp);
void ReadMatrix( int format, char *fileName,
		 vector <double> *sim, int nodeMin,	
		 int netSize, double *average, double *variance);
void UpdatePartition
(int box,int neig,int width, vector< vector <int> >*partition,
 vector<double> *mean,vector<double> meanch, vector<int> *count,
 vector<int> kernelRows, int netSize);
double ComputeDeltaLeastSq
(int netSize,int width,int nbox, int neig,vector<double> meanChange,
 vector<double> mean, vector<int> count,int boxSize1,int boxSize2);
void ComputeChanges
(vector<double> sim,vector<int> kernelChange,vector<int> kernelStay1,
 vector<int> kernelStay2,int netSize,vector<double> *mean, vector<int> kernelRows);
double ComputeLeastSquares(vector<double> sim, vector< vector <int > > partition,vector<double> *mean,
			   vector<int> *ncount,int n, vector<int> kernelRows);
double ComputeLeastSquares(vector<double> sim, vector< vector<int> >partition,int n, vector<int> kernelRows);
void InitializePartition(int n, int nbox, vector<vector<int> > *partition);
void FindOptimalPartition(vector<vector<int> > *partition,vector<int>*label,
			  vector<double> sim,int n,double *chiold,vector < vector< vector <int> > > *,
			  vector< double > *, vector<int> kernelRows);
void FindOptimalPartition(vector<vector<int> > *partition,vector<int>*label,vector<double> sim,int n,double *chiold
			  ,vector < vector< vector <int> > > *, vector< double > *, int, vector<int> kernelRows);
void SplitPartition(vector<vector<int> > *partition,vector<int>*label,vector<double> sim,int n,double *chiold, vector<int> kernelRows);
void SplitPartitionAlways(vector<vector<int> > *partition,vector<int>*label,vector<double> sim,int n,double *chiold,vector<int> kernelRows );
vector<int >  GetBestSplit(vector<vector<int> > partition,int npart,vector<double> sim,int n,double chiold,double *chimin, vector<int> kernelRows);
vector<int >  GetAnyBestSplit(vector<vector<int> > partition,int npart,vector<double> sim,int n,double chiold,double *chimin, vector<int> kernelRows);
vector<vector<int> > AddPartition(vector<vector<int> > partition,int npart,int nodepos);
vector<vector<int> > GetBestPartition(double *,vector < vector< vector <int> > >, vector< double > );
vector<vector<int> > GetBestPartitionAIC(double *,vector < vector< vector <int> > >, vector< double >, int );
vector<vector<int> > GetBestPartitionBIC(double *,vector < vector< vector <int> > >, vector< double >, int);
int GetCount(int kernelIni, int kernelEnd, vector<int> kernelRows, int netSize);
vector<int>  GetKernelRows (vector<double> similarityMatrix, int netSize);
