//==============================================================
// decodeDDBMP.cpp
// 
// This program simulates the DD-BMP agorithm on an AWGN
// channel.
//==============================================================


//--- Standard C++ headers ---//
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <vector>
#include <cmath>
#include <sstream>
#include <time.h>
using namespace std;

//--- Borrowed from Radford Neal's source code ---//
#include "alist.h"
#include "rand.h"


//============ GLOBAL PARAMETERS ============//
int    num_iterations; // Maximum number of iterations for MLSBM (an additional phase of Gallager-A follows after this)

//============ DECODING ALGORITHM PREDEFINES ===============//
void checkNodeUpdates(alist_struct &H, vector<vector<double> > & sym_to_check, vector<vector<double> > & check_to_sym);
void symNodeUpdates(alist_struct &H, vector<double> & y, vector<int> & d, vector<vector<double> > & sym_to_check, vector<vector<double> > & check_to_sym, vector<vector<double> > & sym_memories);
bool checkStoppingCondition(alist_struct &H, vector<int>  & d);


//============= SUPPORTING FUNCTION PREDEFINES =================//
void setupSymMessages(alist_struct & H, vector<vector<double> > & sym_to_check, vector<vector<double> > & sym_memories);
void setupCheckMessages(alist_struct & H, vector<vector<double> > & check_to_sym);
void initializeSymMessages(alist_struct & H,  vector<vector<double> > & sym_to_check, vector<vector<double> > & sym_memories, vector<double> & y);
double sgn(double x);
int find(int symNodes[], int len, int snode);
int countDecisionErrors(vector<int> d, vector<int> c);
double quantize(double x, double Ymax, double Nq);

//============= I/O PREDEFINES ============================//
void printHistogram(vector<int> & h);
void printVector(vector<int> v);
void printVector(vector<double> v);
void printErroneousMessages(vector<vector<int> > v);
void writeErroneousMessagesToFile(alist_struct & H, vector<vector<int> > & check_to_sym, vector<vector<int> > & sym_to_check, vector<int> & c, vector<int> & d, int fid);
void writeErroneousMessagesToFile(alist_struct & H, vector<vector<int> > & check_to_sym, vector<vector<int> > & sym_to_check, vector<int> & c, vector<int> & d, vector<double> & y, vector<int> & yq, vector<int> & r, int fid, int it);



///////////////////////////////////////////////////////////////////////////////
// -----===== MAIN BODY ======------
///////////////////////////////////////////////////////////////////////////////
int main(int argc, char * argv[])
{
  vector<string> command_arguments(0);
  command_arguments.push_back("alist");
  command_arguments.push_back("R");
  command_arguments.push_back("SNR");
  command_arguments.push_back("T");
  command_arguments.push_back("Ymax");
  command_arguments.push_back("Q");
  command_arguments.push_back("logfilename");
  command_arguments.push_back("[codeword filename]");

  // Check arguments and print usage statements:
  if ((argc != command_arguments.size()) && (argc != command_arguments.size()+1))
    {
      cout << "Usage: " << argv[0];
      for (int i=0; i<command_arguments.size(); i++)
	cout << " " << command_arguments[i];
      cout << "\n";
      return 0;
    }

  // Parse command arguments:
  int idx=1;
  alist_struct H = loadFile(argv[idx++]);
  cout << "PARAMETERS: \n alist = \t" << argv[1] << endl;
  double R = atof(argv[idx++]);
  cout << " R = \t" << R << endl;
  double SNR = atof(argv[idx++]);
  cout << " SNR = \t" << SNR << endl;
  num_iterations = atoi(argv[idx++]);
  cout << " T = \t" << num_iterations << endl;
  double Ymax = atof(argv[idx++]);
  cout << " Ymax = \t" << Ymax << endl;
  int Q = atoi(argv[idx++]);
  cout << "Q = \t" << Q << endl;
  string logfilename(argv[idx++]);
  cout << " log = \t" << logfilename << endl;

  ifstream codewordFile;
  if (argc == command_arguments.size()+1)
    {
      cout << "\nUsing codewords from " << argv[idx] << endl;
      codewordFile.open(argv[idx],ios::in);
    }
  else
    cout << "\nUsing all-zero sequence.\n";

  // Compute channel parameters:
  double N0 = pow(10.0,-SNR/10.0)/R;
  double sigma = sqrt(N0/2.0);
  double Nq = pow(2.0,Q);

  // Get code parameters:
  int dv = H.biggest_num_n;
  int dc = H.biggest_num_m;

  // Report initial status messages:
  cout << "Simulating DD-BMP decoding on code with N=" << H.N << ", M=" << H.M << ", R=" << R << ", dv=" << dv << ", dc=" << dc << endl;
  cout << "\nParameters are:\n\tSNR\t" << SNR << "\n\tN0\t" << N0 << "\n\tsigma\t" << sigma << endl;


  // Declare top-level variables:
  vector<int>    c(H.N,1);     // Bipolar codeword (all +1 in this simulation)
  vector<double> x(H.N,1);     // Modulated codeword (all +1 in this simulation)
  vector<double> y(x);          // Channel samples
  vector<double> yq(H.N);       // Quantized channel samples
  vector<int>    d(H.N,0);      // Decoder outputs (+1 or -1 after decoding)
  vector<int>    r(H.N,0);      // Received hard decision

  // Declare and initialize statistics variables:
  long errors = 0;            // Total bit errors
  long uncodedErrors = 0;     // Bit errors in r before decoding
  long totalBits = 0;         // Total number of bits observed
  long totalWords = 0;        // Total number of frames observed
  long wordErrors = 0;        // Number of word errors observed
  long totalIterations = 0;   // Total number of iterations accumulated over all frames.
  vector<int> error_weight_hist(H.N,0);  // Vector to serve as histogram of error-pattern weights (1 up to H.N)

  // Declare and initialize message memories:

  vector<vector<double> > check_to_sym;
  vector<vector<double> > sym_to_check;
  vector<vector<double> > sym_memories;

  // Allocate appropriate sizes for the message memories:                                                                                                                                                         
  setupSymMessages(H,sym_to_check,sym_memories);
  setupCheckMessages(H,check_to_sym);

  /////////////////////////////////////////////////////////////////
  // ------===== MAIN TEST LOOP =====-------
  /////////////////////////////////////////////////////////////////
  ran_seed(time(0)); //(134159);
  int i,j;
   while ((errors < 200) || (wordErrors < 40))
    {
      string s;
      // If a codeword file is specified, load codewords from the file:
      if (argc == command_arguments.size()+1)
	{
	  getline(codewordFile, s);
	  if (codewordFile.eof())
	  {
	    codewordFile.clear();
	    codewordFile.seekg(0);
	    getline(codewordFile, s);
	  }
	  for (i=0; i<H.N; i++)
	  {	    
	    if (s[i] == '1')	      
	      c[i] = -1;
	    else if (s[i] == '0')
	      c[i] = +1;
	    else
	      cout << "Got an invalid symbol at index " << i << endl;
	    x[i] = c[i];
	  }
	}
      // Emulate AWGN transmission      
      for (i=0; i<H.N; i++)
	{
	  y[i] = x[i]*(1.0+sigma*rann());
	  yq[i] = quantize(y[i],Ymax,Nq);
	  if (yq[i] > 0)
	    r[i] = 1;
	  else
	    r[i] = -1;
	  d[i] = r[i];
	  if (r[i]*c[i] < 0)
	    uncodedErrors++;
	}

      initializeSymMessages(H, sym_to_check, sym_memories, yq);


      // Perform decoding iterations:      
      int it;

      
      for (it=0; it<num_iterations; it++)
	{      
	  // First update the check nodes:
	  checkNodeUpdates(H,sym_to_check,check_to_sym);
	  
	  // Then perform Symbol node updates:
	  symNodeUpdates(H, yq, d, sym_to_check, check_to_sym, sym_memories);	  

	  // Check stopping condition:
	  if (checkStoppingCondition(H,d))
	    break;
	}
      
      // --- End of iteration --------------------------------------
      // -------------------------------------------------------------


      // Count remaining errors after decoding:
      int newErrors = countDecisionErrors(d,c);

      if (newErrors > 0)
	{
	  // Report the frame error to the console:
	  cout << "Ferr with " << newErrors << " errors.";
	  cout << endl;

	  // Update statistical information
	  errors += newErrors;
	  error_weight_hist[newErrors-1]++;
	  wordErrors++;
	  
	}

      // Increment frame and bit counters:
      totalWords++;
      totalBits += H.N;
      totalIterations += it;
      
      // ------------------------------------------------
      // Give a status message every 100 frames
      if ((totalWords % 5) == 0)
	{
	  cout << "\nIncremental result: " << errors << " bit errs in " << totalWords << " words, BER=" << (double)errors/totalBits 
	       << ". Average iterations = " << (double) totalIterations/totalWords << ". Word error=" << wordErrors << ". Uncoded errors = " << uncodedErrors << ", uncBER=" << (double)uncodedErrors/totalBits
	       << "\nError weights:\n";
	  printHistogram(error_weight_hist);
	}
      // ------------------------------------------------

    }
  /////////////////////////////////////////////////////////////////
  // ------===== END OF MAIN TEST LOOP =====-------
  /////////////////////////////////////////////////////////////////
  
  // ------------------------------------------------
  // REPORT FINAL RESULTS:
  cout << "\nFinal result: " << errors << " bit errs in " 
       << totalWords << " words, BER=" << (double)errors/totalBits << ". Average iterations = " << (double) totalIterations/totalWords 
       << ". Uncoded errors = " << uncodedErrors << ", uncBER=" 
       << (double)uncodedErrors/totalBits << endl;      

  ofstream of(logfilename.c_str(),ios::app);
  char tab = '\t';
  of << SNR << tab << (double)errors/totalBits << tab << (double) totalIterations/totalWords << tab
     << (double) wordErrors/totalWords << tab
     << num_iterations << tab
     << Ymax << tab 
     << Q << tab;

  of << argv[1]
     << endl;
  of.close();

  freeAlist(H);

  return 0;
}
/////////////////////////////////////////////////////////////////
// ------===== END OF MAIN BODY =====-------
/////////////////////////////////////////////////////////////////


//============================================================//
// Functions follow in no particular order, and without
// adequate comments...
//============================================================//

void  setupSymMessages(alist_struct & H, vector<vector<double> > & sym_to_check, vector<vector<double> > & sym_memories)
{
  for (int i=0; i<H.N; i++)
    {
      vector<double> initial_msgs(H.num_nlist[i],0.0);
      sym_memories.push_back(initial_msgs);
      sym_to_check.push_back(initial_msgs);
    }
}

void setupCheckMessages(alist_struct & H,vector<vector<double> > & check_to_sym)
{
  for (int i=0; i<H.M; i++)
    {
      vector<double> initial_msgs(H.num_mlist[i],0.0);
      check_to_sym.push_back(initial_msgs);
    }
}


void initializeSymMessages(alist_struct & H,  vector<vector<double> > & sym_to_check, vector<vector<double> > & sym_memories, vector<double> & y)
{
  int i,j;
  for (i=0; i<H.N; i++)
    for (j=0; j<H.num_nlist[i]; j++)
      {
	sym_to_check[i][j] = sgn(y[i]);
	sym_memories[i][j] = y[i];
      }
}


void printHistogram(vector<int> & h)
{
  for (int i=0; i<h.size(); i++)
    {
      if (h[i] > 0)
	cout << i+1 << ":\t" << h[i] << endl;
    }
}

int countDecisionErrors(vector<int> d, vector<int> c)
{
  int errs =0 ;
  for (int i=0; i<d.size(); i++)
    {
      if ((isnan(d[i])) || (d[i] == 0))
	cout << "Problem decision at index " << i << " = " << d[i] << endl;
      if (d[i] != c[i])	
	  errs++;
    }
  return errs;
}




void printErroneousMessages(vector<vector<int> > v)
{
  for (int i=0; i<v.size(); i++)
    {
      for (int j=0; j<v[i].size(); j++)
	{
	  if (v[i][j] < 0)
	    cout << " " << i << "." << j;
	}
    }
}

void checkNodeUpdates(alist_struct &H, vector<vector<double> > & sym_to_check, vector<vector<double> > & check_to_sym)
{
  double msg;
  double prod;
  for (int i=0; i<H.M; i++)
    {
      prod = 1.0;
      for (int j=0; j<H.num_mlist[i]; j++)
	{
	  int snode = H.mlist[i][j]-1;
	  int midx = find(H.nlist[snode], H.num_nlist[snode], i);
	  msg = sym_to_check[snode][midx];
	  prod *= sgn(msg);
	}
      for (int j=0; j<H.num_mlist[i]; j++)
	{
	  int snode = H.mlist[i][j]-1;
	  int midx = find(H.nlist[snode], H.num_nlist[snode], i);
	  msg = sym_to_check[snode][midx];
	  check_to_sym[i][j] = prod*sgn(msg);	
	}
    }
}


bool checkStoppingCondition(alist_struct &H, vector<int> & d)
{
  double msg;
  double prod=1.0;
  bool satisfied = true;
  for (int i=0; i<H.M; i++)
    {
      prod=1.0;
      for (int j=0; j<H.num_mlist[i]; j++)
	{
	  int snode = H.mlist[i][j]-1;
	  msg = sgn(d[snode]);
	  prod *= msg;
	}
      if (prod < 0)
	satisfied = false;
    }
    return satisfied;
}


void symNodeUpdates(alist_struct &H, vector<double> & y, vector<int> & d, vector<vector<double> > & sym_to_check, vector<vector<double> > & check_to_sym, vector<vector<double> > & sym_memories)
{ 
  for (int i=0; i<H.N; i++)
    {
      double sum = y[i];
      double dsum = sgn(y[i]);
      for (int j=0; j<H.num_nlist[i]; j++)
	{
	  int cnode = H.nlist[i][j]-1;
	  int midx = find(H.mlist[cnode],H.num_mlist[cnode],i);
	  double msg = check_to_sym[cnode][midx];
	  sum += msg;
	}      
      for (int j=0; j<H.num_nlist[i]; j++) 
	{
	  int cnode = H.nlist[i][j]-1;  
	  int midx = find(H.mlist[cnode],H.num_mlist[cnode],i); 
	  double msg = check_to_sym[cnode][midx]; 
	  sym_memories[i][j] += sum - msg;
	  sym_to_check[i][j] = sgn(sym_memories[i][j]);
	  dsum += sym_to_check[i][j];
	}
      if (dsum > 0)
	d[i] = 1;
      else
	d[i] = -1;	
    }
}


double sgn(double x)
{
  if (x >= 0.0)
    return 1.0;
  return -1.0;
}


double quantize(double x, double Ymax, double Nq)
{
  if (abs(x) > Ymax)
    return sgn(x)*Ymax;
  
  double qval =  sgn(x)*(floor(abs(x)*(Nq-1)/(2.0*Ymax)) + 0.0)*(2*Ymax/(Nq-1));
  if (qval == 0.0)
    qval = sgn(x)*2.0*Ymax/(Nq-1);
  return qval;
}


int find(int symNodes[], int len, int snode)
{
  int result = -1;
  for (int i=0; i<len; i++)
    {
      if (snode == (symNodes[i]-1))
	result = i;      
    }
  return result;
}

void printVector(vector<int> v)
{
  for (int i=0; i<v.size(); i++)
    {
      cout << v[i] << "\t";
    }
}                                                                                                                                                                                                   
void printVector(vector<double> v)
{
  for (int i=0; i<v.size(); i++)
    {
      cout << v[i] << "\t";
    }
}                                           



void writeErroneousMessagesToFile(alist_struct & H, vector<vector<int> > & check_to_sym, vector<vector<int> > & sym_to_check, vector<int> & c, vector<int> & d, int fid)
{
  #ifdef erroneousMessageFile
  stringstream ss;
  ss << erroneousMessageFile << "_" << fid << ".dat";
  ofstream f(ss.str().c_str(),ios::out);

  int i,j,k;
  for (i=0; i<H.N; i++)
    {
      int correct_msg = c[i];
      if (d[i] != correct_msg)
	f << "Dec err at node " << i << endl;
      bool anyErr = false;
      for (j=0; j<H.num_nlist[i]; j++)
	{
	  int msg = sym_to_check[i][j];
	  if (msg != correct_msg)
	    {	      
	      int cnode = H.nlist[i][j]-1;
	      int midx = find(H.mlist[cnode],H.num_mlist[cnode],i);
	      f << "\t StoC err, S " << i << " -> C " << cnode << endl;
	      for (k=0; k<H.num_mlist[cnode]; k++)
		{
		  if (k != midx)
		    {
		      int cmsg = check_to_sym[cnode][k];
		      int snode = H.mlist[cnode][k]-1;
		      int midx2 = find(H.nlist[snode], H.num_nlist[snode], cnode); 
		      int msg2 = sym_to_check[snode][midx2];
		      if (msg2 != c[snode])
			f << "\t\t Ck node " << cnode << " got err from Sym node " << snode << endl;
		    }
		}
	    }	  
	}      
    }
  f.close();
  #endif
}


void writeErroneousMessagesToFile(alist_struct & H, vector<vector<int> > & check_to_sym, vector<vector<int> > & sym_to_check, vector<int> & c, vector<int> & d, vector<double> & y, vector<int> & yq, vector<int> & r, int fid, int it)
{
  #ifdef erroneousMessageFile
  stringstream ss;
  ss << "detail_" << erroneousMessageFile << "_" << fid << ".dat";
  ofstream f;
  if (it == 0)
    f.open(ss.str().c_str(),ios::out);
  else
    f.open(ss.str().c_str(),ios::app);
  int i,j,k;

  f << "\n\nITERATION " << it << endl;
  for (i=0; i<H.N; i++)
    {
      int correct_msg = c[i];
      if (d[i] != correct_msg)
	f << "Dec err at node " << i << " y=" << y[i] << ", yq=" << yq[i] << ", r=" << r[i] << ", d=" << d[i] << ", c=" << c[i] << endl;
      bool anyErr = false;
      for (j=0; j<H.num_nlist[i]; j++)
	{
	  int msg = sym_to_check[i][j];
	  if (msg != correct_msg)
	    {	      
	      int cnode = H.nlist[i][j]-1;
	      int midx = find(H.mlist[cnode],H.num_mlist[cnode],i);
	      f << "\t StoC err, S " << i << " -> C " << cnode << endl;
	      for (k=0; k<H.num_mlist[cnode]; k++)
		{
		  if (k != midx)
		    {
		      int cmsg = check_to_sym[cnode][k];
		      int snode = H.mlist[cnode][k]-1;
		      int midx2 = find(H.nlist[snode], H.num_nlist[snode], cnode); 
		      int msg2 = sym_to_check[snode][midx2];
		      if (msg2 != c[snode])
			f << "\t\t Ck node " << cnode << " got err from Sym node " << snode << endl;
		    }
		}
	    }	  
	}      
    }
  f.close();
  #endif
}
