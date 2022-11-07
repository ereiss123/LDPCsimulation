//==============================================================
// decodeBP.cpp
// 
// This is a reference implementation of the belief propagation
// LDPC decoding algorithm for benchmarking comparisons against
// new decoding algorithms and hardware implementations.
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
int    num_iterations; // Maximum number of iterations 
double MAXLLR;         // Maximum magnitude of LLR messages

//============ DECODING ALGORITHM PREDEFINES ===============//
void checkNodeUpdates(alist_struct &H, vector<vector<double> > & sym_to_check, vector<vector<double> > & check_to_sym);
void symNodeUpdates(alist_struct &H, vector<double> & y, vector<int> & d, vector<vector<double> > & sym_to_check, vector<vector<double> > & check_to_sym);


//============= SUPPORTING FUNCTION PREDEFINES =================//
void setupSymMessages(alist_struct & H, vector<vector<double> > & sym_to_check);
void setupCheckMessages(alist_struct & H, vector<vector<double> > & check_to_sym);
void initializeSymMessages(alist_struct & H,  vector<vector<double> > & sym_to_check, vector<double> & y);
double sgn(double x);
int find(int symNodes[], int len, int snode);
int countDecisionErrors(vector<int> d, vector<int> c);

//============= I/O PREDEFINES ============================//
void printHistogram(vector<int> & h);
void printVector(vector<int> v);
void printVector(vector<double> v);
void printErroneousMessages(vector<vector<int> > v);
void writeErroneousMessagesToFile(alist_struct & H, vector<vector<int> > & check_to_sym, vector<vector<int> > & sym_to_check, vector<int> & c, vector<int> & d, int fid);
void writeErroneousMessagesToFile(alist_struct & H, vector<vector<int> > & check_to_sym, vector<vector<int> > & sym_to_check, vector<int> & c, vector<int> & d, vector<double> & y, vector<int> & yq, vector<int> & r, int fid, int it);
void sendMQTT(int bitErrors, int frameErrors, int totalBits, int totalWords);



///////////////////////////////////////////////////////////////////////////////
// -----===== MAIN BODY ======------
///////////////////////////////////////////////////////////////////////////////
int main(int argc, char * argv[])
{
  MAXLLR = 20;

  vector<string> command_arguments(0);
  command_arguments.push_back("alist");
  command_arguments.push_back("R");
  //  command_arguments.push_back("pchan");
  command_arguments.push_back("SNR");
  command_arguments.push_back("T");
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
  //double pchan = atof(argv[idx++]);
  //cout << " pchan = \t" << pchan << endl;
  double SNR = atof(argv[idx++]);
  cout << " SNR = \t" << SNR << endl;
  num_iterations = atoi(argv[idx++]);
  cout << " T = \t" << num_iterations << endl;
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
  double N0 = pow(10.0,-SNR/10.0)/R;  // power spectral density of noise
  double sigma = sqrt(N0/2.0);  // Standard deviation of channel noise

  // Get code parameters:
  int dv = H.biggest_num_n;
  int dc = H.biggest_num_m;

  // Report initial status messages:
  cout << "Simulating Min-Sum decoding on code with N=" << H.N << ", M=" << H.M << ", R=" << R << ", dv=" << dv << ", dc=" << dc << endl;
  //cout << "\nParameters are:\n\tpchan\t" << pchan << endl; 
  cout << "\nParameters are:\n\tSNR\t" << SNR << endl; 

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

  // Allocate appropriate sizes for the message memories:                                                                                                                                                         
  setupSymMessages(H,sym_to_check);
  setupCheckMessages(H,check_to_sym);

  /////////////////////////////////////////////////////////////////
  // ------===== MAIN TEST LOOP =====-------
  /////////////////////////////////////////////////////////////////
  int minWordErrors = 20;
  if (H.N > 10000) minWordErrors = 10;
  if (H.N > 50000) minWordErrors = 5;
  ran_seed(time(0)); //(134159);
  int i,j;
  while ((errors < 200) || (wordErrors < minWordErrors))
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
      // Emulate Additive White Gaussian Noise (AWGN) transmission      
      for (i=0; i<H.N; i++)
	{
	  /* BSC:
	  y[i] = x[i]; //*(1.0+sigma*rann());
	  double rnum = ranu();
	  if (rnum < pchan) {
	    y[i] = 1.0 - y[i];
	  }
	  */
	  y[i] = x[i]*(1.0+sigma*rann());
	  
	  //yq[i] = log(pchan)/log(1.0-pchan); // y[i]; //
	  
	  yq[i] = 4.0*y[i]/N0;
	  
	  if (abs(yq[i]) > MAXLLR) 
	    yq[i] = sgn(yq[i])*MAXLLR;
	  
	  r[i] = sgn(yq[i]);
	  d[i] = r[i];
	  if (r[i]*c[i] < 0)
	    uncodedErrors++;
	}

      initializeSymMessages(H, sym_to_check, yq);


      // Perform decoding iterations:      
      int it;

      
      for (it=0; it<num_iterations; it++)
	{      
	  // First update the check nodes:
	  checkNodeUpdates(H,sym_to_check,check_to_sym);
	  
	  // Then perform Symbol node updates:
	  symNodeUpdates(H, yq, d, sym_to_check, check_to_sym);	  
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
	  
	  //	  sendMQTT(errors,wordErrors,totalBits+H.N,totalWords+1);
	}

      // Increment frame and bit counters:
      totalWords++;
      totalBits += H.N;
      totalIterations += it;
      
      // ------------------------------------------------
      // Give a status message every 5 frames
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
     << num_iterations << tab;

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

void  setupSymMessages(alist_struct & H, vector<vector<double> > & sym_to_check)
{
  for (int i=0; i<H.N; i++)
    {
      vector<double> initial_msgs(H.num_nlist[i],0.0);
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


void initializeSymMessages(alist_struct & H,  vector<vector<double> > & sym_to_check, vector<double> & y)
{
  int i,j;
  for (i=0; i<H.N; i++)
    for (j=0; j<H.num_nlist[i]; j++)
      sym_to_check[i][j] = y[i];
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
  double outmsg;
  for (int i=0; i<H.M; i++)
    {
      for (int j=0; j<H.num_mlist[i]; j++)
	{
	  prod=1.0;
	  for (int k=0; k<H.num_mlist[i]; k++)
	    {
	      if (j != k)
		{
		  int snode = H.mlist[i][k]-1;
		  int midx = find(H.nlist[snode], H.num_nlist[snode], i);
		  msg = sym_to_check[snode][midx];
		  prod *= tanh(msg/2.0);
		}
	    }
	  outmsg = log((1.0+prod)/(1.0-prod));	      
	  check_to_sym[i][j] = outmsg;	  
	}
    }
}

void symNodeUpdates(alist_struct &H, vector<double> & y, vector<int> & d, vector<vector<double> > & sym_to_check, vector<vector<double> > & check_to_sym)
{ 
  int cnode, midx;
  double msg, outmsg, sum;

  for (int i=0; i<H.N; i++)
    {
      sum = y[i];
      for (int j=0; j<H.num_nlist[i]; j++)
	{
	  cnode = H.nlist[i][j]-1;
	  midx = find(H.mlist[cnode],H.num_mlist[cnode],i);
	  msg = check_to_sym[cnode][midx];
	  sum += msg;
	}      
      for (int j=0; j<H.num_nlist[i]; j++) 
	{
	  cnode = H.nlist[i][j]-1;  
	  midx = find(H.mlist[cnode],H.num_mlist[cnode],i); 
	  msg = check_to_sym[cnode][midx]; 
	  outmsg = sum - msg;
	  if (abs(outmsg) > MAXLLR)
	    outmsg = MAXLLR*sgn(outmsg);
	  sym_to_check[i][j] = outmsg;
	}
      if (sum > 0)
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

/*
void sendMQTT(int bitErrors, int wordErrors, int totalBits, int totalWords) {
  stringstream ss;
  ss << "./send_message.sh \'" 
     << bitErrors << '\t' << wordErrors << '\t'
     << totalBits << '\t' << totalWords << '\t'
     << (double)bitErrors/totalBits << '\t'
     << (double)wordErrors/totalWords << '\'';
  cout << ss.str() << endl;

  system(ss.str().c_str());
}
*/


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
