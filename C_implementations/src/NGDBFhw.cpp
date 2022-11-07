//==============================================================
// GDBF.cpp
// 
// This program performs NGDBF decoding on AWGN channel.
// It performs global threshold adaptation in order to
// regulate flip activity. 
//
// Repeated decoding is also simulated. The received frame
// is decoded with maxPhase repetitions. Each repetition is
// restarted from the same initial condition. 
// 
// When accounting for iteration statistics, the phase with the
// *least* iterations is recorded. This is to simulate an 
// architecture with parallel decoders instead of sequential 
// phases. Once any decoder has completed, all decoders are 
// stopped, so only the minimum number of iterations matters.
// 
// One concern is how errors are counted if all decoders fail.
// What's implemented here is the minimum number of bit errors
// are counted. That may be too optimistic since there is no
// way to know which of the parallel decoders was closest. This
// should be replaced with a systematic procedure for choosing
// the output frame; otherwise only the FER result should be
// relied upon.
//==============================================================


//--- Standard C++ headers ---//
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <vector>
#include <cmath>
#include <sstream>
#include <time.h>
#include <bitset>

using namespace std;

//--- Borrowed from Radford Neal's source code ---//
#include "alist.h"
#include "rand.h"


//============ GLOBAL PARAMETERS ============//
// Example values are shown. These are set
// via the command line.
int    num_iterations = 600;   // Maximum number of iterations per phase
double R              = 0.8413;// Code rate
double w              = 0.185;  // Syndrome weight parameter 
double Ymax           = 1.625;   // Maximum channel sample magnitude
double noiseScale     = 0.95;  // Proportionality between channel noise and perturbation noise
int    maxPhases      = 1;     // Number of decoding repetitions
long   numFrames      = 10000; // Number of frames to simulate
long   seed           = 1234;  // Random number generator seed
const int NQ          = 5;     // Number of bits for quantization
double theta0         = -0.525;

alist_struct H;                // Code definition
string logfilename;            // Filename for output data


//============ GLOBAL VARIABLES ===============//
double SNR            = 3.5;   // Eb/N0 in decibels
double theta          = 8;     // Threshold 
double numFlips       = 0;     // Number of flips in most recent iteration
int    Smult          = 10;    // Syndrome multiplier to account for quantization

//============ DECODING ALGORITHM PREDEFINES ===============//
void checkNodeUpdates(vector<int> & d, vector<int> & syndrome, bool & satisfied);
void symNodeUpdates(vector<double> & yprime, vector<int> & d, vector<int> & syndrome, vector<int> & E, vector<double> & qprime, int qpointer, vector<int> & flip);

//============= SUPPORTING FUNCTION PREDEFINES =================//
void quantize(vector<double> & y, vector<double> & yq);
void quantizebig(vector<double> & y, vector<double> & yq);
int quantize(double y);
unsigned long pack(int sample, int sign);
int unpack(unsigned long sample);
unsigned long packbig(int sample, int sign);
int unpackbig(unsigned long sample);
int find(int symNodes[], int len, int snode);
int countDecisionErrors(vector<int> d, vector<int> c);
double sgn(double y);
vector<string> setupUsage();
void   parseArguments(int argc, char * argv[]);

//============= I/O PREDEFINES ============================//
void printHistogram(vector<int> & h);
void printVector(vector<int> v);
void printVector(vector<double> v);



///////////////////////////////////////////////////////////////////////////////
// -----===== MAIN BODY ======------
///////////////////////////////////////////////////////////////////////////////
int main(int argc, char * argv[])
{

  //=========== Handle command line arguments ============//
  vector<string> command_arguments = setupUsage();

  if ((argc != command_arguments.size()) && (argc != command_arguments.size()+1))
    {
      cout << "Usage: " << argv[0];
      for (int i=0; i<command_arguments.size(); i++)
	cout << " " << command_arguments[i];
      cout << "\n";
      return 0;
    }

  parseArguments(argc, argv);

  ran_seed(seed); 

  ifstream codewordFile;
  if (argc == command_arguments.size()+1)
    {
      cout << "\nUsing codewords from " << argv[argc] << endl;
      codewordFile.open(argv[argc],ios::in);
    }
  else
    cout << "\nUsing all-zero sequence.\n";

  //===========   Initialize Simulation   ============//
  // Compute channel parameters:
  double N0 = pow(10.0,-SNR/10.0)/R;
  double sigma = sqrt(N0/2.0);
  double noiseSigma = sigma*noiseScale;

  // Get code parameters:
  int dv = H.biggest_num_n;
  int dc = H.biggest_num_m;

  // Report initial status messages:
  cout << "Simulating GDBF decoding on code with N=" << H.N << ", M=" << H.M << ", R=" << R << ", dv=" << dv << ", dc=" << dc << endl;
  cout << "\nParameters are:\n\tSNR\t" << SNR << "\n\tN0\t" << N0 << "\n\tsigma\t" << sigma << endl;


  // Declare top-level variables:
  vector<int>    c(H.N,0);      // Bipolar codeword (all +1 in this simulation)
  vector<double> x(H.N,1);      // Modulated codeword (all +1 in this simulation)
  vector<double> y(x);          // Channel samples
  vector<double> ymodified(x);  // Modified channel samples
  vector<double>    yprime(H.N,0); // Modified and quantized channel samples
  vector<int>    r(H.N);        // Received bipolar decisions (+1 or -1)
  vector<int>    d(H.N,0);      // Decoder outputs (0 or 1 after decoding)
  vector<int>    E(H.N,0);      // Flip function
  vector<int>    flip(H.N,0);      // Flip activity

  vector<double> qmodified(2648,0.0);
  vector<double> qprime(2648,0.0);
  int qpointer=0;

  // Declare and initialize statistics variables:
  long errors = 0;            // Total bit errors
  long uncodedErrors = 0;     // Bit errors in r before decoding
  long totalBits = 0;         // Total number of bits observed
  long totalWords = 0;        // Total number of frames observed
  long wordErrors = 0;        // Number of word errors observed
  long totalIterations = 0;   // Total number of iterations accumulated over all frames.

  vector<int> error_weight_hist(H.N,0);       // Vector to serve as histogram of error-pattern weights (1 up to H.N)
  vector<double> itdist(num_iterations,0.0);  // Cumulative distribution of completion times

  // Declare and initialize message memories:
  vector<int> syndrome(H.M,0);

  /////////////////////////////////////////////////////////////////
  // ------===== MAIN TEST LOOP =====-------
  /////////////////////////////////////////////////////////////////
  int minWordErrors = 1;
  int i,j;
  double qmax=pow(2,(NQ));
  double lmax=Ymax/(2.0*w);
  double NL=qmax-1;

  theta = unpack(pack(quantize(2),1));//*(lmax/NL);
  Smult = round(NL/lmax);
  #ifdef LOG_PROCESSING
      stringstream ss1,ss2,ss3;
      ss1 << logfilename << "_" << SNR << "_msgs.dat";
      ss2 << logfilename << "_" << SNR << "_chanin.dat";
      ss3 << logfilename << "_" << SNR << "_noise.dat";
      ofstream ofmsgs(ss1.str().c_str(),ios::trunc);
      ofstream ofchanin(ss2.str().c_str(),ios::trunc);
      ofstream ofnoise(ss3.str().c_str(),ios::trunc);
      std::bitset<NQ+1> bth(theta);
      ofmsgs << "GLOBALS:\n\ttheta = " << theta << "(" << bth << ")" << endl;
      ofmsgs << "\tSmult = " << Smult << endl;
  #endif

  while (totalWords < numFrames)
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
	      c[i] = 1;
	    else if (s[i] == '0')
	      c[i] = 0;
	    else
	      cout << "Got an invalid symbol at index " << i << endl;
	    x[i] = 1-2*c[i];
	  }
	}
      // Emulate AWGN or BSC transmission      
      for (i=0; i<H.N; i++)
	{
	  y[i] = x[i]*(1.0+sigma*rann());

	  if (abs(y[i])>Ymax)
	    y[i] *= Ymax/abs(y[i]);
	  if (y[i] > 0)
	    r[i] = 1;
	  else
	    {
	      r[i] = -1;	      
	    }
	  if (r[i]*c[i] < 0)
	    uncodedErrors++;
	  d[i] = (1-r[i])/2;
	  ymodified[i] = y[i]/(2.0*w);
	  //yprime[i] = ymodified[i];
	}

      quantize(ymodified, yprime);

      for (i=0; i<qprime.size(); i++)
	{
	  double q = noiseSigma*rann();
	  //if (abs(q)>Ymax)
	  //  q = q*Ymax/abs(q);
	  qmodified[i] = ((q-theta0)/(2.0*w) - 1.0);
	  if (qmodified[i]>lmax)
	    qmodified[i] = lmax;
	  else if (qmodified[i] < -lmax)
	    qmodified[i] = -lmax;

	  //qprime[i]=round(128.0*qmodified[i])/128.0;
	}
      quantize(qmodified, qprime);

      bool satisfied;
      int it;


      //-------------- Out Multi-Phase Loop -----------------//
      int leastIterations=num_iterations;
      int leastErrors=H.N;
      #ifdef LOG_PROCESSING
      if (totalWords==0) {
      for (int idx=0; idx<H.N; idx++) {
	unsigned long yul = yprime[idx];
	std::bitset<NQ> by(yul);
	ofchanin << by << endl;
	unsigned long qul = qprime[idx];
	std::bitset<NQ> bn(qul);
	ofnoise << bn << endl;
      }
      for (int idx=H.N; idx<qprime.size(); idx++)
	{
	unsigned long qul = qprime[idx];
	std::bitset<NQ> bn(qul);
	ofnoise << bn << endl;
	}
    }
      #endif

      for (int phase=0; phase<maxPhases; phase++)
	{
	  //	  theta = theta0;
	  for (int idx=0; idx<H.N; idx++)
	    {
	      d[idx] = (1-r[idx])/2;
	    }
     
     
	  //------------ Inner Loop: NGDBF Decoder --------------//
	  for (it=0; it<num_iterations; it++)
	    {      
	      satisfied = true;
	      numFlips = 0;	  
	  
	      //'''''''''''''''''''''''''''''''''''''''''''''''
	      // First update the check nodes:
	      checkNodeUpdates(d,syndrome,satisfied);
	      if (satisfied)
		break;
	  
	      // Then perform Symbol node updates:
	      symNodeUpdates(yprime, d, syndrome, E, qprime,qpointer,flip);

	      #ifdef LOG_PROCESSING
	      if (totalWords==0) {
	      ofmsgs << "IT " << it << endl;
	      for (int idx=0; idx<H.N; idx++)
		{
		  ofmsgs << "S" << idx << ":\n";
		  unsigned long yul = yprime[idx];
		  std::bitset<NQ> by(abs(yul));
		  ofmsgs << "\tchan_msg, x: " << y[idx] << " " << ymodified[idx] << " " <<  yul << " (" << by << ") [" << unpack(yprime[idx]) << "], " << d[idx] << endl;
		  
		  ofmsgs << "\tin_messages: ";
		  int SSum = 0;
		  for (int jdx=0; jdx<H.num_nlist[idx]; jdx++) {
		    int msg = syndrome[H.nlist[idx][jdx]-1];
		    ofmsgs << msg  << " ";
		    SSum += 1-msg;
		  }
		  unsigned long Sul = SSum*Smult;
		  std::bitset<NQ+1> bS(Sul);
		  ofmsgs << "\n\tS: " << SSum << " " << " (" << Sul << "," << bS << ")";
		  unsigned long uq = round(qprime[idx+qpointer]);
		  std::bitset<NQ+1> b(uq);
		  //if (uq>16)
		  //  uq = -(uq-16);
		  ofmsgs << "\n\tq: " << qmodified[idx+qpointer] << " " << qprime[idx+qpointer] << " (" <<  b.to_string<char,std::string::traits_type,std::string::allocator_type>() << ")";
		  ofmsgs << " [" << unpack(qprime[idx+qpointer]) << "]";
		  ofmsgs << "\n\tE: " << E[idx] << endl;
		  ofmsgs << "\ttheta: " << theta << endl;
		  ofmsgs << "\tflip: " << flip[idx] << endl;
		}
	    }
	      #endif


	      //	      symNodeUpdates(y, d, syndrome, noiseSigma); 
	      //,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,
	  
	      // ..............................................
	      // Do threshold adaptation using throttle method:
	      /*
	      double df = f0-numFlips;
	      if (df > f0)
		df = f0;
	      if (df < -f0)
		df = -f0;
	      theta = theta + thetaAdj*df;

	      if (theta > thetaMax)
		theta = thetaMax;
	      */
	      // ..............................................
	      
	      qpointer++;
	      if (qpointer >= (qprime.size()-H.N))
		qpointer=0;
	}
      
      // --- End of iteration --------------------------------------
      // -------------------------------------------------------------


      // Count remaining errors after decoding:
      int newErrors = countDecisionErrors(d,c);
      
      // Count least iterations from repeated decoding:
      if (newErrors < leastErrors)
	leastErrors = newErrors;
      if (it < leastIterations)
	leastIterations = it;
     }
	 
	  

      //==================  ACCOUNTING  ==================//
      if (leastErrors > 0)
	{
	  // Report the frame error to the console:
	  cout << "Ferr with " << leastErrors << " errors.";
	  if (satisfied)
	    cout << " All checks satisfied.\n";
	  else
	    cout << endl;

	  // Update statistical information
	  errors += leastErrors;
	  error_weight_hist[leastErrors-1]++;
	  wordErrors++;
	  
	  cout << " BER=" << (double)errors/(totalBits+H.N) << ", WER=" << (double) wordErrors/(totalWords+1) << endl;
	  // ------------------------------------------------
	  // WRITE ERROR PATTNERS TO FILE
	  // ------------------------------------------------
	  #ifdef writeErrorPatterns
	  stringstream ss1, ss2;
	  ss1 << logfilename << "_" << SNR << "_errpat.dat";
	  ss2 << logfilename << "_" << SNR << "_dec.dat";
	  ofstream oferrpat(ss1.str().c_str(),ios::app);
	  ofstream ofdec(ss2.str().c_str(),ios::app);
	  for (int idx=0; idx<H.N; idx++)
	    {
	      oferrpat << y[idx] << "\t";
	      ofdec << d[idx] << "\t";
	    }
	  oferrpat << endl;
	  ofdec << endl;
	  oferrpat.close();
	  ofdec.close();
	  #endif
	}

      // Increment frame and bit counters:
      totalWords++;
      totalBits += H.N;
      totalIterations += leastIterations;

      // Update cumulative distribution of completion times:
      for (int idx=0; idx<=leastIterations; idx++)
	itdist[idx] = (double)((totalWords-1.0)/totalWords)*itdist[idx] + (double)(1.0/totalWords); 

      // ------------------------------------------------
      // Give a status message every 100 frames
      // ------------------------------------------------
      int reportInterval = 100;
      if ((totalWords % reportInterval) == 0)
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
  // APPEND FINAL RESULTS TO LOG FILE:
  // ------------------------------------------------
  cout << "\nFinal result: " << errors << " bit errs in " 
       << totalWords << " words, BER=" << (double)errors/totalBits << ". Average iterations = " << (double) totalIterations/totalWords 
       << ". Uncoded errors = " << uncodedErrors << ", uncBER=" 
       << (double)uncodedErrors/totalBits << endl;      

  ofstream of(logfilename.c_str(),ios::app);
  char tab = '\t';
  of << SNR << tab << errors << tab << wordErrors << tab << (double)errors/totalBits << tab << (double) totalIterations/totalWords << tab
     << (double) wordErrors/totalWords << tab
     << totalBits << tab << totalWords << tab
     << num_iterations << tab << theta0 << tab;
  of << noiseScale << tab; 
  of << w << tab;
  of << Ymax << tab << NQ << tab;
  of << maxPhases << tab << seed;
  of << endl;

  // ------------------------------------------------
  // WRITE COMPLETION TIME DISTRIBUTION TO FILE
  // ------------------------------------------------
  stringstream ss;
  ss << logfilename << "_" << SNR << "_itdist.dat";
  ofstream ofitdist(ss.str().c_str(),ios::trunc);
  for (int idx=0; idx<itdist.size(); idx++)
    ofitdist << idx << "\t" << itdist[idx] << "\n";
  ofitdist.close();

  return 0;
}
/////////////////////////////////////////////////////////////////
// ------===== END OF MAIN BODY =====-------
/////////////////////////////////////////////////////////////////


vector<string> setupUsage()
{
  vector<string> command_arguments(0);

  command_arguments.push_back("alist");
  command_arguments.push_back("SNR");
  command_arguments.push_back("numFrames");
  command_arguments.push_back("seed");
  command_arguments.push_back("logfilename");

  command_arguments.push_back("[codeword filename]");

  return command_arguments;
}


void parseArguments(int argc, char * argv[])
{
  // Parse command arguments:
  int idx=1;
  H = loadFile(argv[idx++]);
  cout << "PARAMETERS: \n alist = \t" << argv[1] << endl;
  SNR = atof(argv[idx++]);
  cout << " SNR = \t" << SNR << endl;

  numFrames = atoi(argv[idx++]);
  cout << "Simulating for " << numFrames << " frames." << endl;

  seed = atoi(argv[idx++]);
  cout << "Using random seed " << seed << endl;

  logfilename.assign(argv[idx++]);
  cout << " log = \t" << logfilename << endl;
}




//============================================================//
// Functions follow in no particular order, and without
// adequate comments...
//============================================================//

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
      if (isnan(d[i]))
	cout << "Problem decision at index " << i << " = " << d[i] << endl;
      if (d[i] != c[i])	
	  errs++;
    }
  return errs;
}




void checkNodeUpdates(vector<int> & d, vector<int> & syndrome, bool & satisfied)
{
  int msg;
  satisfied = true;
  for (int i=0; i<H.M; i++)
    {
      int prod = 1;
      for (int j=0; j<H.num_mlist[i]; j++)
	{
	  int snode = H.mlist[i][j]-1;
	  msg = d[snode];
	  prod *= (1-2*msg);
	}
      if (prod < 0)
	satisfied = false;
      syndrome[i] = (1-prod)/2;	
    }
}

void symNodeUpdates(vector<double> & yprime, vector<int> & d, vector<int> & syndrome, vector<int> & E, vector<double> & qprime, int qpointer, vector<int> & flip)
{
  double qmax=pow(2,(NQ));
  double lmax=Ymax/(2.0*w);
  double NL=qmax-1;
  
  for (int i=0; i<H.N; i++)
    {
      E[i] = (1-2*d[i])*unpack(yprime[i]);//*(lmax/NL);

      int dv = H.num_nlist[i];
      double SSum=0;
      for (int j=0; j<H.num_nlist[i]; j++)
	{
	  int cnode = H.nlist[i][j]-1;
	  int msg = syndrome[cnode];
	  SSum += 1-msg;	  
	}      
      E[i] += SSum*Smult+unpack(qprime[i+qpointer]);//*(lmax/NL);
      if (E[i] <= theta)
	{
	  flip[i] = 1;
	  d[i] = 1-d[i];      	    
          numFlips++;
	}
      else
	flip[i] = 0;
    }
}


/*
void quantize(vector<double> & y, vector<double> & yq)
{
  int i;
  double qmax=pow(2,(NQ-1));
  double lmax=Ymax/(2.0*w);
  for (i=0; i<y.size(); i++)
    {
      yq[i] = sgn(y[i])*floor((abs(y[i])*qmax)/(2*lmax)+0.5)*(2.0*lmax/qmax);
    }
  
}
*/

 // Alternative quantization:
void quantize(vector<double> & y, vector<double> & yq)
{
  int i;
  double qmax=pow(2,(NQ));
  double lmax=Ymax/(2.0*w);
  double NL=qmax-1;
  for (i=0; i<y.size(); i++)
    {
      //yq[i] = sgn(y[i])*(lmax/NL)*(1.0+2.0*floor((abs(y[i])*NL)/(2*lmax)));
      yq[i] = pack(quantize(y[i]),sgn(y[i])); //sgn(y[i])*round(floor(abs(y[i])*NL/(2*lmax))));
    }
  
}

void quantizebig(vector<double> & y, vector<double> & yq)
{
  int i;
  double qmax=pow(2,(NQ));
  double lmax=Ymax/(2.0*w);
  double NL=qmax-1;
  for (i=0; i<y.size(); i++)
    {
      //yq[i] = sgn(y[i])*(lmax/NL)*(1.0+2.0*floor((abs(y[i])*NL)/(2*lmax)));
      yq[i] = packbig(quantize(y[i]),sgn(y[i])); //sgn(y[i])*round(floor(abs(y[i])*NL/(2*lmax))));
    }
  
}


int quantize(double y)
{
  int i;
  double qmax=pow(2,(NQ));
  double lmax=Ymax/(2.0*w);
  double NL=qmax-1;
  int yq;
  yq = sgn(y)*round(floor(abs(y)*NL/(2*lmax)));
    
  return yq;
}


unsigned long pack(int sample,int sign)
{
  bitset<NQ> b1(abs(sample));
  //b1 = b1 >> 1;
  if (sign<0)
    b1.set(NQ-1);
  unsigned long msg = b1.to_ulong();
  return msg;
}

int unpack(unsigned long sample)
{
  int j;
  std::bitset<NQ+1> b(sample);
  b = b << 1;
  b.set(0);
  if (b.test(NQ))
    {
      b.reset(NQ);
      j=-b.to_ulong();
    }
  else
    j=b.to_ulong();
  return j;
}

unsigned long packbig(int sample,int sign)
{
  bitset<NQ+1> b1(abs(sample));
  //b1 = b1 >> 1;
  if (sign<0)
    b1.set(NQ);
  unsigned long msg = b1.to_ulong();
  return msg;
}

int unpackbig(unsigned long sample)
{
  int j;
  std::bitset<NQ+2> b(sample);
  b = b << 1;
  b.set(0);
  if (b.test(NQ+1))
    {
      b.reset(NQ+1);
      j=-b.to_ulong();
    }
  else
    j=b.to_ulong();
  return j;
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




double sgn(double y)
{
  if (y>0)
    return 1.0;
  else 
    return -1.0;
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



