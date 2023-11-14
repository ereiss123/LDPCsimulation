/*
 *  LDPC_testbench.h
 *  LDPC
 *
 *  Created by Chris Winstead on 5/20/11.
 *  Copyright 2011 Utah State University. All rights reserved.
 *
 */

/*
 *  block_testbench.h
 *
 *  Created by Chris Winstead on 4/12/11.
 *  Copyright 2011 Utah State University. All rights reserved.
 *
 */


#ifndef BLOCK_TB_H
#define BLOCK_TB_H

#include <iostream>
#include <string>
#include <cmath>
#include <fstream>
#include <sstream>
using namespace std;

#include "systemc.h"
#include "sc_vector.h"

#include "ldpcsim.h"
#include "nodes.h"
#include "rand.h"


///////////////////////////////////////////////////////////////
// LDPC Testbench Module
///////////////////////////////////////////////////////////////
class LDPC_testbench : sc_module
{
public:
  //===========================================
  // Ports:
  //===========================================

  sc_vector<sc_out<double > >	y;     // Inputs to decoder
  sc_vector<sc_in<bool> >	d;     // Outputs from decoder
  sc_in<bool>           	finished;
  sc_in<bool>                   ready;
  sc_in<bool>			clk;
  sc_out<bool>			rst;


  //===========================================
  // Constructor:
  //===========================================
  SC_HAS_PROCESS(LDPC_testbench);
 LDPC_testbench(sc_module_name name) : sc_module(name), y("y",p.N), d("d",p.N), finished("stop"), ready("ready"),
    c(p.N), x(p.N), rx(p.N), codewordFile(p.stimfilename,ios::in)
    {
      init_object();
      SC_METHOD(update);
      sensitive << clk.pos();
    }
	
  //===========================================
  // Member Variables:
  //===========================================
 private:	
  bool			initialize;
  bool 			reset;
  ifstream		codewordFile;
  vector<int>		c;   // Codeword bits
  vector<int>		x;   // Antipodal modulated symbols in {+1,-1}
  vector<double>	rx;  // Received signal (modulated plus AWGN)	

  int errors;
  int word_errors;
  int totalbits;
  int total_words;
  int totalIterations;

  int numclocks;
  int numiterations;

	
	
  //=============================================
  // Update Method: Executed every clock cycle
  //=============================================
  void update()
  {
    numclocks++;

    int i;
    if (initialize)
      {			
	initialize = false;
	reset = true;
	numiterations = 0;
			
	// Read a line from the data file:
	string s;
	getline(codewordFile, s);

	if ((codewordFile.eof()) || (numclocks > p.total_clock_cycles))
	  {
	    numclocks = 0;
	    cout << "\n\nIncremental result at SNR=" << p.SNR << ", errors=" << errors << ", BER=";
	    cout << (double)errors/totalbits << " FER=" << (double)word_errors/total_words 
		 << ", iterations = " << (double) totalIterations/total_words << "\n\n";
	  }

	if (codewordFile.eof())
	  {
	    codewordFile.clear();
	    codewordFile.seekg(0);
	  }

	if ((word_errors > 30)&&(errors > 250))
	  {
	    cout << "At SNR=" << p.SNR << ", BER=";
	    cout << (double)errors/totalbits << " FER=" << (double)word_errors/total_words;
	    cout << ", iterations = " << (double) totalIterations/total_words << endl;
	    //----------------------------------------------
	    // DONE SIMULATING. WRITE OUT RESULTS AND STOP.
	    //----------------------------------------------
	    append_result_to_data_file(errors, totalbits, word_errors, total_words, totalIterations);
	    sc_stop();				
	  }

			
	for (i=0; i<p.N; i++)
	  {	    
	    if (s[i] == '1')
	      c[i] = 1;
	    else
	      c[i] = 0;
	  }

			
	// Transmit, add noise and receive:
	for (i=0; i<p.N; i++)
	  {
	    x[i] = (1 - 2.0*c[i]);
	    rx[i] = x[i] + p.sigma*rann();
	    y[i].write(rx[i]);
	  }					      		
      }
    // Burn one clock cycle to ensure the data is ready before asserting rst:
    else if (reset)
      {
	rst.write(true);
	if (ready.read())
	  reset = false;	
      }
    // Begin decoding:
    else
      {
	rst.write(false);				
	numiterations++;		

	// Is decoding finished?
	if (finished.read())
	  {
	    initialize = true;		       
	    int new_errors = 0;				
	    for(int i=0;i<p.N;i++)
	      {
		if (c[i] != d[i].read())
		  new_errors++;				
	      }
				
	    total_words++;
	    totalbits += p.N;
	    totalIterations += numiterations;				
	    if (new_errors > 0)
	      {
		word_errors++;
		errors += new_errors;
	      }				
	    rst.write(true);
	  }
			
      }
  }
	
	
	
  //===========================================
  // Helper Methods:
  //===========================================
	
	
  void init_object()
  {
    int i;
    numclocks = 0;
    numiterations = 0;    

    initialize = true;
    reset = false;
		
    errors = 0;
    word_errors = 0;
    totalbits = 0;
    total_words = 0;
    totalIterations = 0;
		
    // SET PARAMETERS:
    cout << "\nInitializing simulation for SNR=" << p.SNR 
	 << ", N=" << p.N << ", M=" << p.M << ", R=" << p.Rate 
	 << "\n";
  }

};


#endif
