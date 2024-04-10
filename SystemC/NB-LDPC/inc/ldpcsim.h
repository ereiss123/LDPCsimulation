/*=========================================================================
** ldpcsim.h
**
** Date: March, 2024
**
** Authors: Eric Reiss, Chris Winstead
**          Utah State University
**
** Based on original source from
** Sept, 2011
** By Chris Winstead
      Department of Electrical and Computer Engineering
      Utah State University
      chris.winstead@usu.edu
** Description:
   Top header file for ldpcsim tool.
==========================================================================*/

#ifndef LDPCSIM_H
#define LDPCSIM_H

#include <iostream>
#include <string>
#include <cmath>
#include <vector>
#include <itpp/comm/galois.h>
using namespace std;

#include "alist.h"

// ===================================================================
// Type for internal messages between check nodes and symbol nodes:
// ===================================================================
typedef vector<double> message_type;


//================================
// Struct simparams for storing
// global parameters:
//================================
typedef struct
{
  // Stopping and reporting parameters:
  int iterations_per_frame;  // Maximum iterations per frame
  int total_clock_cycles;    // Maximum iterations between incremental progress reports

  // Structural parameters:
  alist_struct  alist;  // Alist specification for the code's parity-check matrix
  int           q;      // Number of bits per symbol
  int           N;      // Number of symbol nodes (symbols per frame)
  int           M;      // Number of check nodes
  int           dc;     // Maximum dc
  int           dv;     // Maximum dv
  double        Rate;   // Code rate

  // Testbench I/O files:
  char *        stimfilename;
  string        fname;

  // Channel parameters:
  double        SNR;          // Signal to Noise Ratio (Eb/N0) in dB
  double        sigma;        // Channel noise standard deviation
  double        N0;

  unsigned int  nEdges;       // Total number of edges in the code graph

} simparams;

// Declare a global parameter struct available to all modules:
extern simparams p;

//=============================
// Function Predefines:
//=============================
void get_arguments(int argc, char * argv[]);
void append_result_to_data_file(unsigned int errors,
				unsigned int totalbits,
				unsigned int word_errors,
				unsigned int total_words,
				unsigned long totalIterations);


#endif
