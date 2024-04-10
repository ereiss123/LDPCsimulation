/*==========================================================================================
** ldpcsim.cpp
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

==============================================================================================*/


#include "alist.h"
#include "ldpcsim.h"
#include "LDPC_testbench.h"
#include "nodes.h"
#include "decoder.h"

#include <iostream>
#include <fstream>
#include <vector>

// Global simulation parameters:
simparams p;


//=====================================
// SC_MAIN BODY
//=====================================
int sc_main(int argc, char * argv[])
{
  get_arguments(argc, argv);

  // -----------------------------------------
  // Declare Top-Level Signals:
  // -----------------------------------------
  sc_vector<sc_signal<double> >	         y("y",p.N)*p.q;  // binary bits
  sc_vector<sc_signal<int>  >	     	 d("d",p.N);      // GF(q) decisions
  sc_signal<bool>	     		 rst;
  sc_signal<bool>                        ready;
  sc_signal<bool>                        finished;
  sc_clock				 clk("clock", 10, SC_NS);

  // Testbench modules:
  LDPC_testbench  tb("tb");
  decoder         dec("dec");


  // -----------------------------------------
  // Connect global clock and reset signals,
  // and the main input/output signals:
  // -----------------------------------------

  for (int i=0; i<p.N*p.q; i++)
  {
    tb.y[i](y[i]);
    dec.y[i](y[i]);
  }
  
  for (int i=0; i<p.N; i++)
  {
    tb.d[i](d[i]);
    dec.d[i](d[i]);
  }

  tb.clk(clk);
  tb.rst(rst);
  tb.ready(ready);
  tb.finished(finished);
  dec.clk(clk);
  dec.rst(rst);
  dec.ready(ready);
  dec.finished(finished);

  // Start Simulation:
  sc_start();

  // Done.
  return 0;
}


//===================================
// READ PARAMETERS FROM COMMAND LINE:
//===================================
void get_arguments(int argc, char * argv[])
{
  if (argc != 13)
  {
    cout << "Usage: \n\t"
	 << argv[0]
	 << " <alist_fname> "
	 << "<stim_fname> "
	 << "<code_rate> "
	 << "<SNR> "
	 << "<iterations> "
	 << "<report_interval> "
	 << "<log_fname>\n";
    exit(0);
  }
  else
  {
    p.alist = loadFile(argv[1]);
    p.stimfilename = argv[2];
    p.Rate = atof(argv[3]);
    p.SNR = atof(argv[4]);
    p.iterations_per_frame = atoi(argv[5]);
    p.total_clock_cycles = atoi(argv[6]);

    stringstream ss;
    ss << argv[12] << ".dat";
    p.fname = ss.str();

    p.N = p.alist.M;
    p.M = p.alist.N;
    p.dv = p.alist.biggest_num_m;
    p.dc = p.alist.biggest_num_n;

    p.nEdges = 0;
    for (int i=0; i<p.N; i++)
    {
      p.nEdges += p.alist.num_mlist[i];
    }
    p.N0 = (1/p.Rate)*pow(10.0, -p.SNR/10.0);
    p.sigma = sqrt(p.N0/2.0);

  }
}


//===================================
// REPORT RESULTS TO FILE:
//===================================
void append_result_to_data_file(unsigned int errors, unsigned int totalbits, unsigned int word_errors, unsigned int total_words, unsigned long totalIterations)
{
  ofstream result_file(p.fname.c_str(),ios::app);
  result_file << endl
	      << (double)errors/totalbits << "\t"
	      << (double)word_errors/total_words << "\t"
	      << (double)totalIterations/total_words << "\t"
	      << p.SNR << "\t"
	      << errors << "\t"
	      << totalbits << "\t"
	      << word_errors << "\t"
	      << total_words;
  result_file.close();
}

//==================================
// GENERATE GF(q) LOOKUP TABLE
//==================================

void generate_LUT()
{
  itpp::GF z = itpp::GF(p.alist.q); // return 0th element
  // iterate over every symbol combination
  for(int i = 0; i < p.alist.q ; i++)
  {
    itpp::GF GF_i(p.alist.q,i);
    vector<int> idxs;

    itpp::GF i_gf = itpp::GF(p.alist.q,i);
    for(int j = 0; j < p.alist.q; j++)
    {
      itpp::GF GF_j(p.alist.q,j);
      if(GF_j + GF_i == z)
      {
        pair<int,int> elements(i,j);
        idxs.emplace_back(elements);
      }
    }
    p.LUT.emplace_back(idxs);
  }
}


