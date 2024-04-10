/*=========================================================================
** decoder.h
**
** Date: March, 2024
**
** Authors: Eric Reiss, Chris Winstead
**          Utah State University
**
** Based on original source from
** Sept, 2011
** Copyright Chris Winstead
** Utah State University
**
** DESCRIPTION:
**
** This header defines the top-level interface and control behavior
** for a decoder. The decoder module should connect directly to the
** testbench module, and contains the entire modular and functional
** heirarchy of the decoder under test.
==========================================================================*/
#ifndef DECODER_H
#define DECODER_H

#include "systemc.h"
#include "sc_vector.h"
#include <cmath>

#include "ldpcsim.h"
#include "nodes.h"



// ===================================================================
// DECODER CLASS
// ===================================================================
class decoder : public sc_module
{
 public:
//=============================
// Allocator for symbol nodes:
//=============================
struct symnode_creator
{
  symnode_creator(alist_struct a)
  {alist = a;}
  
  symnode * operator() (const char * name, size_t idx)
  {
    return new symnode(name, alist.num_mlist[idx]);
  }
  
  alist_struct alist;
};

//=============================
// Allocator for check nodes:
//=============================
struct checknode_creator
{
  checknode_creator(alist_struct a)   {alist = a;}
  checknode * operator() (const char * name, size_t idx)
  {
    return new checknode(name, alist.num_nlist[idx]);
  }
  alist_struct alist;
};


  //-------------------------------------------
  // TOP-LEVEL IO for the DECODER:
  //-------------------------------------------
  sc_vector<sc_in<double> >             y;        // Input from channel
  sc_vector<sc_out<bool> >               d;        // Decisions
  sc_in<bool>                           clk;      // Clock
  sc_in<bool>                           rst;      // Reset (new frame available)
  sc_out<bool>                          ready;    // Ready to decode
  sc_out<bool>                          finished; // Finished decoding

  //-------------------------------------------
  // Constructor for decoder class
  //-------------------------------------------
  SC_HAS_PROCESS(decoder);
 decoder(sc_module_name name) : sc_module(name),
    // SC_VECTOR INTIALIZATIONS HERE:
    y("y",p.N),
    d("d",p.N),
    stoc("stoc",p.nEdges),
    ctos("ctos",p.nEdges),
    stop("stop",p.M),
    r("r",p.N/p.q),
    S("S"),
    C("C"),
    Sd("Sd",p.N/p.q)
    {
      // -----------------------------------------------
      // Initialize Symbol and Check node module arrays:
      // -----------------------------------------------
      // See ALLOCATOR definitions near the bottom of this file
      S.init(p.N/p.q,symnode_creator(p.alist));
      C.init(p.M,checknode_creator(p.alist));

      // -----------------------------------------
      // Establish interleaver connections:
      // -----------------------------------------
      // The connections are defined by
      // a matrix called "mlist" which is
      // part of the alist file. The
      // connections are specified as follows:
      // If mlist[i][j] = k, then symnode i
      // connects to checknode k on its jth edge.

      vector<unsigned int>	symbol_edges(p.N/p.q,0);
      vector<unsigned int>	check_edges(p.M,0);

      int signal_index = 0;
      
      for (int i=0; i<p.N; i++){
        S[i].r(r[i]);
        
        S[i].d(Sd[i]);
        S[i].rst(rst);
        S[i].clk(clk);
	
        for (int j=0; j<p.alist.num_mlist[i]; j++){
          int check_node_index = p.alist.mlist[i][j] - 1;

          S[i].to_check[j](stoc[signal_index]);
          C[check_node_index].to_check[check_edges[check_node_index] ](stoc[signal_index]);

          S[i].from_check[j](ctos[signal_index]);
          C[check_node_index].from_check[check_edges[check_node_index] ](ctos[signal_index]);

          check_edges[check_node_index]++;
          symbol_edges[i]++;
          signal_index++;
        }
      }

      for (int i=0; i<p.M; i++){
        C[i].stop(stop[i]);
        C[i].clk(clk);
        C[i].rst(rst);
      }

      // Declare top-level control behavior:
      SC_METHOD(behavior);
      sensitive << clk.pos();

      // Setup initial state:
      ready.initialize(false);
      first_frame = true;
      clock_count = 0;
    }

  //----------------------------------
  // Declare Internal Signals:
  //----------------------------------
  sc_vector<sc_signal<message_type > >	 stoc;    // Symbol-to-check messages
  sc_vector<sc_signal<message_type > >	 ctos;    // Check-to-symbol messages
  sc_vector<sc_signal<double > >         r;       // Received input at decoder (quantized)
  sc_vector<sc_signal<bool> >            stop;    // Stopping signals from check nodes
  sc_vector<sc_signal<bool> >            Sd;      // Decision signals from symbol nodes

  //----------------------------------
  // Declare Submodules:
  //----------------------------------
  sc_vector<symnode>	S;
  sc_vector<checknode>	C;

  //----------------------------------
  // Declare State Variables and Parameters:
  //----------------------------------
  bool         first_frame;
  unsigned int clock_count;
  unsigned int iteration_count;


  //===========================================
  // BEHAVIOR DEFINITION
  //===========================================
  void behavior()
  {

    // --------------------------
    // RESET behavior:
    // --------------------------
    if (rst.read()){
      iteration_count = 0;
      finished.write(false);

      if (first_frame)
        if (clock_count > p.N){
          first_frame = false;
          ready.write(true);
        }
        else
          clock_count++;
      if (!first_frame){
        // Deliver quantized inputs into the symbol nodes:
        for (int i=0; i<p.N; i++)
          r[i].write(quantize(y[i].read()));
      }
    }

    // --------------------------
    // Termination behavior:
    // --------------------------
    else if ((!finished) && (!first_frame)){
      iteration_count++;

      bool all_stop = true;
      for (int i=0; i<p.M; i++){
        if (!stop[i].read())
          all_stop = false;
      }
      for (int i=0; i<p.N; i++)
	d[i].write(Sd[i].read());

      if (all_stop || (iteration_count > p.iterations_per_frame))
        finished.write(true);
    }
  }
};






#endif

