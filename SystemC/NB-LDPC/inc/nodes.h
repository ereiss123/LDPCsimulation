/*=========================================================================
** nodes.h
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
   Contains module definitions for symnode and checknode.
   These models use the min-sum algorithm. Symnodes and
   checknodes are latched on opposite clock phases. Latches
   are used on both node types in order to prevent glitch
   propagation on the interleaver.
==========================================================================*/

#ifndef NODES_H
#define NODES_H
#include "rand.h"
#include "systemc.h"
#include "sc_vector.h"
#include <iostream>
#include <cmath>
#include <vector>
#include <algorithm>

#include "ldpcsim.h"
#include <itpp/comm/galois.h>

using namespace std;

//======Debug Messages ========//
#define VERBOSE false
#define dblog(s) if (VERBOSE) {cout << s;}

//=================================================
// Symbol Node implementation
//=================================================
class symnode : public sc_module
{
 public:
  sc_in<vector<message_type> >		r;           // Channel information
  sc_vector<sc_in<message_type> >	from_check;  // Vector of symbol probability vectors
  sc_vector<sc_out<message_type> >	to_check;
  sc_out<int> 		          	d;
  sc_in<bool>				clk;
  sc_in<bool>				rst;


  SC_HAS_PROCESS(symnode);

  symnode(sc_module_name name,
	  int    _dv,
	  int    _q) : sc_module(name),
		       dv(_dv),
		       q(_q),
		       from_check("from_check",_dv),
		       to_check("to_check",_dv)
    {

      b = int(log2(q));

      // Declare node behavior method:
      SC_METHOD(behavior);
      sensitive << clk.pos();
    }

 private:
  unsigned int dv;

  int b;    // number of bits in field
  int q;    // field size, 2^b


  //************* SYMNODE BEHAVIOR *********************//
  void behavior()
  {
    //------------------------
    // Reset behavior:
    //------------------------
    if (rst.read())
    {

      message_type soft_symbol;
      soft_symbol = r.read();


      /* move this out to another module, so the symnode's "soft_symbol"
	 input is already a message_type containing probabilities, so the
	 nodes only ever work with message_type inputs and outputs.

      // calculate the likelihood for each symbol
      message_type probabilities;
      double likelihood;

      for(int sym = 0; sym < q; sym++){
        int temp_sym = sym;
        int bit = 0;
        double prob = 1.0;
        // calculate the likelihood for each bit in the symbol
        for(int i=0; i<b; i++){
          bit = temp_sym%2;
          temp_sym /= 2;
          likelihood = 1/(1+exp(2*abs(soft_symbol[i])/sigma_sq));
          prob = (bit==0)?(prob*(1-likelihood)):(prob*likelihood);
        }
        probabilities.emplace_back(prob);
      }
      */


      // Write initial probabilities to all adjacent check nodes
      for (int i=0; i<dv; i++){
        to_check[i].write(soft_symbol);
      }

      // Get initial decision
      // find the index of the max probability
      auto iter = find(soft_symbol.begin(),soft_symbol.end(),max_element(soft_symbol.begin(),soft_symbol.end()));
      // 'iter is an address, subtract from the beginning address to get the index
      int symbol = iter-soft_symbol.begin();
      dblog(iter << "\t" << soft_symbol.begin() << "\t" << symbol << endl);
      d.write(symbol);
    }

    //------------------------
    // Normal behavior:
    //------------------------
    else
    {
      // TODO: Figure out permutation

      // update symbol probability for each symbol for each vector
      message_type prob;
      for(int i=0; i<dv; i++){
        message_type in = from_check[i].read();
        for(int sym=0; sym < q; sym++){ // loop through each element of GF(q)
          if(prob.size() < q){ // initial pass
            prob.emplace_back(in[sym]);
          }else{
            prob[sym] *= in[sym];
          }
        }
      }

      // Write out probabilities
      for(int i=0; i<dv; i++){
        to_check[i].write(prob);
      }

      // Update symbol node decision
      // find the index of the max probability
      auto it = find(prob.begin(),prob.end(),max_element(prob.begin(),prob.end()));
      // 'it' is an address, subtract from the beginning address to get the index
      int symbol = it-probabilities.begin();
      dblog(it << "\t" << probabilities.begin() << "\t" << symbol << endl);
      d.write(symbol);
    }

  } // end behavior
}; // end symnode class


//=================================================
// Check Node implementation
//=================================================
class checknode : public sc_module
{
 public:
  sc_vector<sc_out<message_type > > 	from_check; // check node -> symbol node
  sc_vector<sc_in<message_type > > 	to_check; // symbol_node -> check node
  sc_out<bool>				stop;
  sc_in<bool>				clk;
  sc_in<bool>				rst;

  SC_HAS_PROCESS(checknode);
  checknode(sc_module_name name,
	    int _dc) : sc_module(name),
		       dc(_dc),
		       from_check("from_check",_dc),
		       to_check("to_check",_dc)
  {
    SC_METHOD(behavior);
    for (int i=0; i<dc; i++)
      sensitive << to_check[i];

    stop.initialize(false);

    // Generate LUT
    // for each element in GF(q)
    for(int sym=0; sym<p.alist.q; sym++)
    {
      itpp::GF gf_element(p.alist.q,sym);
      // find combinatations of dc length
      std::vector<int> combos(dc, 0);

      int idx = 0;
      while(!((idx == dc-1) && (combos[idx] == p.alist.q)))
      {
        itpp::GF gf_sum(p.alist.q); // returns the 0 element
        // Sum the combo
        for(int x : combos)
        {
          itpp::GF gf_temp(p.alist.q,x)
          gf_sum += gf_temp;
        }
        // Check if valid combo was found
        if(gf_sum == gf_element)
        {
          // save valid combos for each GF(q) element
          node_LUT[sym].emplace_back(combos);
        }
        // update combo list
        idx = 0;
        // Treat the combo as a base q number and add 1 to go through every option
        while((combos[idx]+1)%p.alist.q == 0)
        {
          combos[idx] = 0;
          idx++;
        }
        combos[idx]++;
      }
    }
  }

  private:
  int dc;
  // Lookup table to calculate probabilities
  // The first index is the element of GF(q)
  std::vector<std::vector<int>> node_LUT;

  //************* CHECKNODE BEHAVIOR *********************//
  void behavior()
  {
    //------------------------
    // Reset behavior:
    //------------------------
    /*
    if (rst.read())
    {
      None Defined...
    }
    */


    //------------------------
    // Normal behavior:
    //------------------------

    // Collect all neighboring messages into a matrix
    std::vector<message_type> to_check_mat;
    for (int i=0; i<dc; i++)
    {
      to_check_mat.emplace_back(to_check[i].read());
    }

    // Iterate through the lookup table to calculate the outgoing message
    message_type result(p.alist.q);
    for(int i=0; i<p.alist.q; i++)
    {
      // i is the element of GF(q) being calculated for
      // p.LUT is a vector of vectors
      // p.LUT[i] is a vector of pairs
      result[i] = 0;
      int idx = 0;
      for(std::pair<int,int> symbol_pair : p.LUT[i])
      {

      }
    }

    for (int i=0; i<dc; i++)
      from_check[i].write(prod);

    if (prod==1)
      stop.write(true);
    else
      stop.write(false);
  }

};


#endif
