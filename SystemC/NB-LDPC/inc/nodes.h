/*=========================================================================
** nodes.h
** Date: Sept, 2011
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
  sc_vector<sc_in<double> >			        r;
  sc_vector<sc_vector<sc_in<message_type> > >	from_check;  // Vector of symbol probability vectors
  sc_vector<sc_vector<sc_out<message_type> > >	to_check;
  sc_out<bool> 		          	d;
  sc_in<bool>				clk;
  sc_in<bool>				rst;
  // sc_in<double>                         rnd_in;
  // sc_out<double>                        rnd_out;

  SC_HAS_PROCESS(symnode);
  // NOTE: Here the theta parameter is passed as a pointer so that it can be globally adapted.
  symnode(sc_module_name name, int _dv, message_type _theta, message_type _lambda, double _sigma) : sc_module(name),
    dv(_dv), theta(_theta), lambda(_lambda), sigma(_sigma), q(_q),
    from_check("from_check",_dv), to_check("to_check",_dv)
    {
      E = 0;
      x = 0;
      w = p.alpha*p.Ymax/dv;
      b = int(log2(q));
      sigma_sq = pow(sigma,2);

      // Declare node behavior method:
      SC_METHOD(behavior);
      sensitive << clk.pos();

      // rnd_out.initialize(0.0);
    }

 private:
  unsigned int dv;
  int x;
  int b; // number of bits in field
  int q; // field size, 2^b
  double E;
  double sigma;
  double sigma_sq;
  double lambda;
  double theta;
  double local_theta;
  double w;

  //************* SYMNODE BEHAVIOR *********************//
  void behavior()
  {
    //------------------------
    // Reset behavior:
    //------------------------
    if (rst.read())
    {
      local_theta = theta;
      // calculate symbol likelihood, read in q bits
      vector<double> soft_symbol;
      for(int i = 0; i < b; i++){
        soft_symbol.emplace_back(r.read()); //get soft symbol from channel, LSB at index 0
      }

      // calculate the likelihood for each symbol

      sc_vector<message_type> probabilities;
      double likelihood;
      for(int sym = 0; sym < q; sym++){
        int temp_sym = sym;
        int bit = 0;
        double prob = 1.0;
        for(int i=0; i<b; i++){
          bit = temp_sym%2;
          temp_sym/=2;
          likelihood = 1/(1+exp(2*abs(soft_symbol[i])/sigma_sq));
          prob = (bit==0)?(prob*(1-likelihood)):(prob*likelihood);
        }
        probabilities.emplace_back(prob);
      }

    // Write initial probabilities to all adjacent check nodes
    for (int i=0; i<dv; i++){
	    to_check[i].write(probabilities);
    }

    // Get initial decision
    auto it = find(probabilites.begin(),probabilities.end(),max(probabilities));
    int symbol = it-probabilities.begin();
    if(symbol > 0){
      d.write(false);
    }else{
      d.write(true);
    }

    }
    //------------------------
    // Normal behavior:
    //------------------------
    else{
      // Need to figure out permutation


      // update symbol probability for each symbol for each vector
      vector<message_type> prob;
      for(int i=0; i<dv; i++){
        // Update vector i probabilities
        prob.clear();
        for(int j=0; j<dv; j++){
          if(i != j){
            for(int sym=0; sym<q; sym++){
              if(prob.size() < 1){
                prob.append(from_check[j][sym]);
              }else{
                prob[i][sym] *= from_check[j][sym]
              }
            }
          }
        }
        // Write out probabilities
        to_check[i].write(prob);
      }

    }
  }
};


//=================================================
// Check Node implementation
//=================================================
class checknode : public sc_module
{
 public:
  sc_vector<sc_vector<sc_out<message_type > > >	from_check;
  sc_vector<sc_vector<sc_in<message_type > > >	to_check;
  sc_out<bool>				stop;
  sc_in<bool>				clk;
  sc_in<bool>				rst;

  SC_HAS_PROCESS(checknode);
  checknode(sc_module_name name, int _dc) : sc_module(name), dc(_dc),from_check("from_check",_dc), to_check("to_check",_dc)
  {
    SC_METHOD(behavior);
    for (int i=0; i<dc; i++)
      sensitive << to_check[i];

    stop.initialize(false);
  }

  private:
  int dc;

  //************* CHECKNODE BEHAVIOR *********************//
  void behavior()
  {
    //------------------------
    // Reset behavior:
    //------------------------
    /*
    if (rst.read())
    {
	    None defined...
    }
    */


    //------------------------
    // Normal behavior:
    //------------------------

  //   dblog("ChkNode: In=");
  //   int prod = 1;

  //   for (int i=0; i<dc; i++){
  //     dblog(to_check[i].read() << ", ");
  //     prod *= to_check[i].read();
  //   }

  //   dblog("prod=" << prod << endl);

  //   for (int i=0; i<dc; i++)
  //     from_check[i].write(prod);

  //   if (prod==1)
  //     stop.write(true);
  //   else
  //     stop.write(false);
  }

};


#endif
