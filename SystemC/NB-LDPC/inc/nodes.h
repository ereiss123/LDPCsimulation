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
    dv(_dv), theta(_theta), lambda(_lambda), sigma(_sigma),
    from_check("from_check",_dv), to_check("to_check",_dv)
    {
      E = 0;
      x = 0;
      w = p.alpha*p.Ymax/dv;

      // Declare node behavior method:
      SC_METHOD(behavior);
      sensitive << clk.pos();

      // rnd_out.initialize(0.0);
    }

 private:
  unsigned int dv;
  int x;
  double E;
  double sigma;
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


    }
    //------------------------
    // Normal behavior:
    //------------------------
    else{

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
