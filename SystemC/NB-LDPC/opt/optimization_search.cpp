/*========================================================================
  optimization_search.cpp
  Chris Winstead, 2013

  This is a simple program for sweeping across decoder parameters to
  identify best performance.
  ========================================================================*/

#include <iostream>
#include <fstream>
#include <sstream>
#include <stdlib.h>
#include <stdio.h>
#include <cmath>

using namespace std;

int main()
{
  int maxIterations = 300;
  int maxClocks = 10000;
  
  int startPrecision = 5;
  int endPrecision = 5;
  
  double startYmax = 1.75;
  double endYmax = 1.75;
  double stepYmax = 0.01;

  double minSNR = 2.00;
  double maxSNR = 3.5;
  double SNRstep = 0.25;
  
  double startLambda = 0.985;
  double endLambda = 0.985;
  double stepLambda = 0.005;

  double initial_theta = -0.8;
  double final_theta = -0.8;
  double theta_step = 0.1;

  double alpha = 1.5;

  string simID("M-NGDBF_300_local_adaptation_sample_reuse_smoothed_delayed12");
  string result_filename("M-NGDBF_300_local_adaptation_sample_reuse_smoothed_delayed12.dat");
  string alist_filename("../codes/PegReg/PEGReg504x1008.alist");
  string data_filename("../codes/PegReg/data.enc");
  double Rate = 0.5;

  string command_name("./ldpcsim.x");
  

  for (int precision=startPrecision; precision <= endPrecision; precision++)
    for (double SNR=minSNR; SNR<=maxSNR; SNR+=SNRstep)
      {     
	for (double Ymax=startYmax; Ymax <= endYmax; Ymax += stepYmax)
	  for (double theta=initial_theta; theta<=final_theta; theta+=theta_step) {
	    for (double lambda=startLambda; lambda<=endLambda; lambda+=stepLambda)
	      {

		stringstream command_string;
		command_string << command_name << " " << alist_filename << " " << data_filename << " "
			       << Rate << " " << SNR << " " << maxIterations << " " << maxClocks << " " 
			       << lambda << " " << theta << " " << precision << " "  << Ymax << " " 
			       << alpha << " " << simID << "; sleep 1";
		cout << "Executing " << command_string.str() << endl;
      
		system(command_string.str().c_str());      	  
	      }
      
	  }
      }
  return 0;
}

