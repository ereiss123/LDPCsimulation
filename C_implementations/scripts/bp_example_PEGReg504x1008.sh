#!/bin/bash

###################################
### ALIST FILE AND ENCODED DATA ###
###################################
# The simulator uses alist files in the same
# format used by Radford Neal's LDPC tools. 
ALIST=./codes/PEGReg504x1008/PEGReg504x1008.alist
datafile=./codes/PEGReg504x1008/data.enc


###################################
### DEFAULT PARAMETERS          ###
###################################
ITER=100
LOGNAME=results/results_BP_${ITER}_PEGReg504x1008
SNR=3.75
RATE=0.5

###################################
### SIMULATION COMMANDS         ###
###################################
for SNR in 1.6 1.8 2.0 2.1 2.2 2.3 2.4 2.5 2.6 
do
 echo Running ./bin/decodeBP $ALIST $RATE $SNR $ITER $LOGNAME  $datafile \> tmp/nohup_bp${SNR}.out
 nohup ./bin/decodeBP $ALIST $RATE $SNR $ITER  $LOGNAME  $datafile > tmp/nohup_bp${SNR}.out &
done

