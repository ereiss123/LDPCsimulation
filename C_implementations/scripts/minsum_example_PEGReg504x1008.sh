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
ITER=8
LOGNAME=results/results_MinSum_${ITER}_PEGReg504x1008
SNR=3.75
RATE=0.5

###################################
### SIMULATION COMMANDS         ###
###################################
for SNR in 1.6 1.8 2.0 2.2 2.4 2.6 2.8 3.0 3.2 3.4 3.6 3.8
do
 echo Running ./bin/decodeMinSum $ALIST $RATE $SNR $ITER $LOGNAME  $datafile \> tmp/nohup_minsum${SNR}.out
 nohup ./bin/decodeMinSum $ALIST $RATE $SNR $ITER  $LOGNAME  $datafile > tmp/nohup_minsum${SNR}.out &
done

