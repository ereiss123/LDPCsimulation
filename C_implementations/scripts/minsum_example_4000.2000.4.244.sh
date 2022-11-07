#!/bin/bash

###################################
### ALIST FILE AND ENCODED DATA ###
###################################
# The simulator uses alist files in the same
# format used by Radford Neal's LDPC tools. 
ALIST=./codes/4000.2000.4.244/4000.2000.4.244.alist
datafile=./codes/4000.2000.4.244/data.enc


###################################
### DEFAULT PARAMETERS          ###
###################################
ITER=15
LOGNAME=results/results_MinSum_${ITER}_4000.2000.4.244
SNR=3.75
RATE=0.5

###################################
### SIMULATION COMMANDS         ###
###################################
for ITER in 10 15
do
for SNR in 3.05 3.1
do
 echo Running ./bin/decodeMinSum $ALIST $RATE $SNR $ITER $LOGNAME  $datafile \> tmp/nohup_minsum${SNR}.out
 nohup ./bin/decodeMinSum $ALIST $RATE $SNR $ITER  $LOGNAME  $datafile > tmp/nohup_minsum${SNR}.out &
done
done
