#!/bin/bash

###################################
### ALIST FILE AND ENCODED DATA ###
###################################
# The simulator uses alist files in the same
# format used by Radford Neal's LDPC tools. 
ALIST=./codes/4376.282.4.9598/4376.282.4.9598.alist

# NO ENCODED DATAFILE FOR THIS CODE. WILL USE 
# ALL-ZERO MESSAGES


###################################
### DEFAULT PARAMETERS          ###
###################################
THETA=-0.7
ITER=300
LOGNAME=results/results_SM-NGDBF_4376.282.4.9598
NOISESCALE=0.65
LAMBDA=0.993
ALPHA=0.75
WINDOWSIZE=64
RATE=0.9356
YMAX=2.5

###################################
### SIMULATION COMMANDS         ###
###################################
# NOTE: These commands are executed in parallel. 
# To run the simulations sequentially, remove the
# "nohup" and the "&" from the command below.
for SNR in 4 4.25 4.75 5.0 5.25 5.5
do
 echo Running ./bin/decodeSMNGDBF $ALIST $RATE $SNR $ITER $THETA $LOGNAME $NOISESCALE $LAMBDA $ALPHA $WINDOWSIZE $YMAX $datafile \> tmp/nohup_ngdbf${SNR}.out
 nohup ./bin/decodeSMNGDBF $ALIST $RATE $SNR $ITER $THETA $LOGNAME $NOISESCALE $LAMBDA $ALPHA $WINDOWSIZE $YMAX $datafile > tmp/nohup_ngdbf${SNR}.out &
done

