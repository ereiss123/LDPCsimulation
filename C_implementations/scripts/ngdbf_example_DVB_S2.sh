#!/bin/bash

###################################
### ALIST FILE AND ENCODED DATA ###
###################################
# The simulator uses alist files in the same
# format used by Radford Neal's LDPC tools. 
ALIST=./codes/dvbs2_1_2/dvbs2_1_2.alist
datafile=


###################################
### DEFAULT PARAMETERS          ###
###################################
THETA=-1.1
ITER=700
LOGNAME=results/results_SM-NGDBF_dvbs2_1_2
NOISESCALE=0.775
LAMBDA=0.987
ALPHA=2.5
WINDOWSIZE=64
RATE=0.5
SNR=2.8
YMAX=2.5

###################################
### SIMULATION COMMANDS         ###
###################################
# NOTE: These commands are executed in parallel. 
# To run the simulations sequentially, remove the
# "nohup" and the "&" from the command below.
for SNR in 3 3.1 3.2 3.3 3.4
do
 echo Running ./bin/decodeSMNGDBF $ALIST $RATE $SNR $ITER $THETA $LOGNAME $NOISESCALE $LAMBDA $ALPHA $WINDOWSIZE $YMAX $datafile \> tmp/nohup_ngdbf${SNR}.out
 nohup ./bin/decodeSMNGDBF $ALIST $RATE $SNR $ITER $THETA $LOGNAME $NOISESCALE $LAMBDA $ALPHA $WINDOWSIZE $YMAX $datafile > tmp/nohup_ngdbf${SNR}.out &
done


