#!/bin/bash

###################################
### ALIST FILE AND ENCODED DATA ###
###################################
# The simulator uses alist files in the same
# format used by Radford Neal's LDPC tools. 
Reg48=./codes/4000.2000.4.244/4000.2000.4.244.alist
datafile=./codes/4000.2000.4.244/data.enc


###################################
### DEFAULT PARAMETERS          ###
###################################
THETA=-0.7
ITER=200
LOGNAME=results/results_SM-NGDBF_4000.2000.4.244
NOISESCALE=0.75
LAMBDA=0.99
ALPHA=2.2
WINDOWSIZE=64
RATE=0.5
YMAX=2.5

###################################
### SIMULATION COMMANDS         ###
###################################
# NOTE: These commands are executed in parallel. 
# To run the simulations sequentially, remove the
# "nohup" and the "&" from the command below.
for SNR in 2 2.2 2.4 2.6 2.8 3.0 3.05 
do
 echo Running ./bin/decodeSMNGDBF $Reg48 $RATE $SNR $ITER $THETA $LOGNAME $NOISESCALE $LAMBDA $ALPHA $WINDOWSIZE $datafile $YMAX \> tmp/nohup_ngdbf${SNR}.out
 nohup ./bin/decodeSMNGDBF $Reg48 $RATE $SNR $ITER $THETA $LOGNAME $NOISESCALE $LAMBDA $ALPHA $WINDOWSIZE $datafile $YMAX > tmp/nohup_ngdbf${SNR}.out &
done


