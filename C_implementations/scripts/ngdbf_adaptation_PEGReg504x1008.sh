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
THETA=-0.9
ITER=100
LOGNAME=results/results_adaptation_M-NGDBF_PEGReg504x1008
NOISESCALE=0.96
LAMBDA=0.95
ALPHA=2.25
RATE=0.5
SNR=4.0

###################################
### SIMULATION COMMANDS         ###
###################################
# NOTE: These commands are executed in parallel. 
# To run the simulations sequentially, remove the
# "nohup" and the "&" from the command below.
for SNR in 3.0 3.5 4.0
do
for LAMBDA in 0.92 0.93 0.94 0.95 0.96 0.965 0.97 0.975 0.98 0.985 0.99 0.995 0.9975
do
 echo Running ./bin/decodeMNGDBF $ALIST $RATE $SNR $ITER $THETA $LOGNAME $NOISESCALE $LAMBDA $ALPHA $datafile \> tmp/nohup_ngdbf${SNR}.out
 nohup ./bin/decodeMNGDBF $ALIST $RATE $SNR $ITER $THETA $LOGNAME $NOISESCALE $LAMBDA $ALPHA $datafile > tmp/nohup_ngdbf${SNR}.out &
done
done
########################################################
# DIFFERENT PARAMETERS WORK BETTER AT LOWER SNR VALUES:
########################################################
#NOISESCALE=0.75
#LAMBDA=0.98

#for SNR in 2 2.25 2.75 3.0 3.25
#do
# echo Running ./bin/decodeSMNGDBF $ALIST $RATE $SNR $ITER $THETA $LOGNAME $NOISESCALE $LAMBDA $ALPHA $WINDOWSIZE $datafile \> tmp/nohup_ngdbf${SNR}.out
# nohup ./bin/decodeSMNGDBF $ALIST $RATE $SNR $ITER $THETA $LOGNAME $NOISESCALE $LAMBDA $ALPHA $WINDOWSIZE $datafile > tmp/nohup_ngdbf${SNR}.out &
#done
