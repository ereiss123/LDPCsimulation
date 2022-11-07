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
ITER=300
LOGNAME=results/results_uniform_SM-NGDBF_PEGReg504x1008
NOISESCALE=0.96
LAMBDA=0.95
ALPHA=2.25
WINDOWSIZE=16
RATE=0.5
SNR=3.75

###################################
### SIMULATION COMMANDS         ###
###################################
# NOTE: These commands are executed in parallel. 
# To run the simulations sequentially, remove the
# "nohup" and the "&" from the command below.
#for SNR in 2.5 2.75 3.0 3.25 3.5 3.75 4.0
for NOISESCALE in 0.8 0.9 1.0 1.1 1.2
do
 echo Running ./bin/decodeSMNGDBF $ALIST $RATE $SNR $ITER $THETA $LOGNAME $NOISESCALE $LAMBDA $ALPHA $WINDOWSIZE $datafile \> tmp/nohup_ngdbf${SNR}.out
 nohup ./bin/decodeSMNGDBF $ALIST $RATE $SNR $ITER $THETA $LOGNAME $NOISESCALE $LAMBDA $ALPHA $WINDOWSIZE $datafile > tmp/nohup_ngdbf${SNR}.out &
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
