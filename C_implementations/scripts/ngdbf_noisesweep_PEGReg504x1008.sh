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
LOGNAME=results/results_noisesweep_SM-NGDBF_PEGReg504x1008
NOISESCALE=0.96
LAMBDA=0.98
ALPHA=2.25
RATE=0.5
SNR=3.0
WINDOWSIZE=64

###################################
### SIMULATION COMMANDS         ###
###################################
# NOTE: These commands are executed in parallel. 
# To run the simulations sequentially, remove the
# "nohup" and the "&" from the command below.
for SNR in 2.5 3.0
do
for NOISESCALE in 0.6 0.65 0.7 0.725 0.75 0.775 0.8 0.85 0.9 0.95 1.05 1.1 1.15 1.2
do
#echo .
 echo Running ./bin/decodeSMNGDBF $ALIST $RATE $SNR $ITER $THETA $LOGNAME $NOISESCALE $LAMBDA $ALPHA $WINDOWSIZE $datafile \> tmp/nohup_ngdbf${SNR}.out
 nohup ./bin/decodeSMNGDBF $ALIST $RATE $SNR $ITER $THETA $LOGNAME $NOISESCALE $LAMBDA $ALPHA $WINDOWSIZE $datafile > tmp/nohup_ngdbf${SNR}.out &
done
done

########################################################
# DIFFERENT PARAMETERS WORK BETTER AT HIGHER SNR VALUES:
########################################################
LAMBDA=0.988
ALPHA=2.3
SNR=3.75

for NOISESCALE in 1.05 1.1 1.15 1.2
do
 echo Running ./bin/decodeSMNGDBF $ALIST $RATE $SNR $ITER $THETA $LOGNAME $NOISESCALE $LAMBDA $ALPHA $WINDOWSIZE $datafile \> tmp/nohup_ngdbf${SNR}.out
 nohup ./bin/decodeSMNGDBF $ALIST $RATE $SNR $ITER $THETA $LOGNAME $NOISESCALE $LAMBDA $ALPHA $WINDOWSIZE $datafile > tmp/nohup_ngdbf${SNR}.out &
done
