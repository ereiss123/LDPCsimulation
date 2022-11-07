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
LOGNAME=results/results_M-NGDBF_PEGReg504x1008
NOISESCALE=0.95
LAMBDA=0.99
ALPHA=2.25
RATE=0.5
SNR=3.75
Ymax=2.5

###################################
### SIMULATION COMMANDS         ###
###################################
# NOTE: These commands are executed in parallel. 
# To run the simulations sequentially, remove the
# "nohup" and the "&" from the command below.
#for SNR in 3.5 3.75 4.0 4.25 4.5
#do
#for LAMBDA in 0.9 0.91 0.92 0.93 0.94 0.95 0.96 0.97 0.98 0.99 0.991 0.992 0.993 0.994 0.995
#do
# echo Running ./bin/decodeMNGDBF $ALIST $RATE $SNR $ITER $THETA $LOGNAME $NOISESCALE $LAMBDA $ALPHA $Ymax $datafile \> tmp/nohup_ngdbf${SNR}.out
# nohup ./bin/decodeMNGDBF $ALIST $RATE $SNR $ITER $THETA $LOGNAME $NOISESCALE $LAMBDA $ALPHA $Ymax $datafile > tmp/nohup_ngdbf${SNR}.out &
#done
#done

########################################################
# DIFFERENT PARAMETERS WORK BETTER AT LOWER SNR VALUES:
########################################################
#NOISESCALE=0.75
#LAMBDA=0.98

for SNR in 3.0 3.25
do
for THETA in -0.6 -0.7 -0.8 -1.0
do
for NOISESCALE in 0.65 0.7 0.75 0.8 0.85 0.9 
do
for Ymax in 2.0 2.2 2.4 2.6 2.8 3.0 3.2 3.4 3.6
do
for ALPHA in 2.0 2.25 2.5
do
 echo Running ./bin/decodeMNGDBF $ALIST $RATE $SNR $ITER $THETA $LOGNAME $NOISESCALE $LAMBDA $ALPHA $Ymax $datafile \> tmp/nohup_ngdbf${SNR}.out
 nohup ./bin/decodeMNGDBF $ALIST $RATE $SNR $ITER $THETA $LOGNAME $NOISESCALE $LAMBDA $ALPHA $Ymax $datafile > tmp/nohup_ngdbf${SNR}.out &
done
done
done
done
done