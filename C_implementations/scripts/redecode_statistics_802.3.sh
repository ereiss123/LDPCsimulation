#!/bin/bash

###################################
### ALIST FILE AND ENCODED DATA ###
###################################
# The simulator uses alist files in the same
# format used by Radford Neal's LDPC tools. 
CNAME=802_3
RATE=0.8125
ALIST=./codes/${CNAME}/${CNAME}_H.alist
datafile=
ALGNAME=redecodeStatistics

###################################
### DEFAULT PARAMETERS          ###
###################################
NR=100
NF=200
ITER=1000
LOGNAME=results/statistics_SMNGDBF_802.3_15_Jan_2015
WINDOWSIZE=64
SNR=4.0
YMAX=2.5

THETA=-0.525
ALPHA=0.5
LAMBDA=1
NOISESCALE=0.92

line="./bin/$ALGNAME $ALIST $RATE $SNR $ITER $NR $NF $THETA ${LOGNAME}_R${RATE}_SNR${SNR} $NOISESCALE $LAMBDA $ALPHA $WINDOWSIZE $YMAX $datafile"

echo Running $line

nohup $line &
