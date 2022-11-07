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
ALGNAME=replayGDBF
STATENUM=52
STATEDATE=15-01-2015

###################################
### DEFAULT PARAMETERS          ###
###################################
NR=100
NF=200
ITER=1000
LOGNAME=replay_SMNGDBF_802.3_${STATENUM}
WINDOWSIZE=64
SNR=4.0
YMAX=2.5

THETA=-0.525
ALPHA=0.5
LAMBDA=1
NOISESCALE=0.92

line="./bin/$ALGNAME $ALIST $RATE $SNR $ITER $NR $THETA tmp/${STATEDATE}_RanState_${STATENUM}.state ${LOGNAME} $NOISESCALE $LAMBDA $ALPHA $WINDOWSIZE $YMAX $datafile"

echo Running $line

nohup $line &
