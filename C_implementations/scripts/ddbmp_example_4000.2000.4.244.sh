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
ITER=100
DATE=`date +%Y-%m-%d`
LOGNAME=results/results_DDBMP_4000.2000.4.244_${DATE}.tab
RATE=0.5
YMAX=1.0
Q=3

###################################
### SIMULATION COMMANDS         ###
###################################
# NOTE: These commands are executed in parallel. 
# To run the simulations sequentially, remove the
# "nohup" and the "&" from the command below.
for SNR in 3.8 3.9 4.0
do
for YMAX in 0.9 1.0 1.2 1.4 1.6 1.8 2.0 2.2
do
for Q in 2 3 4
do
 echo Running ./bin/decodeDDBMP $Reg48 $RATE $SNR $ITER  $YMAX $Q $LOGNAME $datafile \> tmp/nohup_ddbmp${Q}_${YMAX}_${SNR}.out
 nohup ./bin/decodeDDBMP $Reg48 $RATE $SNR $ITER $YMAX $Q $LOGNAME $datafile > tmp/nohup_ddbmp${Q}_${YMAX}_${SNR}.out &
done
done
done


