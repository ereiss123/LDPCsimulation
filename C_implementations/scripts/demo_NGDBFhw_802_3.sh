#!/bin/bash

SNR=4.5
T=600
NF=10000
seed=`date +"%s"`

logfile=test_msgs

# Usage: ./bin/NGDBFhw alist SNR numFrames seed logfilename [codeword filename]
# Test example for debug: ./codes/802_3/802_3_H.alist 4 100 1234 debughw

# Define a timestamp function
timestamp() {
  date +"%m-%d.%H.%M"
}

# Test if log file exists and write header if not
if [ ! -f $logfile ]; then
    echo -e "SNR\tNberr\tNwerr\tBER\tTavg\tFER\tNbit\tNw\tT\ttheta\tnoiseScale\tw\tYmax\tnumPhases\tNQ\tseed" >> $logfile 
fi


./bin/NGDBFhw ./codes/802_3/802_3_H.alist $SNR $NF $seed $logfile > loghw_$(timestamp).tmp 
