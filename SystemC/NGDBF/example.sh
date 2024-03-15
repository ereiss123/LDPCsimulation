#!/bin/bash
export LD_LIBRARY_PATH=/home/reiss/systemc/lib-linux64:$LD_LIBRARY_PATH

ALIST=codes/PegReg/PEGReg504x1008.alist
DATA=codes/PegReg/data.enc
RATE=0.5
SNR=8.0
ITERATIONS=1000
INTERVAL=400
LAMBDA=0.975
THETA=-0.5
YMAX=3
ALPHA=0.95
LOGFILE=example.log
BITS=4

#./ldpcsim.x codes/PegReg/PEGReg504x1008.alist codes/PegReg/data.enc 0.5 4 1000 400 0.975 -0.5 3 3 0.95 test
./ldpcsim.x $ALIST $DATA $RATE $SNR $ITERATIONS $INTERVAL $LAMBDA $THETA $BITS $YMAX $ALPHA $LOGFILE

