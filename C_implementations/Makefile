INC = ./inc
SRC = ./src
OBJ = ./obj
BIN = ./bin
CC = g++
CFLAGS = -g -I$(INC) 
LIBFLAGS = -L/usr/local/lib
LIBS= -lm -lgsl -lgslcblas

all: nrutil r alist decodeStochasticNGDBF decodeMGDBF decodeSGDBF decodeSMGDBF decodeMNGDBF decodeSMNGDBF decodeSATGDBF decodeATGDBF decodeMinSum decodeOffsetMinSum decodeNormalizedMinSum decodeBP decodeDDBMP redecodeStatistics decodeRSMNGDBF replayGDBF errtopng

nrutil:$(SRC)/nrutil.cpp
	$(CC) $(CFLAGS) -c -o $(OBJ)/$@.o $(SRC)/$@.cpp

r:$(SRC)/r.cpp
	$(CC) $(CFLAGS) -c -o $(OBJ)/$@.o $(SRC)/$@.cpp

alist:$(SRC)/alist.cpp
	$(CC) $(CFLAGS) -c -o $(OBJ)/$@.o $(SRC)/$@.cpp

errtopng: $(SRC)/errtopng.cpp
	$(CC) $(CFLAGS) -o bin/$@ $(SRC)/errtopng.cpp -lm -lpng

decodeMGDBF: $(SRC)/decodeGDBF.cpp 
	$(CC) $(CFLAGS) -lm -o bin/$@ -D modeswitching $(OBJ)/*.o $(SRC)/decodeGDBF.cpp 

decodeSGDBF: $(SRC)/decodeGDBF.cpp 
	$(CC) $(CFLAGS) -lm -o bin/$@ -D sequentialmode $(OBJ)/*.o $(SRC)/decodeGDBF.cpp 

decodeStochasticNGDBF: $(SRC)/decodeGDBF.cpp
	$(CC) $(CFLAGS) -lm -o bin/$@ -D quantizeSamples -D quantizeProbabilities -D weightSyndromes -D saturateSamples  $(OBJ)/*.o $(SRC)/decodeGDBF.cpp

decodeMNGDBF: $(SRC)/decodeGDBF.cpp 
	$(CC) $(CFLAGS) -lm -o bin/$@ -D addNoise -D thresholdAdaptation -D weightSyndromes -D saturateSamples  $(OBJ)/*.o $(SRC)/decodeGDBF.cpp 

decodeSMNGDBF: $(SRC)/decodeGDBF.cpp 
	$(CC) $(CFLAGS) -lm -o bin/$@ -D addNoise -D thresholdAdaptation -D weightSyndromes -D outputSmoothing -D saturateSamples  $(OBJ)/*.o $(SRC)/decodeGDBF.cpp 

redecodeStatistics: $(SRC)/newstat.cpp
	$(CC) $(CFLAGS) $(LIBFLAGS) -o bin/$@ -D addNoise -D thresholdAdaptation -D weightSyndromes -D outputSmoothing -D saturateSamples  $(OBJ)/*.o $(SRC)/newstat.cpp $(LIBS)

replayGDBF: $(SRC)/replayGDBF.cpp
	$(CC) $(CFLAGS) $(LIBFLAGS) -o bin/$@ -D addNoise -D thresholdAdaptation -D weightSyndromes -D outputSmoothing -D saturateSamples  $(OBJ)/*.o $(SRC)/replayGDBF.cpp $(LIBS)

decodeRSMNGDBF: $(SRC)/RNGDBF.cpp 
	$(CC) $(CFLAGS) -lm -o bin/$@ -D redecode -D addNoise -D thresholdAdaptation -D weightSyndromes -D outputSmoothing -D saturateSamples  $(OBJ)/*.o $(SRC)/RNGDBF.cpp 


decodeSMGDBF: $(SRC)/decodeGDBF.cpp 
	$(CC) $(CFLAGS) -lm -o bin/$@ -D outputSmoothing  $(OBJ)/*.o $(SRC)/decodeGDBF.cpp 

decodeSATGDBF: $(SRC)/decodeGDBF.cpp 
	$(CC) $(CFLAGS) -lm -o bin/$@ -D thresholdAdaptation -D outputSmoothing  $(OBJ)/*.o $(SRC)/decodeGDBF.cpp 

decodeATGDBF: $(SRC)/decodeGDBF.cpp 
	$(CC) $(CFLAGS) -lm -o bin/$@ -D thresholdAdaptation  $(OBJ)/*.o $(SRC)/decodeGDBF.cpp 

decodeMinSum: $(SRC)/decodeMinSum.cpp
	$(CC) $(CFLAGS) -lm -o bin/$@ $(OBJ)/*.o $(SRC)/decodeMinSum.cpp

decodeOffsetMinSum: $(SRC)/decodeMinSum.cpp
	$(CC) $(CFLAGS) -lm -o bin/$@ -D quantizeSamples -D offsetMS $(OBJ)/*.o $(SRC)/decodeMinSum.cpp

decodeNormalizedMinSum: $(SRC)/decodeMinSum.cpp
	$(CC) $(CFLAGS) -lm -o bin/$@ -D quantizeSamples -D normalizedMS $(OBJ)/*.o $(SRC)/decodeMinSum.cpp

decodeBP: $(SRC)/decodeBP.cpp
	$(CC) $(CFLAGS) -lm -o bin/$@ $(OBJ)/*.o $(SRC)/decodeBP.cpp

decodeDDBMP: $(SRC)/decodeDDBMP.cpp
	$(CC) $(CFLAGS) -lm -o bin/$@ $(OBJ)/*.o $(SRC)/decodeDDBMP.cpp

clean:
	-rm obj/* bin/* *~ core src/*~ inc/*~ 
