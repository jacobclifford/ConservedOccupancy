CC = g++-4.4

GSL_DIR = /usr/local

INC = -I$(GSL_DIR)
CFLAGS = -O3 $(INC)/include 

LIB_DIRS = -L$(GSL_DIR)/lib
LIBS = -lgsl -lgslcblas -lm
LFLAGS = -O3 $(LIB_DIRS) $(LIBS)

BIN =  seq2exp
all: $(BIN)

clean:
	rm -f $(BIN)

Tools.o : Tools.h Tools.cpp
	$(CC) $(CFLAGS) -c Tools.cpp
ExprPredictor.o : Tools.h ExprPredictor.h  ExprPredictor.cpp
	$(CC) $(CFLAGS) -c ExprPredictor.cpp
seq2exp.o : Tools.h ExprPredictor.h seq2exp.cpp
	$(CC) $(CFLAGS) -c seq2exp.cpp

seq2exp : Tools.o ExprPredictor.o seq2exp.o 
	$(CC) $(LIB_DIRS) -o $@ Tools.o ExprPredictor.o seq2exp.o $(LIBS)


