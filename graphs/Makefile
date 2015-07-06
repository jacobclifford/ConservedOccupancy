CC = g++-4.4

GSL_DIR = /usr/local

INC = -I$(GSL_DIR)
CFLAGS = -O3 $(INC)/include 

LIB_DIRS = -L$(GSL_DIR)/lib
LIBS = -lgsl -lgslcblas -lm
LFLAGS = -O3 $(LIB_DIRS) $(LIBS)

BIN = plot
all: $(BIN)

clean:
	rm -f $(BIN)

plot1.o : plot1.cpp
	$(CC) $(CFLAGS) -c plot1.cpp


plot : plot1.o
	$(CC) $(LIB_DIRS) -o $@ plot1.o $(LIBS)


