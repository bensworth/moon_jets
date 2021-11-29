# Makefile for Enceladus simulations

# unix-tools
RM = rm -f

LIB = -L/usr/local/lib
INC = -I/usr/local/include
LDFLAGS = -lgsl -lgslcblas
CFLAGS  = $(INC) $(LIB) -g3 -fopenmp -Ofast

# compiler
HDF5_CXX = /usr/local/bin/g++
CXX = h5c++ -std=c++11

SRCS  = genFunctions.cpp Solver.cpp Jet.cpp 
SRCS1 = LowAlttiudeSim.cpp

OBJS  = ${SRCS:.cpp=.o}
OBJS1 = ${SRCS1:.cpp=.o}
EXEC1 = LowAlttiudeSim

#.SUFFIXES :
.SUFFIXES : .o .cpp .h

all: $(EXEC1)

.h.o:
	$(CXX) $(CFLAGS) -c $<

.cpp.o:
	$(CXX) $(CFLAGS) -c $<

$(EXEC1): $(OBJS) $(OBJS1) 
	$(CXX) $(LDFLAGS) $(OBJS) $(OBJS1) $(CFLAGS) -o $(EXEC1)

	
Jet.o: genFunctions.o
LowAlttiudeSim.o: Jet.o Solver.o

clean: 
	$(RM) $(OBJS) $(OBJS1) $(EXEC1)
	
install:
