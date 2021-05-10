# Makefile for Enceladus simulations

# unix-tools
RM = rm -f

LIB = -L/usr/local/lib
INC = -I/usr/local/include
LDFLAGS = -lgsl -lgslcblas
CFLAGS  = $(INC) $(LIB) -g3 -fopenmp -Wall -Ofast

# compiler
CXX = mpicxx -std=c++11

SRCS  = genFunctions.cpp Solver.cpp Jet.cpp 
SRCS1 = EncParallel.cpp
SRCS2 = EurParallel.cpp

OBJS  = ${SRCS:.cpp=.o}
OBJS1 = ${SRCS1:.cpp=.o}
EXEC1 = EncPar
OBJS2 = ${SRCS2:.cpp=.o}
EXEC2 = EurPar

#.SUFFIXES :
.SUFFIXES : .o .cpp .h

all: $(EXEC1) $(EXEC2)

.h.o:
	$(CXX) $(CFLAGS) -c $<

.cpp.o:
	$(CXX) $(CFLAGS) -c $<

$(EXEC1): $(OBJS) $(OBJS1) 
	$(CXX) $(LDFLAGS) $(OBJS) $(OBJS1) $(CFLAGS) -o $(EXEC1)

$(EXEC2): $(OBJS) $(OBJS2)
	$(CXX) $(LDFLAGS) $(OBJS) $(OBJS2) $(CFLAGS) -o $(EXEC2)
	
Jet.o: genFunctions.o
EncParallel.o: Jet.o Solver.o
EurParallel.o: Jet.o Solver.o

clean: 
	$(RM) $(OBJS) $(OBJS1) $(OBJS2) $(EXEC1) $(EXEC2)
	make -f parallel
	
install:
