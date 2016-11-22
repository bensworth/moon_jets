# Makefile for Enceladus simulations

# unix-tools
RM = rm -f

LIB = -L/usr/local/lib
INC = -I/usr/local/include

# compiler
#CC = g++ -std=c++11 -fopenmp -Wall -lgsl -lgslcblas -Ofast
CC = g++ -std=c++11 -Wall -lgsl -lgslcblas -Ofast
CP = gcc 

SRCS  = genFunctions.cpp Solver.cpp Jet.cpp 
SRCS1 = EncSimulation.cpp
SRCS2 = EurSimulation.cpp

OBJS  = ${SRCS:.cpp=.o}
OBJS1 = ${SRCS1:.cpp=.o}
EXEC1 = EncSim
OBJS2 = ${SRCS2:.cpp=.o}
EXEC2 = EurSim

INCLUDES  = 
CFLAGS    = $(INCLUDES) -g3
   
#.SUFFIXES :
.SUFFIXES : .o .cpp .h

all: $(EXEC1) $(EXEC2)

.h.o:
	$(CC) $(INCLUDES) $(CFLAGS) -c $<

.cpp.o:
	$(CC) $(INCLUDES) $(CFLAGS) -c $<

$(EXEC1): $(OBJS) $(OBJS1) 
	$(CC) $(OBJS) $(OBJS1) $(CFLAGS) -o $(EXEC1)

$(EXEC2): $(OBJS) $(OBJS2)
	$(CC) $(OBJS) $(OBJS2) $(CFLAGS) -o $(EXEC2)

Jet.o: genFunctions.o
EncSimulation.o: Jet.o Solver.o
EurSimulation.o: Jet.o Solver.o

clean: 
	$(RM) $(OBJS) $(OBJS1) $(OBJS2) $(EXEC1) $(EXEC2)
	make
	
install:
