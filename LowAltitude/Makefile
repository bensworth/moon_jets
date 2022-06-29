# Makefile for Enceladus simulations

# unix-tools
RM = rm -f

#LIB = -L/usr/local/lib
#INC = -I/usr/local/include
### GSL is installed in -L/home/southworth/c_libs/install/
### -Wl,-rpath=/home/southworth/c_libs/install/lib is necessary
### to not only comile against GSL, but link the library during runtime
LIB = -L/home/southworth/c_libs/hdf5/hdf5/lib -L/home/southworth/c_libs/install/lib -Wl,-rpath=/home/southworth/c_libs/install/lib
INC = -I/home/southworth/c_libs/hdf5/hdf5/include -I/home/southworth/c_libs/install/include
LDFLAGS = -lgsl -lgslcblas
CFLAGS  = $(INC) $(LIB) -g3 -fopenmp -Ofast

# compiler
# HDF5_CXX = /usr/local/bin/g++
# CXX = g++ -std=c++14
CXX = /home/southworth/c_libs/hdf5/hdf5/bin/h5c++ -std=c++14

SRCS  = genFunctions.cpp Solver.cpp Jet.cpp 
SRCS1 = LowAlttiudeSim.cpp
SRCS2 = AlttiudeDensitySim.cpp

OBJS  = ${SRCS:.cpp=.o}
OBJS1 = ${SRCS1:.cpp=.o}
OBJS2 = ${SRCS2:.cpp=.o}
EXEC1 = LowAlttiudeSim
EXEC2 = AlttiudeDensitySim

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
LowAlttiudeSim.o: Jet.o Solver.o

clean: 
	$(RM) $(OBJS) $(OBJS1) $(EXEC1)
	
install:
