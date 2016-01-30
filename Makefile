CXX = icc
CXXFLAGS = -g -std=c99 -qopenmp -O3 -fp-model fast=2
NOPAR = -no-vec -no-simd

CPUFLAGS = $(CXXFLAGS) -xhost
MICFLAGS = $(CXXFLAGS) -mmic
OPTFLAGS = -qopt-report -qopt-report-file=$@.optrpt

OBJECTS = cpu_serial.o

.SUFFIXES: .o .c

.c.o:
	$(CXX) -c $(CPUFLAGS) $(OPTFLAGS) -o "$@" "$<"

all: cpu_bootstrap cpu_bootstrap2

cpu_bootstrap: $(OBJECTS)
	$(CXX) $(CXXFLAGS) -o cpu_bootstrap $(OBJECTS)


cpu_bootstrap2: $(OBJECTS)
	$(CXX) $(CXXFLAGS) $(NOPAR) -o cpu_bootstrap_no $(OBJECTS)


clean: 
	rm -f *.o cpu_bootstrap
