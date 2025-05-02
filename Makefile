CXX = mpic++
CXXFLAGS = -std=c++14 -Wall -O3 -fopenmp
LDFLAGS = -lmetis

all: mpi_butterfly

mpi_butterfly: openmpi.cpp
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LDFLAGS)

clean:
	rm -f mpi_butterfly

.PHONY: all clean