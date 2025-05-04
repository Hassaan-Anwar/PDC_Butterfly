# Makefile for Butterfly Algorithm implementations

# Compiler settings
CXX = g++
MPICXX = mpic++

# Compiler flags
CXXFLAGS = -Wall -O3
OMPFLAGS = -fopenmp
MPIFLAGS = -lmetis

# Targets
all: butterfly quick mpi_butterfly hybrid

# Normal Sequential version
butterfly: attempt1.cpp
	$(CXX) -Wall -O2 -o butterfly attempt1.cpp

# OpenMP only version
quick: openmp.cpp
	$(CXX) -std=c++11 $(OMPFLAGS) -O3 openmp.cpp -o quick

# OpenMPI only version
mpi_butterfly: openmpi.cpp
	$(MPICXX) -std=c++14 $(CXXFLAGS) $(OMPFLAGS) -o mpi_butterfly openmpi.cpp $(MPIFLAGS)

# Hybrid OpenMPI and OpenMP version
hybrid: openmpi_openmp.cpp
	$(MPICXX) $(OMPFLAGS) $(CXXFLAGS) -std=c++17 -o hybrid openmpi_openmp.cpp $(MPIFLAGS)

# Run targets
run-sequential:
	./butterfly

run-openmp:
	./quick

run-mpi:
	mpiexec -n 4 ./mpi_butterfly

run-hybrid:
	mpiexec -n 4 ./hybrid

# Clean up
clean:
	rm -f butterfly quick mpi_butterfly hybrid

.PHONY: all clean run-sequential run-openmp run-mpi run-hybrid