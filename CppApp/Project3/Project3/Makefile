# WILL CLEAN UP MAKEFILE LATER
CC = g++
CFLAGS = -std=c++11 -Wall -g -DNDEBUG -Ofast -march=native 
CXXFLAGS = -I$(MKLROOT)/include 
#-DEIGEN_USE_MKL_ALL -DCMAKE_BUILD_TYPE=Release -DEIGEN_TEST_SSE2=ON 
# LDLIBS += -L$(MKLROOT)/lib
# LDLIBS += -lmkl_intel_lp64 -lmkl_core -lmkl_intel_thread -liomp5 -lpthread -lm -ldl -m64 
OTHER_FLAGS = $(CXXFLAGS) $(LDLIBS) -framework Accelerate -mfma -ffp-contract=fast -ffast-math -flto -D_GLIBCXX_PARALLEL -fno-math-errno -I/Applications/MATLAB_R2020b.app/extern/include -L/Applications/MATLAB_R2020b.app/bin/maci64 -lmat -lmx -lmex -Xpreprocessor -fopenmp -lomp -I/Users/dinula/Documents/repos/ADMM-Code/LBFGSpp/include -I/usr/local/include/eigen3/
#OTHER_FLAGS = -O3 -march=native -ffp-contract=fast -I/Users/dinula/Documents/repos/optim/header_only_version -I/usr/local/include/eigen3/ -I/Applications/MATLAB_R2020b.app/extern/include -L/Applications/MATLAB_R2020b.app/bin/maci64 -lmat -lmx -lmex -Xpreprocessor -fopenmp -lomp
# DINTEL_MKL -DINTEL_MKL_ML -DEIGEN_USE_MKL_ALL -DMKL_DIRECT_CALL -mtune=native -L$(MKLROOT)opt/intel/compilers_and_libraries/mac/mkl/lib -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5 -lpthread -lm -ldl
default: main

test: optimizeexample.o admm_solver.o Triangle_Energy_Term.o
	$(CC) $(CFLAGS) -o test optimizeexample.o admm_solver.o Triangle_Energy_Term.o $(OTHER_FLAGS)

main: main.o admm_solver.o Triangle_Energy_Term.o
	$(CC) $(CFLAGS) -o main main.o admm_solver.o Triangle_Energy_Term.o $(OTHER_FLAGS)

optimizeexample.o: optimizeexample.cpp
	$(CC) $(CFLAGS) -c optimizeexample.cpp $(OTHER_FLAGS)

main.o: main.cpp
	$(CC) $(CFLAGS) -c main.cpp $(OTHER_FLAGS)

admm_solver.o: admm_solver.cpp admm_solver.h
	$(CC) $(CFLAGS) -c admm_solver.cpp $(OTHER_FLAGS)

Triangle_Energy_Term.o: Triangle_Energy_Term.cpp Triangle_Energy_Term.h
	$(CC) $(CFLAGS) -c Triangle_Energy_Term.cpp $(OTHER_FLAGS)

clean: 
	rm test main optimizeexample.o main.o admm_solver.o Triangle_Energy_Term.o
