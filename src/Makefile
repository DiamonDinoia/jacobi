ffdir=~/fastflow
cores=64
iterations=100
tolerance=-1.f

CXX=icc
CXXFLAGS=-O3 -pthread -mmic -std=c++14 -openmp -openmp-link=static -I $(ffdir)

.PHONY: main clean cleanall test benchmark run run_extended test_no_omp benchmark_no_omp

main: main.cpp jacobi_thread.o
	$(CXX) $(CXXFLAGS) $? -o $@

jacobi_thread.o: jacobi_thread.cpp
	$(CXX) $(CXXFLAGS) -c $<

clean:
	rm -f *.o

cleanall:
	rm -f *.o
	rm -f main
	rm -f *txt
	rm -f *png

	
