ffdir=~/fastflow
cores=8
iterations=500
tolerance=-1.f

CXX=g++
CXXFLAGS=-Ofast -pthread -std=c++14 -I $(ffdir)
OMPFLAGS=-fopenmp -DWITHOMP

.PHONY: main clean cleanall test benchmark run run_extended test_no_omp benchmark_no_omp

main: main.cpp jacobi_thread.o jacobi_map.o
	$(CXX) $(CXXFLAGS) $? -o $@

jacobi_thread.o: thread_jacobi.cpp
	$(CXX) $(CXXFLAGS) -c $<

jacobi_map.o: jacobi_map.cpp
	$(CXX) $(CXXFLAGS) -c $<

clean:
	rm -f *.o

cleanall:
	rm -f *.o
	rm -f main
	rm -f *txt
	rm -f *png

run:
	python generate_matrix.py
	./main $(cores) $(iterations) $(tolerance) matrix0.txt matrix1.txt matrix2.txt matrix3.txt > results.txt
	python parse_results.py -i results.txt -o metrics.txt

run_extended:
	python generate_matrix.py
	./main $(cores) $(iterations) $(tolerance) matrix0.txt matrix1.txt matrix2.txt matrix3.txt matrix4.txt matrix5.txt > results.txt
	python parse_results.py -i results.txt -o metrics.txt

test_no_omp: main run

benchmark_no_omp: main run_extended

test: CXXFLAGS += $(OMPFLAGS)
test: main run

benchmark: CXXFLAGS += $(OMPFLAGS)
benchmark: main run_extended

	
	
