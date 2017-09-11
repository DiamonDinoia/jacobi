ffdir=~/fastflow
cores=8
iterations=500
tolerance=-1.f

CXX=g++
CXXFLAGS= -O3 -pthread -std=c++14 -I $(ffdir)

.PHONY: main clean cleanall test benchmark


main: main.cpp jacobi.hpp utils.hpp jacobi_map.hpp jacobi_par_for.hpp \
jacobi_thread.hpp jacobi_thread.o jacobi_map.o
	$(CXX) $(CXXFLAGS) $? -o $@
	
jacobi_thread.o: jacobi_thread.cpp jacobi_thread.hpp utils.hpp \
 barrier.hpp
	$(CXX) $(CXXFLAGS) -c $<

jacobi_map.o: jacobi_map.cpp jacobi_map.hpp utils.hpp
	$(CXX) $(CXXFLAGS) -c $<
	
clean:
	rm -f *.o

cleanall:
	rm -f *.o
	rm -f main
	rm -f *txt
	rm -f *png

test:
	python generate_matrix.py
	make main
	./main $(cores) $(iterations) $(tolerance) matrix0.txt matrix1.txt matrix2.txt matrix3.txt > results.txt
	python parse_results.py -i results.txt -o metrics.txt

benchmark:
	python generate_matrix.py
	make main
	./main $(cores) $(iterations) $(tolerance) matrix0.txt matrix1.txt matrix2.txt matrix3.txt matrix4.txt matrix5.txt > results.txt
	python parse_results.py -i results.txt -o metrics.txt
	
	
