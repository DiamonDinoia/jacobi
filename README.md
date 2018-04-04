# Jacobi 
This is a parallel implementation of the Jacobi algorithm.
## Jacobi method 
 The **Jacobi method** (or **Jacobi iterative method** is an algorithm for determining the solutions of a [diagonally dominant](https://en.wikipedia.org/wiki/Diagonally_dominant_matrix "Diagonally dominant matrix") [system of linear equations](https://en.wikipedia.org/wiki/System_of_linear_equations "System of linear equations"). Each diagonal element is solved for, and an approximate value is plugged in.
## Implementation
In the source code is possible to find different implementation using different parallel framework. The framework supported are:  
  1. OpenMP  
  2. FastFlow  
  3. Low level thread implementation  

## Usage
Just include the source in your code allocate the data structures and call the correct solving function.
## Tests
It is also included a main containing several tests.

### Usage  
**main** <algorithm> [-w <workers>] [-s <size>] « [-p <filename>][-i <iterations>] [-t <tolerance> [-h] [-d] [-c <seed>] » 

The required argument is algorithm that indicates the algorithm executed from the following list.  
  * **sequential**: sequential jacobi algorithm.  
  * **omp**: OpenMP multi-thread implementation of jacobi algorithm.   
  * **thread:** plain thread implementation of jacobi algorithm.   
  * **fastflow**: fastflowimplementation of jacobi algorithm.  

The optional arguments are:  
  * **[-w]** number of threads used, default 8.  
  * **[-s]** size of the matrix, default 1024.  
  * **[-i]** iteration performed, default 50.  
  * **[-t]** error tolerated, default -1.  
  * **[-p]** filename in case of csv exporting, default null.   
  * **[-h]** prints the helper.  
  * **[-d]** enable debug prints, solution and error.  
  * **[-c]** seed used to generate the matrix, default. 42.  

## Future works
Currently I'm developing an high performance libraryin c++ called [ParallelIterativeMethods](https://github.com/DiamonDinoia/parallelIterativeMethods) that achieve better performance.
