//
// Created by marco on 9/2/17.
//

#ifndef JACOBI_UTILS_H
#define JACOBI_UTILS_H

#include <vector>
#include <chrono>

typedef unsigned long ulong;

typedef std::chrono::high_resolution_clock Time;
typedef std::chrono::duration<double> dsec;


template<typename T>
inline T abs(T a) { return (a < 0.) ? -a : a; }

/**
 * utility function implementing the jacobi method in order to find one solution
 * @param row coeffiient row
 * @param solutions vector solution
 * @param term right term vector
 * @param index index of the solution
 * @return solution component
 */
template<typename T>
inline  __attribute__((always_inline))
T solution_find(const std::vector<T> row, const std::vector<T> solutions, T term, const ulong index) {
#pragma simd
    for (int j = 0; j < row.size(); ++j) {
        if (j == index) continue;
        term -= (solutions[j] * row[j]);
    }
    return term / row[index];
}


#endif //JACOBI_UTILS_H