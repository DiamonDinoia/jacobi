//
// Created by marco on 9/2/17.
//

#ifndef JACOBI_UTILS_H
#define JACOBI_UTILS_H

#include <vector>
#include <chrono>
#include <atomic>

typedef unsigned long ulong;

template<typename T>
inline T max(T a, T b) { return (a > b) ? a : b; }

template<typename T>
inline T abs(T a) { return (a < 0.) ? -a : a; }

template<typename T>
inline T solution_find(const std::vector<T> row, const std::vector<T> solutions, T term, const ulong index) {
#pragma ivdep
    for (int j = 0; j < row.size(); ++j) {
        if (j == index) continue;
        term -= (solutions[j] * row[j]);
    }
    return term / row[index];
}


typedef std::chrono::high_resolution_clock Time;
typedef std::chrono::milliseconds ms;
typedef std::chrono::duration<double> dsec;


template<typename T>
inline T my_atomic_fetch_add(std::atomic<T> &obj, T arg) {
    T expected = obj.load();
    while (!std::atomic_compare_exchange_weak(&obj, &expected, expected + arg));
    return expected;
}


#endif //JACOBI_UTILS_H