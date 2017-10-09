//
// Created by marco on 9/2/17.
//

#ifndef JACOBI_UTILS_H
#define JACOBI_UTILS_H

#include <vector>
#include <chrono>
#include <string>
#include <algorithm>

namespace util {

    typedef unsigned long ulong;

    typedef std::chrono::high_resolution_clock Time;
    typedef std::chrono::duration<double> dsec;

    const std::string separator = " | ";
    const std::string iterations_computed = "iterations computed: ";
    const std::string error_s = "error: ";
    const std::string computation_time_s = "computation time: ";

    template<typename T>
    void printmetrics(std::string name, T start_time, T init_time, T total_time) {

    }

    template<typename T>
    void printmetrics(std::string name, T start_time, T init_time, T total_time, T comp_time, T error_time) {

    }


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
            term -= (solutions[j] * row[j]);
        }
        return (term + (solutions[index] * row[index])) / row[index];
    }

    template<typename T>
    T rand_t(T min, T max) {
        T range = max - min;
        T r = (T) rand() / (T) RAND_MAX;
        return (r * range) + min;
    }

    template<typename T>
    void generate_vector(const ulong size, std::vector<T> &v, const T min, const T max) {
        for (int i = 0; i < size; ++i) {
            v.emplace_back(rand_t(min, max));
        }
    }


    template<typename T>
    void generate_diagonal_dominant_matrix(const ulong size, std::vector<std::vector<T>> &matrix,
                                           const T min, const T max) {
        for (ulong i = 0; i < size; ++i) {
            std::vector<T> tmp;
            generate_vector(size, tmp, min, max);
            matrix.emplace_back(tmp);
        }
        for (ulong i = 0; i < size; ++i) {
            T sum = T(0);
            for (auto &val: matrix[i]) sum += abs(val);
            sum -= abs(matrix[i][i]);
            matrix[i][i] = abs(matrix[i][i]) + sum;
        }
    }
}

#endif //JACOBI_UTILS_H
