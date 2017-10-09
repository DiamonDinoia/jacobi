//
// Created by marco on 9/2/17.
//

#ifndef JACOBI_UTILS_H
#define JACOBI_UTILS_H

#include <vector>
#include <chrono>
#include <string>
#include <algorithm>


typedef unsigned long ulong;

typedef std::chrono::high_resolution_clock Time;
typedef std::chrono::duration<double> dsec;


auto static start_time = Time::now();
auto static init_time = Time::now();
auto static total_time = Time::now();

const auto iterations_computed = "iterations_computed: ";
const auto error_s = "error: ";
const auto computation_time_s = "computation_time: ";
const auto initi_time_s = "initialization_time: ";


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
void print_solution(const std::vector<T> &solution) {
    for (auto &sol: solution) {
        std::cout << sol << " ";
    }
    std::cout << std::endl;
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


template<typename T>
void print_metrics(const ulong iterations, const T error) {
    std::cout << iterations_computed << iterations + 1 << std::endl;
    std::cout << error_s << error << std::endl;
    std::cout << initi_time_s << dsec(init_time - start_time).count() << std::endl;
    std::cout << computation_time_s << dsec(total_time - init_time).count() << std::endl;
}

#endif //JACOBI_UTILS_H
