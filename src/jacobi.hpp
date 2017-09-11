//
// Created by marco on 9/1/17.
//

#ifndef JACOBI_JACOBI_H
#define JACOBI_JACOBI_H

#include <vector>
#include "utils.hpp"
#include <iostream>

/**
 * Serial implementation of the Jacobi method simply iterates until reach the convergence or reach the max number of
 * iterations
 * @tparam T template, type of the values
 * @param coefficients coefficients matrix of the linear system
 * @param terms right vector of the system
 * @param iterations max number of iteration
 * @param tolerance error tolerated
 * @return solution vector
 */
template<typename T>
std::vector<T> serial_jacobi(const std::vector<std::vector<T>> coefficients, const std::vector<T> terms,
                             const ulong iterations, const T tolerance) {

    std::vector<T> old_solutions __attribute__((aligned(64)));
    std::vector<T> solutions __attribute__((aligned(64)));
#pragma ivdep
    for (int i = 0; i < coefficients.size(); ++i) {
        old_solutions.emplace_back(tolerance - tolerance);
        solutions.emplace_back(tolerance - tolerance);
    }
    T error;
    auto start = Time::now();
    for (ulong iteration = 0; iteration < iterations; ++iteration) {
        //calculate solutions
        error = tolerance - tolerance;
#pragma ivdep
        for (ulong i = 0; i < solutions.size(); ++i) {
            solutions[i] = solution_find(coefficients[i], old_solutions, terms[i], i);
        }
        //compute the error
#pragma simd
        for (ulong i = 0; i < solutions.size(); ++i) {
            error += abs(solutions[i] - old_solutions[i]);
            old_solutions[i] = solutions[i];
        }
        // check the error
        error /= solutions.size();
        if (error <= tolerance) {
            std::cout << "serial jacobi | iterations computed: " << iteration << " error: " << error << std::endl;
            return solutions;
        }
    }
    auto end = Time::now();
    std::cout << "serial jacobi | iterations computed: " << iterations << " error: " << error << std::endl;
    std::cout << "serial jacobi | computation time: " << dsec(end - start).count() << std::endl;
    return solutions;
}

#endif //JACOBI_JACOBI_H