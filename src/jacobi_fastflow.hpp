//
// Created by marco on 9/5/17.
//

#ifndef JACOBI_JACOBI_PAR_FOR_H
#define JACOBI_JACOBI_PAR_FOR_H

#include <ff/parallel_for.hpp>
#include <vector>
#include "utils.hpp"
#include <iostream>


/**
 * Parallel implementation using the fastflow parallel for skeleton
 * @tparam T template, type of the values
 * @param coefficients coefficients matrix of the linear system
 * @param terms right vector of the system
 * @param iterations max number of iteration
 * @param tolerance error tolerated
 * @return solution vector
 * @param workers number of threads
 * @return solution vector
 */
template<typename T>
std::vector<T> fastflow_jacobi(const std::vector<std::vector<T>> coefficients, const std::vector<T> terms,
                               const ulong iterations, const T tolerance, const ulong workers) {


    std::vector<float> old_solutions __attribute__((aligned(64)));
    std::vector<float> solutions __attribute__((aligned(64)));
    //vectorize the loop
    for (int i = 0; i < coefficients.size(); ++i) {
        old_solutions.emplace_back(tolerance - tolerance);
        solutions.emplace_back(tolerance - tolerance);
    }
    auto error = tolerance - tolerance;
    ff::ParallelFor pf(workers);
    auto zero = tolerance - tolerance;

    for (ulong iteration = 0; iteration < iterations; ++iteration) {
        error = zero;
        //calculate solutions using a parallel for
        pf.parallel_for_static(0, solutions.size(), 1, 0, [&](const ulong i) {
            solutions[i] = solution_find(coefficients[i], old_solutions, terms[i], i);
        });

#pragma simd
        for (ulong i = 0; i < solutions.size(); ++i)
            error += std::abs(solutions[i] - old_solutions[i]);

        //check the error and terminate in case
        if (error / solutions.size() <= tolerance) break;
        std::swap(solutions, old_solutions);
    }
    return solutions;
}

#endif //JACOBI_JACOBI_PAR_FOR_H
