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
                               const ulong iterations, const T tolerance, const ulong workers, std::ofstream &out) {

    start_time = Time::now();
    std::vector<float> old_solutions __attribute__((aligned(64)));
    std::vector<float> solutions __attribute__((aligned(64)));
    auto zero = tolerance - tolerance;
    auto error = zero;

    //vectorize the loop
    for (int i = 0; i < coefficients.size(); ++i) {
        old_solutions.emplace_back(zero);
        solutions.emplace_back(zero);
    }
    ff::ParallelFor pf(workers);

    init_time = Time::now();

    ulong iteration;
    for (iteration = 0; iteration < iterations; ++iteration) {
        error = zero;
        //calculate solutions using a parallel for
        pf.parallel_for_static(0, solutions.size(), 1, 0, [&](const ulong i) {
            solutions[i] = solution_find(coefficients[i], old_solutions, terms[i], i);
        }, workers);

#pragma simd
        for (ulong i = 0; i < solutions.size(); ++i)
            error += std::abs(solutions[i] - old_solutions[i]);

        //check the error and terminate in case
        if (error / solutions.size() <= tolerance) break;
        std::swap(solutions, old_solutions);
    }
    total_time = Time::now();
    print_metrics(iteration, error);
    write_csv(iteration, error, out);
    return solutions;
}

#endif //JACOBI_JACOBI_PAR_FOR_H
