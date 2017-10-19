//
// Created by marco on 9/15/17.
//

#ifndef JACOBI_JACOBI_OMP_HPP
#define JACOBI_JACOBI_OMP_HPP

#include <omp.h>
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
std::vector<T> omp_jacobi(const std::vector<std::vector<T>> coefficients, const std::vector<T> terms,
                          const ulong iterations, const T tolerance, const ulong workers, std::ofstream &out) {

    start_time = Time::now();

    std::vector<T> old_solutions __attribute__((aligned(64)));
    std::vector<T> solutions __attribute__((aligned(64)));

    omp_set_num_threads((int) workers);

    for (int i = 0; i < coefficients.size(); ++i) {
        old_solutions.emplace_back(tolerance - tolerance);
        solutions.emplace_back(tolerance - tolerance);
    }
    T error;
    init_time = Time::now();
    ulong iteration;
    for (iteration = 0; iteration < iterations; ++iteration) {

        //calculate solutions
        error = tolerance - tolerance;

#pragma omp parallel for schedule(static)
        for (ulong i = 0; i < solutions.size(); ++i) {
            solutions[i] = solution_find(coefficients[i], old_solutions, terms[i], i);
        }
        //compute the error
        for (ulong i = 0; i < solutions.size(); ++i) {
            error += std::abs(solutions[i] - old_solutions[i]);
        }
        // check the error
        error /= solutions.size();
        if (error <= tolerance) break;
    }
    total_time = Time::now();
    std::swap(solutions, old_solutions);

    print_metrics(iteration, error);
    write_csv(iteration, error, out);

    return solutions;
}

#endif //JACOBI_JACOBI_OMP_H
