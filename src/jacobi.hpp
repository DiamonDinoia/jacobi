//
// Created by marco on 9/1/17.
//

#ifndef JACOBI_JACOBI_H
#define JACOBI_JACOBI_H

#include <vector>
#include "utils.hpp"
#include <iostream>

using namespace util;

namespace jacobi {

    static const std::string name = "serial jacobi";


    auto start_time = Time::now();
    auto init_time = Time::now();
#ifdef PROBE
    auto comp_time = Time::now();
    auto error_time = Time::now();
#endif
    auto total_time = Time::now();

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

        start_time = Time::now();
        //allocate solution vectors
        std::vector<T> old_solutions __attribute__((aligned(64)));
        std::vector<T> solutions __attribute__((aligned(64)));

        T zero = tolerance - tolerance;
        T error;

        //initialize solution vectors
        for (int i = 0; i < coefficients.size(); ++i) {
            old_solutions.emplace_back(zero);
            solutions.emplace_back(zero);
        }
        //Starting iterations
        init_time = Time::now();
        for (ulong iteration = 0; iteration < iterations; ++iteration) {
            std::swap(solutions, old_solutions);
            //calculate solutions
            error = zero;
            for (ulong i = 0; i < solutions.size(); ++i) {
                solutions[i] = solution_find(coefficients[i], old_solutions, terms[i], i);
            }
            //compute the error
#pragma simd
            for (ulong i = 0; i < solutions.size(); ++i) {
                error += abs(solutions[i] - old_solutions[i]);
            }

            // check the error
            error /= solutions.size();
            if (error <= tolerance) break;
        }
        total_time = Time::now();
        std::cout << name << separator << iterations_computed << iterations << ' ' << error_s << error << std::endl;
        std::cout << name << separator << computation_time_s << dsec(total_time - start_time).count() << std::endl;
        return solutions;
    }
}

#endif //JACOBI_JACOBI_H