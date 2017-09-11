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
std::vector<T> jacobi_par_for(const std::vector<std::vector<T>> coefficients, const std::vector<T> terms,
                              const ulong iterations, const T tolerance, const ulong workers) {

    std::vector<float> old_solutions __attribute__((aligned(64)));
    std::vector<float> solutions __attribute__((aligned(64)));
    //vectorize the loop
#pragma ivdep
    for (int i = 0; i < coefficients.size(); ++i) {
        old_solutions.emplace_back(tolerance - tolerance);
        solutions.emplace_back(tolerance - tolerance);
    }
    auto error = tolerance - tolerance;
    ff::ParallelFor pf;

    // code used by the reduce not needed now
//    auto reduce = [&](float &var, const ulong i) {
//        var += abs(solutions[i] - old_solutions[i]);
//    };
//    auto reduce2 = [&](const ulong i, float &var) {
//        var += abs(solutions[i] - old_solutions[i]);
//    };
    // synchronization flag, used in order to use lock-free mechanisms
    std::atomic_flag flag;
    flag.clear();

    auto start = Time::now();

    for (ulong iteration = 0; iteration < iterations; ++iteration) {
        error = tolerance - tolerance;
        //calculate solutions using a parallel for
        pf.parallel_for(0, solutions.size(), [&](const ulong i) {
            solutions[i] = solution_find(coefficients[i], old_solutions, terms[i], i);
        }, workers);

        //compute the error using a parallel for
        pf.parallel_for(0, solutions.size(), [&](const ulong i) {
            auto val = abs(solutions[i] - old_solutions[i]);
            old_solutions[i] = solutions[i];
            //busy waiting for the lock (spin-lock)
            //relaxed memory order beacuse only atomicity is needed
            while (!flag.test_and_set(std::memory_order_relaxed)) {}
            error += val;
            flag.clear();
        }, workers);

//            Does not works, generates a segmentation fault
//            pf.parallel_reduce(error, 0.f, 0, solutions.size(), reduce2, reduce, workers);


        //check the error and terminate in case
        if (error / solutions.size() <= tolerance) {
            auto end = Time::now();
            std::cout << "parallel for | iterations computed: " << iteration << " error: " << error << std::endl;
            std::cout << "parallel for | computation time: " << dsec(end - start).count() << std::endl;
            return solutions;
        }

    }
    auto end = Time::now();
    std::cout << "parallel for | iterations computed: " << iterations << " error: " << error << std::endl;
    std::cout << "parallel for | computation time: " << dsec(end - start).count() << std::endl;
    return solutions;
}

#endif //JACOBI_JACOBI_PAR_FOR_H
