//
// Created by marco on 9/5/17.
//

#ifndef JACOBI_JACOBI_PAR_FOR_H
#define JACOBI_JACOBI_PAR_FOR_H

#include <ff/parallel_for.hpp>
#include <vector>
#include "utils.hpp"
#include <iostream>


template<typename T>
std::vector<T> jacobi_par_for(const std::vector<std::vector<T>> coefficients, const std::vector<T> terms,
                              const ulong iterations, const T tolerance, const ulong workers) {

    std::vector<T> old_solutions(coefficients.size(), (tolerance - tolerance));
    std::vector<T> solutions(coefficients.size(), (tolerance - tolerance));
    auto error = tolerance - tolerance;
    ff::ParallelFor pf;

//    auto reduce = [&](float &var, const ulong i) {
//        var += abs(solutions[i] - old_solutions[i]);
//    };
//    auto reduce2 = [&](const ulong i, float &var) {
//        var += abs(solutions[i] - old_solutions[i]);
//    };

    std::atomic_flag flag;
    flag.clear();

    auto start = Time::now();
    for (ulong iteration = 0; iteration < iterations; ++iteration) {
        //calculate solutions
        error = tolerance - tolerance;
        pf.parallel_for(0, solutions.size(), [&](const ulong i) {
            solutions[i] = solution_find(coefficients[i], old_solutions, terms[i], i);
//            my_atomic_fetch_add<float>(error, abs(solutions[i] - old_solutions[i]));
            auto val = abs(solutions[i] - old_solutions[i]);
            old_solutions[i] = solutions[i];
            while (!flag.test_and_set(std::memory_order_relaxed)) {}
            error += val;
            flag.clear();
        }, workers);
//            Does not works, generates a segmentation fault
//            pf.parallel_reduce(error, 0.f, 0, solutions.size(), reduce2, reduce, workers);

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
