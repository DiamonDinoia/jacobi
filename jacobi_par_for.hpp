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
    T error;
    ff::ParallelFor pf(workers);
//    auto reduce = [&](float &var, const ulong i) {
//        var += abs(solutions[i] - old_solutions[i]);
//    };
//    auto reduce2 = [&](const ulong i, float &var) {
//        var += abs(solutions[i] - old_solutions[i]);
//    };

    for (ulong iteration = 0; iteration < iterations; ++iteration) {
        //calculate solutions
        pf.parallel_for(0, solutions.size(), [&](const ulong i) {
            solutions[i] = solution_find(coefficients[i], old_solutions, terms[i], i);
        }, workers);

        // check the error
        error = 0.f;
#pragma ivdep
        for (int i = 0; i < solutions.size(); ++i) {
            error += abs(solutions[i] - old_solutions[i]);
            old_solutions[i] = solutions[i];
        }
//        again it generates a segmentation fault
//        pf.parallel_reduce(error, 0.f, 0, solutions.size(), reduce2, reduce, workers);

        if (error / solutions.size() <= tolerance) {
            std::cout << "iteration computed: " << iteration << " error: " << error << std::endl;
            return solutions;
        }

    }

    std::cout << "iteration computed: " << iterations << " error: " << error << std::endl;
    return solutions;
}

#endif //JACOBI_JACOBI_PAR_FOR_H
