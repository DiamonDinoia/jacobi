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

    for (ulong iteration = 0; iteration < iterations; ++iteration) {
        //calculate solutions
        pf.parallel_for(0, solutions.size(), [&](const ulong i) {
            solutions[i] = solution_find(coefficients[i], old_solutions, terms[i], i);
        }, workers);

        // check the error
        error = abs(solutions[0] - old_solutions[0]);
        old_solutions[0] = solutions[0];
        for (int i = 1; i < solutions.size(); ++i) {
            error = max(abs(solutions[i] - old_solutions[i]), error);
            old_solutions[i] = solutions[i];
        }
        if (error <= tolerance) {
            std::cout << "iteration computed: " << iteration << " error: " << error << std::endl;
            return solutions;
        }
    }

    std::cout << "iteration computed: " << iterations << " error: " << error << std::endl;
    return solutions;
}

#endif //JACOBI_JACOBI_PAR_FOR_H
