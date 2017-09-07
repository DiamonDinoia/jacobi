//
// Created by marco on 9/1/17.
//

#ifndef JACOBI_JACOBI_H
#define JACOBI_JACOBI_H

#include <vector>
#include "utils.hpp"
#include <iostream>

template<typename T>
std::vector<T> serial_jacobi(const std::vector<std::vector<T>> coefficients, const std::vector<T> terms,
                             const ulong iterations, const T tolerance) {

    std::vector<T> old_solutions(coefficients.size(), (tolerance - tolerance));
    std::vector<T> solutions(coefficients.size(), (tolerance - tolerance));
    T error;

    for (ulong iteration = 0; iteration < iterations; ++iteration) {
        //calculate solutions
#pragma ivdep
        for (ulong i = 0; i < solutions.size(); ++i)
            solutions[i] = solution_find(coefficients[i], old_solutions, terms[i], i);

        // check the error
        error = abs(solutions[0] - old_solutions[0]);
        old_solutions[0] = solutions[0];
#pragma ivdep
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

#endif //JACOBI_JACOBI_H