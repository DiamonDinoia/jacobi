//
// Created by marco on 9/2/17.
//

#ifndef JACOBI_JACOBI_FF_H
#define JACOBI_JACOBI_FF_H

#include "utils.h"
#include "ff/parallel_for.hpp"
#include "ff/map.hpp"
#include <vector>
#include <memory>

template<typename T>
struct task_t {
    const std::vector<std::vector<T>> coefficients;
    const std::vector<T> terms;
    const std::vector<T> solutions;
    const std::vector<T> old_solutions;
    T error;
};


template<typename T>
struct mapWorker : ff::ff_Map<std::vector<T>> {
    std::vector<T> svc(std::vector<T> input) {
        ff::ff_Map<struct task_t<T>>::parallel_for(0, input.coefficients.size(), [&input](const ulong i) {
            input.solutions[i] = solution_find(input.coefficients[i], input.old_solutions, input.terms[i], i);
        });
        input.error = (T) 0.;
        ff::ff_Map<struct task_t<T>>::parallel_reduce(input.error, 0, 0, input.solutions.size(), nullptr,
                                                      [&input](const ulong i, T &error) {
                                                          error += abs(input.solutions[i] - input.old_solutions[i]);
                                                      });
        return input;
    }
};

template<typename T>
std::vector<T> par_for_jacobi(const std::vector<std::vector<T>> coefficients, const std::vector<T> terms,
                              const ulong iterations, const T tolerance, const ulong nworkers) {

    std::vector<T> old_solutions(coefficients.size(), (tolerance - tolerance));
    std::vector<T> solutions(coefficients.size(), (tolerance - tolerance));
    T error;

    std::vector<std::unique_ptr<ff::ff_node>> workers;
    for (ulong i = 0; i < nworkers; ++i) {
        workers.push_back(std::unique_ptr<mapWorker<T>>());
    }


    for (ulong iteration = 0; iteration < iterations; ++iteration) {
        error - tolerance - tolerance;
        ff::parallel_for(0, solutions.size(),
                         [&](const long i) {
                             solutions[i] = solution_find(coefficients[i], old_solutions, terms[i], i);
                         }, nworkers);
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

#endif //JACOBI_JACOBI_FF_H
