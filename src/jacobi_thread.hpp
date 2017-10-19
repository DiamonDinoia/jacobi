//
// Created by marco on 9/5/17.
//

#ifndef JACOBI_JACOBI_THREAD_H
#define JACOBI_JACOBI_THREAD_H

#include "utils.hpp"
#include <vector>

/**
 * Parallel implementation using a manual multi-thread implementation
 * @tparam T template, type of the values
 * @param coefficients coefficients matrix of the linear system
 * @param terms right vector of the system
 * @param iterations max number of iteration
 * @param tolerance error tolerated
 * @return solution vector
 * @param workers number of threads
 * @return solution vector
 */
std::vector<float> thread_jacobi(const std::vector<std::vector<float>> &_coefficients, const std::vector<float> &_terms,
                                 ulong _iterations, float _tolerance, ulong nWorkers, std::ofstream &out);


#endif //JACOBI_JACOBI_THREAD_H
