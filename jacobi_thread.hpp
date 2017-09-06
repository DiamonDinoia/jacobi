//
// Created by marco on 9/5/17.
//

#ifndef JACOBI_JACOBI_THREAD_H
#define JACOBI_JACOBI_THREAD_H

#include "utils.hpp"
#include <vector>

std::vector<float> jacobi_thread(const std::vector<std::vector<float>> &_coefficients, const std::vector<float> &_terms,
                                 const ulong _iterations, const float _tolerance, const ulong nWorkers);


#endif //JACOBI_JACOBI_THREAD_H
