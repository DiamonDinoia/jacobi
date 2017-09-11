//
// Created by marco on 9/5/17.
//

#ifndef JACOBI_JACOBI_MAP_H
#define JACOBI_JACOBI_MAP_H

#include "utils.hpp"

/**
 * Parallel implementation of the Jacobi's method using the fastflow Map skeleton
 * @param _coefficients matrix of coefficient of the linear system
 * @param _terms vector of right terms of the system
 * @param _size size of the system
 * @param _iterations max number of iterations
 * @param _tolerance error tolerated
 * @param _nworkers number of threads
 * @return solution vector
 */
float *jacobi_map(float **_coefficients, float *_terms, const ulong _size, const ulong _iterations,
                  const float _tolerance,
                  const ulong _nworkers);


#endif //JACOBI_JACOBI_MAP_H
