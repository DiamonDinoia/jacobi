//
// Created by marco on 9/5/17.
//

#ifndef JACOBI_JACOBI_MAP_H
#define JACOBI_JACOBI_MAP_H

#include "utils.h"

float *jacobi_map(float **_coefficients, float *_terms, const ulong _size, const ulong _iterations,
                  const float _tolerance,
                  const ulong nworkers);


#endif //JACOBI_JACOBI_MAP_H
