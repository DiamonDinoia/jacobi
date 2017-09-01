//
// Created by marco on 9/1/17.
//

#ifndef JACOBI_JACOBI_H
#define JACOBI_JACOBI_H

#include <vector>
#include <cmath>

/*

n = no of equations
a[n][n] = coefficient matrix
b[n]= right hand side vector

x[n] - solution vector

*/

template<typename T>
inline T max(T a, T b) { return (a > b) ? a : b; }

template<typename T>
inline T solution_find(const std::vector<T> row, const std::vector<T> solutions, T term, const ulong index) {
    for (int j = 0; j < row.size(); ++j) {
        if (j == index) continue;
        term -= solutions[j] * row[j];
    }
    return term / row[index];
}

template<typename T>
std::vector<T> serial_jacobi(const std::vector<std::vector<T>> coefficients, const std::vector<T> terms,
                             const ulong iterations, const T tolerance) {

    std::vector<T> old_solutions(coefficients.size(), (tolerance - tolerance));
    std::vector<T> solutions(coefficients.size(), (tolerance - tolerance));
    T error;
    for (ulong i = 0; i < solutions.size(); ++i)
        solutions[i] = solution_find(coefficients[i], old_solutions, terms[i], i);

    for (int iteration = 1; iteration < iterations; ++iteration) {
        for (ulong i = 0; i < solutions.size(); ++i)
            solutions[i] = solution_find(coefficients[i], old_solutions, terms[i], i);

        for (auto &sol: solutions) {
            std::cout << sol << " ";
        }
        std::cout << std::endl;

        error = abs(solutions[0] - old_solutions[0]);
        for (int i = 1; i < solutions.size(); ++i) {
            error = max<T>(abs(solutions[i] - old_solutions[i]), error);
            old_solutions[i] = solutions[i];
        }
        if (error <= tolerance) return solutions;
    }

    return solutions;
}

#endif //JACOBI_JACOBI_H
