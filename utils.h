//
// Created by marco on 9/2/17.
//

#ifndef JACOBI_UTILS_H
#define JACOBI_UTILS_H

template<typename T>
inline T max(T a, T b) { return (a > b) ? a : b; }

template<typename T>
inline T abs(T a) { return (a < 0.) ? -a : a; }

template<typename T>
inline T solution_find(const std::vector<T> row, const std::vector<T> solutions, T term, const ulong index) {
    for (int j = 0; j < row.size(); ++j) {
        if (j == index) continue;
        term -= (solutions[j] * row[j]);
    }
    return term / row[index];
}


#endif //JACOBI_UTILS_H
