//#include "jacobi_ff.h"

#include <iostream>

#include "jacobi.h"
#include "jacobi_map.h"
#include "jacobi_par_for.h"
#include "jacobi_thread.h"

using namespace std;

int main() {

    const ulong size = 4;

    auto matrix = new float *[4];
    matrix[0] = new float[4]{10., -1., 2., 0.};
    matrix[1] = new float[4]{-1., 11., -1., 3.};
    matrix[2] = new float[4]{2., -1., 10., -1.};
    matrix[3] = new float[4]{0.0, 3., -1., 8.};

    auto term = new float[4]{6., 25., -11., 15.};


    vector<vector<float>> coefficients(4, vector<float>(4, 0.));
    for (int i = 0; i < 4; ++i) {
        for (int j = 0; j < 4; ++j) {
            coefficients[i][j] = matrix[i][j];
        }
    }
    vector<float> terms;

    for (int i = 0; i < size; ++i) {
        terms.emplace_back(term[i]);
    }

    vector<float> solutions = serial_jacobi(coefficients, terms, 1000, 0.0f);

    cout << "solution: ";
    for (auto &sol: solutions) {
        cout << sol << " ";
    }
    cout << endl;

    float *solution2 = jacobi_map(matrix, term, size, 1000, 0.0, 8);

    cout << "solution2: ";
    for (int i = 0; i < size; ++i) {
        cout << solution2[i] << " ";

    }
    cout << endl;
    vector<float> solution3 = jacobi_par_for(coefficients, terms, 1000, 0.f, 8);

    cout << "solution3: ";
    for (int i = 0; i < size; ++i) {
        cout << solution3[i] << " ";

    }
    cout << endl;


    auto solution4 = jacobi_thread(coefficients, terms, 1000, 0.f, 8);
    cout << "solution4: ";
    for (auto &sol: solution4) {
        cout << sol << " ";
    }
    cout << endl;
    return 0;
}