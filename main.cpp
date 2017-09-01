#include <iostream>
#include "jacobi.h"

using namespace std;

int main() {

    double matrix[4][4] = {
            {10., -1., 2.,  0.},
            {-1., 11., -1., 3.},
            {2.,  -1., 10., -1.},
            {0.0, 3.,  -1., 8.}
    };
    double term[] = {6., 25., -11., 15.};

    vector<vector<double>> coefficients(4, vector<double>(4, 0.));
    for (int i = 0; i < 4; ++i) {
        for (int j = 0; j < 4; ++j) {
            coefficients[i][j] = matrix[i][j];
        }
    }
    vector<double> terms;
    for (int i = 0; i < 4; ++i) {
        terms.emplace_back(term[i]);
    }
    vector<double> solutions = serial_jacobi(coefficients, terms, 1000, 0.0);

    cout << "solution: ";
    for (auto &sol: solutions) {
        cout << sol << " ";
    }
    cout << endl;

    return 0;
}