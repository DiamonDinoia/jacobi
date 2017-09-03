#include <iostream>
#include "jacobi.h"
#include "jacobi_ff.h"

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
    for (double &i : term) {
        terms.emplace_back(i);
    }
    vector<double> solutions = serial_jacobi(coefficients, terms, 1000, 0.0);
    vector<double> solution2 = par_for_jacobi(coefficients, terms, 1000, 0.0, 8);

    cout << "solution: ";
    for (auto &sol: solutions) {
        cout << sol << " ";
    }
    cout << endl;

    cout << "solution2: ";
    for (auto &sol: solution2) {
        cout << sol << " ";
    }
    cout << endl;

    return 0;
}