
#include <iostream>

#include "jacobi.hpp"
#include "jacobi_map.hpp"
#include "jacobi_par_for.hpp"
#include "jacobi_thread.hpp"
#include "jacobi_omp.hpp"

using namespace std;

const auto size_string = "size:";
const auto matrix_string = "matrix:";
const auto solution_string = "solution:";
const auto terms_string = "terms:";

ulong size;
vector<vector<float>> matrix __attribute__((aligned(64)));
vector<float> terms __attribute__((aligned(64)));
vector<float> solution __attribute__((aligned(64)));


void parse_input(const string &filename) {
    ifstream file(filename);
    string str;
    file >> str >> size;
    if (str != size_string) {
        cout << "Size unknown aborting...";
        exit(1);
    }
    file >> str;
    if (str != matrix_string) {
        cout << "Matrix not present aborting...";
        exit(1);
    }
    float value;
    for (ulong i = 0; i < size; ++i) {
        matrix.emplace_back(vector<float>());
        for (ulong j = 0; j < size; ++j) {
            file >> value;
            matrix[i].emplace_back(value);
        }
    }
    file >> str;
    if (str != terms_string) {
        cout << "Terms not present aborting...";
        exit(1);
    }
    for (ulong j = 0; j < size; ++j) {
        file >> value;
        terms.emplace_back(value);
    }
    file >> str;
    if (str != solution_string) {
        cout << "Solution not present aborting...";
        exit(1);
    }
    for (ulong j = 0; j < size; ++j) {
        file >> value;
        solution.emplace_back(value);
    }
    file.close();
}


void print() {
    cout << matrix_string << ' ' << endl;
    for (int i = 0; i < size; ++i) {
        for (int j = 0; j < size; ++j) {
            cout << matrix[i][j] << ' ';
        }
        cout << endl;
    }
    cout << terms_string << ' ' << endl;
    for (int i = 0; i < size; ++i) {
        cout << terms[i] << ' ';
    }
    cout << endl;
    cout << solution_string << ' ' << endl;
    for (int i = 0; i < size; ++i) {
        cout << solution[i] << ' ';
    }
    cout << endl;
}

void print_solution(const vector<float> &solution, const string &name) {
    cout << name << ' ';
    for (auto &sol: solution) {
        cout << sol << " ";
    }
    cout << endl;
}

ulong workers = 8;
ulong max_iterations = 1000;
float tolerance = 0.f;

int main(const int argc, const char *argv[]) {
    std::ofstream out("results.txt");
    std::cout.rdbuf(out.rdbuf());
    if (argc < 3) {
        cout << "Please insert at least one file name, the number of workers and the number of max_iterations" << endl;
        exit(1);
    }
    workers = (ulong) strtol(argv[1], nullptr, 10);
    max_iterations = (ulong) strtol(argv[2], nullptr, 10);
    tolerance = strtof(argv[3], nullptr);

    for (int arg = 4; arg < argc; ++arg) {
        cout << "RUN: -----------------> " << argv[arg] << endl;
        parse_input(argv[arg]);
        auto test = new float *[size];
        for (ulong i = 0; i < size; ++i) {
            test[i] = &matrix[i][0];
        }
        cout << "workers " << workers << endl;

        auto start = Time::now();
        auto serial_solution = serial_jacobi(matrix, terms, max_iterations, tolerance);
        auto end = Time::now();

        dsec serial_solution_time = end - start;
        cout << "serial jacobi | total time: " << serial_solution_time.count() << endl;
        print_solution(serial_solution, "serial jacobi | solution: ");


        start = Time::now();
        auto omp_solution = jacobi_omp(matrix, terms, max_iterations, tolerance, workers);
        end = Time::now();

        dsec omp_solution_time = end - start;
        cout << "openmp jacobi | total time: " << omp_solution_time.count() << endl;
        print_solution(omp_solution, "openmp jacobi | solution: ");


        start = Time::now();
        auto par_for_solution = jacobi_par_for(matrix, terms, max_iterations, tolerance, workers);
        end = Time::now();

        dsec parallel_for_solution_time = end - start;
        cout << "parallel for | total time: " << parallel_for_solution_time.count() << endl;
        print_solution(par_for_solution, "parallel for | solution: ");

        start = Time::now();
        auto thread_solution = jacobi_thread(matrix, terms, max_iterations, tolerance, workers);
        end = Time::now();

        dsec thread_solution_time = end - start;
        cout << "thread jacobi | total time: " << thread_solution_time.count() << endl;
        print_solution(thread_solution, "thread jacobi | solution: ");

        start = Time::now();
        auto map_solution = jacobi_map(test, &terms[0], size, max_iterations, tolerance, workers);
        end = Time::now();

        dsec map_solution_time = end - start;
        cout << "map jacobi | total time: " << map_solution_time.count() << endl;
        print_solution(vector<float>(map_solution, map_solution + size), "map jacobi | solution: ");

        matrix.clear();
        terms.clear();
        solution.clear();
    }
    return 0;
}
