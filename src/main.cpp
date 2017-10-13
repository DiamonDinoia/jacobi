
#include <iostream>

#include "jacobi.hpp"
#include "jacobi_fastflow.hpp"
#include "jacobi_thread.hpp"
#include "jacobi_omp.hpp"

#ifdef WITHOMP

#endif

using namespace std;

const auto size_string = "size:";
const auto matrix_string = "matrix:";
const auto solution_string = "solution:";
const auto terms_string = "terms:";

const auto sequential_string = "sequential";
const auto omp_string = "omp";
const auto thread_string = "thread";
const auto fastflow_string = "fastflow";

static const float range = 10000.f;

enum METHODS {
    SEQUENTIAL,
    OPENMP,
    THREADS,
    FASTFLOW
};

auto method = SEQUENTIAL;



vector<vector<float>> matrix __attribute__((aligned(64)));
vector<float> terms __attribute__((aligned(64)));

ulong size = 1024;
ulong workers = 8;
ulong iterations = 50;
float tolerance = -1.f;


void print_helper() {
    cout << "Usage: " << "main " << "<algorithm> " << "[-w <workers>] " << "[-s <size>] "
         << "[-i <iterations>] " << "[-t <tolerance>]" << endl << endl;
    cout << "The required arguments are:" << endl;
    cout << '\t' << "algorithm" << '\t' << "indicate the algorithm executed, possible values:" << endl;
    cout << "\t\t" << sequential_string << ": sequential jacobi algorithm" << endl;
    cout << "\t\t" << omp_string << ": OpenMP multi-thread implementation of jacobi algorithm" << endl;
    cout << "\t\t" << thread_string << ": plain thread implementation of jacobi algorithm" << endl;
    cout << "\t\t" << fastflow_string << ": fastflow implementation of jacobi algorithm" << endl;
    cout << "The optional arguments are:" << endl;
    cout << '\t' << "[-w]" << '\t' << "number of threads used, default 8" << endl;
    cout << '\t' << "[-s]" << '\t' << "size of the matrix, default 1024" << endl;
    cout << '\t' << "[-i]" << '\t' << "iteration performed, default 50" << endl;
    cout << '\t' << "[-t]" << '\t' << "error tolerated, default -1" << endl;
    flush(cout);
    exit(EINVAL);
}

void parse_args(const int argc, char *const argv[]) {
    if (argc < 2) {
        print_helper();
    }

    string arg = std::string(argv[1]);

    if (arg == sequential_string) method = SEQUENTIAL;
    else if (arg == omp_string) method = OPENMP;
    else if (arg == thread_string) method = THREADS;
    else if (arg == fastflow_string) method = FASTFLOW;
    else print_helper();

    errno = 0;
    int c;
    while ((c = getopt(argc, argv, "w:s:i:t:")) != -1) {
        switch (c) {
            case 'w':
                workers = static_cast<ulong> (strtol(optarg, nullptr, 10));
                break;
            case 's':
                size = static_cast<ulong> (strtol(optarg, nullptr, 10));
                break;
            case 'i':
                iterations = static_cast<ulong> (strtol(optarg, nullptr, 10)) - 1;
                break;
            case 't':
                tolerance = stof(optarg);
                break;
            default:
                print_helper();
        }
    }
    if (errno) {
        print_helper();
    }
}

int main(const int argc, char *const argv[]) {
//    std::ofstream out("results.txt");
//    std::cout.rdbuf(out.rdbuf());

    parse_args(argc, argv);

    srand(42);

    vector<float> solution;

    generate_diagonal_dominant_matrix(size, matrix, -range, range);
    generate_vector(size, terms, -range, range);

    cout << "matrix_size: " << size << endl;

    cout << "algorithm: ";
    switch (method) {
        case SEQUENTIAL:
            cout << sequential_string << endl;
            solution = serial_jacobi(matrix, terms, iterations, tolerance);
            break;
        case THREADS:
            cout << thread_string << endl;
            solution = thread_jacobi(matrix, terms, iterations, tolerance, workers);
            break;
        case OPENMP:
            cout << omp_string << endl;
            solution = omp_jacobi(matrix, terms, iterations, tolerance, workers);
            break;
        case FASTFLOW:
            cout << fastflow_string << endl;
            solution = fastflow_jacobi(matrix, terms, iterations, tolerance, workers);
            break;
    }


    return 0;
}
