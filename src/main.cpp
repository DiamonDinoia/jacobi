
#include <iostream>

#include "jacobi.hpp"
#include "jacobi_fastflow.hpp"
#include "jacobi_thread.hpp"

#ifdef WITHOMP

#include "jacobi_omp.hpp"

#endif

using namespace std;


vector<vector<float>> matrix __attribute__((aligned(64)));
vector<float> terms __attribute__((aligned(64)));


ulong size = 1024;
ulong workers = 8;
ulong iterations = 50;
float tolerance = 0;

char *filename;

ofstream outfile;
auto to_csv = false;
auto debug = false;
auto flag = false;
unsigned int seed = 42;

void print_helper() {
    cout << "Usage: " << "main " << "<algorithm> " << "[-w <workers>] " << "[-s <size>]  << [-p <filename>]"
         << "[-i <iterations>] " << "[-t <tolerance> " << "[-c <seed>] >>" << endl << endl;
    cout << "The required arguments are:" << endl;
    cout << '\t' << "algorithm" << '\t' << "indicate the algorithm executed, possible values:" << endl;
    cout << "\t\t" << sequential_string << ": sequential jacobi algorithm" << endl;
#ifdef WITHOMP
    cout << "\t\t" << omp_string << ": OpenMP multi-thread implementation of jacobi algorithm" << endl;
#endif
    cout << "\t\t" << thread_string << ": plain thread implementation of jacobi algorithm" << endl;
    cout << "\t\t" << fastflow_string << ": fastflow implementation of jacobi algorithm" << endl;
    cout << "The optional arguments are:" << endl;
    cout << '\t' << "[-w]" << '\t' << "number of threads used, default 8" << endl;
    cout << '\t' << "[-s]" << '\t' << "size of the matrix, default 1024" << endl;
    cout << '\t' << "[-i]" << '\t' << "iteration performed, default 50" << endl;
    cout << '\t' << "[-t]" << '\t' << "error tolerated, default -1" << endl;
    cout << '\t' << "[-c]" << '\t' << "seed used to generate the matrix, default 42" << endl;
    flush(cout);
    if (flag) exit(0);
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
    while ((c = getopt(argc, argv, "w:s:i:t:p:c:dh")) != -1) {
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
            case 'p':
                to_csv = true;
                filename = new char[strlen(optarg)];
                strcpy(filename, optarg);
                break;
            case 'd':
                debug = true;
                break;
            case 'c':
                seed = static_cast<unsigned int> (strtol(optarg, nullptr, 10));
                break;
            case 'h':
                flag = true;
                print_helper();
            default:
                print_helper();
        }
    }
    if (errno) {
        print_helper();
    }
}

int main(const int argc, char *const argv[]) {

    parse_args(argc, argv);

    srand((int) seed);

    vector<float> solution;

    if (to_csv) {
        if (!ifstream(filename)) {
            outfile.open(filename);
            outfile << algorithm << ',' << matrix_size << ',' << nworkers << ',' << init_time_s << ','
                    << computation_time_s
                    << ',' << iterations_computed << ',' << error_s << std::endl;
        } else outfile.open(filename, ios::app);
    }

    generate_diagonal_dominant_matrix(size, matrix, -range, range);
    generate_vector(size, terms, -range, range);

    cout << matrix_size << ' ' << size << endl;
    cout << nworkers << ' ' << workers << endl;
    cout << algorithm << ' ';

    switch (method) {
        case SEQUENTIAL:
            cout << sequential_string << endl;
            if (to_csv)
                outfile << sequential_string << ',' << size << ',' << 1 << ',';
            solution = serial_jacobi(matrix, terms, iterations, tolerance, outfile);
            break;
        case THREADS:
            cout << thread_string << endl;
            if (to_csv)
                outfile << thread_string << ',' << size << ',' << workers << ',';
            solution = thread_jacobi(matrix, terms, iterations, tolerance, workers, outfile);
            break;
#ifdef WITHOMP
        case OPENMP:
            cout << omp_string << endl;
            if (to_csv)
                outfile << omp_string << ',' << size << ',' << workers << ',';
            solution = omp_jacobi(matrix, terms, iterations, tolerance, workers, outfile);
            break;
#endif
        case FASTFLOW:
            cout << fastflow_string << endl;
            if (to_csv)
                outfile << fastflow_string << ',' << size << ',' << workers << ',';
            solution = fastflow_jacobi(matrix, terms, iterations, tolerance, workers, outfile);
            break;
    }

    if (debug) {
        cout << "solution: ";
        print_solution(solution);
        float error = check_error(matrix, terms, solution);
        cout << "error: " << error << endl;
    }

    outfile.flush();
    outfile.close();
    delete (filename);
    return 0;
}
