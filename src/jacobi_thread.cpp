//
// Created by marco on 9/5/17.
//

#include <atomic>
#include <iostream>
#include "jacobi_thread.hpp"
#include <mutex>
#include <condition_variable>
#include <thread>
#include "barrier.hpp"

using namespace std;

auto static start_time = Time::now();
auto static init_time = Time::now();
auto static total_time = Time::now();

namespace {

    vector<std::vector<float>> coefficients;
    vector<float> terms;
    ulong iterations;
    float tolerance;

    typedef struct job {
        ulong start;
        ulong stop;
        job(ulong start, ulong stop) : start(start), stop(stop) {};
    } job;

    vector<job> works;
    ulong nWorkers;

    /**
     * utility function just to divide the rows and perform a static work balancing
     */
    void scheduler() {
        works.clear();
        nWorkers = min(nWorkers, terms.size());
        ulong slice = terms.size() / nWorkers;
        long residual = terms.size() - (nWorkers * slice);
        for (ulong i = 0; i < nWorkers; ++i) {
            works.emplace_back(job(i * slice, (i + 1) * slice));
        }

        if (residual != 0) {
            for (ulong i = 0; i < nWorkers; ++i) {
                if (i != 0) works[i].start = works[i - 1].stop;
                if (residual-- > 0) works[i].stop = works[i].start + slice + 1;
                else works[i].stop = works[i].start + slice;
            }
        }

//        for (ulong i = 0; i < nWorkers; ++i) {
//            cerr << "start: " << works[i].start << endl;
//            cerr << "stop: " << works[i].stop << endl;
//        }
    }


    bool termination = false;

    atomic_flag flag;
    spinning_barrier *barrier = nullptr;


    float error = 0.f;

    ulong iteration_computed;

    /**
     * function execute by each thread.
     */
    inline void task(const ulong id, vector<float> &solutions, vector<float> &old_solutions) {
        flag.clear();
        ulong iteration;
        float myerror;
        for (iteration = 0; iteration < iterations; ++iteration) {
            //calculate solutions
            myerror = 0.f;
            error = 0.f;
            barrier->wait();
#pragma ivdep
            for (ulong i = works[id].start; i < works[id].stop; ++i) {
                solutions[i] = solution_find(coefficients[i], old_solutions, terms[i], i);
            }
            // wait the others
            barrier->wait();
            //calculate the error and update the solution vectors
#pragma simd
            for (ulong i = works[id].start; i < works[id].stop; ++i) {
                myerror += abs(solutions[i] - old_solutions[i]);
                old_solutions[i] = solutions[i];
            }
            //save the error, spin-lock on the global variable
            while (!flag.test_and_set(std::memory_order_relaxed)) {}
            error += myerror;
            flag.clear(std::memory_order_relaxed);
            barrier->wait();
            // similar to #pragma omp once, execute it only one time
            if (!flag.test_and_set()) {
                error /= (float) solutions.size();
                termination = error <= tolerance;
                flag.clear();
            }
            // wait the others and check if terminate
            barrier->wait();
            if (termination)break;
        }
        // similar to #pragma omp once, execute it only one time
        if (!flag.test_and_set()) iteration_computed = iteration;
    }


}

vector<float> thread_jacobi(const std::vector<std::vector<float>> &_coefficients, const std::vector<float> &_terms,
                            const ulong _iterations, const float _tolerance, const ulong _nWorkers) {

    start_time = Time::now();
    // setting up shared data strucutures
    coefficients = _coefficients;
    terms = _terms;
    iterations = _iterations;
    iteration_computed = _iterations;
    tolerance = _tolerance;
    nWorkers = _nWorkers;
    termination = false;
    scheduler();
    barrier = new spinning_barrier(nWorkers);
    // alllocate the solution vectors
    std::vector<float> old_solutions __attribute__((aligned(64)));
    std::vector<float> solutions __attribute__((aligned(64)));
    // initialize the solution vector
    for (int i = 0; i < coefficients.size(); ++i) {
        old_solutions.emplace_back(tolerance - tolerance);
        solutions.emplace_back(tolerance - tolerance);
    }
    vector<thread> threads;
    // create the threads
    for (ulong i = 0; i < nWorkers; ++i)
        threads.emplace_back(thread(task, i, ref(solutions), ref(old_solutions)));

    init_time = Time::now();
    // wait for the termination
    for (int i = 0; i < nWorkers; ++i) threads[i].join();
    total_time = Time::now();
    std::cout << iterations_computed << iterations << ' ' << error_s << error << std::endl;
    std::cout << initi_time_s << dsec(init_time - start_time).count() << std::endl;
    std::cout << computation_time_s << dsec(total_time - init_time).count() << std::endl;
    flag.clear();
    threads.clear();
    delete (barrier);
    return solutions;
}
