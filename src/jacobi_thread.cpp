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

    }


    bool termination = false;

    atomic_flag flag;
    spinning_barrier *barrier = nullptr;


    float error = 0.f;

    ulong iteration;

    /**
     * function execute by each thread.
     */
    inline void task(const ulong id, vector<float> &solutions, vector<float> &old_solutions) {
        flag.clear();
        ulong _iteration;
        float myerror;
        for (_iteration = 0; _iteration < iterations; ++_iteration) {
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

#ifdef PARALLEL_REDUCE
#pragma simd
            for (ulong i = works[id].start; i < works[id].stop; ++i) {
                myerror += std::abs(solutions[i] - old_solutions[i]);
            }
            //save the error, spin-lock on the global variable
            while (!flag.test_and_set(std::memory_order_relaxed)) {}
            error += myerror;
            flag.clear(std::memory_order_relaxed);
            barrier->wait();
#endif

            // similar to #pragma omp once, execute it only one time
            if (!flag.test_and_set()) {

#ifndef PARALLEL_REDUCE
                for (ulong i = 0; i < solutions.size(); ++i) {
                    error += std::abs(solutions[i] - old_solutions[i]);
                }
#endif

                error /= (float) solutions.size();
                termination = error <= tolerance;
                std::swap(solutions, old_solutions);
                iteration = _iteration + 1;
                flag.clear();
            }
            // wait the others and check if terminate
            barrier->wait();
            if (termination)break;
        }
        // similar to #pragma omp once, execute it only one time
        if (!flag.test_and_set()) _iteration = _iteration;
    }


}

vector<float> thread_jacobi(const std::vector<std::vector<float>> &_coefficients, const std::vector<float> &_terms,
                            const ulong _iterations, const float _tolerance, const ulong _nWorkers) {

    start_time = Time::now();
    // setting up shared data strucutures
    coefficients = _coefficients;
    terms = _terms;
    iterations = _iterations;
    iteration = _iterations;
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
    for (auto &th: threads) {
        th.join();
    }
    total_time = Time::now();
    print_metrics(iteration, error);


    flag.clear();
    threads.clear();
    delete (barrier);

    return solutions;
}
