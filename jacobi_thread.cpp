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

    void scheduler() {
        nWorkers = min(nWorkers, terms.size());
        ulong slice = terms.size() / nWorkers;
        ulong residual = terms.size() % nWorkers;
//        cout << "size: " << terms.size() << " workers: " << nWorkers << endl;
//        cout << "slice: " << slice << " residual: " << residual << endl;
        for (ulong i = 0; i < nWorkers; ++i) {
            works.emplace_back(job(i * slice, (i + 1) * slice));
        }

        if (residual == 0) return;

        for (ulong i = 0; i < nWorkers; ++i) {
            if (residual-- > 0) {
                if (i != 0) {
                    works[i].start = works[i - 1].stop;
                    works[i].stop++;
                }
                works[i].stop++;
            } else {
                works[i].start = works[i - 1].stop;
                works[i].stop = works[i].start + slice;
            }
        }
//        for (auto &w : works) {
//            cout << "first:" << w.start << " last: " << w.stop << endl;
//        }
    }


    bool termination = false;

    atomic_flag flag;
    spinning_barrier *barrier = nullptr;


    vector<float> errors;
    float error_computed;


    inline void task(const ulong id, vector<float> &solutions, vector<float> &old_solutions) {
        flag.clear();
        for (ulong iteration = 0; iteration < iterations; ++iteration) {
            //calculate solutions
            errors[id] = 0.f;
            for (ulong i = works[id].start; i < works[id].stop; ++i) {
                solutions[i] = solution_find(coefficients[i], old_solutions, terms[i], i);
                errors[id] += abs(solutions[i] - old_solutions[i]);
                old_solutions[i] = solutions[i];
            }
            barrier->wait();
            if (!flag.test_and_set()) {
                error_computed = 0.f;
                for (ulong i = 0; i < nWorkers; ++i) {
                    error_computed += errors[i];
                }
                error_computed /= (float) solutions.size();
                termination = error_computed <= tolerance;
                flag.clear();
            }
            barrier->wait();
            if (termination) {
                if (!flag.test_and_set())iterations = iteration;
                break;
            }
        }
    }


}

vector<float> jacobi_thread(const std::vector<std::vector<float>> &_coefficients, const std::vector<float> &_terms,
                            const ulong _iterations, const float _tolerance, const ulong _nWorkers) {

    coefficients = _coefficients;
    terms = _terms;
    iterations = _iterations;
    tolerance = _tolerance;
    nWorkers = _nWorkers;

    scheduler();
    barrier = new spinning_barrier(nWorkers);
    std::vector<float> solutions(terms.size(), (tolerance - tolerance));
    std::vector<float> old_solutions(terms.size(), (tolerance - tolerance));
    vector<thread> threads;
    for (int i = 0; i < nWorkers; ++i) errors.emplace_back(0.f);
    for (ulong i = 0; i < nWorkers; ++i) {
        threads.emplace_back(thread(task, i, ref(solutions), ref(old_solutions)));
    }
    for (int i = 0; i < nWorkers; ++i) threads[i].join();
    std::cout << "iterations computed: " << iterations << " error: " << error_computed << std::endl;
    delete (barrier);
    return solutions;
}
