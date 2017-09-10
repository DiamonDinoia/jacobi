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
    ulong max_iterations;
    float tolerance;

    typedef struct job {
        ulong start;
        ulong stop;
        job(ulong start, ulong stop) : start(start), stop(stop) {};
    } job;

    vector<job> works;
    ulong nWorkers;

    void scheduler() {
        works.clear();
        nWorkers = min(nWorkers, terms.size());
//        cerr << "workers " << nWorkers << endl;
//        cerr << "size " << terms.size() << endl;
        ulong slice = terms.size() / nWorkers;
        long residual = terms.size() - (nWorkers * slice);
//        cerr << "slice " << slice << endl;
//        cerr << "residual " << residual << endl;
        for (ulong i = 0; i < nWorkers; ++i) {
            works.emplace_back(job(i * slice, (i + 1) * slice));
        }

        if (residual == 0) {
//            for (int i = 0; i < nWorkers; ++i) {
//                cerr << "start " << works[i].start << endl;
//                cerr << "stop " << works[i].stop << endl;
//            }
            return;
        }
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
//        for (int i = 0; i < nWorkers; ++i) {
//            cerr << "start " << works[i].start << endl;
//            cerr << "stop " << works[i].stop << endl;
//        }

    }


    bool termination = false;

    atomic_flag flag;
    spinning_barrier *barrier = nullptr;


    float errors = 0.f;

    ulong iteration_computed;

    inline void task(const ulong id, vector<float> &solutions, vector<float> &old_solutions) {
        flag.clear();
        ulong iteration;
        float myerror;
        for (iteration = 0; iteration < max_iterations; ++iteration) {
            //calculate solutions
            myerror = 0.f;
            errors = 0.f;
            barrier->wait();
#pragma ivdep
            for (ulong i = works[id].start; i < works[id].stop; ++i) {
                solutions[i] = solution_find(coefficients[i], old_solutions, terms[i], i);
                myerror += abs(solutions[i] - old_solutions[i]);
                old_solutions[i] = solutions[i];
            }
            while (!flag.test_and_set(std::memory_order_relaxed)) {}
            errors += myerror;
            flag.clear();
            barrier->wait();
            if (!flag.test_and_set()) {
                errors /= (float) solutions.size();
                termination = errors <= tolerance;
                flag.clear();
            }
            barrier->wait();
            if (termination)break;
        }
        if (!flag.test_and_set()) iteration_computed = iteration;
    }


}

vector<float> jacobi_thread(const std::vector<std::vector<float>> &_coefficients, const std::vector<float> &_terms,
                            const ulong _iterations, const float _tolerance, const ulong _nWorkers) {

    coefficients = _coefficients;
    terms = _terms;
    max_iterations = _iterations;
    iteration_computed = _iterations;
    tolerance = _tolerance;
    nWorkers = _nWorkers;
    termination = false;

    scheduler();
    barrier = new spinning_barrier(nWorkers);
    std::vector<float> solutions(terms.size(), (tolerance - tolerance));
    std::vector<float> old_solutions(terms.size(), (tolerance - tolerance));
    vector<thread> threads;
    for (ulong i = 0; i < nWorkers; ++i)
        threads.emplace_back(thread(task, i, ref(solutions), ref(old_solutions)));
    auto start = Time::now();
    for (int i = 0; i < nWorkers; ++i) threads[i].join();
    auto end = Time::now();
    std::cout << "thread jacobi | iterations computed: " << iteration_computed << " error: " << errors << std::endl;
    std::cout << "thread jacobi | computation time: " << dsec(end - start).count() << std::endl;
    flag.clear();
    threads.clear();
    delete (barrier);
    return solutions;
}
